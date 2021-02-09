### 29 August 2020 ###






# Import packages and source code -----------------------------------------

library("Rsamtools")
library("readxl")
library("ShortRead")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory  <- file.path(file_directory, "1) Input")
R_objects_directory   <- file.path(file_directory, "2) R objects")

first_plate_directory <- file.path(file_input_directory, "Raw data", "384-well-1_2_B01")
subreads_bam_file     <- file.path(first_plate_directory, "m54073_190912_185203.subreads.bam")
subreads_stats_file   <- file.path(first_plate_directory, "m54073_190912_185203.subreads.stats.csv")

metadata_directory    <- file.path(file_input_directory, "Metadata")
sg_sequences_file     <- file.path(metadata_directory, "1-384 sgRNA summary.xlsm")
barcodes_file         <- file.path(metadata_directory, "Barcoded primers.xlsx")



# Read in data ------------------------------------------------------------

sg_sequences_df <- as.data.frame(read_excel(sg_sequences_file, sheet = 2),
                                 stringsAsFactors = FALSE, check.names = FALSE
                                 )
barcodes_df <- as.data.frame(read_excel(barcodes_file, col_names = FALSE),
                             stringsAsFactors = FALSE, check.names = FALSE
                             )


ccs3_bam <- scanBam(file.path(file_input_directory,
                              "Delete this",
                              "CCS3_m54073_190912_185203.ccs.bam"
                              ))

ccs5_bam <- scanBam(file.path(file_input_directory,
                              "Delete this",
                              "CCS5_m54073_190912_185203.ccs.bam"
                              ))


ccs3_fastq <- readFastq(dirPath = file.path(file_input_directory, "Delete this"),
                        pattern = "CCS3_m54073_190912_185203.Q20.fastq",
                        qualityType = "FastqQuality"
                        )
ccs5_fastq <- readFastq(dirPath = file.path(file_input_directory, "Delete this"),
                        pattern = "CCS5_m54073_190912_185203.Q20.fastq"
                        )


ccs3_stats_df <- read.csv(file.path(file_input_directory, "Raw data", "CCS3",
                                    "ccs_statistics.csv"
                                    ),
                          quote = "", stringsAsFactors = FALSE
                          )



subreads_stats_df <- read.csv(subreads_stats_file, stringsAsFactors = FALSE)

subreads_bam <- scanBam(subreads_bam_file)






# Try stuff ---------------------------------------------------------------

ccs3_sequences_vec <- as.character(ccs3_bam[[1]][["seq"]])
ccs5_sequences_vec <- as.character(ccs5_bam[[1]][["seq"]])


ccs3_fastq_sequences_vec <- as.character(sread(ccs3_fastq))
ccs3_fastq_IDs_vec <- as.character(id(ccs3_fastq))
ccs3_fastq_qualities_object <- quality(ccs3_fastq)

ccs3_score_per_base <- alphabetScore(ccs3_fastq_qualities_object) / nchar(ccs3_fastq_sequences_vec)




ccs5_fastq_sequences_vec <- as.character(sread(ccs5_fastq))
ccs5_fastq_IDs_vec <- as.character(id(ccs5_fastq))
ccs5_fastq_qualities_object <- quality(ccs5_fastq)

ccs5_fastq_holes <- sapply(strsplit(ccs5_fastq_IDs_vec, "/", fixed = TRUE), "[", 2)


ccs5_score_per_base <- alphabetScore(ccs5_fastq_qualities_object) / nchar(ccs5_fastq_sequences_vec)



ccs5_to_ccs3_matches <- match(ccs5_fastq_IDs_vec, ccs3_fastq_IDs_vec)
all.equal(ccs3_score_per_base[ccs5_to_ccs3_matches], ccs5_score_per_base)



are_identical <- mapply(identical, ccs3_fastq_sequences_vec, ccs3_sequences_vec)

identical(sort(ccs3_sequences_vec),
          sort(ccs3_fastq_sequences_vec)
          )


holes_vec_ccs3 <- sapply(strsplit(ccs3_bam[[1]][["qname"]], "/", fixed = TRUE), "[", 2)
holes_vec_ccs5 <- sapply(strsplit(ccs5_bam[[1]][["qname"]], "/", fixed = TRUE), "[", 2)


holes_vec_ccs3stats <- sapply(strsplit(ccs3_stats_df[["qname"]], "/", fixed = TRUE), "[", 2)

ccs3_stats_scores <- ccs3_stats_df[["readscore"]][match(holes_vec_ccs3, holes_vec_ccs3stats)]
cor.test(ccs3_stats_scores, ccs3_score_per_base)


stopifnot(identical(sort(holes_vec_ccs3stats), sort(holes_vec_ccs3)))
hist(ccs3_stats_df[["readscore"]])


table(ccs3_sequences_vec %in% ccs5_sequences_vec)
table(ccs5_sequences_vec %in% ccs3_sequences_vec)

table(holes_vec_ccs5 %in% holes_vec_ccs3)
table(holes_vec_ccs3 %in% holes_vec_ccs5)


table(holes_vec_ccs3 %in% subreads_stats_df[["hole"]])
table(holes_vec_ccs5 %in% subreads_stats_df[["hole"]])


ccs3_sequences_vec[holes_vec_ccs3 == "4194470"]
ccs5_sequences_vec[holes_vec_ccs5 == "4194470"]



holes_vec <- c(
  "69207006",
  "4850529",
  "57344543",
  "72352113",
  "4981175",
  "60621418",
  "51511917",
  "7602285",
  "23200053",
  "71368998"
)


ccs5_score_per_base[ccs5_fastq_holes %in% c("23200053", "71368998")]
ccs5_sequences_vec[ccs5_fastq_holes %in% holes_vec]







ccs5_example_aligned_bam <- scanBam(file.path(file_input_directory,
                              "Delete this",
                              "CCS5_lima.0_16.sorted.bam"
                              ))







# Define functions --------------------------------------------------------

SearchForBarcodes <- function(barcode,
                              start_with = FALSE,
                              end_with = FALSE,
                              reverse_complement = FALSE
                              ) {

  stopifnot(exists("sequences_vec")) # requires sequences_vec in the global environment

  barcode <- toupper(barcode)
  assign("delete_barcode", barcode, envir = globalenv())
  if (reverse_complement) {
    barcode <- as.character(reverseComplement(DNAStringSet(barcode)))
  }
  if (start_with) {
    barcode <- paste0("^", barcode)
  }
  if (end_with) {
    barcode <- paste0(barcode, "$")
  }
  results_vec <- grepl(barcode,
                       sequences_vec,
                       fixed = !(start_with || end_with)
                       )
  return(results_vec)
}





# Explore data ------------------------------------------------------------

names(subreads_bam[[1]])
length(unique(subreads_bam[[1]][["qname"]]))
length(unique(subreads_bam[[1]][["flag"]]))
stopifnot(all(is.na(subreads_bam[[1]][["rname"]])))
stopifnot(all(is.na(subreads_bam[[1]][["strand"]])))
stopifnot(all(is.na(subreads_bam[[1]][["pos"]])))
stopifnot(all(is.na(subreads_bam[[1]][["qwidth"]])))
stopifnot(all(is.na(subreads_bam[[1]][["mapq"]])))
stopifnot(all(is.na(subreads_bam[[1]][["cigar"]])))
stopifnot(all(is.na(subreads_bam[[1]][["mrnm"]])))
stopifnot(all(is.na(subreads_bam[[1]][["mpos"]])))
stopifnot(all(is.na(subreads_bam[[1]][["isize"]])))
length(unique(subreads_bam[[1]][["seq"]]))
length(unique(subreads_bam[[1]][["qual"]]))




# Look for barcodes in the reads ------------------------------------------

sequences_vec <- as.character(subreads_bam[[1]][["seq"]])

goo

are_row_barcodes    <- grepl("2PF_", barcodes_df[[1]], fixed = TRUE)
are_column_barcodes <- grepl("2PR_", barcodes_df[[1]], fixed = TRUE)

row_barcodes    <- barcodes_df[[2]][are_row_barcodes]
column_barcodes <- barcodes_df[[2]][are_column_barcodes]

start_with_row_list   <- lapply(row_barcodes, SearchForBarcodes, start_with = TRUE)
end_with_rev_row_list <- lapply(row_barcodes, SearchForBarcodes, end_with = TRUE, reverse_complement = TRUE)

start_with_column_list   <- lapply(column_barcodes, SearchForBarcodes, start_with = TRUE)
end_with_rev_column_list <- lapply(column_barcodes, SearchForBarcodes, end_with = TRUE, reverse_complement = TRUE)

within_row_list <- lapply(row_barcodes, SearchForBarcodes)
within_rev_row_list <- lapply(row_barcodes, SearchForBarcodes, reverse_complement = TRUE)

within_column_list <- lapply(column_barcodes, SearchForBarcodes)
within_rev_column_list <- lapply(column_barcodes, SearchForBarcodes, reverse_complement = TRUE)










cell_number <- 1L
cell_vec <- rep(NA_integer_, length(sequences_vec))
row_barcode_is_reversed <- rep(NA, length(sequences_vec))
for (row_index in seq_along(row_barcodes)) {
  for (column_index in seq_along(column_barcodes)) {
    are_this_cell_forward  <- start_with_row_list[[row_index]] & end_with_rev_column_list[[column_index]]
    are_this_cell_reversed <- end_with_rev_row_list[[row_index]] & start_with_column_list[[column_index]]
    row_barcode_is_reversed[are_this_cell_reversed] <- TRUE
    row_barcode_is_reversed[are_this_cell_forward] <- FALSE
    cell_vec[are_this_cell_forward | are_this_cell_reversed] <- cell_number
    cell_number <- cell_number + 1L
  }
}



contain_row_barcode_list <- lapply(seq_along(row_barcodes),
                                   function(x) start_with_row_list[[x]] | end_with_rev_row_list[[x]]
                                   )
weird_row_barcode_list <- lapply(seq_along(row_barcodes),
                                 function(x) start_with_row_list[[x]] & end_with_rev_row_list[[x]]
                                 )
vapply(contain_row_barcode_list, sum, integer(1))
vapply(weird_row_barcode_list, sum, integer(1))



contain_column_barcode_list <- lapply(seq_along(column_barcodes),
                                      function(x) start_with_column_list[[x]] | end_with_rev_column_list[[x]]
                                      )

weird_column_barcode_list <- lapply(seq_along(column_barcodes),
                                    function(x) start_with_column_list[[x]] & end_with_rev_column_list[[x]]
                                    )

vapply(contain_column_barcode_list, sum, integer(1))
vapply(weird_column_barcode_list, sum, integer(1))






contain_within_row_barcode_list <- lapply(seq_along(row_barcodes),
                                          function(x) within_row_list[[x]] | within_rev_row_list[[x]]
                                          )
contain_within_column_barcode_list <- lapply(seq_along(column_barcodes),
                                             function(x) within_column_list[[x]] | within_rev_column_list[[x]]
                                             )





ReverseString <- function(string) {
  intToUtf8(rev(utf8ToInt((string))))
}


FourVariants <- function(DNA_string) {
  complement_string <- as.character(reverseComplement(DNAStringSet(DNA_string)))
  c("Original"                    = DNA_string,
    "Reverse complement"          = complement_string,
    "Reversed original"           = ReverseString(DNA_string),
    "Reversed reverse complement" = ReverseString(complement_string)
    )
}


set.seed(1)
sort(sample(seq_along(sequences_vec), 10))


sequences_vec[[2221565]]









# Save data ---------------------------------------------------------------

save(list = c("start_with_row_list", "end_with_rev_row_list",
              "start_with_column_list", "end_with_rev_column_list",
              "within_row_list", "within_rev_row_list",
              "within_column_list", "within_rev_column_list"
              ),
     file = file.path(R_objects_directory,
                      "1) Process PacBio data - barcode searches.RData"
                      )
     )










names(delete_this_bam[[1]])
length(unique(delete_this_bam[[1]][["qname"]]))
length(unique(subreads_bam[[1]][["flag"]]))
stopifnot(all(is.na(subreads_bam[[1]][["rname"]])))
stopifnot(all(is.na(subreads_bam[[1]][["strand"]])))
stopifnot(all(is.na(subreads_bam[[1]][["pos"]])))
stopifnot(all(is.na(subreads_bam[[1]][["qwidth"]])))
stopifnot(all(is.na(subreads_bam[[1]][["mapq"]])))
stopifnot(all(is.na(subreads_bam[[1]][["cigar"]])))
stopifnot(all(is.na(subreads_bam[[1]][["mrnm"]])))
stopifnot(all(is.na(subreads_bam[[1]][["mpos"]])))
stopifnot(all(is.na(subreads_bam[[1]][["isize"]])))
length(unique(delete_this_bam[[1]][["seq"]]))
length(unique(subreads_bam[[1]][["qual"]]))







