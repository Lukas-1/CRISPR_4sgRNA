### 9th April 2021 ###





# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "03) Creating reference plasmid sequences.R"))




# Define folder paths -----------------------------------------------------

sql2_directory              <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")

p1_R_objects_directory      <- file.path(plate1_directory, "3) R objects")
sql2_R_objects_directory    <- file.path(sql2_directory, "3) R objects")

plate1_metadata_directory   <- file.path(plate1_directory, "2) Input", "Metadata")
tracRNAs_file               <- file.path(plate1_metadata_directory, "tracRNAs_4sg.txt")
promoters_file              <- file.path(plate1_metadata_directory, "promoters_4sg.txt")
full_plasmid_file           <- file.path(plate1_metadata_directory, "reference with 4sg and barcode example.fa")

annotated_plasmid_directory <- file.path(sql2_directory, "2) Input", "Annotated plasmid")

file_output_directory       <- file.path(sql2_directory, "5) Output")
reference_output_directory  <- file.path(file_output_directory, "Reference sequences")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))




# Read in data ------------------------------------------------------------

tracRNAs_vec <- read.table(tracRNAs_file, header = FALSE,
                           stringsAsFactors = FALSE
                           )[[1]]
promoters_vec <- read.table(promoters_file, header = FALSE,
                            stringsAsFactors = FALSE
                            )[[1]]
plasmid_lines_vec <- read.table(full_plasmid_file, header = FALSE, skip = 1,
                                stringsAsFactors = FALSE
                                )[[1]]



# Define functions --------------------------------------------------------

PlasmidForWell <- function(well_number) {
  stopifnot("plasmid_list" %in% ls(envir = globalenv()))
  search_string <- paste0("_well", well_number, ".gbk")
  plasmid_name <- grep(search_string, names(plasmid_list), fixed = TRUE)
  if (length(plasmid_name) == 0) {
    plasmid_name <- "pYJA5_4sg_nested_barcodes.gbk"
  }
  return(plasmid_name)
}




# Read in annotated plasmid sequences -------------------------------------

plasmid_files <- list.files(annotated_plasmid_directory)

plasmid_list <- sapply(plasmid_files, function(x) {
  ReadInPlasmids(file.path(annotated_plasmid_directory, x))
}, simplify = FALSE)





# Define the reference sequences ------------------------------------------

rev_tracRNAs_vec <- as.character(reverseComplement(DNAStringSet(tracRNAs_vec)))

guides_ref_list <- lapply(1:4, function(x) {
  sg_column <- paste0("Sequence_sg", x)
  paste0(library_df[[sg_column]], rev_tracRNAs_vec[[x]], "TTTTT")
})

guides_with_promoters_list <- lapply(1:4, function(x) {
  paste0(toupper(promoters_vec[[x]]), guides_ref_list[[x]])
})

plasmid_string <- toupper(paste0(plasmid_lines_vec, collapse = ""))
plasmid_string <- substr(plasmid_string, 11, nchar(plasmid_string) - 10)
library_seq <- seq_len(nrow(library_df))
plasmids_vec <- vapply(library_seq, function(x) {
  sg_N <- paste0(rep("N", 20), collapse = "")
  sg_sequences <- vapply(1:4, function(y) {
    library_df[[paste0("Sequence_sg", y)]][[x]]
  }, "")
  result_string <- sub(sg_N, sg_sequences[[1]], plasmid_string, fixed = TRUE)
  result_string <- sub(sg_N, sg_sequences[[2]], result_string,  fixed = TRUE)
  result_string <- sub(sg_N, sg_sequences[[3]], result_string,  fixed = TRUE)
  result_string <- sub(sg_N, sg_sequences[[4]], result_string,  fixed = TRUE)
  return(result_string)
}, "")





# Insert barcodes and guides into the annotated plasmids ------------------

plate_matches <- match(library_df[["Plate_name"]], plates_df[["Plate_name"]])
plate_barcodes <- plates_df[["Barcode_sequence"]][plate_matches]
rev_plate_barcodes <- as.character(reverseComplement(DNAStringSet(plate_barcodes)))

column_barcodes <- column_bc_vec[library_df[["Well"]]]
row_barcodes <- row_bc_vec[library_df[["Well"]]]

plasmid_lines_list <- lapply(library_seq, function(x) {
  use_plasmid <- PlasmidForWell(x)
  use_lines <- plasmid_list[[use_plasmid]][["sequence"]]
  sg_seqs <- as.character(library_df[x, paste0("Sequence_sg", 1:4)])
  replace_seqs <- c(plate_barcodes[[x]],
                    column_barcodes[[x]],
                    sg_seqs,
                    row_barcodes[[x]],
                    rev_plate_barcodes[[x]]
                    )
  ReplaceNNNLines(use_lines, replace_seqs)
})





# Export the plain plasmid sequences --------------------------------------

barcoded_plasmids <- paste0(plate_barcodes,
                            "T",
                            column_barcodes,
                            plasmids_vec,
                            row_barcodes,
                            "A",
                            rev_plate_barcodes
                            )
names(barcoded_plasmids) <- library_df[["Combined_ID"]]

fasta_titles <- paste0(">", library_df[["Combined_ID"]])
fasta_list <- lapply(library_seq,
                     function(x) c(fasta_titles[[x]], barcoded_plasmids[[x]], "")
                     )

write.table(unlist(fasta_list),
            file = file.path(file_output_directory, "reference_sequences.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            )





# Export the annotated plasmids -------------------------------------------

for (i in library_seq) {
  file_name <- paste0(library_df[["Combined_ID"]][[i]], "_barcoded_cassette.gbk")
  use_plasmid <- PlasmidForWell(i)
  export_lines <- c(plasmid_list[[use_plasmid]][["preceding"]],
                    plasmid_lines_list[[i]],
                    plasmid_list[[use_plasmid]][["following"]]
                    )
  write.table(x         = export_lines,
              file      = file.path(reference_output_directory, file_name),
              col.names = FALSE,
              row.names = FALSE,
              quote     = FALSE
              )
}





# Save data ---------------------------------------------------------------

save(list = c("guides_ref_list", "guides_with_promoters_list",
              "plasmid_string", "plasmids_vec", "barcoded_plasmids"
              ),
     file = file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData")
     )

save(list = "plasmid_list",
     file = file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - annotated plasmid.RData")
     )




