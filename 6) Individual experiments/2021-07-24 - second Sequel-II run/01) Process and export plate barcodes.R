### 3rd April 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("Biostrings") # For taking the reverse complement




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
experiments_directory  <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory         <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
file_input_directory   <- file.path(file_directory, "2) Input")
intermediate_directory <- file.path(file_directory, "4) Intermediate files")
R_objects_directory    <- file.path(file_directory, "3) R objects")
metadata_directory     <- file.path(file_input_directory, "Metadata")

s2r1_directory         <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r1_meta_directory    <- file.path(s2r1_directory, "2) Input", "Metadata")

barcodes_pb_file       <- file.path(s2r1_meta_directory, "Sequel_16_barcodes_v3.fasta")
samples_file           <- file.path(metadata_directory, "EtOH_beads_library_purification_plate_layout_AB.xlsx")
barcoded_adapters_file <- file.path(metadata_directory, "Oligo-Ordering-Sheet-for-Microbial-Mplex_v2.xlsx")





# Load data ---------------------------------------------------------------

load(file.path(s2r1_directory, "3) R objects", "01) Process and export plate barcodes.RData"))

run1_plates_df <- plates_df
rm(plates_df)




# Read in data ------------------------------------------------------------

barcodes_pb <- readDNAStringSet(barcodes_pb_file)
samples_df <- data.frame(read_excel(samples_file),
                         stringsAsFactors = FALSE,
                         check.names = FALSE
                         )
adapters_df <- data.frame(read_excel(barcoded_adapters_file, skip = 4),
                         stringsAsFactors = FALSE,
                         check.names = FALSE
                         )



# Tidy samples data -------------------------------------------------------

plate_names_mat <- as.matrix(samples_df[1:8, as.character(1:12)])
numbers_96wp_mat <- as.matrix(samples_df[11:18, as.character(1:12)])
dimnames(numbers_96wp_mat) <- NULL
mode(numbers_96wp_mat) <- "integer"

all_barcodes_df <- samples_df[25:120, 2:4]
names(all_barcodes_df) <- samples_df[24, 2:4]
row.names(all_barcodes_df) <- NULL
all_barcodes_df[["sample"]] <- as.integer(all_barcodes_df[["sample"]])

pool1_mat <- as.matrix(samples_df[126:133, 3:14])
dimnames(pool1_mat) <- NULL
mode(pool1_mat) <- "integer"

pool2_mat <- as.matrix(samples_df[140:147, 3:14])
dimnames(pool2_mat) <- NULL
mode(pool2_mat) <- "integer"




# Perform checks ----------------------------------------------------------

are_same_barcodes <- all_barcodes_df[["barcode"]][1:48] == all_barcodes_df[["barcode"]][49:96]
table(are_same_barcodes)
all_barcodes_df[!(rep(are_same_barcodes, 2)), ]





# Assemble information on each of the plates ------------------------------

stopifnot(identical(as.integer(t(numbers_96wp_mat)), seq_len(96)))

matches_vec <- match(seq_len(96), all_barcodes_df[["sample"]])

plates_df <- data.frame(
  "Run2_pool"        = NA,
  "Plate_name"       = as.character(t(plate_names_mat)),
  "Plate_number"     = NA,
  "Number_96wp"      = seq_len(96),
  "Position_96wp"    = all_barcodes_df[["well position"]][matches_vec],
  "Barcode_ID"       = all_barcodes_df[["barcode"]][matches_vec],
  "Barcode_sequence" = NA,
  stringsAsFactors = FALSE
)

plates_df[["Plate_name"]][plates_df[["Plate_name"]] == "HO_5 AND HO_53"] <- "HO_5 and HO_53"

pool1_df <- plates_df[plates_df[["Number_96wp"]] %in% pool1_mat, ]
pool1_df[["Run2_pool"]] <- 1L
pool1_df[["Plate_name"]][pool1_df[["Plate_name"]] == "Int.Ctl_20cycles"] <- "IntCtl_pool1"

pool2_df <- plates_df[plates_df[["Number_96wp"]] %in% pool2_mat, ]
pool2_df[["Run2_pool"]] <- 2L
pool2_df[["Plate_name"]][pool2_df[["Plate_name"]] == "Int.Ctl_20cycles"] <- "IntCtl_pool2"

plates_df <- rbind.data.frame(pool1_df,
                              pool2_df,
                              stringsAsFactors = FALSE,
                              make.row.names = FALSE
                              )
plates_df <- plates_df[order(plates_df[["Plate_name"]]), ]

plate_name_splits <- strsplit(plates_df[["Plate_name"]], "[ _]")
plate_numbers_vec <- suppressWarnings(as.integer(sapply(plate_name_splits, "[[", 2)))
new_order <- order(sapply(plate_name_splits, "[[", 1), plate_numbers_vec)
plates_df <- plates_df[new_order, ]
row.names(plates_df) <- NULL

plates_df[["Colony_picked"]] <- grepl("IntCtl", plates_df[["Plate_name"]], fixed = TRUE)




# Strip the barcodes ------------------------------------------------------

pacbio_adapters <- sub("/5Phos/", "", adapters_df[[2]], fixed = TRUE)
fwd_pacbio_barcodes <- substr(pacbio_adapters, 1, 16)
rev_pacbio_barcodes <- substr(pacbio_adapters, nchar(pacbio_adapters) - 16, nchar(pacbio_adapters) - 1)
stopifnot(identical(fwd_pacbio_barcodes,
                    as.character(reverseComplement(DNAStringSet(rev_pacbio_barcodes)))
                    ))




# Add the barcode sequences -----------------------------------------------

PB_matches_vec <- match(plates_df[["Barcode_ID"]], adapters_df[["name"]])
full_other_IDs <- paste0(plates_df[["Barcode_ID"]], "_BAK8A_OA")
other_matches_vec <- match(full_other_IDs, names(barcodes_pb))

are_PB <- !(is.na(PB_matches_vec))
are_other <- !(is.na(other_matches_vec))
stopifnot(all(are_PB + are_other) == 1L)

plates_df[["Barcode_sequence"]] <- ifelse(are_PB,
                                          rev_pacbio_barcodes[PB_matches_vec],
                                          substr(as.character(barcodes_pb)[other_matches_vec], 1, 16) # Omit the final T
                                          )
plates_df[["Barcode_ID"]][are_other] <- full_other_IDs[are_other]




# Number the plates -------------------------------------------------------

plates_df[["Plate_number"]] <- seq_len(nrow(plates_df)) + max(run1_plates_df[["Plate_number"]])





# Prepare barcodes for export ---------------------------------------------

barcode_fasta_titles <- paste0(">", plates_df[["Barcode_ID"]])

barcodes_fastas <- lapply(seq_len(nrow(plates_df)), function(x) {
  c(barcode_fasta_titles[[x]], paste0(plates_df[["Barcode_sequence"]][[x]], "T"))
})





# Export barcodes ---------------------------------------------------------

write.table(unlist(barcodes_fastas[plates_df[["Run2_pool"]] == 1]),
            file = file.path(intermediate_directory, "pool1_plate_barcodes.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )
write.table(unlist(barcodes_fastas[plates_df[["Run2_pool"]] == 2]),
            file = file.path(intermediate_directory, "pool2_plate_barcodes.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )




# Save data ---------------------------------------------------------------

save(list = "plates_df",
     file = file.path(R_objects_directory, "01) Process and export plate barcodes.RData")
     )






