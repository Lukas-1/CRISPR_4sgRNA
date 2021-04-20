### 3rd April 2021 ###



# Import packages and source code -----------------------------------------

library("readODS")
library("Biostrings") # For taking the reverse complement





# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
file_input_directory   <- file.path(file_directory, "2) Input")
intermediate_directory <- file.path(file_directory, "4) Intermediate files")
R_objects_directory    <- file.path(file_directory, "3) R objects")
metadata_directory     <- file.path(file_input_directory, "Metadata")

barcodes_ods_file      <- file.path(metadata_directory, "barcodes_p3239_24065.ods")
barcodes_pb_file       <- file.path(metadata_directory, "Sequel_16_barcodes_v3.fasta")
samples_file           <- file.path(metadata_directory, "Project_3239_sample.csv")



# Read in data ------------------------------------------------------------

barcodes_ods_df <- read_ods(barcodes_ods_file, skip = 1)
barcodes_pb <- readDNAStringSet(barcodes_pb_file)
samples_df <- read.csv(samples_file, stringsAsFactors = FALSE,
                       check.names = FALSE
                       )




# Assemble barcodes -------------------------------------------------------

ods_tube_splits <- strsplit(barcodes_ods_df[["Sample"]], "/", fixed = TRUE)
order_ID <- unique(sapply(ods_tube_splits, "[[", 1))
ods_tube_numbers <- as.integer(sapply(ods_tube_splits, "[[", 2))

full_fgcz_barcodes <- sub("/5Phos/", "", barcodes_ods_df[["sequence"]], fixed = TRUE)
full_fgcz_barcodes <- gsub(" ", "", full_fgcz_barcodes, fixed = TRUE)

fwd_fgcz_barcodes <- substr(full_fgcz_barcodes, 1, 16)
rev_fgcz_barcodes <- substr(full_fgcz_barcodes, nchar(full_fgcz_barcodes) - 16, nchar(full_fgcz_barcodes) - 1)
are_fgcz <- !(is.na(fwd_fgcz_barcodes))
stopifnot(identical(fwd_fgcz_barcodes[are_fgcz],
                    as.character(reverseComplement(DNAStringSet(rev_fgcz_barcodes[are_fgcz])))
                    ))

combined_barcodes <- ifelse(are_fgcz,
                            fwd_fgcz_barcodes,
                            as.character(barcodes_pb)[barcodes_ods_df[["barcode"]]]
                            )
combined_barcodes <- substr(combined_barcodes, 1, 16) # Omit the final 'T'

sample_matches <- match(barcodes_ods_df[["Sample"]], samples_df[["Tube Id"]])

plates_df <- data.frame(
  "Plate_name"       = samples_df[["Name"]][sample_matches],
  "Plate_ID"         = samples_df[[" Id"]][sample_matches],
  "Plate_number"     = ods_tube_numbers,
  "Barcode_ID"       = barcodes_ods_df[["barcode"]],
  "Barcode_sequence" = combined_barcodes,
  stringsAsFactors   = FALSE
)




# Prepare barcodes for export ---------------------------------------------

barcode_fasta_titles <- paste0(">", plates_df[["Barcode_ID"]])

barcodes_fastas <- lapply(seq_len(nrow(plates_df)), function(x) {
  c(barcode_fasta_titles[[x]], plates_df[["Barcode_sequence"]][[x]])
})





# Export barcodes ---------------------------------------------------------

write.table(unlist(barcodes_fastas),
            file = file.path(intermediate_directory, "plate_barcodes.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )




# Save data ---------------------------------------------------------------

save(list = "plates_df",
     file = file.path(R_objects_directory, "01) Process and export plate barcodes.RData")
     )






