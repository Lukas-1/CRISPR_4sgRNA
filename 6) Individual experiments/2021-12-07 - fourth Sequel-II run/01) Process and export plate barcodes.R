### 7th December 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2r4_directory           <- file.path(experiments_directory, "2021-12-07 - fourth Sequel-II run")

s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
s2r4_R_objects_directory <- file.path(s2r4_directory, "3) R objects")

intermediate_directory   <- file.path(s2r4_directory, "4) Intermediate files")





# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "30) Calculate plate weightings for the next round of sequencing.RData"))

run2_plates_df <- plates_df




# Create the new plates_df ------------------------------------------------

use_columns <- c("Run4_pool", intersect(setdiff(names(run2_plates_df), "Run2_pool"), names(run4_df)))
plates_df <- run4_df[run4_df[["Run4_pool"]] %in% c(4, 5), use_columns]

matches_vec <- match(plates_df[, "Barcode_ID"], run2_plates_df[, "Barcode_ID"])
plates_df[["Barcode_sequence"]] <- run2_plates_df[matches_vec, "Barcode_sequence"]





# Prepare barcodes for export ---------------------------------------------

pool4_plates_df <- plates_df[plates_df[, "Run4_pool"] %in% 4, ]
pool5_plates_df <- plates_df[plates_df[, "Run4_pool"] %in% 5, ]

pool4_barcode_fasta_titles <- paste0(">", pool4_plates_df[["Barcode_ID"]])
pool4_barcodes_fastas <- lapply(seq_len(nrow(pool4_plates_df)), function(x) {
  c(pool4_barcode_fasta_titles[[x]], paste0(pool4_plates_df[["Barcode_sequence"]][[x]], "T"))
})

pool5_barcode_fasta_titles <- paste0(">", pool5_plates_df[["Barcode_ID"]])
pool5_barcodes_fastas <- lapply(seq_len(nrow(pool5_plates_df)), function(x) {
  c(pool5_barcode_fasta_titles[[x]], paste0(pool5_plates_df[["Barcode_sequence"]][[x]], "T"))
})



# Export barcodes ---------------------------------------------------------

write.table(unlist(pool4_barcodes_fastas),
            file = file.path(intermediate_directory, "pool4_plate_barcodes.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )

write.table(unlist(pool5_barcodes_fastas),
            file = file.path(intermediate_directory, "pool5_plate_barcodes.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )




# Save data ---------------------------------------------------------------

save(list = "plates_df",
     file = file.path(s2r4_R_objects_directory, "01) Process and export plate barcodes.RData")
     )



