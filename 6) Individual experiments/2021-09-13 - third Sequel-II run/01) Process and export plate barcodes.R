### 13th September 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")

s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")

intermediate_directory   <- file.path(s2r3_directory, "4) Intermediate files")



# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "30) Calculate correction factors to account for mean read counts.RData"))

run2_plates_df <- plates_df



# Create the new plates_df ------------------------------------------------

pool3_names <- extended_df[["Plate_name"]][extended_df[["Run3_pool"]] %in% 3]
are_pool3 <- run2_plates_df[["Plate_name"]] %in% pool3_names

exclude_columns <- c("Number_96wp", "Position_96wp")
plates_df <- run2_plates_df[are_pool3, !(names(run2_plates_df) %in% exclude_columns)]

matches_vec <- match(plates_df[["Plate_name"]], extended_df[["Plate_name"]])

plates_df[["Barcode_ID"]] <- extended_df[["Barcode_ID"]][matches_vec]

matches_vec <- match(plates_df[["Barcode_ID"]], run2_plates_df[["Barcode_ID"]])

plates_df[["Barcode_sequence"]] <- run2_plates_df[["Barcode_sequence"]][matches_vec]



# Prepare barcodes for export ---------------------------------------------

barcode_fasta_titles <- paste0(">", plates_df[["Barcode_ID"]])

barcodes_fastas <- lapply(seq_len(nrow(plates_df)), function(x) {
  c(barcode_fasta_titles[[x]], paste0(plates_df[["Barcode_sequence"]][[x]], "T"))
})



# Export barcodes ---------------------------------------------------------

write.table(unlist(barcodes_fastas),
            file = file.path(intermediate_directory, "pool3_plate_barcodes.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )


# Save data ---------------------------------------------------------------

save(list = "plates_df",
     file = file.path(s2r3_R_objects_directory, "01) Process and export plate barcodes.RData")
     )



