### 5th September 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("Biostrings") # For taking the reverse complement



# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory   <- file.path(file_directory, "2) Input")
intermediate_directory <- file.path(file_directory, "4) Intermediate files")
R_objects_directory    <- file.path(file_directory, "3) R objects")

metadata_directory     <- file.path(file_input_directory, "Metadata")
barcodes_file          <- file.path(metadata_directory, "Barcoded primers.xlsx")




# Read in data ------------------------------------------------------------

barcodes_df <- as.data.frame(read_excel(barcodes_file, col_names = FALSE),
                             stringsAsFactors = FALSE, check.names = FALSE
                             )




# Process barcodes --------------------------------------------------------

TidyBarcodes <- function(barcodes_vec) {
  toupper(substr(barcodes_vec, 1, 10))
}

are_row_barcodes    <- grepl("2PF_", barcodes_df[[1]], fixed = TRUE)
are_column_barcodes <- grepl("2PR_", barcodes_df[[1]], fixed = TRUE)

original_row_barcodes    <- TidyBarcodes(barcodes_df[[2]][are_row_barcodes])
original_column_barcodes <- TidyBarcodes(barcodes_df[[2]][are_column_barcodes])

row_constant_region <- unique(substr(original_row_barcodes, 11, nchar(original_row_barcodes)))
column_constant_region <- unique(substr(original_column_barcodes, 11, nchar(original_column_barcodes)))

row_barcodes    <- TidyBarcodes(original_row_barcodes)
column_barcodes <- TidyBarcodes(original_column_barcodes)




# Create all possible barcode combinations --------------------------------

barcode_combos_df <- expand.grid(seq_along(column_barcodes),
                                 seq_along(row_barcodes),
                                 stringsAsFactors = FALSE,
                                 KEEP.OUT.ATTRS = FALSE
                                 )[, 2:1]
colnames(barcode_combos_df) <- c("Row_number", "Column_number")


barcodes_to_wells_map <- seq_len(384)
names(barcodes_to_wells_map) <- paste0("FwdBC",
                                       barcode_combos_df[["Row_number"]],
                                       "--",
                                       "RevBC",
                                       barcode_combos_df[["Column_number"]]
                                       )




# Create sample names -----------------------------------------------------

sample_name_df <- data.frame(
  "Barcode" = paste0("FwdBC", barcode_combos_df[["Row_number"]], "--",
                     "RevBC", barcode_combos_df[["Column_number"]]
                     ),
  "Name"    = paste0("First384_Well",
                     formatC(seq_len(384), width = 3, flag = "0")
                     )
)
colnames(sample_name_df)[[2]] <- paste0(
  "Bio Sample (delete entire rows of  barcodes not used; allowed characters: ",
  "alphanumeric space dot underscore hyphen. Other characters will be ",
  "automatically removed.)"
)






# Create the barcode name map (PacBio ID to well number) ------------------

barcodes_to_wells_map <- seq_len(384)
names(barcodes_to_wells_map) <- sample_name_df[["Barcode"]]







# Create a vector of the barcode sequences for each well ------------------

row_barcodes_reversed <- as.character(reverseComplement(DNAStringSet(row_barcodes)))

column_bc_vec <- column_barcodes[barcode_combos_df[["Column_number"]]]
row_bc_vec <- row_barcodes_reversed[barcode_combos_df[["Row_number"]]]





# Prepare barcodes for export ---------------------------------------------

row_titles <- paste0("FwdBC", seq_along(row_barcodes))
column_titles <- paste0("RevBC", seq_along(column_barcodes))

export_barcode_titles <- c(row_titles, column_titles)

export_barcode_fasta_titles <- paste0(">", export_barcode_titles)

export_only_barcodes <- c(row_barcodes, column_barcodes)
export_with_constant_regions <- toupper(c(original_row_barcodes, original_column_barcodes))

only_barcodes_fastas <- lapply(seq_along(export_only_barcodes), function(x) {
  c(export_barcode_fasta_titles[[x]], export_only_barcodes[[x]])
})
barcode_and_constant_region_fastas <- lapply(seq_along(export_with_constant_regions), function(x) {
  c(export_barcode_fasta_titles[[x]], export_with_constant_regions[[x]])
})






# Export barcodes ---------------------------------------------------------

write.table(unlist(only_barcodes_fastas),
            file = file.path(intermediate_directory, "barcodes_first384.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            )

write.table(unlist(barcode_and_constant_region_fastas),
            file = file.path(intermediate_directory, "barcodes_first384_with_constant_region.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            )

write.csv(sample_name_df,
          file = file.path(intermediate_directory, "sample_names_first384.csv"),
          quote = FALSE, row.names = FALSE
          )




# Save data ---------------------------------------------------------------

save(list = c("barcode_combos_df",
              "row_barcodes", "column_barcodes",
              "row_constant_region", "column_constant_region",
              "column_bc_vec", "row_bc_vec",
              "barcodes_to_wells_map"
              ),
     file = file.path(R_objects_directory, "01) Process and export barcodes.RData")
     )






