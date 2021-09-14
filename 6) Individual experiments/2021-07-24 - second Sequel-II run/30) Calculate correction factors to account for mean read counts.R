### 22nd August 2021 ###




# Import packages and source code -----------------------------------------

library("readxl")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
experiments_directory  <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2r2_directory         <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_objects_directory <- file.path(s2r2_directory, "3) R objects")
plate_selection_path   <- file.path(s2r2_directory, "2) Input/Metadata/After sequencing",
                                    "Plates can be excluded for further sequencing_20210819.xlsx"
                                    )
file_output_directory  <- file.path(s2r2_directory, "5) Output", "Tables", "Correction factors")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_objects_directory, "29) Correlate median read count with DNA concentration.RData"))




# Read in data ------------------------------------------------------------

selection_df <- data.frame(read_excel(plate_selection_path),
                           stringsAsFactors = FALSE, check.names = FALSE
                           )




# Define plate selections -------------------------------------------------

good_plates <- grep("^H", selection_df[, 1], value = TRUE)
good_plates <- sub("-", "_", good_plates, fixed = TRUE)

replace_plates <- c("HA_21", "HO_13", "HO_16")




# Merge the plate annotation and sequencing yield summary data ------------

extended_df <- merged_plates_df
extended_df[["Mean_read_count"]] <- extended_df[, "Sum_counts"] / extended_df[, "Number_of_wells"]

new_order <- order(extended_df[, "Run2_pool"], extended_df[, "Plate_number"])
extended_df <- extended_df[new_order, ]
row.names(extended_df) <- NULL




# Calculate rounded correction factors ------------------------------------

are_remaining <- !(extended_df[, "Plate_name"] %in% c(good_plates, replace_plates))

mean_of_means <- mean(extended_df[, "Mean_read_count"][are_remaining])
corrections <- mean_of_means / extended_df[, "Mean_read_count"][are_remaining]

rounded_corrections <- as.integer(round(corrections))
rounded_corrections <- ifelse(rounded_corrections == 0,
                              1L,
                              ifelse(rounded_corrections > 2,
                                     4L,
                                     rounded_corrections
                                     )
                              )

extended_df[["Correction_factor"]][are_remaining] <- rounded_corrections

extended_df[["Is_to_replace"]] <- extended_df[["Plate_name"]] %in% replace_plates
extended_df[["Correction_factor"]][extended_df[["Is_to_replace"]]] <- 2L




# Identify plates that should be transferred from pool 2 to pool 1 --------

are_good_plates <- extended_df[, "Plate_name"] %in% good_plates
are_pool1 <- extended_df[, "Run2_pool"] %in% 1
are_pool2 <- extended_df[, "Run2_pool"] %in% 2

leave_plates <- "HA_40"
are_exceptions <- extended_df[, "Plate_name"] %in% leave_plates
to_switch_pool1 <- are_good_plates & are_pool1 & !(are_exceptions)

good_barcodes_pool1 <- extended_df[["Barcode_ID"]][to_switch_pool1]
good_barcodes_pool2 <- extended_df[["Barcode_ID"]][are_good_plates & are_pool2]
replace_barcodes <- extended_df[["Barcode_ID"]][extended_df[["Is_to_replace"]]]
exchange_barcodes <- setdiff(good_barcodes_pool1,
                             c(good_barcodes_pool2, replace_barcodes)
                             )

extended_df[["Switch_pools"]] <- (extended_df[, "Barcode_ID"] %in% exchange_barcodes) &
                                  are_pool2

extended_df[["Run3_pool"]] <- ifelse(extended_df[["Switch_pools"]],
                                     3L,
                                     ifelse(are_good_plates,
                                            NA,
                                            extended_df[, "Run2_pool"] + 2L
                                            )
                                     )




# Check if the distribution of samples is balanced ------------------------

sum(extended_df[["Correction_factor"]][extended_df[, "Run2_pool"] %in% 1], na.rm = TRUE)
sum(extended_df[["Correction_factor"]][extended_df[, "Run2_pool"] %in% 2], na.rm = TRUE)

table(extended_df[, "Run3_pool"])

sum(extended_df[["Correction_factor"]][extended_df[, "Run3_pool"] %in% 3], na.rm = TRUE)
sum(extended_df[["Correction_factor"]][extended_df[, "Run3_pool"] %in% 4], na.rm = TRUE)





# Re-order the data frame -------------------------------------------------

selected_columns <- c(
  "Plate_name", "Number_96wp", "Position_FGCZ_plate", "Barcode_ID",
  "Run2_pool", "Run3_pool", "Switch_pools", "Is_to_replace",
  "Correction_factor", "Colony_picked", "No_DNA",
  "Was_corrected", "Sum_counts_in_k", "Thousands_of_reads", "Number_of_wells",
  "Mean_read_count"
)

selected_columns <- union(selected_columns, names(extended_df))

new_order <- order(extended_df[["Run3_pool"]], extended_df[["Run2_pool"]])

extended_df <- extended_df[new_order, selected_columns]




# Export data -------------------------------------------------------------

exclude_columns <- c(
  "Position_96wp", "Plate_number", "Corrected_concentration",
  "Median_count", "Sum_counts"
)

export_df <- extended_df[, !(names(extended_df) %in% exclude_columns)]

export_df[["Is_to_replace"]][is.na(export_df[["Run3_pool"]])] <- NA

use_columns <- c("Is_to_replace", "Switch_pools")
extra_columns <- c("Colony_picked", "Was_corrected", "No_DNA")

for (column_name in c(use_columns, extra_columns)) {
  if (!(column_name %in% extra_columns)) {
    export_df[[column_name]][is.na(export_df[["Run3_pool"]])] <- NA
  }
  export_df[[column_name]] <- ifelse(export_df[[column_name]], "Yes", "No")
}

write.table(export_df,
            file      = file.path(file_output_directory, "All_plates_new_pools.tsv"),
            row.names = FALSE,
            sep       = "\t",
            na        = ""
            )




# Save data ---------------------------------------------------------------

save(list = "extended_df",
     file = file.path(s2r2_objects_directory, "30) Calculate correction factors to account for mean read counts.RData")
     )


