### 22nd August 2021 ###




# Import packages and source code -----------------------------------------

library("readxl")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
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

religate_plates <- c("HA_21", "HO_13", "HO_16")




# Merge the plate annotation and sequencing yield summary data ------------

extended_df <- merged_plates_df

new_order <- order(extended_df[, "Run2_pool"], extended_df[, "Plate_number"])
extended_df <- extended_df[new_order, ]
row.names(extended_df) <- NULL





# Calculate rounded correction factors based on mean read counts ----------

are_remaining <- !(extended_df[, "Plate_name"] %in% c(good_plates, religate_plates))

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

extended_df[["Corrections_mean_counts"]][are_remaining] <- rounded_corrections

extended_df[["Religation"]] <- extended_df[["Plate_name"]] %in% religate_plates
extended_df[["Corrections_mean_counts"]][extended_df[["Religation"]]] <- 2L





# Calculate corrections based on the number of reads < 0 ------------------

barplot(extended_df[, "Expected_fraction_below_0"][!(extended_df[, "Plate_name"] %in% religate_plates)])

use_vec <- extended_df[, "Expected_fraction_below_0"]
corrections_below_0 <- ifelse(use_vec < 0.005,
                              0L,
                              ifelse(use_vec < 0.02,
                                     1L,
                                     ifelse(use_vec < 0.05,
                                            2L,
                                            4L
                                            )
                                     )
                              )
extended_df[["Corrections_below_0"]] <- ifelse(!(extended_df[, "Plate_name"] %in% good_plates) &
                                                (corrections_below_0 == 0),
                                               1L,
                                               corrections_below_0
                                               )




# Identify plates that should be transferred from pool 2 to pool 1 --------

are_good_plates <- extended_df[, "Plate_name"] %in% good_plates
are_pool1 <- extended_df[, "Run2_pool"] %in% 1
are_pool2 <- extended_df[, "Run2_pool"] %in% 2

leave_plates <- c("HA_40", "HO_7")
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




# Check for pool 3 plates that may be repeated in pool 4 ------------------

are_to_repeat <- ((extended_df[["Corrections_below_0"]] == 4) |
                 ((extended_df[["Corrections_below_0"]] == 2) &
                  (extended_df[["Corrections_mean_counts"]] == 1))) &
                  (extended_df[["Run3_pool"]] == 3)
repeat_barcodes <- extended_df[are_to_repeat, "Barcode_ID"]
pool4_barcodes <- extended_df[extended_df[["Run3_pool"]] == 4, "Barcode_ID"]

free_barcodes <- setdiff(extended_df[, "Barcode_ID"], pool4_barcodes)
extended_df[are_to_repeat, ]

manual_selection <- c(
  "HA_16", "HO_8", "HO_10", "HO_30", # high-priority
  "HO_11", "HA_18" # additional (exceptions)
  # "HA_19", "HA_24", "HA_25" # considered for re-inclusion, but not chosen
  # "HO_31", "HA_50", "HO_9", "HO_23", "HO_32", "HO_34", "HO_47"# flagged by the criteria above, but not included
)

extended_df[extended_df[["Barcode_ID"]] %in% free_barcodes, ]




# Choose selected pool 3 plates to repeat in pool 4 -----------------------

are_chosen <- extended_df[["Plate_name"]] %in% manual_selection
extended_df[["Run3_pool"]] <- ifelse(are_chosen, "3 and 4", extended_df[["Run3_pool"]])
barcode_column <- which(names(extended_df) == "Barcode_ID")
names(extended_df)[[barcode_column]] <- "Original_barcode_ID"
extended_df[["New_barcode_ID"]] <- NA
extended_df[["New_barcode_ID"]][are_chosen] <- free_barcodes[seq_along(manual_selection)]
extended_df[["Religation"]][are_chosen] <- TRUE


# column_indices <- c(seq_len(barcode_column), ncol(extended_df),
#                     (barcode_column + 1L):(ncol(extended_df) - 1L)
#                     )
# extended_df <- extended_df[, column_indices]




# Examine cases where correction factors changed --------------------------

are_to_examine <- ((extended_df[["Corrections_below_0"]] < extended_df[["Corrections_mean_counts"]]) %in% TRUE) &
                  (extended_df[, "Run3_pool"] %in% 4)

extended_df[["Corrections_below_0"]] <- ifelse(are_to_examine &
                                               (extended_df[["Corrections_mean_counts"]] == 2) &
                                                !(extended_df[, "Plate_name"] %in% c("HO_17", "HO_24", "HO_26")),
                                               extended_df[["Corrections_mean_counts"]],
                                               extended_df[["Corrections_below_0"]]
                                               )




# Check if the distribution of samples is balanced ------------------------

sum(extended_df[["Corrections_mean_counts"]][extended_df[, "Run2_pool"] %in% 1], na.rm = TRUE)
sum(extended_df[["Corrections_mean_counts"]][extended_df[, "Run2_pool"] %in% 2], na.rm = TRUE)

table(extended_df[, "Run3_pool"])

sum(extended_df[["Corrections_mean_counts"]][extended_df[, "Run3_pool"] %in% 3], na.rm = TRUE)
sum(extended_df[["Corrections_mean_counts"]][extended_df[, "Run3_pool"] %in% 4], na.rm = TRUE)
sum(extended_df[["Corrections_mean_counts"]][extended_df[, "Run3_pool"] %in% c("4", "3 and 4")], na.rm = TRUE)


sum(extended_df[["Corrections_below_0"]][extended_df[, "Run2_pool"] %in% 1], na.rm = TRUE)
sum(extended_df[["Corrections_below_0"]][extended_df[, "Run2_pool"] %in% 2], na.rm = TRUE)
sum(extended_df[["Corrections_below_0"]][extended_df[, "Run3_pool"] %in% 3], na.rm = TRUE)
sum(extended_df[["Corrections_below_0"]][extended_df[, "Run3_pool"] %in% 4], na.rm = TRUE)

sum(extended_df[["Corrections_below_0"]][extended_df[, "Run3_pool"] %in%  c("4", "3 and 4")], na.rm = TRUE)


are_to_examine <- (extended_df[["Corrections_below_0"]] < extended_df[["Corrections_mean_counts"]]) &
                  (extended_df[, "Run3_pool"] %in% 4)
extended_df[are_to_examine, ]




# Re-order the data frame -------------------------------------------------

selected_columns <- c(
  "Plate_name", "Number_96wp", "Position_FGCZ_plate",
  "Original_barcode_ID", "New_barcode_ID",
  "Run2_pool", "Run3_pool", "Switch_pools", "Religation",
  "Corrections_mean_counts", "Corrections_below_0", "Colony_picked", "No_DNA",
  "Was_corrected", "Sum_counts_in_k", "Thousands_of_reads", "Number_of_wells",
  "Mean_read_count"
)

selected_columns <- union(selected_columns, names(extended_df))

new_order <- order(extended_df[["Run3_pool"]], extended_df[["Run2_pool"]])

extended_df <- extended_df[new_order, selected_columns]
row.names(extended_df) <- NULL




# Export data -------------------------------------------------------------

exclude_columns <- c(
  "Position_96wp", "Plate_number", "Corrected_concentration",
  "Median_read_count", "Sum_counts"
)

export_df <- extended_df[, !(names(extended_df) %in% exclude_columns)]

export_df[["Religation"]][is.na(export_df[["Run3_pool"]])] <- NA

use_columns <- c("Religation", "Switch_pools")
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


