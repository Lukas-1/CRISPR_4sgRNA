### 18th August 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")

CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")
source(file.path(s2r1_R_functions_directory, "06) Drawing scatter plots of plate-level summary data.R"))




# Define folder paths -----------------------------------------------------

experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")

file_output_directory    <- file.path(s2rC_directory, "5) Output", "Tables", "Correction factors")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "29) Correlate median read count with DNA concentration.RData"))

load(file.path(s2rC_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rC_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Read in data ------------------------------------------------------------

pool3_df <- data.frame(read_excel(file.path(s2r2_directory, "5) Output", "Tables", "Correction factors", "2021_08_25__Pool_3_plates.xlsx")),
                       check.names = FALSE, stringsAsFactors = FALSE
                       )




# Produce plate-level summary data ----------------------------------------

summaries_df <- SummarizePlates(ccs7_df_list[["filtered_summary_df"]])




# Integrate the selection for pool 3 --------------------------------------

matches_vec <- match(summaries_df[, "Plate_name"], pool3_df[, "Plate name"])

summaries_df[["In_pool3"]] <- !(is.na(matches_vec))
summaries_df[["Pool3_weight"]] <- pool3_df[matches_vec, "Correct-ion factor (relative amount)"]

matches_vec <- match(summaries_df[, "Plate_name"], merged_plates_df[, "Plate_name"])
summaries_df[["Position_FGCZ_plate"]] <- merged_plates_df[matches_vec, "Position_FGCZ_plate"]





# Define the plates that require re-ligation ------------------------------

replace_plates <- c("HA_21", "HO_13", "HO_16", "HO_15", "HO_19")
summaries_df[["Religation"]] <- summaries_df[["Plate_name"]] %in% replace_plates




# Define the plates that do not need re-sequencing ------------------------

are_run1 <- (summaries_df[["Plate_number"]] <= 23)

have_wells_with_no_reads    <- summaries_df[["Fraction_wells_with_no_reads"]] > 0
have_underrepresented_wells <- summaries_df[["Fraction_good_wells_with_few_reads"]] != 0


are_borderline <- have_wells_with_no_reads & !(have_underrepresented_wells)
summaries_df[are_borderline, ] # borderline wells


are_good <- (summaries_df[["Fraction_good_wells_with_few_reads"]] < 0.01) &
            (!(have_wells_with_no_reads)) &
            (!(summaries_df[["Religation"]]))

are_good[are_borderline %in% TRUE] <- TRUE

do_not_resequence <- c(
  "HA_24", "HA_25", "HA_28", "HA_39", "HA_41",
  "HO_18", "HO_23", "HO_31", "HO_32", "HO_34", "HO_45"
)
are_good[summaries_df[["Plate_name"]] %in% do_not_resequence] <- TRUE

are_good[are_run1 | summaries_df[["Colony_picked"]]] <- NA

summaries_df[are_good %in% TRUE, ]

# are_chosen <- !(are_run1 | are_good)





# Calculate rounded correction factors based on mean read counts ----------

mean_of_means <- mean(summaries_df[, "Mean_read_count"][are_good %in% FALSE])
corrections <- mean_of_means / summaries_df[, "Mean_read_count"][are_good %in% FALSE]

rounded_corrections <- as.integer(round(corrections))
rounded_corrections <- ifelse(rounded_corrections == 0,
                              1L,
                              ifelse(rounded_corrections > 2,
                                     4L,
                                     rounded_corrections
                                     )
                              )

summaries_df[["Corrections_mean_counts"]][are_good %in% FALSE] <- rounded_corrections
summaries_df[["Corrections_mean_counts"]][summaries_df[["Religation"]]] <- 2L




# Calculate corrections based on the number of reads < 0 ------------------

barplot(summaries_df[, "Expected_fraction_below_0"][are_good %in% FALSE])
abline(h = c(0.005, 0.02, 0.05))

barplot(summaries_df[, "Expected_fraction_below_0"][are_good %in% FALSE])
abline(h = c(0.01, 0.05))

use_vec <- summaries_df[are_good %in% FALSE, "Expected_fraction_below_0"]
corrections_below_0 <- ifelse(use_vec < 0.01,
                              1L,
                              ifelse(use_vec < 0.05,
                                     2L,
                                     4L
                                     )
                              )
summaries_df[["Corrections_below_0"]] <- ifelse(are_run1, NA, 0L)
summaries_df[["Corrections_below_0"]][are_good %in% FALSE] <- corrections_below_0






# Integrate both types of correction factors ------------------------------

summaries_df[["Corrections_integrated"]] <- pmax(summaries_df[["Corrections_mean_counts"]],
                                                 summaries_df[["Corrections_below_0"]]
                                                 )
are_divergent <- (summaries_df[["Corrections_mean_counts"]] == 4) &
                 (summaries_df[["Corrections_below_0"]] == 1)
summaries_df[["Corrections_integrated"]][are_divergent] <- 2L





# Identify plates that should be transferred from pool 5 to pool 4 --------

are_pool1 <- summaries_df[, "Run2_pool"] %in% 1
are_pool2 <- summaries_df[, "Run2_pool"] %in% 2

are_to_switch <- (are_good %in% TRUE) & !(summaries_df[["Colony_picked"]])

good_barcodes_pool1 <- summaries_df[["Barcode_ID"]][are_pool1 & are_to_switch]
good_barcodes_pool2 <- summaries_df[["Barcode_ID"]][are_pool2 & are_to_switch]

religate_barcodes <- summaries_df[["Barcode_ID"]][summaries_df[["Religation"]]]
exchange_barcodes <- setdiff(good_barcodes_pool1,
                             c(good_barcodes_pool2, religate_barcodes)
                             )

summaries_df[are_pool2 & (summaries_df[["Barcode_ID"]] %in% exchange_barcodes), ]

exception_plates <- c("HA_20", "HA_23", "HA_57", "HA_58", "HO_51", "HO_52")
are_exceptions <- summaries_df[["Plate_name"]] %in% exception_plates

summaries_df[["Switch_pools"]] <- (summaries_df[, "Barcode_ID"] %in% exchange_barcodes) &
                                   are_pool2 & !(are_exceptions)

summaries_df[summaries_df[["Switch_pools"]], ]


are_to_include <- ((are_good %in% FALSE) | summaries_df[["Colony_picked"]]) &
                  (!(are_run1))
summaries_df[["Run4_pool"]] <- ifelse(summaries_df[["Switch_pools"]],
                                      4L,
                                      ifelse(are_to_include,
                                             summaries_df[, "Run2_pool"] + 3L,
                                             NA
                                             )
                                      )



# Perform manual corrections ----------------------------------------------

summaries_df[["Corrections_integrated"]][summaries_df[["In_pool3"]]] <- as.integer(ceiling(summaries_df[["Corrections_integrated"]][summaries_df[["In_pool3"]]] / 2))
summaries_df[["Corrections_integrated"]][summaries_df[["Plate_name"]] %in% "HA_35"] <- 4L

summaries_df[["Corrections_integrated"]][summaries_df[["Run4_pool"]] %in% 5] <- as.integer(ceiling(summaries_df[["Corrections_integrated"]][summaries_df[["Run4_pool"]] %in% 5] / 2))

CF_two_plates <- c("HA_45", "HA_47", "HO_13", "HO_16", "HO_29",
                   paste0("HO_", 37:40), "HA_43", "HO_14"
                   )
summaries_df[["Corrections_integrated"]][summaries_df[["Plate_name"]] %in% CF_two_plates] <- 2L

summaries_df[["Corrections_integrated"]][summaries_df[["Religation"]]] <- 2L

summaries_df[["Corrections_integrated"]][summaries_df[["Plate_name"]] %in% c("HA_55", "HO_19")] <- 4L

summaries_df[["Corrections_integrated"]][summaries_df[["Colony_picked"]] & !(are_run1)] <- 0.3

summaries_df[["Corrections_integrated"]][summaries_df[["Plate_name"]] %in% c("HA_46", paste0("HO_", 24:28), "HO_47")] <- 1L






# Calculate relative weights between the different pools ------------------

are_pool4 <- summaries_df[["Run4_pool"]] %in% 4
are_pool5 <- summaries_df[["Run4_pool"]] %in% 5

summaries_df[["Weight"]] <- NA
pool4_CF_sum <- sum(summaries_df[["Corrections_integrated"]][are_pool4])
pool5_CF_sum <- sum(summaries_df[["Corrections_integrated"]][are_pool5])
weight_ratio <- pool4_CF_sum / pool5_CF_sum
summaries_df[["Weight"]][are_pool4] <- summaries_df[["Corrections_integrated"]][are_pool4]
summaries_df[["Weight"]][are_pool5] <- summaries_df[["Corrections_integrated"]][are_pool5] * weight_ratio

pool3_CF_sum <- sum(summaries_df[["Pool3_weight"]], na.rm = TRUE)




# Check if the distribution of samples is balanced ------------------------

table(summaries_df[, "Run4_pool"])
table(summaries_df[, "Corrections_below_0"])
table(summaries_df[, "Corrections_integrated"])
table(summaries_df[, "Corrections_mean_counts"])

sum(summaries_df[["Corrections_mean_counts"]][are_pool4], na.rm = TRUE)
sum(summaries_df[["Corrections_mean_counts"]][are_pool5], na.rm = TRUE)

sum(summaries_df[["Corrections_below_0"]][are_pool4], na.rm = TRUE)
sum(summaries_df[["Corrections_below_0"]][are_pool5], na.rm = TRUE)

sum(summaries_df[["Corrections_integrated"]][are_pool4], na.rm = TRUE)
sum(summaries_df[["Corrections_integrated"]][are_pool5], na.rm = TRUE)

sum(summaries_df[["Weight"]][are_pool4], na.rm = TRUE)
sum(summaries_df[["Weight"]][are_pool5], na.rm = TRUE)

stopifnot(!(anyDuplicated(summaries_df[are_pool4, "Barcode_ID"])))
stopifnot(!(anyDuplicated(summaries_df[are_pool5, "Barcode_ID"])))




# Re-order the data frame -------------------------------------------------

selected_columns <- c(
  "Plate_name", "Number_96wp", "Position_FGCZ_plate",
  "Barcode_ID", "Run2_pool", "In_pool3", "Pool3_weight", "Pool3_weight",
  "Run4_pool", "Switch_pools", "Religation", "Weight",
  "Corrections_integrated", "Corrections_below_0", "Corrections_mean_counts",
  "Colony_picked", "Sum_counts_in_k", "Number_of_wells", "Mean_read_count"
)

selected_columns <- union(selected_columns, names(summaries_df))

run4_df <- summaries_df[!(summaries_df[["Run2_pool"]]) %in% 0, ]

new_order <- order(run4_df[["Run4_pool"]],
                   run4_df[["Colony_picked"]],
                   run4_df[["Run2_pool"]]
                   )

run4_df <- run4_df[new_order, selected_columns]
row.names(run4_df) <- NULL

switched_in_pool3 <- (run4_df[["Run2_pool"]] %in% 2) & run4_df[["In_pool3"]]
run4_df[switched_in_pool3, ]


run4_df[run4_df[, "Switch_pools"], ]

run4_df[(run4_df[["Pool3_weight"]] > run4_df[["Corrections_integrated"]]) %in% TRUE, ]


run4_df[run4_df[, "Switch_pools"], ]

new_weight1 <- (run4_df[["In_pool3"]] %in% FALSE) &
               (run4_df[["Corrections_integrated"]] %in% 1)

run4_df[new_weight1, ]





# Export data -------------------------------------------------------------

exclude_columns <- c("Position_96wp", "Plate_number", "DNA_isolation",
                     "Median_read_count", "Sum_counts",
                     "Median_count_adjust_by_no_of_wells"
                     )

export_df <- run4_df[, !(names(run4_df) %in% exclude_columns)]

export_df[["Religation"]][is.na(export_df[["Run4_pool"]])] <- NA

use_columns <- c("Religation", "Switch_pools")
extra_columns <- c("Colony_picked", "In_pool3")

export_df[["Corrections_below_0"]][export_df[["Corrections_below_0"]] == 0] <- NA

for (column_name in c(use_columns, extra_columns)) {
  if (!(column_name %in% extra_columns)) {
    export_df[[column_name]][is.na(export_df[["Run4_pool"]])] <- NA
  }
  export_df[[column_name]] <- ifelse(export_df[[column_name]], "Yes", "No")
}

write.table(export_df,
            file      = file.path(file_output_directory, "All_plates_pools_4_and_5.tsv"),
            row.names = FALSE,
            sep       = "\t",
            na        = ""
            )




# Save data ---------------------------------------------------------------

save(list = "run4_df",
     file = file.path(s2rC_R_objects_directory, "30) Calculate plate weightings for the next round of sequencing.RData")
     )


