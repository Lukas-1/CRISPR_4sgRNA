### 8th June 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
plate1_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory      <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions
source(file.path(R_functions_directory, "20) Summarizing data across wells.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output", "Figures", "Summaries across wells")
PNGs_output_directory    <- file.path(sql2_directory, "5) Output", "PNGs", "Summaries across wells")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Define plate selections -------------------------------------------------

are_beads <- grepl("-beads", plates_df[["Plate_name"]], fixed = TRUE)
are_controls <- plates_df[["Colony_picked"]]

non_library_plates <- c("Vac-1", "PD_A", "PD_O")
are_columns <- !(are_beads | are_controls)

are_library <- are_columns & (!(plates_df[["Plate_name"]] %in% non_library_plates))

column_matched_plates <- sub("-beads", "", plates_df[["Plate_name"]][are_beads], fixed = TRUE)

are_early_batch <- plates_df[["Plate_name"]] %in% c(paste0("HA_", 1:5), paste0("HO_", 1:5))
are_late_batch <- !(are_early_batch | are_beads | are_controls)

plate_selection_list <- list(
  "All plates"              = plates_df[["Plate_name"]][are_columns],
  "Colony-picked"           = "Intctrl",
  "Bead-purified"           = plates_df[["Plate_name"]][are_beads],
  "Matched column-purified" = column_matched_plates,
  "Library plates"          = plates_df[["Plate_name"]][are_library],
  "Early batch"             = plates_df[["Plate_name"]][are_early_batch],
  "Later batch"             = plates_df[["Plate_name"]][are_late_batch]
)

plate_selection_titles_list <- list(
  "All plates"              = "All plates (purified using columns)",
  "Bead-purified"           = "Plates HA-11 and HO-1 (purified using beads)",
  "Matched column-purified" = "Plates HA-11 and HO-1 (purified using columns)",
  "Library plates"          = "4sg library plates (purified using columns)",
  "Early batch"             = "Early batch (CRISPRa/o plates 1-5)",
  "Later batch"             = "Plates from later batches"
)




# Create labels for individual plates -------------------------------------

plate_labels <- paste0("Plate #", plates_df[["Plate_number"]], " \u2013 ", plates_df[["Plate_name"]])
names(plate_labels) <- plates_df[["Plate_name"]]
plate_labels <- plate_labels[order(plates_df[["Plate_rank"]])]
plate_labels <- plate_labels[names(plate_labels) != "Intctrl"]
are_single_plates <- c(rep(FALSE, length(plate_selection_titles_list)),
                       rep(TRUE, length(plate_labels))
                       )
plate_selection_titles_list <- c(plate_selection_titles_list, as.list(plate_labels))





# Prepare file name strings for plate selections --------------------------

plate_matches <- match(names(plate_selection_titles_list)[are_single_plates],
                       plates_df[["Plate_name"]]
                       )
plate_numbers <- plates_df[["Plate_number"]][plate_matches]
plate_selection_prefixes <- c(paste0("0", letters[seq_len(sum(!(are_single_plates)))]),
                              formatC(plate_numbers, width = 2, flag = "0")
                              )




# Produce example plots ---------------------------------------------------

use_df <- ccs7_df_list[["filtered_summary_df"]]

SummaryBoxPlot(use_df, "All plates")


LollipopPlot(use_df, "All plates")
LollipopPlot(use_df, "Bead-purified")


LollipopPlot(use_df,
             "Matched column-purified",
             paste0("Correct_sg", 1:4)
             )
LollipopPlot(use_df,
             "Bead-purified",
             paste0("Correct_sg", 1:4)
             )

LollipopPlot(use_df,
             "Bead-purified",
             c("Num_contaminated_reads",
               "Num_contaminated_reads_aligned",
               "Num_cross_plate_contaminated"
               ),
             use_y_limits = c(0, 0.05)
             )

LollipopPlot(use_df,
             "Bead-purified",
             c("Num_contaminated_reads"),
             use_y_limits = c(0, 0.05)
             )

LollipopPlot(use_df,
             "All plates",
             c("Num_reads_with_sgRNA_deletion",
               "Num_reads_with_deletions_exceeding_20bp",
               "Num_reads_with_deletions_spanning_tracrRNAs",
               "Num_reads_with_deletions_spanning_promoters"
               )
             )





# Export lollipop plots and violin/box plots ------------------------------

DrawAllLollipopsAndViolins()





# Save data ---------------------------------------------------------------

save(list = c("plate_selection_list",
              "plate_selection_titles_list",
              "plate_selection_prefixes"
              ),
     file = file.path(sql2_R_objects_directory, "23) Summarize data across wells - plate selections.RData")
     )





# # Try stuff ---------------------------------------------------------------
#
# ExpectedCorrect <- function(p_mutant, at_least_num_correct = 1) {
#   exactly_0_mutations <- (1 - p_mutant)^4
#   exactly_1_mutations <- p_mutant * ((1 - p_mutant)^3) * choose(4, 1)
#   exactly_2_mutations <- (p_mutant^2) * ((1 - p_mutant)^2) * choose(4, 2)
#   exactly_3_mutations <- (p_mutant^3) * (1 - p_mutant) * choose(4, 3)
#   exactly_4_mutations <- (p_mutant^4)
#   if (at_least_num_correct == 4) {
#     result <- exactly_0_mutations
#   } else if (at_least_num_correct == 3) {
#     result <- exactly_0_mutations + exactly_1_mutations
#   } else if (at_least_num_correct == 2) {
#     result <- exactly_0_mutations + exactly_1_mutations + exactly_2_mutations
#   } else if (at_least_num_correct == 1) {
#     result <- exactly_0_mutations + exactly_1_mutations + exactly_2_mutations + exactly_3_mutations
#   } else if (at_least_num_correct == 0) {
#     result <- exactly_4_mutations
#   }
#   return(result)
# }
#
#
#
# GeneralExpectedCorrect <- function(use_p, use_n, use_k) {
#   choose(use_n, use_k) * (use_p^use_k) * ((1 - use_p)^(use_n - use_k))
# }
#
#
#
# error_rate <- (1 - 0.965)
#
#
#
# ExpectedCorrect(error_rate, 4)
# ExpectedCorrect(error_rate, 3)
# ExpectedCorrect(error_rate, 2)
# ExpectedCorrect(error_rate, 1)
#
#
#
# pbinom(0, size = 4, prob = error_rate)
# pbinom(1, size = 4, prob = error_rate)
# pbinom(2, size = 4, prob = error_rate)
# pbinom(3, size = 4, prob = error_rate)
#
#
#
# pbinom(0, size = 3, prob = error_rate)
# pbinom(1, size = 3, prob = error_rate)
# pbinom(2, size = 3, prob = error_rate)
#




