### 17th April 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))





# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
tables_output_directory  <- file.path(file_output_directory, "Tables")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))






# Export individual reads -------------------------------------------------

ccs_numbers <- c(3, 5, 7)
ccs_folders <- c("CCS3_99", "CCS5_999", "CCS7_9999")

for (i in seq_along(ccs_numbers)) {
  this_df <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))[["individual_reads_df"]]
  plate_numbers <- setdiff(this_df[, "Plate_number"], NA)
  plate_names <- paste0("Plate", formatC(plate_numbers, width = 2, flag = "0"))
  this_dir <- file.path(tables_output_directory, "Individual reads", ccs_folders[[i]])
  for (j in seq_along(plate_numbers)) {
    sub_df <- this_df[this_df[, "Plate_number"] %in% j, ]
    row.names(sub_df) <- NULL
    ExportIndivTable(sub_df,
                     paste0(ccs_folders[[i]], "_individual_reads_", plate_names[[j]]),
                     this_dir
                     )

  }
}





# Export tables -----------------------------------------------------------

ExportSummaryTable(ccs3_df_list[["original_summary_df"]],
                   "CCS3_99_summary_per_well_unfiltered",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )
ExportSummaryTable(ccs3_df_list[["filtered_summary_df"]],
                   "CCS3_99_summary_per_well_filtered",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )
ExportSummaryTable(ccs3_df_list[["filtered_gRNAs_df"]],
                   "CCS3_99_summary_per_well_filtered_gRNAs",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )


ExportSummaryTable(ccs5_df_list[["original_summary_df"]],
                   "CCS5_999_summary_per_well_unfiltered",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )
ExportSummaryTable(ccs5_df_list[["filtered_summary_df"]],
                   "CCS5_999_summary_per_well_filtered",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )
ExportSummaryTable(ccs5_df_list[["filtered_gRNAs_df"]],
                   "CCS5_999_summary_per_well_filtered_gRNAs",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )


ExportSummaryTable(ccs7_df_list[["original_summary_df"]],
                   "CCS7_9999_summary_per_well_unfiltered",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )
ExportSummaryTable(ccs7_df_list[["filtered_summary_df"]],
                   "CCS7_9999_summary_per_well_filtered",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )
ExportSummaryTable(ccs7_df_list[["filtered_gRNAs_df"]],
                   "CCS7_9999_summary_per_well_filtered_gRNAs",
                   file_directory = file.path(tables_output_directory, "Summary tables")
                   )





# Define columns to export ------------------------------------------------

use_columns <- c(
  "Combined_ID", "Plate_number", "Well_number",

  "Known_empty",

  "Count_total",

  "Perc_at_least_1", "Perc_at_least_2", "Perc_at_least_3", "Perc_all_4",
  "Perc_all_4_promoters", "Perc_whole_plasmid",
  "Perc_sg1_cr1", "Perc_sg2_cr2", "Perc_sg3_cr3", "Perc_sg4_cr4",

  "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
  "Count_all_4_promoters", "Count_whole_plasmid",
  "Count_sg1_cr1", "Count_sg2_cr2", "Count_sg3_cr3", "Count_sg4_cr4",

  "Num_under_2kb", "Mean_read_quality",

  "Num_contaminating_genes", "Num_contaminated_reads",
  "Num_from_close_wells", "Expected_from_close_wells",
  # "Mean_distance", "Expected_distance", "Distance_p_value",

  "Correct_TpR_DHFR", "Mutation_TpR_DHFR", "Deletion_TpR_DHFR",
  "Correct_sg1", "Mutation_sg1", "Deletion_sg1", "Contamination_sg1",
  "Correct_sg2", "Mutation_sg2", "Deletion_sg2", "Contamination_sg2",
  "Correct_sg3", "Mutation_sg3", "Deletion_sg3", "Contamination_sg3",
  "Correct_sg4", "Mutation_sg4", "Deletion_sg4", "Contamination_sg4",
  "Correct_sg1_cr1", "Mutation_sg1_cr1", "Deletion_sg1_cr1", "Contamination_sg1_cr1",
  "Correct_sg2_cr2", "Mutation_sg2_cr2", "Deletion_sg2_cr2", "Contamination_sg2_cr2",
  "Correct_sg3_cr3", "Mutation_sg3_cr3", "Deletion_sg3_cr3", "Contamination_sg3_cr3",
  "Correct_sg4_cr4", "Mutation_sg4_cr4", "Deletion_sg4_cr4", "Contamination_sg4_cr4"
)





# Identify problematic wells ----------------------------------------------

use_df <- ccs7_df_list[["filtered_summary_df"]]
use_df[["Known_empty"]] <- ifelse(sg_sequences_df[["Empty_well"]], "Yes", NA)

order_by_count   <- order(use_df[, "Count_total"],
                          use_df[, "Perc_at_least_1"],
                          na.last = FALSE
                          )
order_by_correct <- order(use_df[, "Perc_at_least_1"],
                          use_df[, "Count_total"],
                          na.last = FALSE
                          )

ExportSummaryTable(use_df[order_by_count, use_columns],
                   "CCS7_9999_filtered_ordered_by_number_of_reads",
                   file_directory = file.path(tables_output_directory,
                                              "Summary tables",
                                              "Re-ordered (problematic wells first)"
                                              )
                   )

ExportSummaryTable(use_df[order_by_correct, use_columns],
                   "CCS7_9999_filtered_ordered_by_percentage_correct",
                   file_directory = file.path(tables_output_directory,
                                              "Summary tables",
                                              "Re-ordered (problematic wells first)"
                                              )
                   )





