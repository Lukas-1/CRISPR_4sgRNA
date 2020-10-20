### 8th September 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))




# Define folder paths -----------------------------------------------------

file_input_directory      <- file.path(file_directory, "2) Input")
R_objects_directory       <- file.path(file_directory, "3) R objects")

file_output_directory     <- file.path(file_directory, "5) Output")
tables_output_directory   <- file.path(file_output_directory, "Tables")

raw_data_directory        <- file.path(file_input_directory, "Raw data")
reanalysis_directory      <- file.path(raw_data_directory, "LukasAnalysis")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - demultiplexed - ccs3.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - ccs5 ZMWs.RData"))
load(file.path(R_objects_directory, "08) Extract barcode sequences and quality scores.RData"))





# Read in extra data (for checks) -----------------------------------------

ReadTable <- function(file_path, header = TRUE, sep = "\t") {
  read.table(file_path, sep = sep, quote = "", stringsAsFactors = FALSE,
             header = header, row.names = NULL, check.names = FALSE
             )
}

sl7_CCS3_lima_counts_file <- file.path(reanalysis_directory, "SmrtLink7_CCS3/lima", "lukaslimaCCS3SL7output.lima.counts")
sl7_CCS5_lima_counts_file <- file.path(reanalysis_directory, "SmrtLink7_CCS5/lima", "lukaslimaCCS5SL7output.lima.counts")

sl7_ccs3_counts_df <- ReadTable(sl7_CCS3_lima_counts_file)
sl7_ccs5_counts_df <- ReadTable(sl7_CCS5_lima_counts_file)

sl7_ccs3_counts_384_file <- file.path(raw_data_directory, "FGCZ/ccs_3_99/lima_26/lima.lima.384.counts.txt")
sl7_ccs5_counts_384_file <- file.path(raw_data_directory, "FGCZ/ccs_5_999/lima_26/lima.lima.384.counts.txt")

sl7_ccs3_384_counts_df <- ReadTable(sl7_ccs3_counts_384_file)
sl7_ccs5_384_counts_df <- ReadTable(sl7_ccs5_counts_384_file)

ccs3_fgcz_accuracies_df <- ReadTable(file.path(raw_data_directory, "FGCZ/ccs_3_99/lima_26/blastn/identical.perc.txt"))
ccs5_fgcz_accuracies_df <- ReadTable(file.path(raw_data_directory, "FGCZ/ccs_5_999/lima_26/blastn/identical.perc.txt"))





# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)






# Process the data on the level of individual reads -----------------------

sl7_ccs3_analysis_list <- AnalyzeWells(sl7_ccs3_lima,
                                       sl7_ccs3_report_df,
                                       sl7_barcodes_df,
                                       sl7_ccs3_counts_df,
                                       sl7_ccs3_384_counts_df
                                       )

sl9_ccs3_analysis_list <- AnalyzeWells(sl9_ccs3_lima, sl9_ccs3_report_df, sl9_barcodes_df)





# Create the summary data frames ------------------------------------------

sl7_ccs3_df_list <- SummarizeWells(sl7_ccs3_analysis_list)
sl9_ccs3_df_list <- SummarizeWells(sl9_ccs3_analysis_list)

sl7_ccs5_df_list <- SummarizeWells(sl7_ccs3_analysis_list, use_zmws = sl7_ccs5_lima_zmws)
sl9_ccs5_df_list <- SummarizeWells(sl9_ccs3_analysis_list, use_zmws = sl9_ccs5_lima_zmws)







# Check for equivalence with the read counts of the FGCZ ------------------

CheckCountsAreEqual <- function(reanalysis_df, fgcz_df) {
  stopifnot(identical(reanalysis_df[["Count_total"]], fgcz_df[["tot"]]))
  for (i in 1:4) {
    reanalysis_col <- paste0("Count_sg", i, "_cr", i)
    fgcz_col <- paste0("sg", i, "_cr", i, "_identical")
    stopifnot(identical(reanalysis_df[[reanalysis_col]], fgcz_df[[fgcz_col]]))
  }
  return(invisible(NULL))
}

CheckCountsAreEqual(sl7_ccs3_df_list[["original_summary_df"]],
                    ccs3_fgcz_accuracies_df
                    )
CheckCountsAreEqual(sl7_ccs5_df_list[["original_summary_df"]],
                    ccs5_fgcz_accuracies_df
                    )




# Export tables -----------------------------------------------------------

ExportSummaryTable(sl7_ccs3_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_original"
                   )
ExportSummaryTable(sl7_ccs3_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_filtered"
                   )
ExportIndivTable(sl7_ccs3_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS3_99_individual_reads"
                 )
ExportTable(sl7_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS3_99_contaminations"
            )


ExportSummaryTable(sl7_ccs5_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_original"
                   )
ExportSummaryTable(sl7_ccs5_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_filtered"
                   )
ExportIndivTable(sl7_ccs5_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS5_999_individual_reads"
                 )
ExportTable(sl7_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS5_999_contaminations"
            )



ExportSummaryTable(sl9_ccs3_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_original"
                   )
ExportSummaryTable(sl9_ccs3_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_filtered"
                   )
ExportIndivTable(sl9_ccs3_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS3_99_individual_reads"
                 )
ExportTable(sl9_ccs3_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS3_99_contaminations"
            )

ExportSummaryTable(sl9_ccs5_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_original"
                   )
ExportSummaryTable(sl9_ccs5_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_filtered"
                   )
ExportIndivTable(sl9_ccs5_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS5_999_individual_reads"
                 )
ExportTable(sl9_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS5_999_contaminations"
            )






# Save data ---------------------------------------------------------------

save(list = c("sl7_ccs3_df_list", "sl7_ccs5_df_list",
              "sl9_ccs3_df_list", "sl9_ccs5_df_list"
              ),
     file = file.path(R_objects_directory, "10) Process demultiplexed PacBio reads.RData")
     )












