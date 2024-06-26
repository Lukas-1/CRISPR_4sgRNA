### 20th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))




# Define folder paths -----------------------------------------------------

plate2_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory  <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory  <- file.path(plate2_directory, "2) R objects")
file_output_directory   <- file.path(plate2_directory, "3) Output")
tables_output_directory <- file.path(file_output_directory, "Tables")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data.RData"))
load(file.path(p2_R_objects_directory, "05) Extract barcode sequences and quality scores.RData"))
load(file.path(p2_R_objects_directory, "06) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(p2_R_objects_directory, "07) Characterize contaminations (using aligned reads).RData"))
load(file.path(p2_R_objects_directory, "08) Identify and characterize deletions.RData"))





# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)





# Process the data on the level of individual reads -----------------------

sl7_ccs3_analysis_list <- AnalyzeWells(sl7_ccs_df,
                                       sg_sequences_df,
                                       sl7_barcodes_df,
                                       sl7_extracted_df,
                                       )
sl9_ccs3_analysis_list <- AnalyzeWells(sl9_ccs_df,
                                       sg_sequences_df,
                                       sl9_barcodes_df,
                                       sl9_extracted_df,
                                       )




# Create the summary data frames ------------------------------------------

sl7_ccs3_df_list <- SummarizeWells(sl7_ccs3_analysis_list, unique_IDs = sg_sequences_df[["Well_number"]],
                                   deletions_df = sl7_deletions_df, aligned_contam_df = sl7_contam_df
                                   )
sl9_ccs3_df_list <- SummarizeWells(sl9_ccs3_analysis_list, unique_IDs = sg_sequences_df[["Well_number"]],
                                   deletions_df = sl9_deletions_df, aligned_contam_df = sl9_contam_df
                                   )

sl7_ccs5_lima_zmws <- GetCCS5_ZMWs(sl7_ccs_df, wells_vec = sg_sequences_df[["Well_number"]])
sl9_ccs5_lima_zmws <- GetCCS5_ZMWs(sl9_ccs_df, wells_vec = sg_sequences_df[["Well_number"]])
sl7_ccs7_lima_zmws <- GetCCS7_ZMWs(sl7_ccs_df, wells_vec = sg_sequences_df[["Well_number"]])
sl9_ccs7_lima_zmws <- GetCCS7_ZMWs(sl9_ccs_df, wells_vec = sg_sequences_df[["Well_number"]])


sl7_ccs5_df_list <- SummarizeWells(sl7_ccs3_analysis_list, use_zmws = sl7_ccs5_lima_zmws,
                                   unique_IDs = sg_sequences_df[["Well_number"]],
                                   deletions_df = sl7_deletions_df, aligned_contam_df = sl7_contam_df
                                   )
sl7_ccs7_df_list <- SummarizeWells(sl7_ccs3_analysis_list, use_zmws = sl7_ccs7_lima_zmws,
                                   unique_IDs = sg_sequences_df[["Well_number"]],
                                   deletions_df = sl7_deletions_df, aligned_contam_df = sl7_contam_df
                                   )
sl9_ccs5_df_list <- SummarizeWells(sl9_ccs3_analysis_list, use_zmws = sl9_ccs5_lima_zmws,
                                   unique_IDs = sg_sequences_df[["Well_number"]],
                                   deletions_df = sl9_deletions_df, aligned_contam_df = sl9_contam_df
                                   )
sl9_ccs7_df_list <- SummarizeWells(sl9_ccs3_analysis_list, use_zmws = sl9_ccs7_lima_zmws,
                                   unique_IDs = sg_sequences_df[["Well_number"]],
                                   deletions_df = sl9_deletions_df, aligned_contam_df = sl9_contam_df
                                   )





# Export tables -----------------------------------------------------------

ExportSummaryTable(sl7_ccs3_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_unfiltered"
                   )
ExportSummaryTable(sl7_ccs3_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_filtered"
                   )
ExportSummaryTable(sl7_ccs3_df_list[["filtered_gRNAs_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_filtered_gRNAs"
                   )
ExportIndivTable(sl7_ccs3_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS3_99_individual_reads"
                 )
ExportTable(sl7_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS3_99_contaminations"
            )


ExportSummaryTable(sl7_ccs5_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_unfiltered"
                   )
ExportSummaryTable(sl7_ccs5_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_filtered"
                   )
ExportSummaryTable(sl7_ccs5_df_list[["filtered_gRNAs_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_filtered_gRNAs"
                   )
ExportIndivTable(sl7_ccs5_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS5_999_individual_reads"
                 )
ExportTable(sl7_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS5_999_contaminations"
            )


ExportSummaryTable(sl7_ccs7_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS7_9999_summary_per_well_unfiltered"
                   )
ExportSummaryTable(sl7_ccs7_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS7_9999_summary_per_well_filtered"
                   )
ExportSummaryTable(sl7_ccs7_df_list[["filtered_gRNAs_df"]],
                   "SmrtLink7/SmrtLink7_CCS7_9999_summary_per_well_filtered_gRNAs"
                   )
ExportIndivTable(sl7_ccs7_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS7_9999_individual_reads"
                 )
ExportTable(sl7_ccs7_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS7_9999_contaminations"
            )





ExportSummaryTable(sl9_ccs3_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_unfiltered"
                   )
ExportSummaryTable(sl9_ccs3_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_filtered"
                   )
ExportSummaryTable(sl9_ccs3_df_list[["filtered_gRNAs_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_filtered_gRNAs"
                   )
ExportIndivTable(sl9_ccs3_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS3_99_individual_reads"
                 )
ExportTable(sl9_ccs3_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS3_99_contaminations"
            )


ExportSummaryTable(sl9_ccs5_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_unfiltered"
                   )
ExportSummaryTable(sl9_ccs5_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_filtered"
                   )
ExportSummaryTable(sl9_ccs5_df_list[["filtered_gRNAs_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_filtered_gRNAs"
                   )
ExportIndivTable(sl9_ccs5_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS5_999_individual_reads"
                 )
ExportTable(sl9_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS5_999_contaminations"
            )


ExportSummaryTable(sl7_ccs7_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS7_9999_summary_per_well_unfiltered"
                   )
ExportSummaryTable(sl7_ccs7_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS7_9999_summary_per_well_filtered"
                   )
ExportSummaryTable(sl7_ccs7_df_list[["filtered_gRNAs_df"]],
                   "SmrtLink9/SmrtLink9_CCS7_9999_summary_per_well_filtered_gRNAs"
                   )
ExportIndivTable(sl7_ccs7_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS7_9999_individual_reads"
                 )
ExportTable(sl7_ccs7_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS7_9999_contaminations"
            )




# Save data ---------------------------------------------------------------

save(list = c("sl7_ccs3_df_list", "sl7_ccs5_df_list", "sl7_ccs7_df_list",
              "sl9_ccs3_df_list", "sl9_ccs5_df_list", "sl9_ccs7_df_list"
              ),
     file = file.path(p2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData")
     )






