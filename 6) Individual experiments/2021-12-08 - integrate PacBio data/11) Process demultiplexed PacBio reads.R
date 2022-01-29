### 8th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))




# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rI_R_objects_directory, "07.5) Pre-filter reads - ccs_df.RData"))
load(file.path(s2rI_R_objects_directory, "07) Extract barcode sequences and quality scores.RData"))
load(file.path(s2rI_R_objects_directory, "08-b) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(s2rI_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData"))
load(file.path(s2rI_R_objects_directory, "10) Identify and characterize deletions.RData"))




# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)




# Examine wells with no reads ---------------------------------------------

sg_sequences_df[!(sg_sequences_df[["Combined_ID"]] %in% ccs_df[["Combined_ID"]]), 1:9]




# Process the data on the level of individual reads -----------------------

plates_analysis_list <- AnalyzePlates(ccs_df,
                                      sg_sequences_df,
                                      barcodes_df,
                                      extracted_df
                                      )



# Create the summary data frames ------------------------------------------

ccs_df[["Passed_filters"]] <- ccs_df[["Plate_passed_filters"]] &
                              (ccs_df[["Well_passed_filters"]] %in% TRUE)

ccs3_zmws <- GetCCS3_ZMWs(ccs_df)
ccs7_zmws <- GetCCS7_ZMWs(ccs_df)

ccs3_df_list <- SummarizeWells(plates_analysis_list,
                               use_zmws           = ccs3_zmws,
                               ID_column          = "Combined_ID",
                               unique_IDs         = sg_sequences_df[["Combined_ID"]],
                               deletions_df       = deletions_df,
                               aligned_contam_df  = contam_df,
                               filter_cross_plate = TRUE
                               )

ccs7_df_list <- SummarizeWells(plates_analysis_list,
                               use_zmws           = ccs7_zmws,
                               ID_column          = "Combined_ID",
                               unique_IDs         = sg_sequences_df[["Combined_ID"]],
                               deletions_df       = deletions_df,
                               aligned_contam_df  = contam_df,
                               filter_cross_plate = TRUE
                               )



# Remove unnecessary data -------------------------------------------------

ccs3_df_list <- ccs3_df_list[c("original_summary_df", "individual_reads_df")]
stopifnot(identical(ccs7_df_list[["original_summary_df"]], ccs7_df_list[["filtered_summary_df"]]))
ccs7_df_list <- ccs7_df_list[names(ccs7_df_list) != "original_summary_df"]




# Save data ---------------------------------------------------------------

save(list = "plates_analysis_list",
     file = file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - plates_analysis_list.RData")
     )

save(list = c("ccs3_df_list", "ccs7_df_list"),
     file = file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData")
     )



