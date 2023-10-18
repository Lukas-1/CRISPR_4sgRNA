### 20th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
plate1_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "13) Importing PacBio reads.R"))





# Define folder paths -----------------------------------------------------

plate2_directory          <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")

p1_R_objects_directory    <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory    <- file.path(plate2_directory, "2) R objects")

file_input_directory      <- file.path(plate2_directory, "1) Input")
raw_data_directory        <- file.path(file_input_directory, "Raw data", "Analysis")

ccs_sam_name <- "m54073_200523_132654.ccs.sam"
sl7_CCS3_ccs_sam_file <- file.path(raw_data_directory, "SmrtLink7_CCS3/SAM", ccs_sam_name)
sl9_CCS3_ccs_sam_file <- file.path(raw_data_directory, "SmrtLink9_CCS3/SAM", ccs_sam_name)

sl7_CCS3_lima_report_file <- file.path(raw_data_directory, "SmrtLink7_CCS3/lima", "lukas_lima_CCS3_SL7.lima.report")
sl7_CCS5_lima_report_file <- file.path(raw_data_directory, "SmrtLink7_CCS5/lima", "lukas_lima_CCS5_SL7.lima.report")
sl9_CCS3_lima_report_file <- file.path(raw_data_directory, "SmrtLink9_CCS3/lima", "lukas_lima_CCS3_SL9.lima.report")
sl9_CCS5_lima_report_file <- file.path(raw_data_directory, "SmrtLink9_CCS5/lima", "lukas_lima_CCS5_SL9.lima.report")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))





# Read in data ------------------------------------------------------------

sl7_ccs3_sam_df <- ProcessSAM(sl7_CCS3_ccs_sam_file)
sl9_ccs3_sam_df <- ProcessSAM(sl9_CCS3_ccs_sam_file)

sl7_ccs3_report_df <- ReadLimaReport(sl7_CCS3_lima_report_file)
sl7_ccs5_report_df <- ReadLimaReport(sl7_CCS5_lima_report_file)
sl9_ccs3_report_df <- ReadLimaReport(sl9_CCS3_lima_report_file)
sl9_ccs5_report_df <- ReadLimaReport(sl9_CCS5_lima_report_file)




# Process data ------------------------------------------------------------

sl7_ccs_df <- IntegrateReportDf(sl7_ccs3_sam_df, sl7_ccs3_report_df, sl7_ccs5_report_df)
sl9_ccs_df <- IntegrateReportDf(sl9_ccs3_sam_df, sl9_ccs3_report_df, sl9_ccs5_report_df)





# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_ccs_df"),
     file = file.path(p2_R_objects_directory, "03) Read in PacBio data.RData")
     )





