### 8th September 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "13) Importing PacBio reads.R"))




# Define folder paths -----------------------------------------------------

file_input_directory      <- file.path(file_directory, "2) Input")
R_objects_directory       <- file.path(file_directory, "3) R objects")

raw_data_directory        <- file.path(file_input_directory, "Raw data")
reanalysis_directory      <- file.path(raw_data_directory, "LukasAnalysis")

ccs_sam_name <- "m54073_190912_185203.ccs.sam"
sl7_CCS3_ccs_sam_file <- file.path(reanalysis_directory, "SmrtLink7_CCS3/SAM", ccs_sam_name)
sl7_CCS5_ccs_sam_file <- file.path(reanalysis_directory, "SmrtLink7_CCS5/SAM", ccs_sam_name)
sl9_CCS3_ccs_sam_file <- file.path(reanalysis_directory, "SmrtLink9_CCS3/SAM", ccs_sam_name)
sl9_CCS5_ccs_sam_file <- file.path(reanalysis_directory, "SmrtLink9_CCS5/SAM", ccs_sam_name)

sl7_CCS3_lima_report_file <- file.path(reanalysis_directory, "SmrtLink7_CCS3/lima", "lukas_lima_CCS3_SL7.lima.report")
sl7_CCS5_lima_report_file <- file.path(reanalysis_directory, "SmrtLink7_CCS5/lima", "lukas_lima_CCS5_SL7.lima.report")
sl9_CCS3_lima_report_file <- file.path(reanalysis_directory, "SmrtLink9_CCS3/lima", "lukas_lima_CCS3_SL9.lima.report")
sl9_CCS5_lima_report_file <- file.path(reanalysis_directory, "SmrtLink9_CCS5/lima", "lukas_lima_CCS5_SL9.lima.report")


# SmrtLink 9 lima GUI
# sl9_CCS3_lima_GUI_fastq_dir   <- file.path(reanalysis_directory, "SmrtLink9_CCS3/lima_GUI/call-auto_ccs_outputs_barcoded")
# sl9_CCS3_lima_GUI_report_file <- file.path(reanalysis_directory,
#                                            "SmrtLink9_CCS3/lima_GUI/call-demultiplex_barcodes",
#                                            "call-lima/demultiplex.lima.report"
#                                            )
# sl9_CCS5_lima_GUI_fastq_dir   <- file.path(reanalysis_directory, "SmrtLink9_CCS5/lima_GUI/call-auto_ccs_outputs_barcoded")
# sl9_CCS5_lima_GUI_report_file <- file.path(reanalysis_directory,
#                                            "SmrtLink9_CCS5/lima_GUI/call-demultiplex_barcodes",
#                                            "call-lima/demultiplex.lima.report"
#                                            )



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))





# Define import functions -------------------------------------------------

# ReadAndProcessFastq <- function(fastq_directory) {
#   lima_fastq_pattern <- "^demultiplex\\..+\\.fastq$"
#   readFastq_results <- ShortRead::readFastq(dirPath = fastq_directory,
#                                             pattern = lima_fastq_pattern
#                                             )
#   result_list <- list(
#     "qname" = as.character(ShortRead::id(readFastq_results)),
#     "seq"   = ShortRead::sread(readFastq_results),
#     "qual"  = PhredQuality(quality(quality(readFastq_results)))
#   )
#   return(result_list)
# }



# Read in data ------------------------------------------------------------

sl7_ccs3_sam_df <- ProcessSAM(sl7_CCS3_ccs_sam_file)
sl9_ccs3_sam_df <- ProcessSAM(sl9_CCS3_ccs_sam_file)

sl7_ccs3_report_df <- ReadLimaReport(sl7_CCS3_lima_report_file)
sl7_ccs5_report_df <- ReadLimaReport(sl7_CCS5_lima_report_file)
sl9_ccs3_report_df <- ReadLimaReport(sl9_CCS3_lima_report_file)
sl9_ccs5_report_df <- ReadLimaReport(sl9_CCS5_lima_report_file)




# Process data ------------------------------------------------------------

sl7_ccs3_ccs_sam_df <- ProcessSAM(sl7_CCS3_ccs_sam_file)
sl9_ccs3_ccs_sam_df <- ProcessSAM(sl9_CCS3_ccs_sam_file)

sl7_ccs_df <- IntegrateReportDf(sl7_ccs3_ccs_sam_df, sl7_ccs3_report_df, sl7_ccs5_report_df)
sl9_ccs_df <- IntegrateReportDf(sl9_ccs3_ccs_sam_df, sl9_ccs3_report_df, sl9_ccs5_report_df)






# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_ccs_df"),
     file = file.path(R_objects_directory, "05) Read in PacBio data.RData")
     )












