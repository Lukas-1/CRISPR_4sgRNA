### 8th September 2020 ###



# Import packages and source code -----------------------------------------

library("Rsamtools")

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "2) Functions for analyzing reads.R"))





# Define folder paths -----------------------------------------------------

file_input_directory      <- file.path(file_directory, "2) Input")
R_objects_directory       <- file.path(file_directory, "3) R objects")

raw_data_directory        <- file.path(file_input_directory, "Raw data")
reanalysis_directory      <- file.path(raw_data_directory, "LukasAnalysis")

ccs_bam_name <- "m54073_190912_185203.ccs.bam"
ccs_sl7_sub_path <- file.path("ccs/tasks/pbcoretools.tasks.auto_ccs_outputs-0", ccs_bam_name)
sl7_CCS3_ccs_bam_file     <- file.path(reanalysis_directory, "SmrtLink7_CCS3", ccs_sl7_sub_path)
sl7_CCS5_ccs_bam_file     <- file.path(reanalysis_directory, "SmrtLink7_CCS5", ccs_sl7_sub_path)
sl9_CCS3_ccs_bam_file     <- file.path(reanalysis_directory, "SmrtLink9_CCS3/ccs", ccs_bam_name)
sl9_CCS5_ccs_bam_file     <- file.path(reanalysis_directory, "SmrtLink9_CCS5/ccs", ccs_bam_name)

sl7_CCS3_lima_bam_file    <- file.path(reanalysis_directory, "SmrtLink7_CCS3/lima", "lukaslimaCCS3SL7output.bam")
sl7_CCS5_lima_bam_file    <- file.path(reanalysis_directory, "SmrtLink7_CCS5/lima", "lukaslimaCCS5SL7output.bam")
sl9_CCS3_lima_bam_file    <- file.path(reanalysis_directory, "SmrtLink9_CCS3/lima", "lukas_lima_CS3_SL9.bam")
sl9_CCS5_lima_bam_file    <- file.path(reanalysis_directory, "SmrtLink9_CCS5/lima", "lukas_lima_CS5_SL9.bam")

sl7_CCS3_lima_report_file <- file.path(reanalysis_directory, "SmrtLink7_CCS3/lima", "lukaslimaCCS3SL7output.lima.report")
sl7_CCS5_lima_report_file <- file.path(reanalysis_directory, "SmrtLink7_CCS5/lima", "lukaslimaCCS5SL7output.lima.report")
sl9_CCS3_lima_report_file <- file.path(reanalysis_directory, "SmrtLink9_CCS3/lima", "lukas_lima_CS3_SL9.lima.report")
sl9_CCS5_lima_report_file <- file.path(reanalysis_directory, "SmrtLink9_CCS5/lima", "lukas_lima_CS5_SL9.lima.report")


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

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))




# Define functions --------------------------------------------------------

ExtractZMWs <- function(report_df, lima_list) {
  ccs_well_numbers <- GetWellNumbers(report_df)
  report_zmws <- as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]])))
  lima_zmws <- as.integer(substr(lima_list[["qname"]], 22, nchar(lima_list[["qname"]]) - 4))

  result_zmws <- report_zmws[!(is.na(ccs_well_numbers)) & (report_zmws %in% lima_zmws)]
  return(result_zmws)
}




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

ReadLimaReport <- function(file_path) {
  read.table(file_path,
             header = TRUE, sep = "\t",
             quote = "", stringsAsFactors = FALSE
             )
}



# Read in data ------------------------------------------------------------

use_params <- c("qname", "seq", "qual")


sl7_ccs3_ccs <- scanBam(sl7_CCS3_ccs_bam_file,
                        param = ScanBamParam(what = use_params)
                        )[[1]]
sl7_ccs5_ccs <- scanBam(sl7_CCS5_ccs_bam_file,
                        param = ScanBamParam(what = use_params)
                        )[[1]]
sl9_ccs3_ccs <- scanBam(sl9_CCS3_ccs_bam_file,
                        param = ScanBamParam(what = use_params)
                        )[[1]]
sl9_ccs5_ccs <- scanBam(sl9_CCS5_ccs_bam_file,
                        param = ScanBamParam(what = use_params)
                        )[[1]]


sl7_ccs3_lima <- scanBam(sl7_CCS3_lima_bam_file,
                         param = ScanBamParam(what = use_params)
                         )[[1]]
sl7_ccs5_lima <- scanBam(sl7_CCS5_lima_bam_file,
                         param = ScanBamParam(what = use_params)
                         )[[1]]

sl9_ccs3_lima <- scanBam(sl9_CCS3_lima_bam_file,
                         param = ScanBamParam(what = use_params)
                         )[[1]]
sl9_ccs5_lima <- scanBam(sl9_CCS5_lima_bam_file,
                         param = ScanBamParam(what = use_params)
                         )[[1]]

sl7_ccs3_report_df <- ReadLimaReport(sl7_CCS3_lima_report_file)
sl7_ccs5_report_df <- ReadLimaReport(sl7_CCS5_lima_report_file)
sl9_ccs3_report_df <- ReadLimaReport(sl9_CCS3_lima_report_file)
sl9_ccs5_report_df <- ReadLimaReport(sl9_CCS5_lima_report_file)






# Extract the ZMWs that are retained in CCS5 ------------------------------

sl7_ccs5_lima_zmws <- ExtractZMWs(sl7_ccs5_report_df, sl7_ccs5_lima)
sl9_ccs5_lima_zmws <- ExtractZMWs(sl9_ccs5_report_df, sl9_ccs5_lima)







# Save data ---------------------------------------------------------------

ccs3_combos <- paste0("sl", c(7, 9), "_ccs3_")
ccs5_combos <- paste0("sl", c(7, 9), "_ccs5_")

save(list = paste0(ccs3_combos, "ccs"),
     file = file.path(R_objects_directory, "5) Read in PacBio data - consensus reads - ccs3.RData")
     )

save(list = c(paste0(ccs3_combos, "lima"),
              paste0(ccs3_combos, "report_df")
              ),
     file = file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed - ccs3.RData")
     )



save(list = c("sl7_ccs5_lima_zmws",  "sl9_ccs5_lima_zmws"),
     file = file.path(R_objects_directory, "5) Read in PacBio data - ccs5 ZMWs.RData")
     )



save(list = paste0(ccs5_combos, "ccs"),
     file = file.path(R_objects_directory, "5) Read in PacBio data - consensus reads - ccs5.RData")
     )

save(list = c(paste0(ccs5_combos, "lima"),
              paste0(ccs5_combos, "report_df")
              ),
     file = file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed - ccs5.RData")
     )









