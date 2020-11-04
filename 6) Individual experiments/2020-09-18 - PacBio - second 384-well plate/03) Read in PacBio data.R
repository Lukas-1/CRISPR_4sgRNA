### 20th October 2020 ###



# Import packages and source code -----------------------------------------

library("Rsamtools")

CRISPR_root_directory  <- "~/CRISPR"
plate1_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))





# Define folder paths -----------------------------------------------------

plate2_directory          <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")

p1_R_objects_directory    <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory    <- file.path(plate2_directory, "2) R objects")

file_input_directory      <- file.path(plate2_directory, "1) Input")
raw_data_directory        <- file.path(file_input_directory, "Raw data", "Analysis")

ccs_bam_name <- "m54073_200523_132654.ccs.bam"
ccs_sl7_sub_path <- file.path("ccs/tasks/pbcoretools.tasks.auto_ccs_outputs-0", ccs_bam_name)
sl7_CCS3_ccs_bam_file     <- file.path(raw_data_directory, "SmrtLink7_CCS3", ccs_sl7_sub_path)
sl7_CCS5_ccs_bam_file     <- file.path(raw_data_directory, "SmrtLink7_CCS5", ccs_sl7_sub_path)
sl9_CCS3_ccs_bam_file     <- file.path(raw_data_directory, "SmrtLink9_CCS3/ccs", ccs_bam_name)
sl9_CCS5_ccs_bam_file     <- file.path(raw_data_directory, "SmrtLink9_CCS5/ccs", ccs_bam_name)

sl7_CCS3_lima_bam_file    <- file.path(raw_data_directory, "SmrtLink7_CCS3/lima", "lukas_lima_CCS3_SL7_output.bam")
sl7_CCS5_lima_bam_file    <- file.path(raw_data_directory, "SmrtLink7_CCS5/lima", "lukas_lima_CCS5_SL7_output.bam")
sl9_CCS3_lima_bam_file    <- file.path(raw_data_directory, "SmrtLink9_CCS3/lima", "lukas_lima_CCS3_SL9.bam")
sl9_CCS5_lima_bam_file    <- file.path(raw_data_directory, "SmrtLink9_CCS5/lima", "lukas_lima_CCS5_SL9.bam")

sl7_CCS3_lima_report_file <- file.path(raw_data_directory, "SmrtLink7_CCS3/lima", "lukas_lima_CCS3_SL7_output.lima.report")
sl7_CCS5_lima_report_file <- file.path(raw_data_directory, "SmrtLink7_CCS5/lima", "lukas_lima_CCS5_SL7_output.lima.report")
sl9_CCS3_lima_report_file <- file.path(raw_data_directory, "SmrtLink9_CCS3/lima", "lukas_lima_CCS3_SL9.lima.report")
sl9_CCS5_lima_report_file <- file.path(raw_data_directory, "SmrtLink9_CCS5/lima", "lukas_lima_CCS5_SL9.lima.report")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))






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

sl7_ccs5_lima_zmws <- ExtractZMWs(sl7_ccs5_report_df, sl7_ccs5_lima,
                                  wells_vec = sg_sequences_df[["Well_number"]]
                                  )
sl9_ccs5_lima_zmws <- ExtractZMWs(sl9_ccs5_report_df, sl9_ccs5_lima,
                                  wells_vec = sg_sequences_df[["Well_number"]]
                                  )




# Save data ---------------------------------------------------------------

ccs3_combos <- paste0("sl", c(7, 9), "_ccs3_")
ccs5_combos <- paste0("sl", c(7, 9), "_ccs5_")


save(list = paste0(ccs3_combos, "ccs"),
     file = file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs3.RData")
     )

save(list = c(paste0(ccs3_combos, "lima"),
              paste0(ccs3_combos, "report_df")
              ),
     file = file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs3.RData")
     )

save(list = paste0(ccs5_combos, "lima_zmws"),
     file = file.path(p2_R_objects_directory, "03) Read in PacBio data - ccs5 ZMWs.RData")
     )

save(list = paste0(ccs5_combos, "ccs"),
     file = file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs5.RData")
     )

save(list = c(paste0(ccs5_combos, "lima"),
              paste0(ccs5_combos, "report_df")
              ),
     file = file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs5.RData")
     )





