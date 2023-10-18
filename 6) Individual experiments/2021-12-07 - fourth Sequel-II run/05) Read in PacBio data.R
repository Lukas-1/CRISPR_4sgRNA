### 7th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "02) Analyzing reads.R"))
source(file.path(p1_R_functions_directory, "13) Importing PacBio reads.R"))
source(file.path(s2r1_R_functions_directory, "01) Importing PacBio reads.R"))



# Define folder paths -----------------------------------------------------

s2r4_directory           <- file.path(experiments_directory, "2021-12-07 - fourth Sequel-II run")
file_input_directory     <- file.path(s2r4_directory, "2) Input")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r4_R_objects_directory <- file.path(s2r4_directory, "3) R objects")
raw_data_directory       <- file.path(file_input_directory, "Raw data")

sam_directory            <- file.path(raw_data_directory, "SAM")

pool4_ccs_sam_file       <- file.path(sam_directory, "pool4_reads.sam")
pool5_ccs_sam_file       <- file.path(sam_directory, "pool5_reads.sam")

pool4_plates_lima_file   <- file.path(raw_data_directory, "demuxed_plates",
                                      "pool4_plates_demuxed.lima.report"
                                      )
pool5_plates_lima_file   <- file.path(raw_data_directory, "demuxed_plates",
                                      "pool5_plates_demuxed.lima.report"
                                      )

pool4_wells_lima_file    <- file.path(raw_data_directory, "demuxed_wells",
                                      "pool4_wells_demuxed.lima.report"
                                      )
pool5_wells_lima_file    <- file.path(raw_data_directory, "demuxed_wells",
                                      "pool5_wells_demuxed.lima.report"
                                      )


# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2r4_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r4_R_objects_directory, "03) Import and process sgRNA sequences.RData"))




# Prepare subsets of plates_df --------------------------------------------

pool4_plates_df <- plates_df[plates_df[["Run4_pool"]] %in% 4, ]
pool5_plates_df <- plates_df[plates_df[["Run4_pool"]] %in% 5, ]
row.names(pool4_plates_df) <- NULL
row.names(pool5_plates_df) <- NULL




# Read in data ------------------------------------------------------------

pool4_ccs_sam_df <- ProcessSAM(pool4_ccs_sam_file)
pool5_ccs_sam_df <- ProcessSAM(pool5_ccs_sam_file)

pool4_plates_report_df <- ReadLimaReport(pool4_plates_lima_file)
pool5_plates_report_df <- ReadLimaReport(pool5_plates_lima_file)

pool4_wells_report_df <- ReadLimaReport(pool4_wells_lima_file)
pool5_wells_report_df <- ReadLimaReport(pool5_wells_lima_file)




# Process data ------------------------------------------------------------

pool4_ccs_df <- IntegrateReportDfs(pool4_ccs_sam_df,
                                   pool4_plates_report_df,
                                   pool4_wells_report_df,
                                   pool4_plates_df
                                   )
pool5_ccs_df <- IntegrateReportDfs(pool5_ccs_sam_df,
                                   pool5_plates_report_df,
                                   pool5_wells_report_df,
                                   pool5_plates_df
                                   )




# Prepare for merging pools -----------------------------------------------

pool4_ccs_df[["SmrtCell"]] <- "Sequel2_run4_pool4"
pool5_ccs_df[["SmrtCell"]] <- "Sequel2_run4_pool5"

pool4_ccs_df[["Run"]] <- 4L
pool5_ccs_df[["Run"]] <- 4L

pool4_ccs_df[["Pool"]] <- 4L
pool5_ccs_df[["Pool"]] <- 5L

pool4_ccs_df[["Original_ZMW"]] <- pool4_ccs_df[["ZMW"]]
pool5_ccs_df[["Original_ZMW"]] <- pool5_ccs_df[["ZMW"]]

zeros_string <- paste0(rep("0", 7), collapse = "")
pool4_ccs_df[["ZMW"]] <- as.integer(paste0(4, zeros_string)) + seq_len(nrow(pool4_ccs_df))
pool5_ccs_df[["ZMW"]] <- as.integer(paste0(5, zeros_string)) + seq_len(nrow(pool5_ccs_df))




# Merge pools -------------------------------------------------------------

ccs_df <- rbind.data.frame(pool4_ccs_df,
                           pool5_ccs_df,
                           make.row.names = FALSE,
                           stringsAsFactors = FALSE
                           )

new_columns <- c("ZMW", "Original_ZMW", "SmrtCell")
new_columns <- c(new_columns, setdiff(names(ccs_df), new_columns))
ccs_df <- ccs_df[, new_columns]




# Add a 'Well_exists' column ----------------------------------------------

ccs_df <- AddWellExistsColumn(ccs_df, library_df)




# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(s2r4_R_objects_directory, "05) Read in PacBio data.RData")
     )





