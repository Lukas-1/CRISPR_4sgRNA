### 26th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "02) Analyzing reads.R"))
source(file.path(p1_R_functions_directory, "13) Importing PacBio reads.R"))
source(file.path(s2r1_R_functions_directory, "01) Importing PacBio reads.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
file_input_directory     <- file.path(s2r2_directory, "2) Input")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
raw_data_directory       <- file.path(file_input_directory, "Raw data")

sam_directory      <- file.path(raw_data_directory, "SAM")

pool1_ccs_sam_file <- file.path(sam_directory, "sql2_run2_pool1_reads.sam")
pool2_ccs_sam_file <- file.path(sam_directory, "sql2_run2_pool2_reads.sam")

pool1_plates_lima_file <- file.path(raw_data_directory, "demuxed_plates",
                                    "pool1", "sql2_run2_pool1_plates_demuxed.lima.report"
                                    )
pool2_plates_lima_file <- file.path(raw_data_directory, "demuxed_plates",
                                    "pool2", "sql2_run2_pool2_plates_demuxed.lima.report"
                                    )
pool1_wells_lima_file <- file.path(raw_data_directory, "demuxed_wells",
                                   "pool1", "sql2_run2_pool1_wells_demuxed.lima.report"
                                   )
pool2_wells_lima_file <- file.path(raw_data_directory, "demuxed_wells",
                                   "pool2", "sql2_run2_pool2_wells_demuxed.lima.report"
                                   )



# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))





# Prepare subsets of plates_df --------------------------------------------

pool1_plates_df <- plates_df[plates_df[["Run2_pool"]] %in% 1, ]
pool2_plates_df <- plates_df[plates_df[["Run2_pool"]] %in% 2, ]
row.names(pool1_plates_df) <- NULL
row.names(pool2_plates_df) <- NULL




# Read in data ------------------------------------------------------------

pool1_ccs_sam_df <- ProcessSAM(pool1_ccs_sam_file)
pool2_ccs_sam_df <- ProcessSAM(pool2_ccs_sam_file)

pool1_plates_report_df <- ReadLimaReport(pool1_plates_lima_file)
pool2_plates_report_df <- ReadLimaReport(pool2_plates_lima_file)

pool1_wells_report_df <- ReadLimaReport(pool1_wells_lima_file)
pool2_wells_report_df <- ReadLimaReport(pool2_wells_lima_file)




# Process data ------------------------------------------------------------

pool1_ccs_df <- IntegrateReportDfs(pool1_ccs_sam_df,
                                   pool1_plates_report_df,
                                   pool1_wells_report_df,
                                   pool1_plates_df
                                   )
pool2_ccs_df <- IntegrateReportDfs(pool2_ccs_sam_df,
                                   pool2_plates_report_df,
                                   pool2_wells_report_df,
                                   pool2_plates_df
                                   )



# Merge SMRT cells --------------------------------------------------------

pool1_ccs_df[["SmrtCell"]] <- "Sequel2_run2_pool1"
pool2_ccs_df[["SmrtCell"]] <- "Sequel2_run2_pool2"

pool1_ccs_df[["Original_ZMW"]] <- pool1_ccs_df[["ZMW"]]
pool2_ccs_df[["Original_ZMW"]] <- pool2_ccs_df[["ZMW"]]

all_zmws <- c(pool1_ccs_df[["ZMW"]], pool2_ccs_df[["ZMW"]])
num_chars_zmw <- nchar(as.character(max(all_zmws)))
add_to_zmws <- as.integer(paste0("1", paste0(rep("0", num_chars_zmw), collapse = "")))

new_zmws <- seq_len(nrow(pool1_ccs_df) + nrow(pool2_ccs_df)) + add_to_zmws

ccs_df <- rbind.data.frame(pool1_ccs_df,
                           pool2_ccs_df,
                           make.row.names = FALSE,
                           stringsAsFactors = FALSE
                           )
ccs_df[["ZMW"]] <- new_zmws

new_columns <- c("ZMW", "Original_ZMW", "SmrtCell")
new_columns <- c(new_columns, setdiff(names(ccs_df), new_columns))
ccs_df <- ccs_df[, new_columns]




# Add a 'Well_exists' column ----------------------------------------------

ccs_df <- AddWellExistsColumn(ccs_df, library_df)





# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(s2r2_R_objects_directory, "05) Read in PacBio data.RData")
     )





