### 14th September 2021 ###



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

s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")
file_input_directory     <- file.path(s2r3_directory, "2) Input")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")
raw_data_directory       <- file.path(file_input_directory, "Raw data")

sam_directory            <- file.path(raw_data_directory, "SAM")
pool3_ccs_sam_file       <- file.path(sam_directory, "run3_pool3_reads.sam")
pool3_plates_lima_file   <- file.path(raw_data_directory, "demuxed_plates",
                                      "run3_pool3_plates_demuxed.lima.report"
                                      )
pool3_wells_lima_file    <- file.path(raw_data_directory, "demuxed_wells",
                                      "run3_pool3_wells_demuxed.lima.report"
                                      )



# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2r3_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r3_R_objects_directory, "03) Import and process sgRNA sequences.RData"))



# Read in data ------------------------------------------------------------

pool3_ccs_sam_df       <- ProcessSAM(pool3_ccs_sam_file)
pool3_plates_report_df <- ReadLimaReport(pool3_plates_lima_file)
pool3_wells_report_df  <- ReadLimaReport(pool3_wells_lima_file)



# Process data ------------------------------------------------------------

ccs_df <- IntegrateReportDfs(pool3_ccs_sam_df,
                             pool3_plates_report_df,
                             pool3_wells_report_df,
                             plates_df
                             )


# Merge SMRT cells --------------------------------------------------------

ccs_df[["SmrtCell"]] <- "Sequel2_run3_pool3"

new_columns <- c("ZMW", "SmrtCell")
new_columns <- c(new_columns, setdiff(names(ccs_df), new_columns))
ccs_df <- ccs_df[, new_columns]



# Add a 'Well_exists' column ----------------------------------------------

ccs_df <- AddWellExistsColumn(ccs_df, library_df)




# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(s2r3_R_objects_directory, "05) Read in PacBio data.RData")
     )





