### 9th April 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
plate1_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
sql2_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
sql2_R_functions_directory <- file.path(sql2_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "02) Analyzing reads.R"))
source(file.path(p1_R_functions_directory, "13) Importing PacBio reads.R"))
source(file.path(sql2_R_functions_directory, "01) Importing PacBio reads.R"))




# Define folder paths -----------------------------------------------------

file_input_directory     <- file.path(sql2_directory, "2) Input")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
raw_data_directory       <- file.path(file_input_directory, "Raw data")

sam_directory    <- file.path(raw_data_directory, "lima", "SAM")
ccs_sam_file     <- file.path(sam_directory, "m64141e_210227_112845.sam")

plates_lima_report_file <- file.path(raw_data_directory, "lima", "demuxed_plates",
                                     "FirstSQ2_plates_demuxed.lima.report"
                                     )
wells_lima_report_file <- file.path(raw_data_directory, "lima", "demuxed_wells",
                                    "FirstSQ2_wells_demuxed.lima.report"
                                    )




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))





# Read in data ------------------------------------------------------------

ccs_sam_df <- ProcessSAM(ccs_sam_file)
plates_lima_report_df <- ReadLimaReport(plates_lima_report_file)
wells_lima_report_df <- ReadLimaReport(wells_lima_report_file)




# Process data ------------------------------------------------------------

ccs_df <- IntegrateReportDfs(ccs_sam_df,
                             plates_lima_report_df,
                             wells_lima_report_df,
                             plates_df
                             )





# Add a 'Well_exists' column ----------------------------------------------

ccs_df <- AddWellExistsColumn(ccs_df, library_df)





# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(sql2_R_objects_directory, "05) Read in PacBio data.RData")
     )





