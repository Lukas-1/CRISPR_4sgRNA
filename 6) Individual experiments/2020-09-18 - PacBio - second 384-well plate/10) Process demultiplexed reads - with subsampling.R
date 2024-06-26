### 24th January 2021 ###





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





# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)






# Process reads, with subsampling -----------------------------------------

sl7_subsampled_list <- ProcessWithSubsampling(sl7_ccs_df,
                                              sl7_barcodes_df,
                                              sl7_extracted_df,
                                              wells_vec = sg_sequences_df[["Well_number"]]
                                              )
sl9_subsampled_list <- ProcessWithSubsampling(sl9_ccs_df,
                                              sl9_barcodes_df,
                                              sl9_extracted_df,
                                              wells_vec = sg_sequences_df[["Well_number"]]
                                              )




# Save data ---------------------------------------------------------------

save(list = c("sl7_subsampled_list", "sl9_subsampled_list"),
     file = file.path(p2_R_objects_directory, "10) Process demultiplexed reads - with subsampling.RData")
     )






