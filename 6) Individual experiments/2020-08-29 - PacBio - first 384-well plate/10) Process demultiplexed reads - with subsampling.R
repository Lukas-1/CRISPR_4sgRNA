### 24th January 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))
source(file.path(R_functions_directory, "04) Exporting FASTA and FASTQ files.R"))




# Define folder paths -----------------------------------------------------

file_input_directory    <- file.path(file_directory, "2) Input")
R_objects_directory     <- file.path(file_directory, "3) R objects")
file_output_directory   <- file.path(file_directory, "5) Output")
tables_output_directory <- file.path(file_output_directory, "Tables")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(R_objects_directory, "07) Extract barcode sequences and quality scores.RData"))
load(file.path(R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))





# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)







# Process reads, with subsampling -----------------------------------------

sl7_subsampled_list <- ProcessWithSubsampling(sl7_ccs_df,
                                              sl7_barcodes_df,
                                              sl7_extracted_df
                                              )
sl9_subsampled_list <- ProcessWithSubsampling(sl9_ccs_df,
                                              sl9_barcodes_df,
                                              sl9_extracted_df
                                              )




# Save data ---------------------------------------------------------------

save(list = c("sl7_subsampled_list", "sl9_subsampled_list"),
     file = file.path(R_objects_directory, "10) Process demultiplexed reads - with subsampling.RData")
     )









