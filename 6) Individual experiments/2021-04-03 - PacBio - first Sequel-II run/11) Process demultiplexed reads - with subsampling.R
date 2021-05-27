### 26th May 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")





# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(sql2_R_objects_directory, "07) Extract barcode sequences and quality scores.RData"))
load(file.path(sql2_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))





# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)





# Process reads, with subsampling -----------------------------------------

subsampled_list <- ProcessWithSubsampling(ccs_df,
                                          barcodes_df,
                                          extracted_df,
                                          use_fractions = c(1, 0.5, 0.25, 0.1, 0.05, 0.01),
                                          num_repetitions = 3L
                                          )




# Save data ---------------------------------------------------------------

save(list = c("sl7_subsampled_list", "sl9_subsampled_list"),
     file = file.path(R_objects_directory, "11) Process demultiplexed reads - with subsampling.RData")
     )





