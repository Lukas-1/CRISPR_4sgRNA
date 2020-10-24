### 20th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "06) Extracting barcodes.R"))



# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - raw sequences.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs3.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs3.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs5.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs5.RData"))
load(file.path(p2_R_objects_directory, "04) Perform pairwise alignments with the reference sequence.RData"))




# Extract barcodes --------------------------------------------------------

############################
### DELETE THIS LATER!!! ###
############################
sl7_ccs3_ccs <- sl7_ccs5_ccs
sl7_ccs3_lima <- sl7_ccs5_lima
sl7_ccs3_report_df <- sl7_ccs5_report_df
############################
############################
############################

sl7_barcodes_df <- GetBarcodes(use_sl7 = TRUE, wells_vec = sg_sequences_df[["Well_number"]])
# sl9_barcodes_df <- GetBarcodes(use_sl7 = FALSE)





# Save data ---------------------------------------------------------------

# save(list = paste0(c("sl7", "sl9"), "_barcodes_df"),
#      file = file.path(p2_R_objects_directory, "06) Extract barcode sequences and quality scores.RData")
#      )

save(list = paste0(c("sl7"), "_barcodes_df"),
     file = file.path(p2_R_objects_directory, "05) Extract barcode sequences and quality scores.RData")
     )









