### 21st October 2020 ###



# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "06) Extract barcode sequences and quality scores.RData"))
load(file.path(p2_R_objects_directory, "08) Process demultiplexed PacBio reads.RData"))





# Check for outliers ------------------------------------------------------

summary_df <- sl7_ccs5_df_list[["original_summary_df"]]

counts_mat <- as.matrix(summary_df[, paste0("Count_sg", 1:4, "_cr", 1:4)])
have_low_counts <- apply(counts_mat < 20, 1, any)

summary_df[have_low_counts, ]











