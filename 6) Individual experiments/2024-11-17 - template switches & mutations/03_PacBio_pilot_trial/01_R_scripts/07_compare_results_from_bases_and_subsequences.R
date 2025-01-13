## 2024-12-30


# Load packages and source code -------------------------------------------

root_dir     <- "~/CRISPR_4sgRNA"
exper_dir    <- file.path(root_dir, "6) Individual experiments")
project_dir  <- file.path(exper_dir, "2024-11-17 - template switches & mutations")
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R"))
# source(file.path(project_dir, "01_R_functions", "02_computing_error_rates.R"))



# Define paths ------------------------------------------------------------

rdata_dir       <- file.path(project_dir, "03_PacBio_pilot_trial", "02_R_objects")
first_rdata_dir <- file.path(exper_dir, "2022-04-06 - PacBio pooled 4sg - first trial", "03_R_objects")



# Load data ---------------------------------------------------------------

# load(file.path(first_rdata_dir, "07_assign_sgRNAs_to_plasmids.RData"))
load(file.path(rdata_dir, "01_extract_and_categorize_subsequences__features_df.RData"))
load(file.path(rdata_dir, "01_extract_and_categorize_subsequences__categorized_df.RData"))
load(file.path(rdata_dir, "04_categorize_each_base.RData"))



# Take insertions into account --------------------------------------------

is_correct_mat[has_insertion_mat] <- FALSE




# Categorize subsequences -------------------------------------------------

feature_categ_df <- CategorizeFeaturesFromBases(is_correct_mat, is_deleted_mat)



nrow(feature_categ_df)
nrow(categorized_df)

table(feature_categ_df[, "Feature"] == categorized_df[, "Feature"])

df1 <- feature_categ_df[!(feature_categ_df[, "Feature"] %in% paste0("sg", 1:4)), ]
df2 <- categorized_df[!(categorized_df[, "Feature"] %in% paste0("sg", 1:4)), ]

are_different <- df1[, "Is_correct"] != df2[, "Is_correct"]

table(df1[, "Is_correct"])
table(df2[, "Is_correct"])


table(df1[, "Is_correct"], df2[, "Is_correct"])





