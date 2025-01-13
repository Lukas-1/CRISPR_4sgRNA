## 2025-01-01


# Load packages and source code -------------------------------------------

root_dir     <- "~/CRISPR_4sgRNA"
exper_dir    <- file.path(root_dir, "6) Individual experiments")
nanopore_dir <- file.path(exper_dir, "2022-01-05 - first nanopore sequencing run")
project_dir  <- file.path(exper_dir, "2024-11-17 - template switches & mutations")

source(file.path(nanopore_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "02_PacBio_plasmids", "02_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_categorize_each_base.RData"))



# Prepare features_df -----------------------------------------------------

features_df <- TweakFeaturesDf(FeaturesListToDf(features_list))



# Take insertions into account --------------------------------------------

is_correct_mat[has_insertion_mat] <- FALSE



# Categorize subsequences -------------------------------------------------

feature_categ_df <- CategorizeFeaturesFromBases(is_correct_mat, is_deleted_mat)
feature_categ_df <- data.frame("ZMW" = rep(selected_zmws, each = nrow(features_df)),
                               feature_categ_df
                               )



# Save data ---------------------------------------------------------------

save(feature_categ_df,
     file = file.path(rdata_dir, "02_categorize_subsequences_from_bases.RData")
     )


