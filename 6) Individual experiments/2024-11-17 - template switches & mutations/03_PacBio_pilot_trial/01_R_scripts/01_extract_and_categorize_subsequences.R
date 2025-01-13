## 2024-11-17


# Load packages and source code -------------------------------------------

root_dir <- "~/CRISPR_4sgRNA"
exper_dir <- file.path(root_dir, "6) Individual experiments")

first_QC_dir <- file.path(exper_dir, "2020-08-29 - PacBio - first 384-well plate")
nanopore_dir <- file.path(exper_dir, "2022-01-05 - first nanopore sequencing run")
project_dir  <- file.path(exper_dir, "2024-11-17 - template switches & mutations")

source(file.path(first_QC_dir, "1) R functions", "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(nanopore_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "03_PacBio_pilot_trial", "02_R_objects")
first_rdata_dir <- file.path(exper_dir, "2022-04-06 - PacBio pooled 4sg - first trial", "03_R_objects")



# Read in data ------------------------------------------------------------

amplicon_ref <- read.table(file.path(nanopore_dir, "02_input_data", "amplicon_4sg.txt"),
                           quote = "", stringsAsFactors = FALSE
                           )[, 1]



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "02_align_reads.RData"))
load(file.path(first_rdata_dir, "07_assign_sgRNAs_to_plasmids.RData"))



# Prepare features_df -----------------------------------------------------

features_df <- TweakFeaturesDf(FeaturesListToDf(features_list))
features_indices_list <- lapply(seq_len(nrow(features_df)),
                                function(y) features_df[y, "Start"]:features_df[y, "End"]
                                )


# Extract aligned sequences -----------------------------------------------

extracted_df <- ExtractAllFeatures(FilterAlignmentsDf(alignments_df, pb_df))



# Categorize extracted sequences ------------------------------------------

categorized_df <- AddThreeCategories(extracted_df, features_df)



# Incorporate data on sgRNAs ----------------------------------------------

categorized_df <- AddGuideData(categorized_df, pb_df)



# Save data ---------------------------------------------------------------

save(categorized_df, file = file.path(rdata_dir, "01_extract_and_categorize_subsequences__categorized_df.RData"))
save(features_df, file = file.path(rdata_dir, "01_extract_and_categorize_subsequences__features_df.RData"))


