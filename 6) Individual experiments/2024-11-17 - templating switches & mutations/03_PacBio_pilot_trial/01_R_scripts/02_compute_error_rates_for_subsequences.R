## 2024-11-18


# Load packages and source code -------------------------------------------

root_dir    <- "~/CRISPR_4sgRNA"
exper_dir   <- file.path(root_dir, "6) Individual experiments")
project_dir <- file.path(exper_dir, "2024-11-17 - templating switches & mutations")
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R")) # For CheckThatIntegerVectorIsInOrder
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_sequences.R")) # For CheckThatIntegerVectorIsInOrder
source(file.path(project_dir, "01_R_functions", "02_computing_error_rates_for_subsequences.R"))



# Define paths ------------------------------------------------------------

rdata_dir       <- file.path(project_dir, "03_PacBio_pilot_trial", "02_R_objects")
first_rdata_dir <- file.path(exper_dir, "2022-04-06 - PacBio pooled 4sg - first trial", "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_extract_and_categorize_sequences__categorized_df.RData"))
load(file.path(first_rdata_dir, "07_assign_sgRNAs_to_plasmids.RData"))



# Prepare for computing statistics ----------------------------------------

sg_pairs_mat <- GetSgPairsMat()



# Compute statistics for all pairs of sgRNAs ------------------------------

message("Computing statistics using only fully mapped reads...")
full_reads_errors_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing error rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y,
                             include_reads = "fully mapped"
                             )
})
names(full_reads_errors_df_list) <- colnames(sg_pairs_mat)

full_reads_deletions_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing deletion rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y,
                             only_deletions = TRUE,
                             include_reads = "fully mapped"
                             )
})
names(full_reads_deletions_df_list) <- colnames(sg_pairs_mat)


message("\n\nComputing statistics using all reads...")
all_reads_errors_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing error rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y)
})
names(all_reads_errors_df_list) <- colnames(sg_pairs_mat)


all_reads_deletions_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing deletion rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y,
                             only_deletions = TRUE
                             )
})
names(all_reads_deletions_df_list) <- paste0("sg", sg_pairs_mat[1, ], "_sg", sg_pairs_mat[2, ])



# Explore results ---------------------------------------------------------

use_columns <- c("Fraction_incorrect_switched", "Fraction_incorrect_nonswitched")
vapply(full_reads_errors_df_list,    function(x) max(x[, use_columns]), numeric(1))
vapply(full_reads_deletions_df_list, function(x) max(x[, use_columns]), numeric(1))
vapply(all_reads_errors_df_list,     function(x) max(x[, use_columns]), numeric(1))
vapply(all_reads_deletions_df_list,  function(x) max(x[, use_columns]), numeric(1))



# Save data ---------------------------------------------------------------

save(list = c("full_reads_errors_df_list", "full_reads_deletions_df_list",
              "all_reads_errors_df_list", "all_reads_deletions_df_list"
              ),
     file = file.path(rdata_dir, "02_compute_error_rates.RData")
     )

