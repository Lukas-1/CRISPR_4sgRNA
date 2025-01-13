## 2024-12-30


# Load packages and source code -------------------------------------------

root_dir    <- "~/CRISPR_4sgRNA"
exper_dir   <- file.path(root_dir, "6) Individual experiments")
project_dir <- file.path(exper_dir, "2024-11-17 - template switches & mutations")
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R")) # For CheckThatIntegerVectorIsInOrder and UnMatrix
source(file.path(project_dir, "01_R_functions", "02_computing_error_rates.R"))



# Define paths ------------------------------------------------------------

rdata_dir       <- file.path(project_dir, "03_PacBio_pilot_trial", "02_R_objects")
first_rdata_dir <- file.path(exper_dir, "2022-04-06 - PacBio pooled 4sg - first trial", "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "07_assign_sgRNAs_to_plasmids.RData"))
load(file.path(rdata_dir, "01_extract_and_categorize_subsequences__categorized_df.RData"))
load(file.path(rdata_dir, "04_categorize_each_base.RData"))



# Prepare for computing statistics ----------------------------------------

sg_pairs_mat <- GetSgPairsMat()




# Take insertions into account --------------------------------------------

is_correct_mat[has_insertion_mat] <- FALSE



# Replicate categorized_df ------------------------------------------------

all_read_numbers <- unique(categorized_df[, "Read_number"])
num_reads <- length(all_read_numbers)

categorized_df <- data.frame(
  "Read_number"    = rep(all_read_numbers, each = 2225),
  "Feature"        = rep(seq_len(2225), times = num_reads),
  "Mostly_deleted" = UnMatrix(t(is_deleted_mat)), # For compatibility, it's actually 100% deleted since it's only one base
  "Is_correct"     = UnMatrix(t(is_correct_mat))
)



# Perform checks ----------------------------------------------------------

test_mat <- CategorDfToMat(categorized_df, "Is_correct")
dimnames(test_mat) <- NULL
stopifnot(identical(test_mat, is_correct_mat))

test_mat <- CategorDfToMat(categorized_df, "Mostly_deleted")
dimnames(test_mat) <- NULL
stopifnot(identical(test_mat, is_deleted_mat))

rm(test_mat)



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

full_reads_deletions_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing deletion rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y,
                             only_deletions = TRUE,
                             include_reads = "fully mapped"
                             )
})


message("\n\nComputing statistics using all reads...")
all_reads_errors_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing error rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y)
})


all_reads_deletions_df_list <- lapply(1:6, function(x) {
  sg_X <- sg_pairs_mat[1, x]
  sg_Y <- sg_pairs_mat[2, x]
  message("Computing deletion rates for sg", sg_X, " and sg", sg_Y, "...")
  CompareSwitchedNonswitched(categorized_df, pb_df, sg_X, sg_Y,
                             only_deletions = TRUE
                             )
})




# Modify results ----------------------------------------------------------

df_list_names <- c(
  "full_reads_errors_df_list", "full_reads_deletions_df_list",
  "all_reads_errors_df_list", "all_reads_deletions_df_list"
)
for (list_name in df_list_names) {
  new_list <- get(list_name)
  names(new_list) <- colnames(sg_pairs_mat)
  new_list <- lapply(new_list, function(x) {
    names(x)[names(x) == "Feature"] <- "Base_number"
    x[, "Base_number"] <- as.integer(x[, "Base_number"])
    x
  })
  assign(list_name, new_list)
}



# Explore results ---------------------------------------------------------

use_columns <- c("Fraction_incorrect_switched", "Fraction_incorrect_nonswitched")
vapply(full_reads_errors_df_list,    function(x) max(x[, use_columns], na.rm = TRUE), numeric(1))
vapply(full_reads_deletions_df_list, function(x) max(x[, use_columns], na.rm = TRUE), numeric(1))
vapply(all_reads_errors_df_list,     function(x) max(x[, use_columns], na.rm = TRUE), numeric(1))
vapply(all_reads_deletions_df_list,  function(x) max(x[, use_columns], na.rm = TRUE), numeric(1))



# Save data ---------------------------------------------------------------

save(list = df_list_names,
     file = file.path(rdata_dir, "05_compute_error_rates_for_each_base.RData")
     )

