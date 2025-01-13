## 2024-01-01


# Load packages and source code -------------------------------------------

root_dir    <- "~/CRISPR_4sgRNA"
exper_dir   <- file.path(root_dir, "6) Individual experiments")
project_dir <- file.path(exper_dir, "2024-11-17 - template switches & mutations")
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R")) # For CheckThatIntegerVectorIsInOrder and UnMatrix
source(file.path(project_dir, "01_R_functions", "02_computing_error_rates.R"))



# Define paths ------------------------------------------------------------

s2rI_dir  <- file.path(exper_dir, "2021-12-08 - integrate PacBio data", "3) R objects")
rdata_dir <- file.path(project_dir, "02_PacBio_plasmids", "02_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(s2rI_dir, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))
rm(ccs3_df_list)
load(file.path(rdata_dir, "01_categorize_each_base.RData"))
load(file.path(rdata_dir, "02_categorize_subsequences_from_bases.RData"))



# Define functions --------------------------------------------------------

GetSingleBaseStats <- function(are_included) {
  results_mat <- cbind(
    "Fraction_incorrect" = colSums(!(is_correct_mat[are_included, ])) / sum(are_included),
    "Fraction_deleted"   = colSums(is_deleted_mat[are_included, ]) / sum(are_included)
  )
  return(results_mat)
}

GetSubsequenceStats <- function(are_included) {
  results_mat <- cbind(
    "Fraction_incorrect" = colSums(!(features_correct_mat[are_included, ])) / sum(are_included),
    "Fraction_deleted"   = colSums(features_deleted_mat[are_included, ]) / sum(are_included)
  )
  return(results_mat)
}


# Prepare for computing statistics ----------------------------------------

sg_pairs_mat <- GetSgPairsMat()



# Take insertions into account --------------------------------------------

is_correct_mat[has_insertion_mat] <- FALSE



# Identify correct sgRNAs -------------------------------------------------

reads_df <- ccs7_df_list[["individual_reads_df"]]
rm(ccs7_df_list)
matches_vec <- match(selected_zmws, reads_df[["ZMW"]])
stopifnot(!(anyNA(matches_vec)))

sg_categ_mat <- as.matrix(reads_df[, paste0("sg", 1:4, "_category")])[matches_vec, ]
correct_categories <- c("Correct", "Contamination", "Flanking insertion")
correct_categories <- "Correct"
sg_correct_mat <- apply(sg_categ_mat, 2, function(x) x %in% correct_categories)



# Identify reads where both of a pair of sgRNAs are present ---------------

have_both_mat <- do.call(cbind, lapply(1:6, function(x) {
  sg_numbers <- sg_pairs_mat[, x]
  rowSums(sg_correct_mat[, sg_numbers]) == 2
}))
colnames(have_both_mat) <- colnames(sg_pairs_mat)



# Compute error rates for individual bases --------------------------------
## Statistics are computed for "fully mapped" reads and each pairing

have_all <- rowSums(sg_correct_mat) == 4

base_error_mat_list <- c(
  list(GetSingleBaseStats(have_all)),
  lapply(1:6, function(x) GetSingleBaseStats(have_both_mat[, x]))
)
names(base_error_mat_list) <- c("Full", colnames(sg_pairs_mat))



# Compute error rates for subsequences ------------------------------------

features_correct_mat <- CategorDfToMat(feature_categ_df, "Is_correct")
features_deleted_mat <- CategorDfToMat(feature_categ_df, "Mostly_deleted")
subsequence_error_mat_list <- c(
  list(GetSubsequenceStats(have_all)),
  lapply(1:6, function(x) GetSubsequenceStats(have_both_mat[, x]))
)
names(subsequence_error_mat_list) <- c("Full", colnames(sg_pairs_mat))



# Save data ---------------------------------------------------------------

save(base_error_mat_list,
     file = file.path(rdata_dir, "03_compute_error_rates__individual_bases.RData")
     )
save(subsequence_error_mat_list,
     file = file.path(rdata_dir, "03_compute_error_rates__subsequences.RData")
     )

