## 2024-12-30


# Load packages and source code -------------------------------------------

root_dir    <- "~/CRISPR_4sgRNA"
exper_dir   <- file.path(root_dir, "6) Individual experiments")
project_dir <- file.path(exper_dir, "2024-11-17 - template switches & mutations")

source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R"))



# Define paths ------------------------------------------------------------

first_rdata_dir <- file.path(exper_dir, "2022-04-06 - PacBio pooled 4sg - first trial", "03_R_objects")
nanopore_dir    <- file.path(exper_dir, "2022-01-05 - first nanopore sequencing run")
rdata_dir       <- file.path(project_dir, "03_PacBio_pilot_trial", "02_R_objects")



# Read in data ------------------------------------------------------------

amplicon_ref <- read.table(file.path(nanopore_dir, "02_input_data", "amplicon_4sg.txt"),
                           quote = "", stringsAsFactors = FALSE
                           )[, 1]



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "02_align_reads.RData"))
load(file.path(first_rdata_dir, "07_assign_sgRNAs_to_plasmids.RData"))



# Extract aligned bases ---------------------------------------------------

filtered_df <- FilterAlignmentsDf(alignments_df, pb_df)

stopifnot(identical(nchar(filtered_df[, "Aligned_ref"]), nchar(filtered_df[, "Aligned_read"])))

ref_char_list <- strsplit(filtered_df[, "Aligned_ref"], "", fixed = TRUE)
read_char_list <- strsplit(filtered_df[, "Aligned_read"], "", fixed = TRUE)

ref_are_gaps_list <- lapply(ref_char_list, function(x) x == "-")
ref_char_numbers_list <- lapply(ref_are_gaps_list, function(x) cumsum(!(x)))

use_indices <- seq_len(2225)

indices_vec_list <- lapply(use_indices, function(x) {
  message("Finding the aligned index at position #", x, "...")
  vapply(ref_char_numbers_list, function(y) which(y == x)[[1]], integer(1))
})

char_vec_list <- lapply(seq_along(use_indices), function(x) {
  message("Extracting the base at position #", x, "...")
  mapply(function(y, z) y[[z]], read_char_list, indices_vec_list[[x]])
})



# Compare to the reference ------------------------------------------------

amplicon_ref_vec <- strsplit(amplicon_ref, "", fixed = TRUE)[[1]]

is_correct_mat <- do.call(cbind, lapply(seq_along(use_indices), function(x) {
  char_vec_list[[x]] == amplicon_ref_vec[[x]]
}))



# Identify insertions and deletions ---------------------------------------

insertion_vec_list <- lapply(use_indices, function(x) {
  message("Counting the number of inserted bases at position #", x, "...")
  vapply(ref_char_numbers_list, function(y) sum(y == x) - 1L, integer(1))
})

num_insertions_mat <- do.call(cbind, insertion_vec_list)
has_insertion_mat <- num_insertions_mat > 0

char_mat <- do.call(cbind, char_vec_list)
is_deleted_mat <- char_mat == "-"



# Save data ---------------------------------------------------------------

save(list = c("is_deleted_mat", "is_correct_mat", "has_insertion_mat"),
     file = file.path(rdata_dir, "04_categorize_each_base.RData")
     )


