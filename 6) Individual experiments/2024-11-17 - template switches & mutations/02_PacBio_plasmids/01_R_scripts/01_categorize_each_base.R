## 2025-01-01


# Load packages and source code -------------------------------------------

root_dir    <- "~/CRISPR_4sgRNA"
exper_dir   <- file.path(root_dir, "6) Individual experiments")
project_dir <- file.path(exper_dir, "2024-11-17 - template switches & mutations")

source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R"))



# Define paths ------------------------------------------------------------

s2rI_dir     <- file.path(exper_dir, "2021-12-08 - integrate PacBio data", "3) R objects")
nanopore_dir <- file.path(exper_dir, "2022-01-05 - first nanopore sequencing run")
rdata_dir    <- file.path(project_dir, "02_PacBio_plasmids", "02_R_objects")



# Read in data ------------------------------------------------------------

amplicon_ref <- read.table(file.path(nanopore_dir, "02_input_data", "amplicon_4sg.txt"),
                           quote = "", stringsAsFactors = FALSE
                           )[, 1]



# Load data ---------------------------------------------------------------

load(file.path(s2rI_dir, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rI_dir, "01) Process and export plate barcodes.RData"))
load(file.path(s2rI_dir, "07.5) Pre-filter reads - ccs_df.RData"))
load(file.path(s2rI_dir, "07.5) Pre-filter reads - alignments_df.RData"))



# Perform checks ----------------------------------------------------------

stopifnot(identical(ccs_df[, "ZMW"], alignments_df[, "ZMW"]))



# Filter data -------------------------------------------------------------

are_hifi <- (ccs_df[, "Num_full_passes"] >= 7) &
            (ccs_df[, "Read_quality"] >= 0.9999)

HA_plates <- plates_df[, "Plate_number"][grepl("^HA_[0-9]", plates_df[, "Plate_name"])]
are_HA <- ccs_df[, "Plate_number"] %in% HA_plates
pass_filters <- are_hifi & are_HA



# Randomly select a certain number of reads per well ----------------------

max_reads <- 10

wells_vec <- ccs_df[, "Combined_ID"][pass_filters]
wells_fac <- factor(wells_vec, levels = unique(wells_vec))
zmws_split <- split(ccs_df[, "ZMW"][pass_filters], wells_fac)

set.seed(1)
random_zmws_list <- lapply(zmws_split, function(x) {
  num_reads <- length(x)
  if (num_reads > max_reads) {
    use_indices <- sort(sample(seq_len(num_reads), max_reads))
    x[use_indices]
  } else {
    x
  }
})

selected_zmws <- unlist(random_zmws_list, use.names = FALSE)
are_selected <- ccs_df[, "ZMW"] %in% selected_zmws
alignments_df <- alignments_df[are_selected, ]
row.names(alignments_df) <- NULL

ccs_df <- ccs_df[are_selected, ]
row.names(ccs_df) <- NULL



# Extract aligned bases ---------------------------------------------------

stopifnot(identical(nchar(alignments_df[, "Aligned_plasmid"]), nchar(alignments_df[, "Aligned_read"])))

ref_char_list <- strsplit(alignments_df[, "Aligned_plasmid"], "", fixed = TRUE)
read_char_list <- strsplit(alignments_df[, "Aligned_read"], "", fixed = TRUE)

ref_are_gaps_list <- lapply(ref_char_list, function(x) x == "-")
ref_char_numbers_list <- lapply(ref_are_gaps_list, function(x) cumsum(!(x)))

use_indices <- 28:2252

indices_vec_list <- lapply(use_indices, function(x) {
  message("Finding the aligned index at position #", x, "...")
  vapply(ref_char_numbers_list, function(y) which(y == x)[[1]], integer(1))
})

char_vec_list <- lapply(seq_along(use_indices), function(x) {
  message("Extracting base #", x, "...")
  mapply(function(y, z) y[[z]], read_char_list, indices_vec_list[[x]])
})



# Compare to the reference ------------------------------------------------

matches_vec <- match(ccs_df[, "Combined_ID"], sg_sequences_df[, "Combined_ID"])
stopifnot(!(anyNA(matches_vec)))
reference_vec <- sg_sequences_df[, "Whole_plasmid"][matches_vec]

original_vec_list <- lapply(seq_along(use_indices), function(x) {
  substr(reference_vec, x, x)
})

is_correct_mat <- do.call(cbind, lapply(seq_along(use_indices), function(x) {
  char_vec_list[[x]] == original_vec_list[[x]]
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

save(list = c("is_deleted_mat", "is_correct_mat", "has_insertion_mat", "selected_zmws"),
     file = file.path(rdata_dir, "01_categorize_each_base.RData")
     )


