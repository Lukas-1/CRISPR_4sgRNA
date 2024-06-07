## 2024-05-30


# Load packages and source code -------------------------------------------

library("readxl")
library("Biostrings")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
pacbio_seq_functions_directory  <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate", "1) R functions")
source(file.path(pacbio_seq_functions_directory, "15) Examining the Hamming distances of barcodes.R"))



# Define folder paths -----------------------------------------------------

first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")

project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
input_dir <- file.path(project_dir, "02_input")
rdata_dir <- file.path(project_dir, "03_R_objects")

barcodes_path <- file.path(input_dir, "barcodes", "Barcoded primers.xlsx")



# Read in data ------------------------------------------------------------

barcodes_df <- data.frame(read_excel(barcodes_path, range = cell_cols("A:B"), col_names = FALSE))

library_df <- read.delim(file.path(first_nanopore_dir, "02_input_data", "CRISPRa_4sg_ordered_by_well.tsv"),
                         quote = "", stringsAsFactors = FALSE, na.strings = c("NA", "")
                         )


# Define barcodes ---------------------------------------------------------

names(barcodes_df) <- c("Name", "Sequence")

condition_barcode_names <- paste0("2PF_", LETTERS[1:6])[-c(2, 5)]
replicate_barcode_names <- paste0("2PR_", 1:2)

condition_barcodes <- barcodes_df[barcodes_df[, "Name"] %in% condition_barcode_names, "Sequence"]
replicate_barcodes <- barcodes_df[barcodes_df[, "Name"] %in% replicate_barcode_names, "Sequence"]



# Examine the Hamming distance of barcodes --------------------------------

include_flanking <- 6L

all_barcodes <- c(condition_barcodes, replicate_barcodes)
all_barcodes <- toupper(substr(all_barcodes, 1, 10 + include_flanking))
all_barcodes <- c(all_barcodes, as.character(reverseComplement(DNAStringSet(all_barcodes))))
hamming_distances_df <- GetDistances(all_barcodes)



# Check for instances of barcodes in the library --------------------------

matches_mat <- t(vcountPDict(DNAStringSet(all_barcodes),
                             DNAStringSet(toupper(library_df[, "sgRNA_sequence"])),
                             max.mismatch = 0
                             ))
stopifnot(rowSums(matches_mat) == 0)
matches_mat <- t(vcountPDict(DNAStringSet(all_barcodes),
                             DNAStringSet(toupper(library_df[, "sgRNA_sequence"])),
                             max.mismatch = 1
                             ))
table(rowSums(matches_mat))



# Define search strings ---------------------------------------------------

use_condition_barcodes <- substr(condition_barcodes, 1, 10 + include_flanking)
fwd_condition_barcodes <- as.character(reverseComplement(DNAStringSet(use_condition_barcodes)))

use_replicate_barcodes <- substr(replicate_barcodes, 1, 10 + include_flanking)
fwd_replicate_barcodes <- toupper(use_replicate_barcodes)

condition_names <- c(
  "Pre-pooled T0", "Pre-pooled T12",
  "Post-pooled T0", "Post-pooled T12"
)



# Save data ---------------------------------------------------------------

save(list = c("fwd_replicate_barcodes", "fwd_condition_barcodes",
              "condition_names", "include_flanking"
              ),
     file = file.path(rdata_dir, "01_prepare_barcodes.RData")
     )



