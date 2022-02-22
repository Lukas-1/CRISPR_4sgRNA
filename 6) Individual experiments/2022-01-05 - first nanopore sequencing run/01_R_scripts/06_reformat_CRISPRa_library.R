## 2022-02-15


# Load packages and source code -------------------------------------------

library("Biostrings")
library("stringdist")



# Define paths ------------------------------------------------------------

project_dir <- "~/CRISPR/6) Individual experiments/2022-01-05 - first nanopore sequencing run"
input_dir <- file.path(project_dir, "02_input_data")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "05_extract_aligned_sgRNAs.RData"))



# Read in data ------------------------------------------------------------

library_df <- read.delim(file.path(input_dir, "CRISPRa_4sg_ordered_by_well.tsv"),
                         quote = "", stringsAsFactors = FALSE
                         )



# Create a new data frame with one row per CRISPR library plasmid ---------

library_df[, "Plate_ID"] <- sapply(strsplit(library_df[, "Plate_string"], "_"), "[[", 2)
plasmids_vec <- paste0(library_df[, "Entrez_ID"],
                       "__", library_df[, "Plate_ID"],
                       "__", library_df[, "Well_number"]
                       )
library_df[, "Plasmid_ID"] <- plasmids_vec
plasmids_fac <- factor(plasmids_vec, levels = unique(plasmids_vec))
are_plasmid_specific <- vapply(names(library_df), function(x) {
  message("Checking column '", x, "' whether it relates to the plasmid (rather than to sgRNAs)... ")
  all(tapply(library_df[, x], plasmids_fac, function(y) length(unique(y)) == 1))
}, logical(1))

use_columns <- c("sgRNA_sequence", names(library_df)[are_plasmid_specific])
plasmids_df_list <- split(library_df[, use_columns], plasmids_fac)

unique(split(library_df[, "Rank"], plasmids_fac))

plasmids_df_list <- lapply(plasmids_df_list, function(x) {
  results_list <- as.list(x[, "sgRNA_sequence"])
  names(results_list) <- paste0("Sequence_sg", 1:4)
  results_list <- c(as.list(x[1, names(x) != "sgRNA_sequence"]), results_list)
  return(results_list)
})

sg_sequences_df <- do.call(rbind.data.frame,
                           c(plasmids_df_list,
                             stringsAsFactors = FALSE,
                             make.row.names = FALSE
                           ))



# Save data ---------------------------------------------------------------

save(sg_sequences_df,
     file = file.path(rdata_dir, "06_reformat_CRISPRa_library.RData")
     )






