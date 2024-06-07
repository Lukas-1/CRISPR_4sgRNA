## 2023-11-03


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "05_looking_up_aligned_sgRNAs.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2023-10-05 - prepooled vs postpooled - Nanopore")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_extract_aligned_sgRNAs_from_SAM_file.RData"))
load(file.path(rdata_dir, "04_reformat_CRISPRa_library.RData"))



# Look up aligned sgRNAs --------------------------------------------------

extracted_df <- extracted_df[, "Sequence_sg4", drop = FALSE]
names(extracted_df) <- "Aligned_read_sg4"

gc()

matched_sg4_df <- CheckGuides(extracted_df,
                              4,
                              large_chunk_size = 250000,
                              small_chunk_size = 25000
                              )


# Save data ---------------------------------------------------------------

save(matched_sg4_df,
     file = file.path(rdata_dir, "05_look_up_aligned_sgRNA4.RData")
     )

