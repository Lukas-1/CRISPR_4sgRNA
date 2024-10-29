## 2023-06-09


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "05_looking_up_aligned_sgRNAs.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "04_extract_aligned_sgRNAs_from_SAM_file__1st_half.RData"))
load(file.path(rdata_dir, "04_extract_aligned_sgRNAs_from_SAM_file__2nd_half.RData"))
load(file.path(rdata_dir, "05_reformat_CRISPRa_library.RData"))



# Look up aligned sgRNAs --------------------------------------------------

first_half_extracted_df <- first_half_extracted_df[, "Sequence_sg3", drop = FALSE]
names(first_half_extracted_df) <- "Aligned_read_sg3"

second_half_extracted_df <- second_half_extracted_df[, "Sequence_sg3", drop = FALSE]
names(second_half_extracted_df) <- "Aligned_read_sg3"

extracted_df <- rbind.data.frame(first_half_extracted_df,
                                 second_half_extracted_df,
                                 make.row.names = FALSE
                                 )
rm(list = c("first_half_extracted_df", "second_half_extracted_df"))

gc()

matched_sg3_df <- CheckGuides(extracted_df,
                              3,
                              large_chunk_size = 250000,
                              small_chunk_size = 25000,
                              num_cores = 12
                              )


# Save data ---------------------------------------------------------------

save(matched_sg3_df,
     file = file.path(rdata_dir, "06_look_up_aligned_sgRNA_3.RData")
     )

