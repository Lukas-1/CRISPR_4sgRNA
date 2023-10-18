## 2022-04-11


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

project_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))
rm(project_dir)



# Define paths ------------------------------------------------------------

first_Q20_dir <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")
rdata_dir <- file.path(first_Q20_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_align_reads.RData"))



# Extract aligned sequences -----------------------------------------------

extracted_df <- ExtractAlignedSgRNAs(alignments_df)



# Save data ---------------------------------------------------------------

save(extracted_df, file = file.path(rdata_dir, "03_extract_aligned_sgRNAs.RData"))




