### 2024-08-11


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_pacbio_dir      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate", "1) R functions")
project_dir           <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
source(file.path(first_pacbio_dir, "02) Analyzing reads.R")) # For GetMeanQuality
source(file.path(first_pacbio_dir, "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "2_extracting_aligned_sgRNAs_from_SAM_file.R"))



# Define folder paths -----------------------------------------------------

rdata_dir   <- file.path(project_dir, "03_R_objects")
import_dir  <- file.path(project_dir, "04_intermediate_files", "2_input_to_R")

bam_paths <- list.files(import_dir, pattern = "\\.bam$", full.names = TRUE)



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_demuxed_meta_df.RData"))



# Read in data ------------------------------------------------------------

second_half_bam_data <- ReadBamData(bam_paths[4:6])



# Explore metadata --------------------------------------------------------

table(strand(second_half_bam_data))
table(mcols(second_half_bam_data)[, "flag"])
seqlevels(second_half_bam_data)



# Extract sgRNAs ----------------------------------------------------------

second_half_extracted_df <- BamExtractGuides(second_half_bam_data, demuxed_meta_df)



# Save data ---------------------------------------------------------------

save(second_half_extracted_df,
     file = file.path(rdata_dir, "04_extract_aligned_sgRNAs_from_SAM_file__2nd_half.RData")
     )
