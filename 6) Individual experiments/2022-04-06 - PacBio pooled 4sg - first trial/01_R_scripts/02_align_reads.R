## 2022-04-06


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_pacbio_dir <- file.path(experiments_directory, "2022-04-06 - PacBio pooled 4sg - first trial")

project_dir <- first_nanopore_dir
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))
rm(project_dir)



# Define folder paths -----------------------------------------------------

rdata_dir <- file.path(first_pacbio_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_read_in_data.RData"))



# Align reads -------------------------------------------------------------

alignments_df <- ParallelAlignInChunks(ccs_df)



# Save data ---------------------------------------------------------------

save(alignments_df, file = file.path(rdata_dir, "02_align_reads.RData"))



