## 2022-04-10


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_Q20_dir <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")

project_dir <- first_nanopore_dir
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))
rm(project_dir)



# Define folder paths -----------------------------------------------------

input_dir <- file.path(first_Q20_dir, "02_input_data")
rdata_dir <- file.path(first_Q20_dir, "03_R_objects")

fastq_file_name <- "20220330_1138_X3_FAS18220_5478982b.pass.fastq.gz" # list.files(reads_dir)
fastq_reads <- ShortRead::readFastq(file.path(input_dir, fastq_file_name))



# Align reads -------------------------------------------------------------

alignments_df <- ParallelAlignInChunks(fastq_reads)



# Save data ---------------------------------------------------------------

save(alignments_df, file = file.path(rdata_dir, "01_align_reads.RData"))



