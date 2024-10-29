## 2024-07-14


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
pacbio_seq_functions_dir <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate", "1) R functions")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")


source(file.path(pacbio_seq_functions_dir, "02) Analyzing reads.R")) # For GetMeanQuality
project_dir <- first_nanopore_dir
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))
rm(project_dir)
project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "1_demultiplexing_Nanopore_reads.R"))



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
input_dir <- file.path(project_dir, "02_input")
rdata_dir <- file.path(project_dir, "03_R_objects")

fastq_files <- list.files(file.path(input_dir, "Raw_reads"), full.names = TRUE)
stopifnot(all(grepl(".pass.fastq.gz", fastq_files, fixed = TRUE)))



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_prepare_barcodes.RData"))
load(file.path(rdata_dir, "02_count_reads__fastq_df.RData"))



# Demultiplex reads -------------------------------------------------------

DemultiplexFastqFiles(fastq_files[[4]],
                      file_numbers = 4,
                      max_num_chunks = 6,
                      start_at_chunk = 13,
                      fastq_num_reads = fastq_df[, "records"][fastq_df[, "file_name"] == basename(fastq_files[[4]])],
                      readerBlockSize = 1e7
                      )


