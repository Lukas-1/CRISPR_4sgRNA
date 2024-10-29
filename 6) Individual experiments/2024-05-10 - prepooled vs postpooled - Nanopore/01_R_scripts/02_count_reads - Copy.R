## 2024-05-30


# Load packages and source code -------------------------------------------

library("ShortRead")


# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
input_dir <- file.path(project_dir, "02_input")
rdata_dir <- file.path(project_dir, "03_R_objects")

fastq_files <- list.files(file.path(input_dir, "Raw_reads"), full.names = TRUE)
stopifnot(all(grepl(".pass.fastq.gz", fastq_files, fixed = TRUE)))



# Count reads -------------------------------------------------------------

fastq_new <- "C:/R_root/CRISPR_4sgRNA/6) Individual experiments/2024-05-10 - prepooled vs postpooled - Nanopore/04_intermediate_files/1_output_from_R/file_1.fastq.gz"
fastq_old <- "C:/Users/user/Documents/shared_folder/file_1.fastq.gz"

new_count <- countFastq(fastq_new)
old_count <- countFastq(fastq_old)

reads_new <- readFastq(fastq_new)
reads_old <- readFastq(fastq_old)


# Save data ---------------------------------------------------------------

save(fastq_df,
     file = file.path(rdata_dir, "02_count_reads__fastq_df.RData")
     )
