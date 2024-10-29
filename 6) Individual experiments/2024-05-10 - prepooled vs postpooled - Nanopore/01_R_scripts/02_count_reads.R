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



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_prepare_barcodes.RData"))



# Count reads -------------------------------------------------------------

reads_list <- sapply(fastq_files, countFastq, simplify = FALSE)
fastq_mat <- do.call(rbind, unname(reads_list))[, 1:2]
row.names(fastq_mat) <- NULL

fastq_df <- data.frame(
  "file_name" = basename(fastq_files),
  fastq_mat,
  stringsAsFactors = FALSE
)



# Save data ---------------------------------------------------------------

save(fastq_df,
     file = file.path(rdata_dir, "02_count_reads__fastq_df.RData")
     )
