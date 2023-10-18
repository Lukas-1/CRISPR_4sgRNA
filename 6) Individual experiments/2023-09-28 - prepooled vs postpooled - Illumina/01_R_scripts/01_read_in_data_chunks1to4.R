## 2023-09-28


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")
R_functions_dir       <- file.path(project_dir, "01_R_scripts", "R_functions")
source(file.path(R_functions_dir, "01_reading_in_data.R"))



# Define folder paths -----------------------------------------------------

input_dir <- file.path(project_dir, "02_input")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Read in FASTQ files -----------------------------------------------------

all_files <- list.files(path.expand(input_dir), full.names = TRUE)
all_reads_chunks1to4 <- sapply(all_files[1:8],
                               function(x) ShortRead::readFastq(x),
                               simplify = FALSE
                               )



# Assemble fastq data frames ----------------------------------------------

fastq_df <- MakeFastqDf(all_reads_chunks1to4)

chunk1_fastq_df <- fastq_df[1:45000000, ]
chunk2_fastq_df <- fastq_df[45000001:90000000, ]
chunk3_fastq_df <- fastq_df[90000001:135000000, ]
chunk4_fastq_df <- fastq_df[135000001:(nrow(fastq_df)), ]

row.names(chunk1_fastq_df) <- NULL
row.names(chunk2_fastq_df) <- NULL
row.names(chunk3_fastq_df) <- NULL
row.names(chunk4_fastq_df) <- NULL



# Save data ---------------------------------------------------------------

save(chunk1_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk1.RData"))
save(chunk2_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk2.RData"))
save(chunk3_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk3.RData"))
save(chunk4_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk4.RData"))




