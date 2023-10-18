## 2022-09-02


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
R_functions_dir       <- file.path(project_dir, "01_R_scripts", "R_functions")
source(file.path(R_functions_dir, "01_reading_in_data.R"))



# Define folder paths -----------------------------------------------------

input_dir   <- file.path(project_dir, "02_input_data")
rdata_dir   <- file.path(project_dir, "03_R_objects")
reads_dir   <- file.path(input_dir, "Raw reads")



# Read in FASTQ files -----------------------------------------------------

all_files <- list.files(path.expand(reads_dir), full.names = TRUE)
all_reads_chunks5to8 <- sapply(all_files[7:12],
                               function(x) ShortRead::readFastq(x),
                               simplify = FALSE
                               )



# Assemble fastq data frames ----------------------------------------------

fastq_df <- MakeFastqDf(all_reads_chunks5to8)

chunk5_fastq_df <- fastq_df[1:47000000, ]
chunk6_fastq_df <- fastq_df[47000001:94000000, ]
chunk7_fastq_df <- fastq_df[94000001:141000000, ]
chunk8_fastq_df <- fastq_df[141000000:(nrow(fastq_df)), ]

row.names(chunk5_fastq_df) <- NULL
row.names(chunk6_fastq_df) <- NULL
row.names(chunk7_fastq_df) <- NULL
row.names(chunk8_fastq_df) <- NULL



# Save data ---------------------------------------------------------------

save(chunk5_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk5.RData"))
save(chunk6_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk6.RData"))
save(chunk7_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk7.RData"))
save(chunk8_fastq_df, file = file.path(rdata_dir, "01_read_in_data_chunk8.RData"))




