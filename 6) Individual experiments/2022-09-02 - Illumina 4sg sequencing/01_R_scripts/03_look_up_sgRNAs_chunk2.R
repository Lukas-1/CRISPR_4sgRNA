## 2022-09-03


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "05_looking_up_aligned_sgRNAs.R"))



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_read_in_data_chunk2.RData"))
load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))




# Prepare data frames for input to the CheckGuides function ---------------------------------------

extracted_df <- chunk2_fastq_df[, c("Sequence_sg2", "Sequence_sg3")]
names(extracted_df) <- c("Aligned_read_sg2", "Aligned_read_sg3")



# Look up aligned sgRNAs --------------------------------------------------

extract_df_list <- lapply(2:3, function(x) {
  CheckGuides(extracted_df, x, large_chunk_size = 500000, small_chunk_size = 50000)
})

chunk2_matched_df <- do.call(data.frame, c(list(chunk2_fastq_df[, !(names(chunk2_fastq_df) %in% c("Quality_sg2", "Quality_sg3"))]),
                                           extract_df_list,
                                           stringsAsFactors = FALSE
                                           )
                             )



# Save data ---------------------------------------------------------------

save(chunk2_matched_df,
     file = file.path(rdata_dir, "03_look_up_sgRNAs_chunk2.RData")
     )



