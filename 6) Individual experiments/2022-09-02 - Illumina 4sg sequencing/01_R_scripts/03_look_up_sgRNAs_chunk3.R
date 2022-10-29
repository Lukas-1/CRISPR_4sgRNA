## 2022-09-05


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "05_looking_up_aligned_sgRNAs.R"))



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_read_in_data_chunk3.RData"))
load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))




# Prepare data frames for input to the CheckGuides function ---------------------------------------

extracted_df <- chunk3_fastq_df[, c("Sequence_sg2", "Sequence_sg3")]
names(extracted_df) <- c("Aligned_read_sg2", "Aligned_read_sg3")



# Look up aligned sgRNAs --------------------------------------------------

rm(chunk3_fastq_df)
gc(); gc();
extract_df_2 <- CheckGuides(extracted_df, 2, large_chunk_size = 500000, small_chunk_size = 50000)
# gc(); gc();
# extract_df_3 <- CheckGuides(extracted_df, 3, large_chunk_size = 500000, small_chunk_size = 50000)

chunk3_matched_df <- do.call(data.frame, c(list(chunk3_fastq_df[, !(names(chunk3_fastq_df) %in% c("Quality_sg2", "Quality_sg3"))],
                                                extract_df_2, extract_df_3,
                                                stringsAsFactors = FALSE
                                                )
                                           )
                             )



# Save data ---------------------------------------------------------------

save(extract_df_2,
     file = file.path(rdata_dir, "03_look_up_sgRNAs_chunk3_extract_df_2.RData")
     )
stop()
save(extract_df_3,
     file = file.path(rdata_dir, "03_look_up_sgRNAs_chunk3_extract_df_3.RData")
     )

save(chunk3_matched_df,
     file = file.path(rdata_dir, "03_look_up_sgRNAs_chunk3.RData")
     )




