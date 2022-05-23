## 2022-04-23


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "05_looking_up_aligned_sgRNAs.R"))



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 4sg - first trial")
rdata_dir   <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_read_in_data_run2_chunk1.RData"))
load(file.path(rdata_dir, "02_annotate_CRISPRoff_library.RData"))



# Prepare data frames for input to the CheckGuides function ---------------------------------------

extracted_df <- run2_chunk1_fastq_df[, paste0("Sequence_sg", 1:2)]
names(extracted_df) <- paste0("Aligned_read_sg", 1:2)

sg_sequences_df <- CRISPRoff_df[, c("protospacer_A", "protospacer_B")]
names(sg_sequences_df) <- paste0("Sequence_sg", 1:2)
sg_sequences_df[, "Sequence_sg1"] <- substr(sg_sequences_df[, "Sequence_sg1"], 2, 20)
sg_sequences_df[, "Sequence_sg2"] <- substr(sg_sequences_df[, "Sequence_sg2"], 2, 20)




# Look up aligned sgRNAs --------------------------------------------------

extract_df_list <- lapply(1:2, function(x) {
  CheckGuides(extracted_df, x, large_chunk_size = 500000, small_chunk_size = 50000)
})

run2_chunk1_matched_df <- do.call(data.frame, c(list(run2_chunk1_fastq_df[, !(names(run2_chunk1_fastq_df) %in% c("Quality_sg1", "Quality_sg2"))]),
                                                extract_df_list,
                                                stringsAsFactors = FALSE
                                                )
                                  )




# Save data ---------------------------------------------------------------

save(run2_chunk1_matched_df,
     file = file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk1.RData")
     )



