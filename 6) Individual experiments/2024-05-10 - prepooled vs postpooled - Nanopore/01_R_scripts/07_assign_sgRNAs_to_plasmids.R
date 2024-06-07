## 2023-11-03


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
source(file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions", "04_obtaining_counts_for_2sg_plasmids.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2023-10-05 - prepooled vs postpooled - Nanopore")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_extract_aligned_sgRNAs_from_SAM_file.RData"))
load(file.path(rdata_dir, "04_reformat_CRISPRa_library.RData"))
for (i in 1:4) {
  load(file.path(rdata_dir, paste0("05_look_up_aligned_sgRNA", i, ".RData")))
}



# Create a combined matched_df --------------------------------------------

matched_df_chunks <- paste0("matched_sg", 1:4, "_df")
matched_df <- do.call(data.frame,
                      c(extracted_df,
                        lapply(matched_df_chunks, get),
                        stringsAsFactors = FALSE
                        )
                      )
rm(list = c("extracted_df", "matched_df_chunks"))
gc()

total_num_reads <- nrow(matched_df)



# Prepare "matched_df" for the Assign_gRNAs function ----------------------

for (i in 1:4) {
  names(matched_df)[names(matched_df) == paste0("Sequence_sg", i)] <- paste0("Aligned_read_sg", i)
}



# Create a data frame combining all relevant data -------------------------

nano_df <- Assign_gRNAs(sg_sequences_df,
                        matched_df,
                        include_columns = 2:4
                        )
nano_df <- nano_df[, c(1:4, 13:ncol(nano_df))]




# Count the per-sample number of reads with zero mapped sgRNAs ------------

GetSamplesFac <- function(input_df) {
  samples_fac <- interaction(input_df[, "Condition"], input_df[, "Replicate"], sep = "_", lex.order = TRUE)
  levels(samples_fac) <- paste0(rep(c("prepool_T0", "prepool_T12", "postpool_T0", "postpool_T12"), each = 2),
                                rep(paste0("_rep", 1:2), times = 4)
                                )
  return(samples_fac)
}

num_reads_total <- table(GetSamplesFac(matched_df))
num_reads_mapped <- table(GetSamplesFac(nano_df))

num_reads_mat <- cbind(
  "Total_num_reads"     = num_reads_total,
  "Num_mapped_reads"    = num_reads_mapped,
  "Num_0_mapped_sgRNAs" = num_reads_total - num_reads_mapped
)



# Save data ---------------------------------------------------------------

save(num_reads_mat,
     file = file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__num_reads_mat.RData")
     )
save(nano_df,
     file = file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__nano_df.RData")
     )


