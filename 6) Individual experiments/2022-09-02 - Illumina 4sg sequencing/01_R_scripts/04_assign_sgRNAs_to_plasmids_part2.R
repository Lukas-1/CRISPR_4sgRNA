## 2022-09-12


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
source(file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions", "04_obtaining_counts_for_2sg_plasmids.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))
for (i in 1:8) {
  load(file.path(rdata_dir, paste0("03_look_up_sgRNAs_chunk", i, ".RData")))
}



# Create a combined matched_df --------------------------------------------

matched_df_chunks <- paste0("chunk", 1:8, "_matched_df")
matched_df <- do.call(rbind.data.frame,
                      c(lapply(matched_df_chunks, get),
                        stringsAsFactors = FALSE, make.row.names = FALSE
                        ))
rm(list = matched_df_chunks)



# Prepare "matched_df" for the Assign_gRNAs function ----------------------

names(matched_df)[names(matched_df) == "Sequence_sg2"] <- "Aligned_read_sg2"
names(matched_df)[names(matched_df) == "Sequence_sg3"] <- "Aligned_read_sg3"



# Create a data frame combining all relevant data -------------------------

part2a_matched_df <- matched_df[200000001:300000000, ]
part2b_matched_df <- matched_df[300000001:nrow(matched_df), ]
row.names(part2a_matched_df) <- NULL
row.names(part2b_matched_df) <- NULL
rm(matched_df)
gc()

part2a_lumi_df <- Assign_gRNAs(sg_sequences_df,
                               part2a_matched_df,
                               sg_numbers = 2:3, include_columns = 2:6
                               )
part2b_lumi_df <- Assign_gRNAs(sg_sequences_df,
                               part2b_matched_df,
                               sg_numbers = 2:3, include_columns = 2:6
                               )



# Combine all data --------------------------------------------------------

load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_part1.RData"))

lumi_df <- rbind.data.frame(
  part1a_lumi_df, part1b_lumi_df, part2a_lumi_df, part2b_lumi_df,
  stringsAsFactors = FALSE, make.row.names = FALSE
)

lumi_df[, "Num_template_switches"] <- as.logical(lumi_df[, "Num_template_switches"])
names(lumi_df)[names(lumi_df) == "Num_template_switches"] <- "Has_template_switch"




# Save data ---------------------------------------------------------------

save(list = c("part2a_lumi_df", "part2b_lumi_df"),
     file = file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_part2.RData")
     )



