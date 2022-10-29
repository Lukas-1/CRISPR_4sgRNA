## 2022-09-12


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
source(file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions", "04_obtaining_counts_for_2sg_plasmids.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_part1.RData"))
load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_part2.RData"))




# Combine data frames -----------------------------------------------------

lumi_df_names <- c("part1a_lumi_df", "part1b_lumi_df", "part2a_lumi_df", "part2b_lumi_df")

lumi_df <- do.call(rbind.data.frame,
                   c(lapply(lumi_df_names, get),
                     stringsAsFactors = FALSE, make.row.names = FALSE
                     ))

rm(list = lumi_df_names)
gc()



# Modify lumi_df ----------------------------------------------------------

lumi_df[, "Num_template_switches"] <- as.logical(lumi_df[, "Num_template_switches"])
names(lumi_df)[names(lumi_df) == "Num_template_switches"] <- "Has_template_switch"

samples_in_order <- c("Tbefore_R1", "Tbefore_R2", "T0_R1", "T0_R2", "T12_R1", "T12_R2")
lumi_df[, "Sample_number"] <- match(lumi_df[, "Sample_name"],  samples_in_order)

new_order <- order(lumi_df[, "Sample_number"])

lumi_df <- lumi_df[new_order, ]
row.names(lumi_df) <- NULL



# Count the number of mapped reads per sample -----------------------------

num_mapped_df <- NumMappedReads(lumi_df, sg_numbers = 2:3)




# Save data ---------------------------------------------------------------

save(num_mapped_df,
     file = file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_num_mapped_reads.RData")
     )
save(lumi_df,
     file = file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_parts_combined.RData")
     )




