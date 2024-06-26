## 2022-04-12



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(nanopore_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "04_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "05_look_up_aligned_sgRNAs.RData"))



# Create a data frame combining all relevant data -------------------------

nano_df <- Assign_gRNAs(sg_sequences_df, matched_df)



# Obtain read counts per plasmid ------------------------------------------

counts_df <- MakeCountsDf(sg_sequences_df, nano_df)



# Save data ---------------------------------------------------------------

total_num_reads <- nrow(matched_df)
save(list = c("nano_df", "counts_df", "total_num_reads"),
     file = file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids.RData")
     )


