## 2022-02-17


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

project_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "06_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "07_look_up_aligned_sgRNAs.RData"))



# Create a data frame combining all relevant data -------------------------

nano_df <- Assign_gRNAs(sg_sequences_df, matched_df)



# Obtain read counts per plasmid ------------------------------------------

counts_df <- MakeCountsDf(sg_sequences_df, nano_df)



# Save data ---------------------------------------------------------------

save(list = c("nano_df", "counts_df"),
     file = file.path(rdata_dir, "08_assign_sgRNAs_to_plasmids.RData")
     )


