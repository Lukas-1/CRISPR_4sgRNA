## 2022-01-24


# Load packages and source code -------------------------------------------

project_dir <- file.path("~", "NP1_4sg")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_align_reads__file_1.RData"))
load(file.path(rdata_dir, "02_align_reads__all_other_files.RData"))



# Combine data frames -----------------------------------------------------

alignments_df <- rbind.data.frame(alignments_file_1_df,
                                  alignments_other_df,
                                  stringsAsFactors = FALSE
                                  )



# Save data ---------------------------------------------------------------

save(alignments_df, file = file.path(rdata_dir, "03_combine_aligned_reads.RData"))



