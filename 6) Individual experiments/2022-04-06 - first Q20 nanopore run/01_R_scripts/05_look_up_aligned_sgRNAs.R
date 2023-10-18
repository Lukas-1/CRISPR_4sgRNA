## 2022-04-11


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "05_looking_up_aligned_sgRNAs.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_extract_aligned_sgRNAs.RData"))
load(file.path(rdata_dir, "04_reformat_CRISPRa_library.RData"))



# Look up aligned sgRNAs --------------------------------------------------

extract_df_list <- lapply(1:4, function(x) {
  CheckGuides(extracted_df, x)
})

matched_df <- do.call(data.frame, c(list(extracted_df), extract_df_list))




# Save data ---------------------------------------------------------------

save(matched_df,
     file = file.path(rdata_dir, "05_look_up_aligned_sgRNAs.RData")
     )






