## 2022-02-15


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

project_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "04_reformatting_4sg_libraries.R"))



# Define paths ------------------------------------------------------------

input_dir <- file.path(project_dir, "02_input_data")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Read in data ------------------------------------------------------------

library_df <- read.delim(file.path(input_dir, "CRISPRa_4sg_ordered_by_well.tsv"),
                         quote = "", stringsAsFactors = FALSE
                         )



# Create a new data frame with one row per CRISPR library plasmid ---------

sg_sequences_df <- ReformatLibrary(library_df)



# Save data ---------------------------------------------------------------

save(sg_sequences_df,
     file = file.path(rdata_dir, "06_reformat_CRISPRa_library.RData")
     )






