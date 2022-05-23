## 2022-04-11


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "04_reformatting_4sg_libraries.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Read in data ------------------------------------------------------------

library_df <- read.delim(file.path(first_nanopore_dir, "02_input_data", "CRISPRa_4sg_ordered_by_well.tsv"),
                         quote = "", stringsAsFactors = FALSE
                         )



# Create a new data frame with one row per CRISPR library plasmid ---------

sg_sequences_df <- ReformatLibrary(library_df)



# Save data ---------------------------------------------------------------

save(sg_sequences_df,
     file = file.path(rdata_dir, "04_reformat_CRISPRa_library.RData")
     )




