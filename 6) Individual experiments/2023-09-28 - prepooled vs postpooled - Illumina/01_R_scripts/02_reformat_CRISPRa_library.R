## 2022-09-03


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "04_reformatting_4sg_libraries.R"))
source(file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions", "03_disambiguating_and_annotating_guides.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Read in data ------------------------------------------------------------

library_df <- read.delim(file.path(first_nanopore_dir, "02_input_data", "CRISPRa_4sg_ordered_by_well.tsv"),
                         quote = "", stringsAsFactors = FALSE, na.strings = c("NA", "")
                         )



# Create a new data frame with one row per CRISPR library plasmid ---------

sg_sequences_df <- ReformatLibrary(library_df)



# Export input for GuideScan2 ---------------------------------------------

use_columns <- c("Plasmid_ID", "Rank", "sgRNA_sequence", "Chromosome", "Strand", "Cut_location")
input_df <- AddPlasmidIDs(library_df)[, use_columns]
names(input_df)[names(input_df) == "sgRNA_sequence"] <- "Sequence"
input_df[, "Sequence"] <- toupper(input_df[, "Sequence"])
input_df[, "Start"] <- ifelse(input_df[, "Strand"] == "+",
                              input_df[, "Cut_location"] - 17L,
                              input_df[, "Cut_location"] - 3L
                              )

input_df[, "Combined_ID"] <- paste0(input_df[, "Plasmid_ID"], "__sg", input_df[, "Rank"])

guidescan2_df <- FormatForGuideScan2(input_df)
guidescan2_df <- guidescan2_df[!(is.na(guidescan2_df[, "sense"])), ]

write.csv(guidescan2_df,
          file.path(project_dir, "05_intermediate_files", "Tgonfio_for_GuideScan2.csv"),
          quote = FALSE, row.names = FALSE
          )



# Import output from GuideScan2 -------------------------------------------

GuideScan2_output_df <- read.csv(file.path(project_dir, "05_intermediate_files", "Tgonfio_output.csv"),
                                 stringsAsFactors = FALSE, quote = ""
                                 )

sg_sequences_df <- AddGuideScan2Output(sg_sequences_df,
                                       GuideScan2_output_df,
                                       split_string = "__sg",
                                       ID_column = "Plasmid_ID"
                                       )


# Save data ---------------------------------------------------------------

save(sg_sequences_df,
     file = file.path(rdata_dir, "02_reformat_CRISPRa_library.RData")
     )




