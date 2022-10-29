## 2022-09-13


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
load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_parts_combined.RData"))




# Modify lumi_df ----------------------------------------------------------

names(lumi_df) <- sub("_sg2", "_sg1", names(lumi_df), fixed = TRUE)
names(lumi_df) <- sub("_sg3", "_sg2", names(lumi_df), fixed = TRUE)

lumi_df <- lumi_df[, c("Sample_name", "Sample_number",
                       "Plasmid_sg1", "Plasmid_sg2",
                       "Has_template_switch", "Num_MM_sg1", "Num_MM_sg2",
                       "Num_matched_sgRNAs"
                       )]
gc()




# Obtain read counts per plasmid -----------------------------------------

either0or1_including_switch_counts_mat <- AllSamplesCounts(lumi_df,
                                                           only_0MM = FALSE,
                                                           no_template_switch = FALSE
                                                           )
only0MM_including_switch_counts_mat    <- AllSamplesCounts(lumi_df,
                                                           only_0MM = TRUE,
                                                           no_template_switch = FALSE
                                                           )
either0or1_without_switch_counts_mat   <- AllSamplesCounts(lumi_df,
                                                           only_0MM = FALSE,
                                                           no_template_switch = TRUE
                                                           )
only0MM_without_switch_counts_mat      <- AllSamplesCounts(lumi_df,
                                                           only_0MM = TRUE,
                                                           no_template_switch = TRUE
                                                           )



# Obtain read counts while excluding only evident template switches -------

## This is helpful for understanding cases where only one sgRNA out of the two
## is represented in the library. Is the second sgRNA affected by a template
## switch, or is it replaced by some other sequence?

either0or1_no_evident_switch_vec <- GetCounts2sg(lumi_df,
                                                 no_template_switch = TRUE,
                                                 no_template_switch_requires_both_sgRNAs = FALSE
                                                 )


# Combine count data ------------------------------------------------------

sample_names <- colnames(either0or1_including_switch_counts_mat)
AddPrefix <- function(input_mat, add_prefix) {
  stopifnot(identical(colnames(input_mat), sample_names))
  colnames(input_mat) <- paste0(add_prefix, "_", colnames(input_mat))
  return(input_mat)
}

counts_df <- data.frame(
  sg_sequences_df[, c("Plasmid_ID", "Gene_symbol", "Entrez_ID")],
  "Sum_MaySwitch_xMM"       = rowSums(either0or1_including_switch_counts_mat),
  "Sum_MaySwitch_0MM"       = rowSums(only0MM_including_switch_counts_mat),
  "Sum_NoSwitch_xMM"        = rowSums(either0or1_without_switch_counts_mat),
  "Sum_NoSwitch_0MM"        = rowSums(only0MM_without_switch_counts_mat),
  "Sum_NoEvidentSwitch_xMM" = unname(either0or1_no_evident_switch_vec),
  AddPrefix(either0or1_including_switch_counts_mat, "MaySwitch_xMM"),
  AddPrefix(only0MM_including_switch_counts_mat,    "MaySwitch_0MM"),
  AddPrefix(either0or1_without_switch_counts_mat,   "NoSwitch_xMM"),
  AddPrefix(only0MM_without_switch_counts_mat,      "NoSwitch_0MM"),
  stringsAsFactors = FALSE
)



# Save data ---------------------------------------------------------------

save(list = "counts_df",
     file = file.path(rdata_dir, "05_obtain_read_counts.RData")
     )



