### 2022-04-09


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir    <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
project_dir           <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
source(file.path(project_dir, "01_R_scripts", "R_functions", "04_obtaining_counts_for_2sg_plasmids.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run1.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk1.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk2.RData"))



# Create a combined matched_df --------------------------------------------

matched_df <- rbind.data.frame(data.frame("Run" = 1L, run1_matched_df),
                               data.frame("Run" = 2L, run2_chunk1_matched_df),
                               data.frame("Run" = 2L, run2_chunk2_matched_df),
                               stringsAsFactors = FALSE,
                               make.row.names = FALSE
                               )
matched_df <- data.frame("Read_number" = seq_len(nrow(matched_df)),
                         matched_df, stringsAsFactors = FALSE
                         )
rm(list = c("run1_matched_df", "run2_chunk1_matched_df", "run2_chunk2_matched_df"))
gc()



# Prepare "matched_df" for the Assign_gRNAs function ----------------------

names(matched_df)[names(matched_df) == "Sequence_sg1"] <- "Aligned_read_sg1"
names(matched_df)[names(matched_df) == "Sequence_sg2"] <- "Aligned_read_sg2"



# Prepare "sg_sequences_df" for the Assign_gRNAs function -----------------

use_columns <- c("Plasmid_ID", "Entrez_ID", "Gene_symbol", c("protospacer_A", "protospacer_B"))
sg_sequences_df <- CRISPRoff_df[, use_columns]
names(sg_sequences_df)[4:5] <- paste0("Sequence_sg", 1:2)
sg_sequences_df[, "Sequence_sg1"] <- substr(sg_sequences_df[, "Sequence_sg1"], 2, 20)
sg_sequences_df[, "Sequence_sg2"] <- substr(sg_sequences_df[, "Sequence_sg2"], 2, 20)



# Create a data frame combining all relevant data -------------------------

lumi_df <- Assign_gRNAs(sg_sequences_df, matched_df,
                        sg_numbers = 1:2, include_columns = 2:6
                        )
lumi_df[, "Num_template_switches"] <- as.logical(lumi_df[, "Num_template_switches"])
names(lumi_df)[names(lumi_df) == "Num_template_switches"] <- "Has_template_switch"



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



# Check whether sgRNAs from the library are found in the data -------------

all_sg1_vec <- unique(ifelse(matched_df[, "Num_MM_sg1"] == 0,
                             matched_df[, "Aligned_read_sg1"],
                             matched_df[, "Correct_sgRNA_sg1"]
                             ))
all_sg2_vec <- unique(ifelse(matched_df[, "Num_MM_sg2"] == 0,
                             matched_df[, "Aligned_read_sg2"],
                             matched_df[, "Correct_sgRNA_sg2"]
                             ))

found_sg1 <- toupper(sg_sequences_df[, "Sequence_sg1"]) %in% all_sg1_vec
found_sg2 <- toupper(sg_sequences_df[, "Sequence_sg2"]) %in% all_sg2_vec



# Combine count data ------------------------------------------------------

sample_names <- colnames(either0or1_including_switch_counts_mat)
AddPrefix <- function(input_mat, add_prefix) {
  stopifnot(identical(colnames(input_mat), sample_names))
  colnames(input_mat) <- paste0(add_prefix, "_", colnames(input_mat))
  return(input_mat)
}

counts_df <- data.frame(
  sg_sequences_df[, c("Plasmid_ID", "Gene_symbol", "Entrez_ID")],
  "Data_contains_sg1"       = found_sg1,
  "Data_contains_sg2"       = found_sg2,
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
     file = file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData")
     )

save(list = "lumi_df",
     file = file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__lumi_df.RData")
     )


