## 2022-04-12



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")
rdata_dir <- file.path(project_dir, "03_R_objects")
figures_dir <- file.path(project_dir, "04_output_data", "Figures")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "04_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "05_look_up_aligned_sgRNAs.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids.RData"))




# Display histograms of count data ----------------------------------------

ExportCountHistograms(counts_df, title_postfix = " (Q20 Nanopore data)")




# Explore summary statistics of the Nanopore Q20 sequencing results -------

## Prepare for exploring summary statistics
nano4_df <- nano_df[nano_df[, "Num_matched_sgRNAs"] == 4, ]
row.names(nano4_df) <- NULL
perfect_df <- nano4_df[nano4_df[["Num_template_switches"]] == 0, ]
row.names(perfect_df) <- NULL

## 46.7% of aligned sgRNA sequences could be assigned to an sgRNA from the library.
num_assigned <- sum(nano_df[, "Num_matched_sgRNAs"])
num_gRNAs <- nrow(matched_df) * 4
num_assigned / num_gRNAs

## Of these, 90.1% were a perfect match, and 9.9% showed a (unique) match with one mismatched base.
## In absolute terms, this corresponds to 42.1% of aligned sequences with a
## perfect match, and 4.6% with one mismatch.

num_MM_mat <- as.matrix(nano_df[, paste0("Num_MM_sg", 1:4)])
num_MM <- sum(num_MM_mat, na.rm = TRUE)
num_MM / num_assigned
sum(num_MM_mat == 0, na.rm = TRUE) / num_gRNAs
sum(num_MM_mat, na.rm = TRUE) / num_gRNAs

## 73.8% of reads have at least one gRNA that could be assigned to an sgRNA from the library.
nrow(nano_df) / nrow(matched_df)

## Of these, 21.3% had only one matching sgRNA, 26.1% had exactly two matching sgRNAs,
## 30.6% had exactly 3 matching sgRNAs, and 22.0% had 4 matching sgRNAs.
table(nano_df[, "Num_matched_sgRNAs"]) / nrow(nano_df)

## For 16.2% of reads, all 4 sgRNA sequences could be assigned to an sgRNA from the library.
nrow(nano4_df) / nrow(matched_df)

## 65.3% of reads with 4 sgRNAs showed a template switch.
## 23.0% showed two or more template switches!
## However, only 8.2% targeted 3 or more plasmids, indicating that switches to another plasmid,
## and then back to the first plasmid, were quite common. This affected 16.1% of reads.
table(nano4_df[, "Num_template_switches"]) / nrow(nano4_df)
sum(nano4_df[, "Num_template_switches"] >= 1) / nrow(nano4_df)
sum(nano4_df[, "Num_template_switches"] >= 2) / nrow(nano4_df)
sum(nano4_df[, "Num_targeted_plasmids"] >= 3) / nrow(nano4_df)
sum(nano4_df[, "Num_switch_backs"] > 0) / nrow(nano4_df)


## 34.7% of reads targeted just one plasmid, 57.1% targeted exactly two plasmids,
## 7.8% targeted 3 plasmids, and 0.4% targeted 4 different plasmids.
table(nano4_df[, "Num_targeted_plasmids"]) / nrow(nano4_df)


## The likelihood of a template switch between two neighbouring sgRNAs is 41.2%
## between sg1 and sg2, 23.7% between sg2 and sg3, and 26.6% between sg3 and
## sg4. It makes sense that the switch rate between sg1 and sg2 would be higher,
## since the stretch of DNA between these two sgRNAs is longer (and contains the
## trimethoprim resistance element encoding for dihydrofolate reductase).
table(nano4_df[, "Switch_sg1_to_sg2"]) / nrow(nano4_df)
table(nano4_df[, "Switch_sg2_to_sg3"]) / nrow(nano4_df)
table(nano4_df[, "Switch_sg3_to_sg4"]) / nrow(nano4_df)

## Indeed, for 28.1% of plasmids where sg1 and sg4 matched one plasmid,
## sg2 and/or sg3 matched another plasmid. This has implications for strategies
## where only two sgRNAs are sequenced. Even when targeting sg1 and sg4,
## the prevalence of template switching would still have to be reduced.
sg1_matches_sg4 <- (nano4_df[, "Plasmid_sg1"] == nano4_df[, "Plasmid_sg4"])
double_switch <- sg1_matches_sg4 & (nano4_df[, "Num_template_switches"] > 0)
table(double_switch) / nrow(nano4_df)
sum(nano4_df[, "Plasmid_sg1"] == nano4_df[, "Plasmid_sg4"]) / nrow(nano4_df)
sum(double_switch) / sum(sg1_matches_sg4)


## When including only ideal reads (with four recognized sgRNAs that all map
## to the same plasmid(s)), 84.7% of plasmids in the library were represented
## by at least one read.

upper_library_mat <- toupper(as.matrix(sg_sequences_df[, paste0("Sequence_sg", 1:4)]))
sg_combos_vec <- do.call(paste, c(as.list(data.frame(upper_library_mat)), sep = "_"))
combo_plasmids_list <- split(sg_sequences_df[, "Plasmid_ID"], sg_combos_vec)
combo_plasmids_vec <- vapply(combo_plasmids_list, paste0, collapse = ", ", "")
perfect_vec <- do.call(paste, c(as.list(data.frame(nano4_df[, paste0("Sequence_sg", 1:4)])), sep = "_"))
table(names(combo_plasmids_vec) %in% perfect_vec) / length(combo_plasmids_vec)

are_present <- counts_df[, "Count_perfect"] > 0
sum(are_present) / nrow(counts_df)
median(counts_df[, "Count_perfect"])

## When including only reads where a) at least two sequences match sgRNAs in
## library and b) no template switches are observed, 98.6% of plasmids are
## represented at least once, and 65.7% of plasmids are represented
## by at least 10 reads.
sum(counts_df[, "Count_2_filtered"] > 0) / nrow(counts_df)
sum(counts_df[, "Count_2_filtered"] >= 10) / nrow(counts_df)
median(counts_df[, "Count_2_filtered"])

## When including only reads where no template switches are observed, 99.4% of plasmids are
## represented at least once, and 79.8% of plasmids are represented
## by at least 10 reads.
sum(counts_df[, "Count_filtered"] > 0) / nrow(counts_df)
sum(counts_df[, "Count_filtered"] >= 10) / nrow(counts_df)
median(counts_df[, "Count_filtered"])

## When including all reads, 99.8% of plasmids are represented at least once,
## 18.9% of plasmids are represented by at least 100 reads, and
## 53.4% of plasmids are represented by at least 50 reads.
sum(counts_df[, "Count_unfiltered"] > 0) / nrow(counts_df)
sum(counts_df[, "Count_unfiltered"] >= 100) / nrow(counts_df)
sum(counts_df[, "Count_unfiltered"] >= 50) / nrow(counts_df)
median(counts_df[, "Count_unfiltered"])


table(nano_df[, "Plasmid_sg2"] == nano_df[, "Plasmid_sg3"])




## In 76.3% of reads (of all reads for which all four sgRNA sequences are found in the
## library), sg2 and sg3 matched the same plasmid (i.e. there  was no template
## switch between them). Out of these, 54.5% had a template switch in either
## sg1 or sg4, and among those with a template switch, 66.7% had one more sgRNA
## that belonged to the same plasmid as sg2 and sg3.
## Stated differently, among all reads where sg2 and sg3 matched the same
## plasmid, in 45.5% of reads, all 4 sgRNAs matched this plasmid,
## in a further 36.3% of reads, there were 3 sgRNAs matching this plasmid, and
## in 18.2%, there were  only 2 sgRNAs (sg2 and sg3) matching this plasmid.
sg2_matches_sg3 <- (nano4_df[, "Plasmid_sg2"] == nano4_df[, "Plasmid_sg3"])
sum(sg2_matches_sg3) / nrow(nano4_df)
switch_outside <- sg2_matches_sg3 & (nano4_df[, "Num_template_switches"] > 0)
table(switch_outside) / nrow(nano4_df)
table(nano4_df[, "Num_template_switches"][switch_outside])
sum(nano4_df[, "Plasmid_sg2"] == nano4_df[, "Plasmid_sg3"]) / nrow(nano4_df)
sum(switch_outside) / sum(sg2_matches_sg3)
sg1_matches_sg2 <- (nano4_df[, "Plasmid_sg1"] == nano4_df[, "Plasmid_sg4"])
sg3_matches_sg4 <- (nano4_df[, "Plasmid_sg3"] == nano4_df[, "Plasmid_sg4"])
table((sg1_matches_sg2 & sg3_matches_sg4)[switch_outside]) / sum(switch_outside)
table((sg1_matches_sg2 | sg3_matches_sg4)[switch_outside]) / sum(switch_outside)
table((sg1_matches_sg2 & sg3_matches_sg4)[sg2_matches_sg3]) / sum(sg2_matches_sg3)
table((sg1_matches_sg2 | sg3_matches_sg4)[sg2_matches_sg3]) / sum(sg2_matches_sg3)
(sum((sg1_matches_sg2 | sg3_matches_sg4)[sg2_matches_sg3]) -
 sum((sg1_matches_sg2 & sg3_matches_sg4)[sg2_matches_sg3])) / sum(sg2_matches_sg3)


## In 73.4% of reads (of all reads for which all four sgRNA sequences are found in the
## library), sg3 and sg4 matched the same plasmid (i.e. there  was no template
## switch between them). Out of these, 52.7% had a template switch in either
## sg1 or sg2, and among those with a template switch, 69.6% had one more sgRNA
## that belonged to the same plasmid as sg3 and sg4.
## Stated differently, among all reads where sg3 and sg4 matched the same
## plasmid, in 47.3% of reads, all 4 sgRNAs matched this plasmid,
## in a further 36.7% of reads, there were 3 sgRNAs matching this plasmid, and
## in 16.0%, there were  only 2 sgRNAs (sg3 and sg4) matching this plasmid.
sg3_matches_sg4 <- (nano4_df[, "Plasmid_sg3"] == nano4_df[, "Plasmid_sg4"])
table(sg3_matches_sg4) / nrow(nano4_df)
switch_outside <- sg3_matches_sg4 & (nano4_df[, "Num_template_switches"] > 0)
table(switch_outside) / nrow(nano4_df)
sum(nano4_df[, "Plasmid_sg3"] == nano4_df[, "Plasmid_sg4"]) / nrow(nano4_df)
sum(switch_outside) / sum(sg3_matches_sg4)
sg1_matches_sg3 <- (nano4_df[, "Plasmid_sg1"] == nano4_df[, "Plasmid_sg3"])
sg2_matches_sg3 <- (nano4_df[, "Plasmid_sg2"] == nano4_df[, "Plasmid_sg3"])
table((sg1_matches_sg3 & sg2_matches_sg3)[switch_outside]) / sum(switch_outside)
table((sg1_matches_sg3 | sg2_matches_sg3)[switch_outside]) / sum(switch_outside)
table((sg1_matches_sg3 & sg2_matches_sg3)[sg3_matches_sg4]) / sum(sg3_matches_sg4)
table((sg1_matches_sg3 | sg2_matches_sg3)[sg3_matches_sg4]) / sum(sg3_matches_sg4)
(sum((sg1_matches_sg3 | sg2_matches_sg3)[sg3_matches_sg4]) -
 sum((sg1_matches_sg3 & sg2_matches_sg3)[sg3_matches_sg4])) / sum(sg3_matches_sg4)






sum(nano4_df[, "Plasmid_sg1"] != nano4_df[, "Plasmid_sg2"]) / nrow(nano4_df)





