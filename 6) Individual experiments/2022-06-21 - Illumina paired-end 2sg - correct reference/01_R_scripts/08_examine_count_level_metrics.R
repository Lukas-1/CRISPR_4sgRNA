### 2022-06-23



# Define paths ------------------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_trial_rdata_dir    <- file.path(first_illumina_trial_dir, "03_R_objects")
project_dir              <- file.path(experiments_directory, "2022-06-21 - Illumina paired-end 2sg - correct reference")
rdata_dir                <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(first_trial_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_trial_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))



# Examine the distribution of the proportion of template switches ---------

counts_switch_vec <- counts_df[, "Sum_MaySwitch_xMM"] - counts_df[, "Sum_NoEvidentSwitch_xMM"]
fraction_switch_vec <- counts_switch_vec / counts_df[, "Sum_MaySwitch_xMM"]
hist(fraction_switch_vec, breaks = 200)



# Check for plasmids that are not represented in the data -----------------

## Out of a total of 21554 2sg CRISPRoff plasmids, fewer than 0.3% were
## not represented by any reads in the sequencing data (counts of zero
## in all samples). The exact proportions were:
## -- 29 (0.13%) when tolerating template switches and 1 base mismatch
## -- 29 (0.13%) when tolerating only template switches (zero mismatches allowed)
## -- 49 (0.23%) when including only reads without a template switch (and with two sgRNAs from the library), but tolerating 1 base mismatch
## -- 54 (0.25%) when including only reads without a template switch (and with two sgRNAs from the library) and without any mismatched bases

table(counts_df[, "Sum_MaySwitch_xMM"] == 0)
table(counts_df[, "Sum_MaySwitch_0MM"] == 0)
table(counts_df[, "Sum_NoSwitch_xMM"] == 0)
table(counts_df[, "Sum_NoSwitch_0MM"] == 0)

round(table(counts_df[, "Sum_MaySwitch_xMM"] == 0) / nrow(counts_df), digits = 4)
round(table(counts_df[, "Sum_MaySwitch_0MM"] == 0) / nrow(counts_df), digits = 4)
round(table(counts_df[, "Sum_NoSwitch_xMM"] == 0) / nrow(counts_df), digits = 4)
round(table(counts_df[, "Sum_NoSwitch_0MM"] == 0) / nrow(counts_df), digits = 4)



## When considering sg1 and sg2 individually, 38 sg1 sequences (0.18%) and
## 34 sg2 sequences (0.16%) were not found at all in the sequencing data.
## Only 14 plasmids were represented either only by sg1 or only by
## sg2, whereas zero reads were found for the other sgRNA.
## For 5 plasmids, only sg1 was found, but not sg2, whereas for 9 plasmids,
## only sg2 was found, but not sg1.

table(counts_df[, "Data_contains_sg1"])
table(counts_df[, "Data_contains_sg2"])
table(counts_df[, "Data_contains_sg1"] | counts_df[, "Data_contains_sg1"])
only_has_sg1 <- (counts_df[, "Data_contains_sg1"] & (!(counts_df[, "Data_contains_sg2"])))
only_has_sg2 <- (counts_df[, "Data_contains_sg2"] & (!(counts_df[, "Data_contains_sg1"])))
only_1sg_represented <- only_has_sg1 | only_has_sg2

contain_either <- counts_df[, "Data_contains_sg1"] | counts_df[, "Data_contains_sg2"]
stopifnot(contain_either == (counts_df[, "Sum_MaySwitch_xMM"] != 0))

fisher.test(table(
  c(rep(FALSE, sum(contain_either)), rep(TRUE, sum(contain_either))),
  c(only_has_sg1[contain_either], only_has_sg2[contain_either])
))



# Explore the essentiality scores of missing sgRNAs -----------------------

## Examine all genes for which essentiality data is available from DepMap:
#### When examining all genes and timepoints, 6 out of 2022 essential genes
#### were missing (0.30%), whereas 21 out of 17420 non-essential genes were missing (0.12%).
matches_vec <- match(counts_df[, "Entrez_ID"], essential_df[, "Entrez_ID"])
are_essential <- essential_df[, "CRISPR_common"][matches_vec]
fisher.test(table(are_essential, contain_either))

#### When considering just the “Tbefore” timepoint, 6 out of 2022 essential
#### genes were missing (0.30%), whereas 23 out of 17420 non-essential genes
#### were missing (0.13%); this difference did NOT reach statistical
#### significance (odds ratio 2.3, p value = 0.12).
#### For this analysis, template switches and 1 base mismatch were tolerated.
are_present <- rowSums(counts_df[, c("MaySwitch_xMM_Tbefore_R1", "MaySwitch_xMM_Tbefore_R1")]) != 0
fisher.test(table(are_essential, are_present))



## Examine only the genes used in Nunez et al. (intersection of Blomen and Hart):

are_missing <- rowSums(counts_df[, c("NoSwitch_0MM_Tbefore_R1", "NoSwitch_0MM_Tbefore_R2")]) == 0
are_present <- rowSums(counts_df[, c("NoSwitch_0MM_Tbefore_R1", "NoSwitch_0MM_Tbefore_R1")]) != 0
missing_entrezs <- setdiff(counts_df[, "Entrez_ID"][are_missing], NA)
present_entrezs <- setdiff(counts_df[, "Entrez_ID"][are_present], NA)

common_entrezs <- intersect(missing_entrezs, present_entrezs)
length(common_entrezs)
missing_entrezs <- setdiff(missing_entrezs, common_entrezs)
present_entrezs <- setdiff(present_entrezs, common_entrezs)

## When examining the defined lists of essential and non-essential genes
## from DepMap that were used in Nunez et al. (of which 1310 essential genes
## and 766 non-essential genes were represented in the library), only 2
## essential and 2 non-essential genes were missing (4 genes overall).
## For the purpose of this analysis, I only consider reads without template
## switches and without base mismatches, i.e. both sgRNAs are perfect matches
## to a plasmid in the library. Furthermore, I only consider the two replicates
## from the "Tbefore" timepoint before introduction of the CRISPRoff system,
## because we expect sgRNAs from essential genes to drop out at later timepoints.

are_missing <- ifelse(counts_df[, "Entrez_ID"] %in% missing_entrezs,
                      "Missing",
                      ifelse(counts_df[, "Entrez_ID"] %in% present_entrezs, "Present", NA)
                      )
are_essential <- ifelse(counts_df[, "Entrez_ID"] %in% essentials_2020Q2_df[, "Entrez_ID"],
                        "Essential",
                        ifelse(counts_df[, "Entrez_ID"] %in% non_essentials_2020Q2_df[, "Entrez_ID"], "Non-essential", NA)
                        )
table(are_missing, are_essential)


