### 2022-09-13


# Define paths ------------------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_trial_rdata_dir    <- file.path(first_illumina_trial_dir, "03_R_objects")
project_dir              <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir                <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(first_trial_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_trial_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))
load(file.path(rdata_dir, "06_extend_read_counts.RData"))



# Examine the distribution of the proportion of template switches ---------

counts_switch_vec <- counts_df[, "Sum_MaySwitch_xMM"] - counts_df[, "Sum_NoEvidentSwitch_xMM"]
fraction_switch_vec <- counts_switch_vec / counts_df[, "Sum_MaySwitch_xMM"]
hist(fraction_switch_vec, breaks = 200)




# Check for plasmids that are not represented in the data -----------------

## Out of a total of 22573 4sg plasmids, fewer than 0.3% were
## not represented by any reads in the sequencing data (counts of zero
## in all samples). The exact proportions were:
## -- 37 (0.16%) when tolerating template switches and 1 base mismatch
## -- 38 (0.13%) when tolerating only template switches (zero mismatches allowed)
## -- 68 (0.23%) when including only reads without a template switch (and with two sgRNAs from the library), but tolerating 1 base mismatch
## -- 72 (0.25%) when including only reads without a template switch (and with two sgRNAs from the library) and without any mismatched bases

table(counts_df[, "Sum_MaySwitch_xMM"] == 0)
table(counts_df[, "Sum_MaySwitch_0MM"] == 0)
table(counts_df[, "Sum_NoSwitch_xMM"] == 0)
table(counts_df[, "Sum_NoSwitch_0MM"] == 0)

round(table(counts_df[, "Sum_MaySwitch_xMM"] == 0) / nrow(counts_df), digits = 4)
round(table(counts_df[, "Sum_MaySwitch_0MM"] == 0) / nrow(counts_df), digits = 4)
round(table(counts_df[, "Sum_NoSwitch_xMM"] == 0) / nrow(counts_df), digits = 4)
round(table(counts_df[, "Sum_NoSwitch_0MM"] == 0) / nrow(counts_df), digits = 4)



## When considering sg2 and sg3 individually, 41 sg2 sequences (0.18%) and
## 62 sg3 sequences (0.27%) were not found at all in the sequencing data.
## 41 plasmids were represented either only by sg2 or only by
## sg3, whereas zero reads were found for the other sgRNA.
## For 25 plasmids, only sg2 was found, but not sg3, whereas for 4 plasmids,
## only sg3 was found, but not sg2.

table(counts_df[, "Data_contains_sg2"])
table(counts_df[, "Data_contains_sg3"])
table(counts_df[, "Data_contains_sg2"] | counts_df[, "Data_contains_sg2"])
only_has_sg2 <- (counts_df[, "Data_contains_sg2"] & (!(counts_df[, "Data_contains_sg3"])))
only_has_sg3 <- (counts_df[, "Data_contains_sg3"] & (!(counts_df[, "Data_contains_sg2"])))
only_1sg_represented <- only_has_sg2 | only_has_sg3

contain_either <- counts_df[, "Data_contains_sg2"] | counts_df[, "Data_contains_sg3"]
stopifnot(contain_either == (counts_df[, "Sum_MaySwitch_xMM"] != 0))

fisher.test(table(
  c(rep(FALSE, sum(contain_either)), rep(TRUE, sum(contain_either))),
  c(only_has_sg2[contain_either], only_has_sg3[contain_either])
))



# Explore the essentiality scores of missing sgRNAs -----------------------

## Examine all genes for which essentiality data is available from DepMap:
#### When examining all genes and timepoints, 2 out of 2056 essential genes
#### were missing (0.10%), whereas 32 out of 18266 non-essential genes were missing (0.18%).
matches_vec <- match(counts_df[, "Entrez_ID"], essential_df[, "Entrez_ID"])
are_essential <- essential_df[, "CRISPR_common"][matches_vec]
fisher.test(table(are_essential, contain_either))

#### When considering just the “Tbefore” timepoint, 2 out of 2056 essential
#### genes were missing (0.10%), whereas 33 out of 18266 non-essential genes
#### were missing (0.18%); this difference did NOT reach statistical
#### significance (odds ratio 0.54, p value = 0.5).
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
## from DepMap that were used in Nunez et al. (of which 1319 essential genes
## and 841 non-essential genes were represented in the library), only 6
## essential and 0 non-essential genes were missing.
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


