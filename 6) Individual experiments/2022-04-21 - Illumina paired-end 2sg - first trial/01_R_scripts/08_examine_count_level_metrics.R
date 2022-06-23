## 2022-05-23



# Define paths ------------------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
rdata_dir             <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))



# Check for plasmids that are not represented in the data -----------------

## Out of a total of 21378 2sg CRISPRoff plasmids, close to 20% were
## not represented by any reads in the sequencing data (counts of zero
## in all samples). The exact proportions were:
## -- 2114 (9.9%) when tolerating template switches and 1 base mismatch
## -- 2576 (12.0%) when tolerating only template switches (zero mismatches allowed)
## -- 4516 (21.1%) when including only reads without a template switch (and with two sgRNAs from the library), but tolerating 1 base mismatch
## -- 4518 (21.1%) when including only reads without a template switch (and with two sgRNAs from the library) and without any mismatched bases

table(counts_df[, "Sum_MaySwitch_xMM"] == 0)
table(counts_df[, "Sum_MaySwitch_0MM"] == 0)
table(counts_df[, "Sum_NoSwitch_xMM"] == 0)
table(counts_df[, "Sum_NoSwitch_0MM"] == 0)

round(table(counts_df[, "Sum_MaySwitch_xMM"] == 0) / nrow(counts_df), digits = 3)
round(table(counts_df[, "Sum_MaySwitch_0MM"] == 0) / nrow(counts_df), digits = 3)
round(table(counts_df[, "Sum_NoSwitch_xMM"] == 0) / nrow(counts_df), digits = 3)
round(table(counts_df[, "Sum_NoSwitch_0MM"] == 0) / nrow(counts_df), digits = 3)



## When considering sg1 and sg2 individually, 3003 sg1 sequences (24.3%) and
## 3451 sg2 sequences (27.9%) were not found at all in the sequencing data.
## Interestingly, 2226 plasmids were represented either only by sg1 or only by
## sg2, whereas zero reads were found for the other sgRNA.
## For 1337 plasmids, only sg1 was found, but not sg2, whereas for 889
## plasmids, only sg2 was found, but not sg1.
##
## This may indicate that one of the sgRNAs has undergone a mutation or is
## difficult to sequence. Alternatively, the plasmid may be represented only by
## cells infected with a lentivirus where a template switch has occurred, and
## the template-switched plasmids always lose the same sgRNA. This may indicate
## that a low number of cells were infected, or that the template switches
## are correlated, so that many virions are produced with the same switch.
## UPDATE: The true explanation has been found: A different version of the
## library was used, where the sgRNA selection for a minority of plasmids
## has changed.

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

table(only_1sg_represented[contain_either], counts_df[, "Data_contains_sg1"][contain_either])

counts_switch_vec <- counts_df[, "Sum_MaySwitch_xMM"] - counts_df[, "Sum_NoEvidentSwitch_xMM"]
fraction_switch_vec <- counts_switch_vec / counts_df[, "Sum_MaySwitch_xMM"]

hist(fraction_switch_vec[only_has_sg1], breaks = 200)
hist(fraction_switch_vec[only_has_sg2], breaks = 200)
hist(fraction_switch_vec, breaks = 200)




# Explore the essentiality scores of missing sgRNAs -----------------------

## Examine all genes for which essentiality data is available from DepMap:

#### When examining all genes and timepoints, 376 out of 1944 essential genes
#### were missing (19.3%), whereas 1261 out of 17332 non-essential genes were missing (7.3%).
matches_vec <- match(counts_df[, "Entrez_ID"], essential_df[, "Entrez_ID"])
are_essential <- essential_df[, "CRISPR_common"][matches_vec]
fisher.test(table(are_essential, contain_either))

#### When considering just the “Tbefore” timepoint, 447 out of 1944 essential
#### genes were missing (23.0%), whereas 1428 out of 17332 non-essential genes
#### were missing (7.3%); this is highly significant (odds ratio 3.3, p value < 10-74).
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
## from DepMap that were used in Nunez et al. (of which 1257 essential genes
## and 766 non-essential genes were represented in the library), there seems
## to be an over-representation of essential genes among the missing plasmids.
## For the purpose of this analysis, I only consider reads without template
## switches and without base mismatches, i.e. both sgRNAs are perfect matches
## to a plasmid in the library. Furthermore, I only consider the two replicates
## from the "Tbefore" timepoint before introduction of the CRISPRoff system,
## because we expect sgRNAs from essential genes to drop out at later timepoints.
## Of the missing plasmids that were part of above-mentioned DepMap gene sets,
## 672 were essential and 26 were non-essential (i.e, 96.3% were essential).
## Of the plasmids that were present in the sequencing data, 588 were essential
## and 740 were non-essential (i.e., only 43.5% were essential).
## This corresponds to an odds ratio of 33.5 (p < 10^-142, Fisher's exact test).

are_missing <- ifelse(counts_df[, "Entrez_ID"] %in% missing_entrezs,
                      "Missing",
                      ifelse(counts_df[, "Entrez_ID"] %in% present_entrezs, "Present", NA)
                      )
are_essential <- ifelse(counts_df[, "Entrez_ID"] %in% essentials_2020Q2_df[, "Entrez_ID"],
                        "Essential",
                        ifelse(counts_df[, "Entrez_ID"] %in% non_essentials_2020Q2_df[, "Entrez_ID"], "Non-essential", NA)
                        )
table(are_missing, are_essential)


