### 2022-09-13



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_parts_combined.RData"))
load(file.path(rdata_dir, "07_examine_read_level_metrics_1_matched_df.RData"))




# Examine the incidence of template switching -----------------------------

num_reads
## There were 390'204'947 paired-end reads in total.

## 99.6% of paired-end reads contained at least one sgRNA sequence
## from the library (for either sg1 or sg2).
## One mismatch within the sgRNA sequence was allowed, as long as there was an
## unambiguous match to an sgRNA in the library.
nrow(lumi_df) / num_reads


## 87.9% of paired-end reads contained two sgRNA sequences from the library.
## (In relative terms, 88.3% of just those reads that contained one sgRNA from
## the library also featured a second sgRNA from the library.)
sum(lumi_df[, "Num_matched_sgRNAs"] == 2) / num_reads
sum(lumi_df[, "Num_matched_sgRNAs"] == 2) / nrow(lumi_df)


## 59.1% of reads (with two sgRNA sequences that are found in the library)
## had a template switch, i.e., the two sgRNAs belong to different plasmids.
## (For the purpose of this statistic, I only considered reads for which both
## sgRNAs were found in the library.)
have_2sg <- lumi_df[, "Num_matched_sgRNAs"] == 2
sum(lumi_df[, "Has_template_switch"][have_2sg]) / sum(have_2sg)



