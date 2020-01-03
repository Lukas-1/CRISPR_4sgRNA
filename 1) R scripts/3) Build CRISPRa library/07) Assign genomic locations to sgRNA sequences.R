### 25th July 2019 ###




# Legacy mode -------------------------------------------------------------

legacy_mode <- TRUE




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory        <- "~/CRISPR"
CRISPR_input_directory       <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory              <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory      <- file.path(RData_directory, "2) CRISPRa")
file_output_directory        <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRa-v2 - CRISPRa_df.RData"))
load(file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))
load(file.path(CRISPRa_RData_directory, "05) Filter the output from GuideScan for TSS regions.RData"))
load(file.path(CRISPRa_RData_directory, "06) Find matches for sgRNA sequences in the human genome - genome_search_df.RData"))





# Include the results of a genome search for the sgRNA sequences ----------

extended_CRISPRa_df <- ExtendWithGenomeSearch(CRISPRa_df, genome_search_df, allow_5pG = TRUE)

table(extended_CRISPRa_df[extended_CRISPRa_df[, "Num_5G_MM"] > 0, "Num_5G_MM"])




# Merge with GuideScan data -----------------------------------------------

full_merged_CRISPRa_df <- MergeTSSandGuideScan(extended_CRISPRa_df, guidescan_all_genes_df)

mapped_by_GuideScan <- is.na(full_merged_CRISPRa_df[, "Hits_start"]) & !(is.na(full_merged_CRISPRa_df[, "GuideScan_start"]))
were_confirmed <- vapply(which(mapped_by_GuideScan), function(x) {
  locations_vec <- strsplit(full_merged_CRISPRa_df[x, "Locations_0MM"], "; ", fixed = TRUE)[[1]]
  full_merged_CRISPRa_df[x, "GuideScan_start"] %in% LocationStringToDf(locations_vec)[, "Start"]
}, logical(1))
head(full_merged_CRISPRa_df[mapped_by_GuideScan, c("GuideScan_chromosome", "GuideScan_strand", "GuideScan_start", "GuideScan_end", "Locations_0MM")])

merged_CRISPRa_df <- AdjustPositionColumns(full_merged_CRISPRa_df, guidescan_all_genes_df, allow_5pG_MM = TRUE)

if (legacy_mode) {
  source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
  merged_CRISPRa_df <- RankCRISPRDf(merged_CRISPRa_df)
}








# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "07) Assign genomic locations to sgRNA sequences.RData")
     )





















