### 8th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory        <- "~/CRISPR"
CRISPR_input_directory       <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory              <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRi_RData_directory      <- file.path(RData_directory, "4) CRISPRi")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRi-v2 - CRISPRi_df.RData"))
load(file.path(CRISPRi_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))
load(file.path(CRISPRi_RData_directory, "05) Filter the output from GuideScan for TSS regions.RData"))
load(file.path(CRISPRi_RData_directory, "06) Find matches for sgRNA sequences in the human genome - genome_search_df.RData"))





# Include the results of a genome search for the sgRNA sequences ----------

extended_CRISPRi_df <- ExtendWithGenomeSearch(CRISPRi_df, genome_search_df, allow_5pG = TRUE)

table(extended_CRISPRi_df[["Num_5G_MM"]][extended_CRISPRi_df[["Num_5G_MM"]] > 0])




# Merge with GuideScan data -----------------------------------------------

full_merged_CRISPRi_df <- MergeTSSandGuideScan(extended_CRISPRi_df, guidescan_all_genes_df, combined_TSS_CRISPRi_df)

merged_CRISPRi_df <- AdjustPositionColumns(full_merged_CRISPRi_df, guidescan_all_genes_df, combined_TSS_CRISPRi_df,
                                           allow_5pG_MM = TRUE, minimal_version = TRUE
                                           )





# Check for discrepancies with GuideScan locations ------------------------

mapped_by_GuideScan <- is.na(full_merged_CRISPRi_df[["Hits_start"]]) & !(is.na(full_merged_CRISPRi_df[["GuideScan_start"]]))
were_confirmed <- vapply(which(mapped_by_GuideScan), function(x) {
  locations_vec <- strsplit(full_merged_CRISPRi_df[["Locations_0MM"]][[x]], "; ", fixed = TRUE)[[1]]
  full_merged_CRISPRi_df[["GuideScan_start"]][[x]] %in% LocationStringToDf(locations_vec)[["Start"]]
}, logical(1))
head(full_merged_CRISPRi_df[mapped_by_GuideScan, c("GuideScan_chromosome", "GuideScan_strand", "GuideScan_start", "GuideScan_end", "Locations_0MM")])






# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRi_df",
     file = file.path(CRISPRi_RData_directory, "07) Assign genomic locations to sgRNA sequences.RData")
     )





