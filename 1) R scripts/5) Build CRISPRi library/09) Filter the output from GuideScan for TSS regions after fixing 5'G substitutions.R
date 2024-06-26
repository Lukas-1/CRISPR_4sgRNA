### 9th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR_4sgRNA"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRi_RData_directory   <- file.path(RData_directory, "4) CRISPRi")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRi", "GuideScan")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))
load(file.path(CRISPRi_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData"))





# Read in data ------------------------------------------------------------

guidescan_raw_df <- ReadGuideScanOutput("GuideScan_output_CRISPRi_all_TSSs.csv")





# Process GuideScan data --------------------------------------------------

replaced_guidescan_all_genes_df <- BuildGuideScanDf(guidescan_raw_df, combined_TSS_CRISPRi_df, replaced_merged_CRISPRi_df)





# Save data ---------------------------------------------------------------

save(list = "replaced_guidescan_all_genes_df",
     file = file.path(CRISPRi_RData_directory, "09) Filter the output from GuideScan for TSS regions after fixing 5'G substitutions.RData")
     )



