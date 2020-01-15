### 16th September 2019 ###






# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory   <- file.path(RData_directory, "2) CRISPRa")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "GuideScan")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))
load(file.path(CRISPRa_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData"))





# Read in data ------------------------------------------------------------

guidescan_raw_df <- ReadGuideScanOutput("GuideScan_output_CRISPRa_all_genes.csv")





# Process GuideScan data --------------------------------------------------

replaced_guidescan_all_genes_df <- BuildGuideScanDf(guidescan_raw_df, combined_TSS_CRISPRa_df, replaced_merged_CRISPRa_df)





# Save data ---------------------------------------------------------------

save(list = "replaced_guidescan_all_genes_df",
     file = file.path(CRISPRa_RData_directory, "09) Filter the output from GuideScan for TSS regions after fixing 5'G substitutions.RData")
     )



