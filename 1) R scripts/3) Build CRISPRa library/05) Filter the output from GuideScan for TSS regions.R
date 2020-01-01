### 25th July 2019 ###





# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory   <- file.path(RData_directory, "2) CRISPRa")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "GuideScan")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData"))
load(file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))





# Read in data ------------------------------------------------------------

guidescan_candidates_raw_df <- ReadGuideScanOutput("GuideScan_output_CRISPRa_candidate_genes.csv")
guidescan_raw_df            <- ReadGuideScanOutput("GuideScan_output_CRISPRa_all_genes.csv")




# Process the candidate genes ---------------------------------------------

guidescan_candidates_df <- GuideScanOutputToDf(guidescan_candidates_raw_df)




# Process all genes -------------------------------------------------------

guidescan_all_genes_df <- BuildGuideScanDf(guidescan_raw_df, combined_TSS_CRISPRa_df, CRISPRa_df)





# Save data ---------------------------------------------------------------

save(list = c("guidescan_all_genes_df", "guidescan_candidates_df"),
     file = file.path(CRISPRa_RData_directory, "05) Filter the output from GuideScan for TSS regions.RData")
     )








