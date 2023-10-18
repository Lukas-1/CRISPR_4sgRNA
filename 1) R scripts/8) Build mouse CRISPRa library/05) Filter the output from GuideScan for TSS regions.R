### 8th October 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR_4sgRNA"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory   <- file.path(RData_directory, "7) Mouse - CRISPRa")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Mouse - CRISPRa", "GuideScan")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData"))
load(file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))





# Read in data ------------------------------------------------------------

guidescan_raw_df <- ReadGuideScanOutput("GuideScan_output_CRISPRa_all_TSSs.csv")





# Process all genes -------------------------------------------------------

guidescan_all_genes_df <- BuildGuideScanDf(guidescan_raw_df, combined_TSS_CRISPRa_df, CRISPRa_df)





# Save data ---------------------------------------------------------------

save(list = "guidescan_all_genes_df",
     file = file.path(CRISPRa_RData_directory, "05) Filter the output from GuideScan for TSS regions.RData")
     )







