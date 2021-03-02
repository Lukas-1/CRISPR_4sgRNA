### 24th February 2021 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory  <- file.path(RData_directory, "8) Mouse - CRISPRko")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Mouse - CRISPRko", "GuideScan")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData"))






# Read in and process the output from GuideScan ---------------------------

guidescan_sgRNAs_df <- GetCRISPRkoGuideScanOutput()





# Find the matching regions in the output from GuideScan ------------------

sgRNAs_input_vec <- sgRNAStringForGuideScan(extended_CRISPRko_df)

matches_vec <- match(CRISPRStringForMatching(extended_CRISPRko_df),
                     GuideScanStringForMatching(guidescan_sgRNAs_df)
                     )






# Add the data from GuideScan to the CRISPRko library ---------------------

merged_CRISPRko_df <- extended_CRISPRko_df
guidescan_columns <- c("GuideScan_efficiency", "GuideScan_specificity", "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM")

for (column in guidescan_columns) {
  merged_CRISPRko_df[[column]] <- guidescan_sgRNAs_df[[column]][matches_vec]
}
merged_CRISPRko_df[["GuideScan_offtarget_category"]] <- GetOffTargetCategory(merged_CRISPRko_df)






# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "07) Integrate the output from GuideScan.RData")
     )



