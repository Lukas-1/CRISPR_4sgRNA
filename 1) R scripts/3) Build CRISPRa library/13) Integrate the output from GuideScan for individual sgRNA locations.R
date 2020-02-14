### 23rd October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory   <- file.path(RData_directory, "2) CRISPRa")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "GuideScan")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "11) Refine the genomic locations of sgRNA sequences after fixing 5'G substitutions.RData"))





# Read in data ------------------------------------------------------------

guidescan_sgRNAs_raw_df <- ReadGuideScanOutput("GuideScan_output_CRISPRa_individual_locations.csv")





# Process the output from GuideScan ---------------------------------------

guidescan_sgRNAs_df <- GuideScanOutputToDf(guidescan_sgRNAs_raw_df)
tidy_guidescan_sgRNAs_df <- TidyGuideScanColumns(guidescan_sgRNAs_df)





# Find the matching regions in the output from GuideScan ------------------

sgRNAs_input_vec <- sgRNAStringForGuideScan(merged_replaced_CRISPRa_df)

matches_vec <- match(sgRNAs_input_vec, guidescan_sgRNAs_df[["Region"]])
are_to_be_replaced <- (merged_replaced_CRISPRa_df[["TSS_searched_by_GuideScan"]] %in% c("No", "Not this gene")) & !(is.na(matches_vec))






# Update the CRISPRa library with additional data from GuideScan ----------

guidescan_columns <- c("GuideScan_efficiency", "GuideScan_specificity", "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM")

for (column in guidescan_columns) {
  merged_replaced_CRISPRa_df[[column]][are_to_be_replaced] <- tidy_guidescan_sgRNAs_df[[column]][matches_vec[are_to_be_replaced]]
}

merged_replaced_CRISPRa_df[["GuideScan_offtarget_category"]] <- GetOffTargetCategory(merged_replaced_CRISPRa_df)





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "13) Integrate the output from GuideScan for individual sgRNA locations.RData")
     )














