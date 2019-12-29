### 30th October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory  <- file.path(RData_directory, "3) CRISPRko")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "GuideScan")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData"))





# Read in data ------------------------------------------------------------

guidescan_sgRNAs_raw_df <- read.csv(file.path(GuideScan_files_directory, "GuideScan_output_CRISPRko_individual_sgRNAs.csv"),
                                    header = FALSE, row.names = NULL, quote = "\"", stringsAsFactors = FALSE, comment.char = ""
                                    )





# Process the output from GuideScan ---------------------------------------

guidescan_sgRNAs_df <- GuideScanOutputToDf(guidescan_sgRNAs_raw_df)
tidy_guidescan_sgRNAs_df <- TidyGuideScanColumns(guidescan_sgRNAs_df)





# Find the matching regions in the output from GuideScan ------------------

sgRNAs_input_vec <- sgRNAStringForGuideScan(extended_CRISPRko_df)

matches_vec <- match(sgRNAs_input_vec, guidescan_sgRNAs_df[, "Region"])

guidescan_columns <- c("GuideScan_efficiency", "GuideScan_specificity", "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM")




# Update the CRISPR libraries with additional data from GuideScan ---------

merged_CRISPRko_df <- extended_CRISPRko_df

for (column in guidescan_columns) {
  merged_CRISPRko_df[, column] <- tidy_guidescan_sgRNAs_df[matches_vec, column]
}


new_columns_order <- c(
 "Rule_set_2_score", "PAM"
)

merged_CRISPRko_df[, "GuideScan_offtarget_category"] <- GetOffTargetCategory(merged_CRISPRko_df)





# Rank and re-order the sgRNAs --------------------------------------------

original_merged_CRISPRko_df <- merged_CRISPRko_df

merged_CRISPRko_df <- RankCRISPRDf(merged_CRISPRko_df)

RankCRISPRDf(merged_CRISPRko_df[merged_CRISPRko_df[, "Gene_symbol"] %in% "NAT1", ])
merged_CRISPRko_df[merged_CRISPRko_df[, "Gene_symbol"] %in% "NAT1", ]





# Make adjustments --------------------------------------------------------

merged_CRISPRko_df[, "Num_1MM"] <- rowSums(merged_CRISPRko_df[, c("Num_5G_MM", "Num_1MM")])
merged_CRISPRko_df <- merged_CRISPRko_df[, colnames(merged_CRISPRko_df) != "Num_5G_MM"]





# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "07) Integrate the output from GuideScan.RData")
     )














