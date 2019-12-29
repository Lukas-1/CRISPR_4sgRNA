### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
CRISPOR_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRa_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData"))






# Identify genes with incomplete GuideScan data ---------------------------

unique_IDs <- unique(merged_replaced_CRISPRa_df[, "Combined_ID"])

have_incomplete_GuideScan <- vapply(unique_IDs, function(x) {
  are_this_ID <- merged_replaced_CRISPRa_df[, "Combined_ID"] == x
  num_with_GuideScan <- tapply(merged_replaced_CRISPRa_df[are_this_ID, "GuideScan_specificity"],
                               merged_replaced_CRISPRa_df[are_this_ID, "AltTSS_ID"],
                               function(x) sum(!(is.na(x)))
                               )
  num_overall <- table(merged_replaced_CRISPRa_df[are_this_ID, "AltTSS_ID"])
  are_incomplete <- (num_with_GuideScan < 4) & (num_with_GuideScan != num_overall)
  return(any(are_incomplete))
}, logical(1))


incomplete_GuideScan_IDs <- unique_IDs[have_incomplete_GuideScan]






# Select subsets of genes for submission to CRISPOR -----------------------

TF_combined_IDs <- intersect(all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Combined_ID"], merged_replaced_CRISPRa_df[, "Combined_ID"])

TF_noGuideScan_combined_IDs <- intersect(TF_combined_IDs, incomplete_GuideScan_IDs)





# Prepare data frames that can be exported to .bed files ------------------

TF_bed_df             <- MakeBedDf(merged_replaced_CRISPRa_df, TF_combined_IDs)
TF_noGuideScan_bed_df <- MakeBedDf(merged_replaced_CRISPRa_df, TF_noGuideScan_combined_IDs)





# Write input files for CRISPOR to disk -----------------------------------

write.table(TF_bed_df,
            file = file.path(CRISPOR_files_directory, "Input_for_CRISPOR_CRISPRa_TFs.bed"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )

write.table(TF_noGuideScan_bed_df,
            file = file.path(CRISPOR_files_directory, "Input_for_CRISPOR_CRISPRa_TFs_noGuideScan.bed"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )
















