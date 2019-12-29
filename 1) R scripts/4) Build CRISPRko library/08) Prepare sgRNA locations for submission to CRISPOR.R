### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData"))






# Select subsets of genes for submission to CRISPOR -----------------------

merged_CRISPRko_df <- extended_CRISPRko_df

unique_IDs <- unique(merged_CRISPRko_df[, "Combined_ID"])

TF_combined_IDs <- intersect(all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Combined_ID"], merged_CRISPRko_df[, "Combined_ID"])

# have_incomplete_GuideScan <- vapply(unique_IDs, function(x) { # Identify genes with incomplete GuideScan data
#   are_this_ID <- merged_CRISPRko_df[, "Combined_ID"] == x
#   num_with_GuideScan <- sum(!(is.na(merged_CRISPRko_df[are_this_ID, "GuideScan_specificity"])))
#   num_overall <- sum(are_this_ID)
#   is_incomplete <- (num_with_GuideScan < 4) && (num_with_GuideScan != num_overall)
#   return(is_incomplete)
# }, logical(1))
# incomplete_GuideScan_IDs <- unique_IDs[have_incomplete_GuideScan]
# TF_noGuideScan_combined_IDs <- intersect(TF_combined_IDs, incomplete_GuideScan_IDs)





# Prepare data frames that can be exported to .bed files ------------------

TF_bed_df             <- MakeBedDf(merged_CRISPRko_df, TF_combined_IDs)
# TF_noGuideScan_bed_df <- MakeBedDf(merged_CRISPRko_df, TF_noGuideScan_combined_IDs)





# Write input files for CRISPOR to disk -----------------------------------

write.table(TF_bed_df,
            file = file.path(CRISPOR_files_directory, "Input_for_CRISPOR_CRISPRko_TFs.bed"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )

# write.table(TF_noGuideScan_bed_df,
#             file = file.path(CRISPOR_files_directory, "Input_for_CRISPOR_CRISPRko_TFs_noGuideScan.bed"),
#             quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
#             )
















