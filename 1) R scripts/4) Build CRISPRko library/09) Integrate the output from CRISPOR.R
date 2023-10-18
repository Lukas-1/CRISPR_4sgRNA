### 24th December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # For TruncateLongEntries





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRko_RData_directory, "07) Integrate the output from GuideScan.RData"))






# Identify files to read in -----------------------------------------------

output_files_list <- IdentifyCRISPOROutputFiles()





# Read in data ------------------------------------------------------------

CRISPOR_bed_df              <- ReadCRISPOROutputFiles(output_files_list[["bed"]],              is_FASTA = FALSE)
CRISPOR_FASTA_df            <- ReadCRISPOROutputFiles(output_files_list[["FASTA"]],            is_FASTA = TRUE)
CRISPOR_offtargets_bed_df   <- ReadCRISPOROutputFiles(output_files_list[["bed_offtargets"]],   is_FASTA = FALSE, show_messages = TRUE)
CRISPOR_offtargets_FASTA_df <- ReadCRISPOROutputFiles(output_files_list[["FASTA_offtargets"]], is_FASTA = TRUE)






# Add the output from CRISPOR to the data frame ---------------------------

merged_CRISPRko_df <- AddCRISPORBedData(merged_CRISPRko_df,   CRISPOR_bed_df,   CRISPOR_offtargets_bed_df)
merged_CRISPRko_df <- AddCRISPORFASTAData(merged_CRISPRko_df, CRISPOR_FASTA_df, CRISPOR_offtargets_FASTA_df)





# Truncate PAM_0MM entries, now that they are no longer needed ------------

merged_CRISPRko_df[["PAM_0MM"]] <- TruncateLongEntries(merged_CRISPRko_df[["PAM_0MM"]])






# Calculate the proportions of guides that meet specificity cutoffs -------

both_are_available <- !(is.na(merged_CRISPRko_df[["CRISPOR_CFD_specificity"]])) &
                      !(is.na(merged_CRISPRko_df[["GuideScan_specificity"]])) &
                      !(is.na(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]]))

GetProportion(merged_CRISPRko_df[["CRISPOR_CFD_specificity"]][both_are_available] >= 80)
GetProportion(merged_CRISPRko_df[["GuideScan_specificity"]][both_are_available] >= 0.2)
GetProportion(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]][both_are_available] >= 0.2)
GetProportion(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]][both_are_available] >= 0.04)





# Display sgRNAs for which only GuideScan data is available ---------------

# TF_combined_IDs <- all_TF_df[["Combined_ID"]][all_TF_df[["Is_TF"]] == "Yes"]
# CRISPRko_TF_sgRNAs_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% TF_combined_IDs, ]
# only_GuideScan_present <- !(is.na(CRISPRko_TF_sgRNAs_df[["GuideScan_specificity"]])) &
#                             is.na(CRISPRko_TF_sgRNAs_df[["CRISPOR_3MM_specificity"]])

only_GuideScan_present <- !(is.na(merged_CRISPRko_df[["GuideScan_specificity"]])) &
                            is.na(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]])

if (any(only_GuideScan_present)) {
  warning(paste0("CRISPOR specificity scores are missing for ", sum(only_GuideScan_present),
                 " guides for which Guidescan specificity scores are available!"
                 )
          )
  print(head(merged_CRISPRko_df[only_GuideScan_present, specificity_demo_columns]))
}



# Check which guides have no CRISPOR scores -------------------------------

check_columns <- c(
  "Combined_ID", "Entrez_ID", "Gene_symbol", "Source",
  "sgRNA_sequence", "PAM", "Original_PAM", "PAM_0MM",
  "Num_0MM", "Num_1MM",
  "CRISPOR_Doench_efficacy",
  "Cut_location", "Start", "Location_ID",
  "Chromosome", "Entrez_chromosome", "Locations_0MM"
)

have_no_scores <- is.na(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]])
table(have_no_scores)
head(merged_CRISPRko_df[have_no_scores, check_columns])







# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "09) Integrate the output from CRISPOR.RData")
     )






