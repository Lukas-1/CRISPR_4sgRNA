### 21st November 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # For TruncateLongEntries





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "12) Rat - CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "Rat - CRISPRko", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData"))






# Identify files to read in -----------------------------------------------

output_files_list <- IdentifyCRISPOROutputFiles()





# Read in data ------------------------------------------------------------

CRISPOR_bed_df              <- ReadCRISPOROutputFiles(output_files_list[["bed"]],              is_FASTA = FALSE)
# CRISPOR_FASTA_df            <- ReadCRISPOROutputFiles(output_files_list[["FASTA"]],            is_FASTA = TRUE)
CRISPOR_offtargets_bed_df   <- ReadCRISPOROutputFiles(output_files_list[["bed_offtargets"]],   is_FASTA = FALSE, show_messages = TRUE)
# CRISPOR_offtargets_FASTA_df <- ReadCRISPOROutputFiles(output_files_list[["FASTA_offtargets"]], is_FASTA = TRUE)






# Add the output from CRISPOR to the data frame ---------------------------

merged_CRISPRko_df <- AddCRISPORBedData(extended_CRISPRko_df, CRISPOR_bed_df,   CRISPOR_offtargets_bed_df)
# merged_CRISPRko_df <- AddCRISPORFASTAData(merged_CRISPRko_df, CRISPOR_FASTA_df, CRISPOR_offtargets_FASTA_df)





# Truncate PAM_0MM entries, now that they are no longer needed ------------

merged_CRISPRko_df[["PAM_0MM"]] <- TruncateLongEntries(merged_CRISPRko_df[["PAM_0MM"]])








# Calculate the proportions of guides that meet specificity cutoffs -------

both_are_available <- !(is.na(merged_CRISPRko_df[["CRISPOR_CFD_specificity"]])) &
                      !(is.na(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]]))

GetProportion(merged_CRISPRko_df[["CRISPOR_CFD_specificity"]][both_are_available] >= 80)
GetProportion(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]][both_are_available] >= 0.2)
GetProportion(merged_CRISPRko_df[["CRISPOR_3MM_specificity"]][both_are_available] >= 0.04)







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




# Add GuideScan columns ---------------------------------------------------

for (column_name in c("GuideScan_Num_2MM", "GuideScan_Num_3MM")) {
  merged_CRISPRko_df[[column_name]] <- NA_integer_
}
for (column_name in c("GuideScan_efficiency", "GuideScan_specificity")) {
  merged_CRISPRko_df[[column_name]] <- NA_real_
}




# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "09) Integrate the output from CRISPOR.RData")
     )






