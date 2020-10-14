### 10th October 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # For TruncateLongEntries




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")
CRISPOR_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Mouse - CRISPRa", "CRISPOR")





# Identify files to read in -----------------------------------------------

output_files_list <- IdentifyCRISPOROutputFiles()





# Read in data ------------------------------------------------------------

CRISPOR_bed_df              <- ReadCRISPOROutputFiles(output_files_list[["bed"]], is_FASTA = FALSE)
CRISPOR_FASTA_df            <- ReadCRISPOROutputFiles(output_files_list[["FASTA"]], is_FASTA = TRUE)
CRISPOR_offtargets_bed_df   <- ReadCRISPOROutputFiles(output_files_list[["bed_offtargets"]], is_FASTA = FALSE, show_messages = TRUE)
CRISPOR_offtargets_FASTA_df <- ReadCRISPOROutputFiles(output_files_list[["FASTA_offtargets"]], is_FASTA = TRUE)




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData"))






# Add the output from CRISPOR to the data frame ---------------------------
# The 'missing_offtargets_df' data frame is created purely for illustrative purposes
# (check the section 'Display sgRNAs that only have off-target scores,
#  but for which detailed off-target information was unavailable' below)

missing_offtargets_df <- AddCRISPORBedData(merged_replaced_CRISPRa_df, CRISPOR_bed_df, CRISPOR_offtargets_bed_df,
                                           resolve_missing_offtargets = FALSE
                                           )
# missing_offtargets_df <- AddCRISPORFASTAData(missing_offtargets_df, CRISPOR_FASTA_df, CRISPOR_offtargets_FASTA_df,
#                                              resolve_missing_offtargets = FALSE
#                                              )

merged_replaced_CRISPRa_df <- AddCRISPORBedData(merged_replaced_CRISPRa_df,   CRISPOR_bed_df,   CRISPOR_offtargets_bed_df)
# merged_replaced_CRISPRa_df <- AddCRISPORFASTAData(merged_replaced_CRISPRa_df, CRISPOR_FASTA_df, CRISPOR_offtargets_FASTA_df)






# Truncate PAM_0MM entries, now that they are no longer needed ------------

merged_replaced_CRISPRa_df[["PAM_0MM"]] <- TruncateLongEntries(merged_replaced_CRISPRa_df[["PAM_0MM"]])






# Calculate the proportions of guides that meet specificity cutoffs -------

both_are_available <- !(is.na(merged_replaced_CRISPRa_df[["CRISPOR_CFD_specificity"]])) &
                      !(is.na(merged_replaced_CRISPRa_df[["GuideScan_specificity"]]))
both_are_available <- both_are_available & !(is.na(merged_replaced_CRISPRa_df[["CRISPOR_3MM_specificity"]]))

GetProportion(merged_replaced_CRISPRa_df[["CRISPOR_CFD_specificity"]][both_are_available] >= 80)
GetProportion(merged_replaced_CRISPRa_df[["GuideScan_specificity"]][both_are_available] >= 0.2)
GetProportion(merged_replaced_CRISPRa_df[["CRISPOR_3MM_specificity"]][both_are_available] >= 0.2)
GetProportion(merged_replaced_CRISPRa_df[["CRISPOR_4MM_specificity"]][both_are_available] >= 0.04)






# Display sgRNAs that only have off-target scores -------------------------
# ... but for which detailed off-target information was unavailable

lack_detailed_offtargets <- is.na(missing_offtargets_df[["CRISPOR_4MM_specificity"]]) &
                            !(is.na(missing_offtargets_df[["CRISPOR_CFD_specificity"]]))

head(missing_offtargets_df[lack_detailed_offtargets, specificity_demo_columns])

table(missing_offtargets_df[["CRISPOR_CFD_specificity"]][lack_detailed_offtargets])
hist(missing_offtargets_df[["GuideScan_specificity"]][lack_detailed_offtargets],
     xlim = c(0, 1), breaks = 100, col = "black",
     las = 1, mgp = c(2.7, 0.6, 0), tcl = -0.4, cex.main = 1,
     xlab = "GuideScan specificity score",
     main = "GuideScan specificity scores for sgRNAs\nthat lack detailed off-target data from CRISPOR"
     )





# Display sgRNAs for which only GuideScan data is available ---------------

only_GuideScan_present <- !(is.na(merged_replaced_CRISPRa_df[["GuideScan_specificity"]])) &
                            is.na(merged_replaced_CRISPRa_df[["CRISPOR_3MM_specificity"]])

if (any(only_GuideScan_present)) {
  warning(paste0("CRISPOR specificity scores are missing for ", sum(only_GuideScan_present),
                 " guides for which Guidescan specificity scores are available!"
                 )
          )
  print(merged_replaced_CRISPRa_df[only_GuideScan_present, specificity_demo_columns])
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

have_no_scores <- is.na(merged_replaced_CRISPRa_df[["CRISPOR_3MM_specificity"]])
table(have_no_scores)
head(merged_replaced_CRISPRa_df[have_no_scores, check_columns])






# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData")
     )



