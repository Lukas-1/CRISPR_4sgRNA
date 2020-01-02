### 24th December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
CRISPOR_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "CRISPOR")
output_plots_directory  <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa", "Plots")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRa_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData"))





# Read in data ------------------------------------------------------------

TFs_CRISPOR_bed_df              <- ReadCRISPOROutput("Output_from_CRISPOR_CRISPRa_TFs.tsv")
TFs_CRISPOR_FASTA_df            <- ReadCRISPOROutput("Output_from_CRISPOR_FASTA_CRISPRa_TFs.tsv")
TFs_CRISPOR_offtargets_bed_df   <- ReadCRISPOROutput("Output_from_CRISPOR_CRISPRa_TFs_offs.tsv")
TFs_CRISPOR_offtargets_FASTA_df <- ReadCRISPOROutput("Output_from_CRISPOR_FASTA_CRISPRa_TFs_offs.tsv")





# Add the output from CRISPOR to the data frame ---------------------------

# The 'missing_offtargets_df' data frame is created purely for illustrative purposes
# (check the section 'Display sgRNAs that only have off-target scores, but for which detailed off-target information was unavailable' below)
missing_offtargets_df <- AddCRISPORBedData(merged_replaced_CRISPRa_df, TFs_CRISPOR_bed_df, TFs_CRISPOR_offtargets_bed_df,
                                           resolve_missing_offtargets = FALSE
                                           )
missing_offtargets_df <- AddCRISPORFASTAData(missing_offtargets_df, TFs_CRISPOR_FASTA_df, TFs_CRISPOR_offtargets_FASTA_df,
                                             resolve_missing_offtargets = FALSE
                                             )


merged_replaced_CRISPRa_df <- AddCRISPORBedData(merged_replaced_CRISPRa_df, TFs_CRISPOR_bed_df, TFs_CRISPOR_offtargets_bed_df)
merged_replaced_CRISPRa_df <- AddCRISPORFASTAData(merged_replaced_CRISPRa_df, TFs_CRISPOR_FASTA_df, TFs_CRISPOR_offtargets_FASTA_df)






# Calculate the proportions of guides that meet specificity cutoffs -------

both_are_available <- !(is.na(merged_replaced_CRISPRa_df[, "CRISPOR_CFD_specificity"])) &
                      !(is.na(merged_replaced_CRISPRa_df[, "GuideScan_specificity"]))
both_are_available <- both_are_available & !(is.na(merged_replaced_CRISPRa_df[, "CRISPOR_3MM_specificity"]))


GetProportion(merged_replaced_CRISPRa_df[both_are_available, "CRISPOR_CFD_specificity"] >= 80)
GetProportion(merged_replaced_CRISPRa_df[both_are_available, "GuideScan_specificity"] >= 0.2)
GetProportion(ConvertCFDScores(merged_replaced_CRISPRa_df[both_are_available, "CRISPOR_CFD_specificity"]) >= 0.04)

GetProportion(merged_replaced_CRISPRa_df[both_are_available, "CRISPOR_3MM_specificity"] >= 0.2)






# Display sgRNAs that only have off-target scores -------------------------
# ... but for which detailed off-target information was unavailable


lack_detailed_offtargets <- is.na(missing_offtargets_df[, "CRISPOR_4MM_specificity"]) &
                            !(is.na(missing_offtargets_df[, "CRISPOR_CFD_specificity"]))


head(missing_offtargets_df[lack_detailed_offtargets, specificity_demo_columns])

table(missing_offtargets_df[lack_detailed_offtargets, "CRISPOR_CFD_specificity"])
hist(missing_offtargets_df[lack_detailed_offtargets, "GuideScan_specificity"],
     xlim = c(0, 1), breaks = 100, col = "black",
     las = 1, mgp = c(2.7, 0.6, 0), tcl = -0.4, cex.main = 1,
     xlab = "GuideScan specificity score",
     main = "GuideScan specificity scores for sgRNAs\nthat lack detailed off-target data from CRISPOR"
     )






# Display sgRNAs for which only GuideScan data is available ---------------

TF_combined_IDs <- all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Combined_ID"]
CRISPRko_TF_sgRNAs_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Combined_ID"] %in% TF_combined_IDs, ]

only_GuideScan_present <- !(is.na(CRISPRko_TF_sgRNAs_df[, "GuideScan_specificity"])) &
                            is.na(CRISPRko_TF_sgRNAs_df[, "CRISPOR_3MM_specificity"])

stopifnot(!(any(only_GuideScan_present)))








# Plot CRISPOR vs. GuideScan specificity scores ---------------------------

DrawAllSpecificityScatterPlots(merged_replaced_CRISPRa_df, append_to_file_name = "CRISPRa")






# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData")
     )








