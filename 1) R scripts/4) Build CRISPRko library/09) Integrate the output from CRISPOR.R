### 24th December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))






# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "CRISPOR")
output_plots_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko", "Plots")






# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRko_RData_directory, "07) Integrate the output from GuideScan.RData"))





# Read in data ------------------------------------------------------------

TFs_CRISPOR_bed_df <- read.table(file.path(CRISPOR_files_directory, "Output_from_CRISPOR_CRISPRko_TFs.tsv"),
                                 stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = NULL,
                                 sep = "\t", quote = "", comment.char = ""
                                 )

TFs_CRISPOR_offtargets_bed_df <- read.table(file.path(CRISPOR_files_directory, "Output_from_CRISPOR_CRISPRko_TFs_offs.tsv"),
                                            stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = NULL,
                                            sep = "\t", quote = "", comment.char = ""
                                            )

TFs_CRISPOR_FASTA_df <- read.table(file.path(CRISPOR_files_directory, "Output_from_CRISPOR_FASTA_CRISPRko_TFs.tsv"),
                                   stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = NULL,
                                   sep = "\t", quote = "", comment.char = ""
                                   )

TFs_CRISPOR_offtargets_FASTA_df <- read.table(file.path(CRISPOR_files_directory, "Output_from_CRISPOR_FASTA_CRISPRko_TFs_offs.tsv"),
                                              stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = NULL,
                                              sep = "\t", quote = "", comment.char = ""
                                              )





# Add the output from CRISPOR to the data frame ---------------------------

merged_CRISPRko_df <- AddCRISPORBedData(merged_CRISPRko_df, TFs_CRISPOR_bed_df, TFs_CRISPOR_offtargets_bed_df)
merged_CRISPRko_df <- AddCRISPORFASTAData(merged_CRISPRko_df, TFs_CRISPOR_FASTA_df, TFs_CRISPOR_offtargets_FASTA_df)






# Calculate the proportions of guides that meet specificity cutoffs -------

both_are_available <- !(is.na(merged_CRISPRko_df[, "CRISPOR_CFD_specificity"])) &
                      !(is.na(merged_CRISPRko_df[, "GuideScan_specificity"]))
both_are_available <- both_are_available & !(is.na(merged_CRISPRko_df[, "CRISPOR_3MM_specificity"]))


GetProportion(merged_CRISPRko_df[both_are_available, "CRISPOR_CFD_specificity"] >= 80)
GetProportion(merged_CRISPRko_df[both_are_available, "GuideScan_specificity"] >= 0.2)
GetProportion(ConvertCFDScores(merged_CRISPRko_df[both_are_available, "CRISPOR_CFD_specificity"]) >= 0.04)

GetProportion(merged_CRISPRko_df[both_are_available, "CRISPOR_3MM_specificity"] >= 0.2)





# Display sgRNAs for which only GuideScan data is available ---------------

TF_combined_IDs <- all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Combined_ID"]
CRISPRko_TF_sgRNAs_df <- merged_CRISPRko_df[merged_CRISPRko_df[, "Combined_ID"] %in% TF_combined_IDs, ]

only_GuideScan_present <- !(is.na(CRISPRko_TF_sgRNAs_df[, "GuideScan_specificity"])) &
                            is.na(CRISPRko_TF_sgRNAs_df[, "CRISPOR_3MM_specificity"])

CRISPRko_TF_sgRNAs_df[only_GuideScan_present, specificity_demo_columns]






# Plot CRISPOR vs. GuideScan specificity scores ---------------------------

DrawAllSpecificityScatterPlots(merged_CRISPRko_df, append_to_file_name = "CRISPRko")






# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "09) Integrate the output from CRISPOR.RData")
     )








