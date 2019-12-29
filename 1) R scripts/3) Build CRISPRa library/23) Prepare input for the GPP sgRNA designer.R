### 20th November 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R")) # For EntrezIDsToSymbols (for the status message)
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # For FormatFixedWidthInteger
source(file.path(general_functions_directory, "18) Preparing input for the GPP sgRNA designer.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory   <- file.path(RData_directory, "2) CRISPRa")
GPP_input_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "GPP sgRNA designer", "1) Input files")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRa_RData_directory, "20) Summarize the human transcription factor sub-library.RData"))






# Subset data / define sublibraries ---------------------------------------

replaced_TF_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]






# Collect Entrez IDs for submission to the GPP sgRNA designer -------------

submit_entrezs <- FindProblematicEntrezs(replaced_TF_CRISPRa_df, TF_summary_df)





# Build a data frame for submission to the GPP sgRNA designer -------------

submit_df <- BuildDfForGPP(submit_entrezs)





# Write input files for the GPP sgRNA designer to disk --------------------

WriteGPPDf(submit_df, GPP_input_files_directory)












