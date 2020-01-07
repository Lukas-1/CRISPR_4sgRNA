### 20th November 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory  <- file.path(RData_directory, "3) CRISPRko")
GPP_input_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "GPP sgRNA designer", "1) Input files")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))





# Collect Entrez IDs for submission to the GPP sgRNA designer -------------

submit_entrezs <- FindProblematicEntrezs(merged_CRISPRko_df, sgRNAs_overview_df)





# Separate the Entrez IDs into chunks -------------------------------------

submit_entrez_chunks_list   <- lapply(entrez_chunks_list, function(x) intersect(x, submit_entrezs))
optional_entrez_chunks_list <- lapply(entrez_chunks_list, function(x) setdiff(x, submit_entrezs))





# Build data frames for submission to the GPP sgRNA designer --------------

submit_df_list <- lapply(submit_entrez_chunks_list, BuildDfForGPP)
optional_df_list <- lapply(optional_entrez_chunks_list, BuildDfForGPP)





# Write input files for the GPP sgRNA designer to disk --------------------

for (chunk_ID in names(submit_df_list)) {
  WriteGPPInputDf(submit_df_list[[chunk_ID]], chunk_ID, GPP_input_files_directory)
}

for (chunk_ID in names(optional_df_list)) {
  WriteGPPInputDf(optional_df_list[[chunk_ID]], chunk_ID, GPP_input_files_directory, input_prefix = "optional_")
}







