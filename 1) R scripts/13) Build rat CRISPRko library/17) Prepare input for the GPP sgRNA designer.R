### 18th November 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "10) Rat - General")
CRISPRko_RData_directory  <- file.path(RData_directory, ") CRISPRko")
GPP_input_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Rat - CRISPRko", "GPP sgRNA designer", "1) Input files")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "04) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(general_RData_directory, "06) Read in gene lists.RData"))





# Collect Entrez IDs for submission to the GPP sgRNA designer -------------

pharma_entrezs <- pharma_df[["Entrez_ID"]]




# Build data frames for submission to the GPP sgRNA designer --------------

pharma_submit_df <- BuildDfForGPP(pharma_entrezs)




# Write input files for the GPP sgRNA designer to disk --------------------

WriteGPPInputDf(pharma_submit_df,
                "pharmacology",
                file.path(GPP_input_files_directory, "Pharmacology genes")
                )






