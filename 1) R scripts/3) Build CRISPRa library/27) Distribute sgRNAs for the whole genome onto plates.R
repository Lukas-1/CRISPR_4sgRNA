### 11th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # for FormatFixedWidthInteger
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))





# Exclude most genes that are not protein-coding --------------------------

all_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% unlist(sublibraries_all_entrezs_list, use.names = FALSE), ]




# Exclude incomplete transcripts, or those with shared subsequences -------

are_top4_mat <- CRISPRaAreTop4Mat(all_CRISPRa_df)

are_valid_top4 <- are_top4_mat[, "Are_top4"] & are_top4_mat[, "Have_valid_guides"]
are_invalid_top4 <- are_top4_mat[, "Are_top4"] & !(are_top4_mat[, "Have_valid_guides"])










