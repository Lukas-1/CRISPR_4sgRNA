### 29 October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory    <- file.path(RData_directory, "3) CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData"))





# Search the human genome for matches to sgRNAs ---------------------------

unique_sequences <- unique(toupper(CRISPRko_df[, "sgRNA_sequence"]))

CRISPRko_df[is.na(CRISPRko_df[, "sgRNA_sequence"]), ]

are_20mers <- nchar(unique_sequences) == 20

if (!(all(are_20mers))) {
  stop("Use the FindVariableLengthSequences function for sgRNAs that are not 20-mers!")
}

complete_sequences_df <- FindSequences(unique_sequences)






# Save data ---------------------------------------------------------------

save(list = "complete_sequences_df",
     file = file.path(CRISPRko_RData_directory, "02) Map CRISPRko sgRNA sequences to the human genome.RData")
     )





