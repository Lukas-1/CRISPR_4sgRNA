### 9th September 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData"))





# Search the human genome for matches to sgRNAs ---------------------------

replaced_unique_sequences <- unique(toupper(replaced_merged_CRISPRa_df[, "sgRNA_sequence"]))

are_20mers <- nchar(replaced_unique_sequences) == 20

sequences_not_20mers_df <- FindVariableLengthSequences(replaced_unique_sequences[!(are_20mers)])
sequences_20mers_df     <- FindSequences(replaced_unique_sequences[are_20mers])

replaced_complete_sequences_df <- rbind.data.frame(sequences_20mers_df, sequences_not_20mers_df, stringsAsFactors = FALSE, make.row.names = FALSE)




# Extend sequence matches with additional data (e.g. nearby genes) --------

replaced_all_sequences_df <- FindNearestGenes(replaced_complete_sequences_df)
replaced_all_sequences_df[, "PAM"] <- GetNGGPAM(replaced_all_sequences_df)




# Summarize data on all matches for a given sequence ----------------------

replaced_genome_search_df <- SummarizeFoundSequencesDf(replaced_all_sequences_df, all_sequences = replaced_unique_sequences)




# Save data ---------------------------------------------------------------

save(list = "replaced_all_sequences_df",
     file = file.path(CRISPRa_RData_directory, "10) Find matches for sgRNA sequences after fixing 5'G substitutions - replaced_all_sequences_df.RData")
     )
save(list = "replaced_genome_search_df",
     file = file.path(CRISPRa_RData_directory, "10) Find matches for sgRNA sequences after fixing 5'G substitutions - replaced_genome_search_df.RData")
     )







