### 26 July 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRa-v2 - CRISPRa_df.RData"))





# Search the human genome for matches to sgRNA sequences ------------------

unique_sequences <- unique(toupper(CRISPRa_df[, "sgRNA_sequence"]))

are_20mers <- nchar(unique_sequences) == 20

sequences_not_20mers_df <- FindVariableLengthSequences(unique_sequences[!(are_20mers)])
sequences_20mers_df     <- FindSequences(unique_sequences[are_20mers])

complete_sequences_df <- rbind.data.frame(sequences_20mers_df, sequences_not_20mers_df, stringsAsFactors = FALSE, make.row.names = FALSE)





# Extend sequence matches with additional data (e.g. nearby genes) --------

all_sequences_df <- FindNearestGenes(complete_sequences_df)
all_sequences_df[, "PAM"] <- GetNGGPAM(all_sequences_df)




# Summarize data on all matches for a given sequence ----------------------

genome_search_df <- SummarizeFoundSequencesDf(all_sequences_df, all_sequences = unique_sequences)





# Save data ---------------------------------------------------------------

save(list = "all_sequences_df", file = file.path(CRISPRa_RData_directory, "06) Find matches for sgRNA sequences in the human genome - all_sequences_df.RData"))
save(list = "genome_search_df", file = file.path(CRISPRa_RData_directory, "06) Find matches for sgRNA sequences in the human genome - genome_search_df.RData"))




