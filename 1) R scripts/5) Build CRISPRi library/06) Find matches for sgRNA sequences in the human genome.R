### 8th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRi_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRi-v2 - CRISPRi_df.RData"))





# Search the human genome for matches to sgRNA sequences ------------------

unique_sequences <- unique(toupper(CRISPRi_df[["sgRNA_sequence"]]))

are_20mers <- nchar(unique_sequences) == 20

sequences_not_20mers_df <- FindVariableLengthSequences(unique_sequences[!(are_20mers)])
sequences_20mers_df     <- FindSequences(unique_sequences[are_20mers])

all_sequences_df <- rbind.data.frame(sequences_20mers_df, sequences_not_20mers_df, stringsAsFactors = FALSE, make.row.names = FALSE)





# Extend sequence matches with additional data ----------------------------

all_sequences_df[["PAM"]] <- GetNGGPAM(all_sequences_df)




# Summarize data on all matches for a given sequence ----------------------

genome_search_df <- SummarizeFoundSequencesDf(all_sequences_df, all_sequences = unique_sequences)




# Save data ---------------------------------------------------------------

for (object_name in c("genome_search_df", "all_sequences_df")) {
  save(list = object_name,
       file = file.path(CRISPRi_RData_directory,
                        paste0("06) Find matches for sgRNA sequences in the human genome - ", object_name, ".RData")
                        )
       )
}




