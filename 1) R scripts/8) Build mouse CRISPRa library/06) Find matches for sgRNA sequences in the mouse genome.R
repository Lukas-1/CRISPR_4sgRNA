### 8th October 2020 ###



# Import packages and source code -----------------------------------------

library("BSgenome.Mmusculus.UCSC.mm10")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from mCRISPRa-v2 - CRISPRa_df.RData"))




# Load global variables ---------------------------------------------------

chromosome_names <- paste0("chr", c(as.character(1:19), "X", "Y", "M"))
message("Loading the mouse genome into RAM...")
chromosome_sequences_list <- sapply(chromosome_names, function(x) BSgenome.Mmusculus.UCSC.mm10[[x]], simplify = FALSE)





# Search the mouse genome for matches to sgRNA sequences ------------------

unique_sequences <- unique(toupper(CRISPRa_df[["sgRNA_sequence"]]))

all_sequences_df <- FindSequences(unique_sequences)




# Extend sequence matches with additional data (e.g. nearby genes) --------

all_sequences_df[["PAM"]] <- GetNGGPAM(all_sequences_df, BSgenome.Mmusculus.UCSC.mm10)




# Summarize data on all matches for a given sequence ----------------------

genome_search_df <- SummarizeFoundSequencesDf(all_sequences_df, all_sequences = unique_sequences)





# Save data ---------------------------------------------------------------

for (object_name in c("genome_search_df", "all_sequences_df")) {
  save(list = object_name,
       file = file.path(CRISPRa_RData_directory,
                        paste0("06) Find matches for sgRNA sequences in the mouse genome - ", object_name, ".RData")
                        )
       )
}




