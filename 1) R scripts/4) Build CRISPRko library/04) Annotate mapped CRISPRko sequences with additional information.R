### 30th October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "02) Map CRISPRko sgRNA sequences to the human genome.RData"))
load(file.path(CRISPRko_RData_directory, "03) Disambiguate gene IDs by mapping hg19 genomic coordinates to genes - CRISPRko_df.RData"))





# Search the human genome for matches to sgRNAs ---------------------------

all_sequences_df <- FindOverlappingGenes(complete_sequences_df)
all_sequences_df[["PAM"]] <- GetNGGPAM(all_sequences_df)

unique_sequences <- unique(toupper(CRISPRko_df[["sgRNA_sequence"]]))
genome_search_df <- SummarizeFoundSequencesDf(all_sequences_df, all_sequences = unique_sequences)





# Save data ---------------------------------------------------------------

save(list = "genome_search_df", file = file.path(CRISPRko_RData_directory, "04) Annotate mapped CRISPRko sequences with additional information.RData"))





