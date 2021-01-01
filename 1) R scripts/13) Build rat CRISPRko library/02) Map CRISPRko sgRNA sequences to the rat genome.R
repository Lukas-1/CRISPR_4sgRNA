### 18 November 2020 ###



# Import packages and source code -----------------------------------------

library("BSgenome.Rnorvegicus.UCSC.rn6")
library("TxDb.Rnorvegicus.UCSC.rn6.refGene")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "10) Rat - General")
CRISPRko_RData_directory <- file.path(RData_directory, "12) Rat - CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Rn.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData"))





# Load global variables ---------------------------------------------------

chromosome_names <- paste0("chr", c(as.character(1:20), "X", "Y", "M"))
message("Loading the rat genome into RAM...")
chromosome_sequences_list <- sapply(chromosome_names, function(x) BSgenome.Rnorvegicus.UCSC.rn6[[x]], simplify = FALSE)

rat_genes_rn6_GRanges <- genes(TxDb.Rnorvegicus.UCSC.rn6.refGene)
rm(human_genes_GRanges)





# Search the rat genome for matches to sgRNA sequences --------------------

unique_sequences <- unique(toupper(CRISPRko_df[["sgRNA_sequence"]]))

stopifnot(all(nchar(unique_sequences) == 20))

all_sequences_df <- FindSequences(unique_sequences)




# Extract the PAM ---------------------------------------------------------

all_sequences_df <- FindOverlappingGenes(all_sequences_df, gene_models_GRanges = rat_genes_rn6_GRanges, is_rat = TRUE)
all_sequences_df[["PAM"]] <- GetNGGPAM(all_sequences_df, use_genome = BSgenome.Rnorvegicus.UCSC.rn6)




# Summarize data on all matches for a given sequence ----------------------

genome_search_df <- SummarizeFoundSequencesDf(all_sequences_df, all_sequences = unique_sequences)





# Save data ---------------------------------------------------------------

for (object_name in c("genome_search_df", "all_sequences_df")) {
  save(list = object_name,
       file = file.path(CRISPRko_RData_directory,
                        paste0("02) Map CRISPRko sgRNA sequences to the rat genome - ", object_name, ".RData")
                        )
       )
}



