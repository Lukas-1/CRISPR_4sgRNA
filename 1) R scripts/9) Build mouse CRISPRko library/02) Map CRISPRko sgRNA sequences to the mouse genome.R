### 18 November 2020 ###



# Import packages and source code -----------------------------------------

library("BSgenome.Mmusculus.UCSC.mm10")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "6) Mouse - General")
CRISPRko_RData_directory <- file.path(RData_directory, "8) Mouse - CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData"))





# Load global variables ---------------------------------------------------

chromosome_names <- paste0("chr", c(as.character(1:19), "X", "Y", "M"))
message("Loading the mouse genome into RAM...")
chromosome_sequences_list <- sapply(chromosome_names, function(x) BSgenome.Mmusculus.UCSC.mm10[[x]], simplify = FALSE)

mouse_genes_mm10_GRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
rm(human_genes_GRanges)





# Search the mouse genome for matches to sgRNA sequences  -----------------

unique_sequences <- unique(toupper(CRISPRko_df[["sgRNA_sequence"]]))

stopifnot(all(nchar(unique_sequences) == 20))

all_sequences_df <- FindSequences(unique_sequences)




# Extract the PAM ---------------------------------------------------------

all_sequences_df <- FindOverlappingGenes(all_sequences_df, gene_models_GRanges = mouse_genes_mm10_GRanges, is_mouse = TRUE)
all_sequences_df[["PAM"]] <- GetNGGPAM(all_sequences_df, use_genome = BSgenome.Mmusculus.UCSC.mm10)




# Summarize data on all matches for a given sequence ----------------------

genome_search_df <- SummarizeFoundSequencesDf(all_sequences_df, all_sequences = unique_sequences)





# Save data ---------------------------------------------------------------

for (object_name in c("genome_search_df", "all_sequences_df")) {
  save(list = object_name,
       file = file.path(CRISPRko_RData_directory,
                        paste0("02) Map CRISPRko sgRNA sequences to the mouse genome - ", object_name, ".RData")
                        )
       )
}



