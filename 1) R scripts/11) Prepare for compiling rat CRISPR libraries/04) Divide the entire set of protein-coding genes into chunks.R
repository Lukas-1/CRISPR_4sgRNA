### 18th November 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "24) Assigning genes to sublibraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "10) Rat - General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))





# Decide on the number of genes per chunk ---------------------------------

num_genes_per_chunk <- 1550
num_genes <- length(collected_entrez_IDs)
num_chunks_unrounded <- num_genes / num_genes_per_chunk
num_chunks <- ceiling(num_chunks_unrounded)

chunk_sequence <- rep(seq_len(num_chunks), each = num_genes_per_chunk)
chunk_sequence <- chunk_sequence[seq_len(num_genes)]





# Produce the final list of Entrez ID chunks ------------------------------

entrez_chunks_list <- split(collected_entrez_IDs, chunk_sequence)

names(entrez_chunks_list) <- LETTERS[seq_along(entrez_chunks_list)]





# Save data ---------------------------------------------------------------

save(entrez_chunks_list,
     file = file.path(general_RData_directory, "04) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData")
     )








