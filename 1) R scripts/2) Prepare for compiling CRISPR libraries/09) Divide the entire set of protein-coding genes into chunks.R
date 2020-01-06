### 6th January 2020 ###




# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))






# Define the Entrez IDs to divide into chunks -----------------------------

TF_entrez_IDs <- all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Entrez_ID"]
TF_entrez_IDs <- TF_entrez_IDs[!(is.na(TF_entrez_IDs))]
TF_entrez_IDs <- TF_entrez_IDs[order(as.integer(TF_entrez_IDs))]

not_TF_entrez_IDs <- setdiff(collected_entrez_IDs, TF_entrez_IDs)

# Just out of interest:
setdiff(TF_entrez_IDs, collected_entrez_IDs) # The gene with the Entrez ID 729288 (ZNF286B) is considered to be a pseudogene





# Decide on the number of genes per chunk ---------------------------------

num_genes_per_chunk <- 1550
num_genes <- length(not_TF_entrez_IDs)
num_chunks_unrounded <- num_genes / num_genes_per_chunk
num_chunks <- ceiling(num_chunks_unrounded)

chunk_sequence <- rep(seq_len(num_chunks), each = num_genes_per_chunk)
chunk_sequence <- chunk_sequence[seq_len(num_genes)]





# Produce the final list of Entrez ID chunks ------------------------------

entrez_chunks_list <- split(not_TF_entrez_IDs, chunk_sequence)

entrez_chunks_list <- c(list(TF_entrez_IDs), entrez_chunks_list)

names(entrez_chunks_list) <- LETTERS[seq_along(entrez_chunks_list)]
names(entrez_chunks_list)[[1]] <- paste0(names(entrez_chunks_list)[[1]], "_TF")





# Save data ---------------------------------------------------------------

save(entrez_chunks_list,
     file = file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData")
     )








