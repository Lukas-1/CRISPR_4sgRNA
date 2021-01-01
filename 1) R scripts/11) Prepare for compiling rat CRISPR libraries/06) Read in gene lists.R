### 19th April 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(RData_directory, "10) Rat - General")
gene_lists_directory    <- file.path(CRISPR_root_directory, "2) Input data", "Rat gene lists")
pharma_genes_path       <- file.path(gene_lists_directory,
                                     "Pharmacology",
                                     "gProfiler_rnorvegicus_11-11-2020_2-24-55 PM.csv"
                                     )



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Rn.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))




# Read in data ------------------------------------------------------------

pharma_genes_df <- read.csv(pharma_genes_path, stringsAsFactors = FALSE)




# Collect Entrez IDs from the membrane screen -----------------------------

pharma_symbols_vec <- paste0(substr(pharma_genes_df[["initial_alias"]], 1, 1),
                             tolower(substr(pharma_genes_df[["initial_alias"]], 2, nchar(pharma_genes_df[["initial_alias"]])))
                             )

pharma_mapped_df <- MapToEntrezs(symbols_vec = pharma_symbols_vec, is_rat = TRUE)

pharma_df <- data.frame(
  pharma_mapped_df[, c("Entrez_ID", "Gene_symbol")],
  "Ensembl_gene_ID" = pharma_genes_df[["converted_alias"]],
  stringsAsFactors = FALSE
)




# Save data ---------------------------------------------------------------

save(list = "pharma_df",
     file = file.path(general_RData_directory, "06) Read in gene lists.RData")
     )











