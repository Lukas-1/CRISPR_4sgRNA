### 22nd July 2019 ###





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))





# Define functions --------------------------------------------------------

SymbolToEntrezID <- function(my_symbol) {
  results_entrez <- symbol_to_entrez_list[[my_symbol]]
  if (length(results_entrez) >= 2) {
    print(paste0("More than one Entrez ID found for this gene symbol, and only the first was returned: ",
                 my_symbol, " (", paste(results_entrez, collapse = ", "), ")!"
                 )
          )
  }
  return(results_entrez[1])
}


SymbolsToEntrezIDs <- function(my_symbols) {
  result_entrezs <- lapply(my_symbols, SymbolToEntrezID)
  result_entrezs <- unlist(lapply(result_entrezs, function(x) if (is.null(x)) NA else x[1]))
  return(result_entrezs)
}


EntrezIDsToSymbols <- function(my_entrezs) {
  result_symbols <- lapply(my_entrezs, function(x) entrez_to_symbol_list[[x]])
  result_symbols <- unlist(lapply(result_symbols, function(x) if (is.null(x)) NA_character_ else x))
  return(result_symbols)
}

EntrezIDsToGeneNames <- function(my_entrezs) {
  result_gene_names <- lapply(my_entrezs, function(x) entrez_to_genename_list[[x]])
  result_gene_names <- unlist(lapply(result_gene_names, function(x) if (is.null(x)) NA_character_ else x))
  return(result_gene_names)
}






