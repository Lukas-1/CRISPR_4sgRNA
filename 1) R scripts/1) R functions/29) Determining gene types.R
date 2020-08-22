### 18th August 2020




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))






# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "19) Compile the information on gene type.RData"))





# Define vectors ----------------------------------------------------------

gene_type_order <- c(
  "protein-coding", "tRNA", "rRNA", "miRNA", "scaRNA", "snoRNA", "snRNA",
  "vtRNA", "ribozyme", "lncRNA",  "ncRNA",
  "pseudogene", "other", "unknown"
)




# Define functions --------------------------------------------------------

CollapseNonNA <- function(char_vec, order_vec = NULL) {
  if (all(is.na(char_vec))) {
    return(NA_character_)
  } else {
    unique_vec <- unique(char_vec[!(is.na(char_vec))])
    if (!(is.null(order_vec))) {
      unique_vec <- unique_vec[order(match(unique_vec, order_vec))]
    }
    return(paste0(unique_vec, collapse = ", "))
  }
}



CollapsedTypes <- function(gene_ID_vec,
                           type_reference_df,
                           type_column_name,
                           gene_ID_column_name
                           ) {
    ID_splits <- strsplit(gene_ID_vec, ", ", fixed = TRUE)
    ID_df <- ExpandList(ID_splits)
    ID_matches <- match(ID_df[["Value"]],
                        as.character(type_reference_df[[gene_ID_column_name]])
                        )
    ID_df[["Gene_type"]] <- type_reference_df[[type_column_name]][ID_matches]
    results_vec <- tapply(ID_df[["Gene_type"]],
                          ID_df[["List_index"]],
                          function(x) CollapseNonNA(x, gene_type_order)
                          )
    return(results_vec)
}



GetDfGeneTypes <- function(genes_df) {
  stopifnot(all(c("Entrez_ID", "Ensembl_gene_ID") %in% colnames(genes_df)))
  GetGeneTypes(genes_df[["Entrez_ID"]],
               genes_df[["Ensembl_gene_ID"]]
               )
}




GetGeneTypes <- function(entrezs_vec = NULL, ENSG_vec = NULL) {

  num_entries <- length(entrezs_vec)
  stopifnot(num_entries == length(ENSG_vec))

  if (is.null(entrezs_vec) && is.null(ENSG_vec)) {
    stop("Either the 'entrezs_vec' or 'ENSG_vec' argument must be supplied!")
  }

  if (!(is.null(entrezs_vec))) {
    entrez_gene_types <- CollapsedTypes(entrezs_vec,
                                        entrez_to_gene_type_df,
                                        "Consensus_type",
                                        "Entrez_ID"
                                        )
  } else {
    entrez_gene_types <- rep(NA, num_entries)
  }
  if (!(is.null(ENSG_vec))) {
    ENSG_BioMart_gene_types <- CollapsedTypes(ENSG_vec,
                                              BioMart_ENSG_gene_type_df,
                                              "Recoded_type",
                                              "Ensembl_gene_ID"
                                              )
    ENSG_GENCODE_gene_types <- CollapsedTypes(ENSG_vec,
                                              gencode_ENSG_gene_type_df,
                                              "Recoded_type",
                                              "Ensembl_gene_ID"
                                              )
    ENSG_gene_types <- ifelse(is.na(ENSG_BioMart_gene_types),
                              ENSG_GENCODE_gene_types,
                              ENSG_BioMart_gene_types
                              )
  } else {
    ENSG_gene_types <- rep(NA, num_entries)
  }
  results_vec <- ifelse(is.na(entrez_gene_types),
                        ENSG_gene_types,
                        entrez_gene_types
                        )
  return(results_vec)
}









