### 13th February 2020 ###




# Import packages and source code -----------------------------------------

library("biomaRt")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))




# Define functions --------------------------------------------------------

MapEnsemblIDs <- function(input_df, warn = TRUE) {

  stopifnot(all(c("Gene_symbol", "Ensembl_gene_ID") %in% colnames(input_df)))

  ## Map TF gene symbols to Entrez IDs
  symbol_mappings_df <- MapToEntrezs(symbols_vec = input_df[["Gene_symbol"]])
  symbol_mappings_df <- symbol_mappings_df[, names(symbol_mappings_df) != "Original_entrez"]
  # problematic_symbols_df <- symbol_mappings_df[!(symbol_mappings_df[["Entrez_source"]] %in% 1), ]

  ## Create an Ensembl-Entrez mapping using BioMart
  message("Attempting to look up Ensembl-Entrez mappings in Ensembl's BioMart (requires an internet connection...)")
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))# host = "useast.ensembl.org"
  ensembl_to_entrez_df <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    values = input_df[["Ensembl_gene_ID"]],
    mart = mart
  )
  message("The BioMart lookup was successful.")

  ## Map TF Ensembl IDs to Entrez IDs
  ensembl_mappings_df <- data.frame(
    input_df[, c("Ensembl_gene_ID", "Gene_symbol")],
    "OrgHs_entrez"         = vapply(ensembl_to_entrez_list[input_df[["Ensembl_gene_ID"]]],
                                    function(x) if (length(x) == 0) NA_character_ else paste0(x, collapse = ", "),
                                    ""
                                    ),
    "OrgHs_num_mappings"   = lengths(ensembl_to_entrez_list[input_df[["Ensembl_gene_ID"]]]),
    "BioMart_entrez"       = as.character(ensembl_to_entrez_df[["entrezgene_id"]][match(input_df[["Ensembl_gene_ID"]], ensembl_to_entrez_df[["ensembl_gene_id"]])]),
    "Symbol_to_entrez"     = symbol_mappings_df[["Entrez_ID"]],
    "Symbol_entrez_symbol" = symbol_mappings_df[["Gene_symbol"]],
    "Symbol_entrez_source" = symbol_mappings_df[["Entrez_source"]],
    stringsAsFactors       = FALSE,
    row.names              = NULL
  )

  are_not_identical <- !(mapply(identical, ensembl_mappings_df[["OrgHs_entrez"]], ensembl_mappings_df[["BioMart_entrez"]], USE.NAMES = FALSE))
  are_NA <- is.na(ensembl_mappings_df[["BioMart_entrez"]])

  ensembl_mappings_df[["Consensus_entrez"]] <- vapply(seq_len(nrow(ensembl_mappings_df)), function(x) {
    this_vec <- unlist(ensembl_mappings_df[x, c("BioMart_entrez", "OrgHs_entrez", "Symbol_to_entrez")])
    if (!(are_NA[[x]]) && !(are_not_identical[[x]])) {
      this_vec[["BioMart_entrez"]]
    } else if (!(is.na(this_vec[["BioMart_entrez"]]))) {
      if ((!(is.na(this_vec[["OrgHs_entrez"]]))) &&
          identical(this_vec[["OrgHs_entrez"]], this_vec[["Symbol_to_entrez"]])
      ) { # majority vote
        this_vec[["OrgHs_entrez"]]
      } else if ((!(is.na(this_vec[["OrgHs_entrez"]]))) &&
                 any(strsplit(this_vec[["OrgHs_entrez"]], ", ", fixed = TRUE)[[1]] %in% this_vec[["Symbol_to_entrez"]])
      ) { # this had to be introduced to keep the number of transcription factors at 2765
        this_vec[["Symbol_to_entrez"]]
      } else {
        this_vec[["BioMart_entrez"]]
      }
    } else if (!(is.na(this_vec[["OrgHs_entrez"]]))) {
      if (grepl(", ", this_vec[["OrgHs_entrez"]], fixed = TRUE)) {
        if (!(grepl(", ", this_vec[["Symbol_to_entrez"]], fixed = TRUE)) &&
            !(is.na(this_vec[["Symbol_to_entrez"]])) &&
            grepl(this_vec[["Symbol_to_entrez"]], this_vec[["OrgHs_entrez"]], fixed = TRUE)
            ) {
          this_vec[["Symbol_to_entrez"]]
        } else {
          NA_character_
        }
      } else {
        this_vec[["OrgHs_entrez"]]
      }
    } else if (!(is.na(this_vec[["Symbol_to_entrez"]])) && !(grepl(", ", this_vec[["Symbol_to_entrez"]], fixed = TRUE))) {
      this_vec[["Symbol_to_entrez"]]
    } else {
      NA_character_
    }
  }, "")

  if (any(grepl(", ", ensembl_mappings_df[["Consensus_entrez"]]))) {
    stop("Ambiguous Entrez ID mappings were found!")
  }

  # Check for duplicate mappings that would lead to fewer than 2765 entries being found in the final database!
  num_occurrences_vec <- table(ensembl_mappings_df[["Consensus_entrez"]])[ensembl_mappings_df[["Consensus_entrez"]]]
  if (any(num_occurrences_vec > 1, na.rm = TRUE)) {
    message("\nProblematic genes:")
    print(ensembl_mappings_df[(num_occurrences_vec > 1) %in% TRUE, ])
    warning_text <- "Multiple ensembl IDs seemed to map to same gene or Entrez ID!"
    if (warn) {
      warning(warning_text)
    } else {
      message(warning_text)
    }
  }

  final_symbols_df <- MapToEntrezs(ensembl_mappings_df[["Consensus_entrez"]], ensembl_mappings_df[["Gene_symbol"]])

  ensembl_mappings_df[["Consensus_symbol"]] <- ifelse(is.na(ensembl_mappings_df[["Consensus_entrez"]]), NA_character_, final_symbols_df[["Gene_symbol"]])
  ensembl_mappings_df[["Original_symbol"]] <- final_symbols_df[["Original_symbol"]]
  ensembl_mappings_df[["Combined_ID"]] <- ifelse(is.na(ensembl_mappings_df[["Consensus_entrez"]]),
                                                 toupper(ensembl_mappings_df[["Original_symbol"]]),
                                                 ensembl_mappings_df[["Consensus_entrez"]]
                                                 )

  ensembl_mappings_df[["Are_conflicting"]] <- are_not_identical
  ensembl_mappings_df[["Are_problematic"]] <- are_NA | are_not_identical

  return(ensembl_mappings_df)
}
