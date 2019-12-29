### 24th July 2019 ###




# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "2) Map gene symbols to Entrez IDs.RData"))





# Define functions --------------------------------------------------------

ExpandList <- function(my_list) {
  group_names <- as.character(seq_along(my_list))
  names(my_list) <- paste0(group_names, ".")
  multiples_unlisted <- unlist(my_list, recursive = FALSE)
  groups_unlisted <- strsplit(names(multiples_unlisted), ".", fixed = TRUE)
  groups_unlisted <- sapply(groups_unlisted, "[", 1)
  results_df <- data.frame("Value"      = multiples_unlisted,
                           "List_index" = as.integer(groups_unlisted),
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  return(results_df)
}


GetMinEntrez <- function(entrez_IDs_vec) {
  vapply(strsplit(entrez_IDs_vec, ", ", fixed = TRUE), function(x) min(as.integer(x)), integer(1))
}




MapToEntrezs <- function(entrez_IDs_vec = NULL, symbols_vec = NULL, entrez_IDs_separator = ", ") {

  if (is.null(entrez_IDs_vec) && is.null(symbols_vec)) {
    stop("Either the entrez_IDs_vec or the symbols_vec parameter must be supplied!")
  }
  if (is.null(entrez_IDs_vec)) {
    entrez_IDs_vec <- rep.int(NA_character_, length(symbols_vec))
  }
  if (is.null(symbols_vec)) {
    symbols_vec <- rep.int(NA_character_, length(entrez_IDs_vec))
  }
  stopifnot(length(entrez_IDs_vec) == length(symbols_vec))

  are_split_vec <- grepl(entrez_IDs_separator, entrez_IDs_vec, fixed = TRUE)

  if (any(are_split_vec)) {
    entrez_IDs_split <- strsplit(entrez_IDs_vec, entrez_IDs_separator, fixed = TRUE)
    result_entrezs_vec <- vapply(entrez_IDs_split, function(x) {
      if (all(is.na(x))) {
        return(NA_character_)
      } else {
        are_found <- x %in% entrez_to_symbol_df[, "Entrez_ID"]
        if (!(any(are_found))) {
          return(NA_character_)
        } else {
          return(paste0(x[are_found], collapse = entrez_IDs_separator))
        }
      }
    }, "")
  } else {
    result_entrezs_vec <- ifelse(is.na(entrez_IDs_vec), NA_character_, ifelse(entrez_IDs_vec %in% entrez_to_symbol_df[, "Entrez_ID"], entrez_IDs_vec, NA_character_))
  }

  symbols_to_entrezs_vec <- rep.int(NA_character_, length(symbols_vec))
  translation_method_vec <- rep.int(NA_integer_, length(symbols_vec))

  to_translate <- !(is.na(symbols_vec)) & (toupper(symbols_vec) %in% toupper(symbol_to_entrez_df[, "Symbol"]))

  my_matches <- match(toupper(symbols_vec[to_translate]), toupper(symbol_to_entrez_df[, "Symbol"]))
  are_there <- !(is.na(my_matches))
  my_matches <- my_matches[are_there]

  translated_entrezs_vec <- rep.int(NA_character_, sum(to_translate))
  translated_methods_vec <- rep.int(NA_integer_, sum(to_translate))
  translated_entrezs_vec[are_there] <- ifelse(as.integer(symbol_to_entrez_df[my_matches, "Chosen_via"]) %in% 1:6,
                                              symbol_to_entrez_df[my_matches, "Chosen_entrez"],
                                              symbol_to_entrez_df[my_matches, "Entrez_IDs_all"]
                                              )

  translated_methods_vec[are_there] <- as.integer(symbol_to_entrez_df[my_matches, "Chosen_via"])

  symbols_to_entrezs_vec[to_translate] <- translated_entrezs_vec
  translation_method_vec[to_translate] <- translated_methods_vec

  translation_method_vec[!(is.na(result_entrezs_vec))] <- 0L
  result_entrezs_vec <- ifelse(is.na(result_entrezs_vec), symbols_to_entrezs_vec, result_entrezs_vec)

  expanded_entrezs_df <- ExpandList(strsplit(result_entrezs_vec, entrez_IDs_separator, fixed = TRUE))
  entrez_matches <- match(expanded_entrezs_df[, "Value"], entrez_to_symbol_df[, "Entrez_ID"])

  backtranslated_symbols_vec_expanded <- ifelse(is.na(entrez_to_symbol_df[entrez_matches, "Symbol_Org_Hs_eg_db"]),
                                                entrez_to_symbol_df[entrez_matches, "Symbol_NCBI_Hs_info"],
                                                entrez_to_symbol_df[entrez_matches, "Symbol_Org_Hs_eg_db"]
                                                )
  backtranslated_symbols_vec <- sapply(split(backtranslated_symbols_vec_expanded, expanded_entrezs_df[, "List_index"]),
                                       function(x) if (all(is.na(x))) NA_character_ else paste0(x, collapse = entrez_IDs_separator)
                                       )

  are_identical_entrezs <- vapply(seq_along(entrez_IDs_vec), function(x) identical(entrez_IDs_vec[[x]], result_entrezs_vec[[x]]), logical(1))
  original_entrezs_vec <- ifelse(is.na(entrez_IDs_vec) | are_identical_entrezs, "", entrez_IDs_vec)

  are_identical_symbols <- vapply(seq_along(symbols_vec), function(x) identical(symbols_vec[[x]], backtranslated_symbols_vec[[x]]), logical(1))
  original_symbols_vec <- ifelse(is.na(symbols_vec) | are_identical_symbols, "", symbols_vec)

  results_df <- data.frame(
    "Entrez_ID"       = result_entrezs_vec,
    "Gene_symbol"     = backtranslated_symbols_vec,
    "Original_entrez" = original_entrezs_vec,
    "Original_symbol" = original_symbols_vec,
    "Entrez_source"   = translation_method_vec,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}


SymbolToEntrezDf <- function(symbols_vec) {
  results_df <- MapToEntrezs(symbols_vec = symbols_vec)
  results_df <- results_df[, colnames(results_df) != "Original_entrez"]
  return(results_df)
}



FilterDfEntrezSymbol <- function(gene_df, entrez_IDs, symbols_as_backup) {

  stopifnot(all(c("Entrez_ID", "Gene_symbol") %in% colnames(gene_df)))

  if (grepl(",", entrez_IDs, fixed = TRUE)) {
    entrez_IDs_splits <- strsplit(entrez_IDs, ", ", fixed = TRUE)[[1]]
    are_this_entrez_list <- lapply(entrez_IDs_splits, function(x) grepl(paste0("(^|, )", x, "($|, )"), gene_df[, "Entrez_ID"]))
    are_this_gene <- apply(do.call(cbind, are_this_entrez_list), 1, any)
  } else {
    are_this_gene <- grepl(paste0("(^|, )", entrez_IDs, "($|, )"), gene_df[, "Entrez_ID"])
  }

  if (!(any(are_this_gene))) {
    if (grepl(",", symbols_as_backup, fixed = TRUE)) {
      symbols_splits <- strsplit(symbols_as_backup, ", ", fixed = TRUE)[[1]]
      are_this_symbol_list <- lapply(symbols_splits, function(x) grepl(paste0("(^|, )", x, "($|, )"), gene_df[, "Gene_symbol"]))
      are_this_gene <- apply(do.call(cbind, are_this_symbol_list), 1, any)
    } else {
      are_this_gene <- grepl(paste0("(^|, )", symbols_as_backup, "($|, )"), gene_df[, "Gene_symbol"])
    }
  }

  if (!(any(are_this_gene))) {
    return(NULL)
  } else {
    if (any(duplicated(gene_df[are_this_gene, ]))) {
      assign("delete_gene_df", gene_df, envir = globalenv())
      assign("delete_gene_sub_df", gene_df[are_this_gene, , drop = FALSE], envir = globalenv())
      assign("delete_are_this_gene", are_this_gene, envir = globalenv())
      assign("delete_entrez_IDs", entrez_IDs, envir = globalenv())
      assign("delete_symbols_as_backup", symbols_as_backup, envir = globalenv())
      stop("Duplicated rows were found by the FilterDfEntrezSymbol function!")
    }
    return(gene_df[are_this_gene, , drop = FALSE])
  }
}



FilterAndCombineEntrezSymbol <- function(gene_df, entrez_IDs_vec, symbols_as_backup_vec, IDs_vec) {
  stopifnot(length(unique(c(length(entrez_IDs_vec), length(symbols_as_backup_vec), length(IDs_vec)))) == 1)
  results_list <- lapply(seq_along(entrez_IDs_vec), function(x) FilterDfEntrezSymbol(gene_df, entrez_IDs_vec[[x]], symbols_as_backup_vec[[x]]))
  new_results_list <- lapply(seq_along(IDs_vec), function(x) {
    if (is.null(results_list[[x]])) {
      return(NULL)
    } else {
      return(data.frame("Combined_ID" = IDs_vec[[x]], results_list[[x]], stringsAsFactors = FALSE, row.names = NULL))
    }
  })
  results_df <- do.call(rbind.data.frame, c(new_results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}



GetGenes <- function(symbols_vec, CRISPR_df = CRISPRa_df) {

  symbol_entrez_df <- SymbolToEntrezDf(symbols_vec)

  dummy_df <- CRISPR_df[, c("Entrez_ID", "Original_symbol")]
  colnames(CRISPR_df)[colnames(CRISPR_df) == "Gene_symbol"] <- "New_symbol"
  colnames(CRISPR_df)[colnames(CRISPR_df) == "Original_symbol"] <- "Gene_symbol"

  combined_IDs_vec <- ifelse(is.na(symbol_entrez_df[, "Entrez_ID"]), symbol_entrez_df[, "Gene_symbol"], symbol_entrez_df[, "Entrez_ID"])
  filtered_df <- FilterAndCombineEntrezSymbol(CRISPR_df, symbol_entrez_df[, "Entrez_ID"], symbol_entrez_df[, "Original_symbol"], combined_IDs_vec)
  were_not_found <- !(combined_IDs_vec %in% filtered_df[, "Combined_ID"])

  if (any(were_not_found)) {
    message(paste0("The following ",
                   if (sum(were_not_found) > 1) paste0(sum(were_not_found), " genes were") else "gene was",
                   " not found: ", paste0(symbols_vec[were_not_found], collapse = ", ")
                   ))
  }

  colnames(filtered_df)[colnames(filtered_df) == "Gene_symbol"] <- "Original_symbol"
  colnames(filtered_df)[colnames(filtered_df) == "New_symbol"] <- "Gene_symbol"
  filtered_df <- filtered_df[, 2:ncol(filtered_df)]
  colnames(filtered_df)[[1]] <- "Combined_ID"

  return(filtered_df)
}












