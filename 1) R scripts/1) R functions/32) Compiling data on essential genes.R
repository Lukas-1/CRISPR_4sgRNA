### 16th August 2021 ###




# Define functions --------------------------------------------------------

ProcessAchillesGenesDf <- function(genes_df) {
  # genes_df is a single column in the format A1BG (1) {1 being the Entrez ID}

  gene_splits <- strsplit(genes_df[, "gene"], " ", fixed = TRUE)
  entrezs_vec <- sapply(gene_splits, "[[", 2)
  stopifnot(all(substr(entrezs_vec, 1, 1) == "("))
  stopifnot(all(substr(entrezs_vec, nchar(entrezs_vec), nchar(entrezs_vec)) == ")"))

  stripped_entrezs <- as.integer(substr(entrezs_vec, 2, nchar(entrezs_vec) - 1))
  stripped_entrezs <- as.integer(stripped_entrezs)
  original_symbols <- sapply(gene_splits, "[[", 1)
  results_df <- data.frame(
    "Entrez_ID"       = stripped_entrezs,
    "Original_symbol" = original_symbols,
    stringsAsFactors = FALSE
  )
  CheckEntrezs(stripped_entrezs, results_df)
  return(results_df)
}



CheckEntrezs <- function(entrezs_vec, results_df = NULL) {
  are_present <- as.character(entrezs_vec) %in% names(entrez_to_symbol_list)
  if (!(all(are_present))) {
    message(paste0("The following Entrez IDs were not found in ",
                   "org.Hs.egSYMBOL: ",
                   paste0(entrezs_vec[!(are_present)], collapse = ", ")
                   ))
    if (!(is.null(results_df))) {
      print(results_df[!(are_present), ])
    }
    message("")
  }
  return(invisible(NULL))
}


ProcessAchillesDataDf <- function(data_df, cell_lines_df) {

  data_mat <- t(as.matrix(data_df[, 2:ncol(data_df)]))
  cell_line_matches <- match(data_df[[1]],#[, "DepMap_ID"],
                             cell_lines_df[["DepMap_ID"]]
                             )
  colnames(data_mat) <- cell_lines_df[["stripped_cell_line_name"]][cell_line_matches]
  genes_df <- data.frame("gene" = row.names(data_mat),
                         stringsAsFactors = FALSE
                         )
  genes_df <- ProcessAchillesGenesDf(genes_df)

  results_df <- data.frame(
    genes_df,
    "Mean" = rowMeans(data_mat, na.rm = TRUE),
    data_mat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



ProcessDEMETERDataDf <- function(data_df) {
  stopifnot(all(c("depmap_samples_df", "DEMETER2_samples_df") %in% ls(envir = globalenv())))

  new_order <- order(as.integer(sub("ACH-", "", data_df[[1]])))
  data_df <- data_df[new_order, ]
  row.names(data_df) <- NULL

  data_mat <- t(as.matrix(data_df[, 2:ncol(data_df)]))
  cell_line_matches <- match(data_df[[1]], depmap_samples_df[["DepMap_ID"]])

  DEMETER_matches <- match(data_df[[1]], DEMETER2_samples_df[[1]])
  DEMETER_names <- DEMETER2_samples_df[["Marcotte_name"]][DEMETER_matches]
  have_no_match <- is.na(cell_line_matches)
  replaced_names <- DEMETER_names[have_no_match]
  stopifnot(!(anyNA(replaced_names)))
  replaced_names <- toupper(replaced_names)
  stopifnot(!(any(replaced_names %in% depmap_samples_df[["stripped_cell_line_name"]])))

  colnames(data_mat) <- ifelse(have_no_match,
                               replaced_names,
                               depmap_samples_df[["stripped_cell_line_name"]][cell_line_matches]
                               )

  gene_splits <- strsplit(rownames(data_mat), " ", fixed = TRUE)
  gene_symbols <- sapply(gene_splits, "[[", 1)
  entrez_IDs <- sapply(gene_splits, "[[", 2)
  entrez_IDs <- substr(entrez_IDs, 2, nchar(entrez_IDs) - 1)
  symbol_splits <- strsplit(gene_symbols, "&", fixed = TRUE)
  entrez_splits <- strsplit(entrez_IDs, "&", fixed = TRUE)
  stopifnot(identical(lengths(symbol_splits), lengths(entrez_splits)))
  stopifnot(!(any(unlist(entrez_IDs[lengths(entrez_IDs) == 1]) %in% unlist(entrez_IDs[lengths(entrez_IDs) > 1]))))

  symbols_vec <- unlist(symbol_splits, use.names = FALSE)
  entrezs_vec <- as.integer(unlist(entrez_splits, use.names = FALSE))

  CheckEntrezs(entrezs_vec)
  expanded_indices <- rep(seq_along(entrez_IDs), lengths(entrez_splits))

  expanded_df <- data.frame(
    "Entrez_ID" = entrezs_vec,
    "Gene_symbol" = symbols_vec,
    "Num_targeted_genes" = lengths(entrez_splits)[expanded_indices],
    "Mean" = rowMeans(data_mat[expanded_indices, ], na.rm = TRUE),
    data_mat[expanded_indices, ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  expanded_df <- expanded_df[order(expanded_df[["Entrez_ID"]]), ]
  row.names(expanded_df) <- NULL

  return(expanded_df)
}





GetGeneEssentiality <- function(entrezs_vec, datasets_list) {

  required_objects <- c("achilles_depend_df", "CRISPR_depend_df",
                        "CRISPR_effects_df", "DEMETER2_combined_depend_df"
                        )
  stopifnot(all(required_objects %in% ls(envir = globalenv())))

  stopifnot(identical("integer", typeof(entrezs_vec)))

  AreEssential <- function(essential_genes, all_genes) {
    stopifnot(all(essential_genes %in% all_genes))
    ifelse(!(entrezs_vec %in% all_genes),
           NA,
           entrezs_vec %in% essential_genes
           )
  }

  MetricsFromDependDf <- function(use_depend_df, column_prefix = NULL) {
    matches_vec <- match(entrezs_vec, use_depend_df[, "Entrez_ID"])
    use_depend_mat <- as.matrix(use_depend_df[, 4:ncol(use_depend_df)])
    use_depend_mat <- use_depend_mat[matches_vec, ]
    results_df <- data.frame(
      "Mean_probability" = use_depend_df[["Mean"]][matches_vec],
      "Num_essential"    = ifelse(is.na(matches_vec),
                                  NA,
                                  rowSums(use_depend_mat > 0.5, na.rm = TRUE)
                                  ),
      "Num_cell_lines"   = ifelse(is.na(matches_vec),
                                  NA,
                                  rowSums(!(is.na(use_depend_mat)))
                                  ),
      stringsAsFactors = FALSE
    )
    if (!(is.null(column_prefix))) {
      names(results_df) <- paste0(column_prefix, tolower(names(results_df)))
    }
    return(results_df)
  }

  symbols_vec <- MapToEntrezs(entrez_IDs_vec = as.character(all_entrezs))[["Gene_symbol"]]

  CRISPR_matches <- match(entrezs_vec, CRISPR_effects_df[, "Entrez_ID"])

  genes_list <- list(
    "Entrez_ID"   = entrezs_vec,
    "Gene_symbol" = symbols_vec,
    MetricsFromDependDf(CRISPR_depend_df, "CRISPR_"),
    "CRISPR_mean_effect" = CRISPR_effects_df[["Mean"]][CRISPR_matches],
    MetricsFromDependDf(achilles_depend_df, "Achilles_"),
    MetricsFromDependDf(DEMETER2_combined_depend_df, "DEMETER2_")
  )

  genes_list <- c(
    genes_list,
    lapply(datasets_list, function(x) AreEssential(x[["essential"]], x[["all"]]))
  )

  genes_df <- do.call(data.frame, c(genes_list, stringsAsFactors = FALSE))

  for (i in seq_along(genes_df)) {
    if (is.logical(genes_df[[i]])) {
      genes_df[[i]] <- ifelse(genes_df[[i]], "Essential", "Non-essential")
    }
  }

  return(genes_df)
}



ReplaceEssentialNAs <- function(input_df) {
  NA_columns <- c("Achilles_common", "CRISPR_common", "Hart_3_or_more_lines",
                  "Hart_HeLa", "Blomen_HAP1_KBM7_intersect"
                  )
  no_entrez <- is.na(input_df[["Entrez_ID"]])
  for (column_name in NA_columns) {
    input_df[, column_name] <- ifelse(is.na(input_df[, column_name]),
                                      ifelse(no_entrez, "", "N/A"),
                                      as.character(input_df[, column_name])
                                      )
  }

  for (column_name in setdiff(names(input_df), NA_columns)) {
    input_df[, column_name] <- ifelse(is.na(input_df[, column_name]),
                                      "",
                                      as.character(input_df[, column_name])
                                      )
  }
  return(input_df)
}


