### 16th August 2021 ###




# Functions for processing the datasets from DepMap -----------------------

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
  # stopifnot(!(anyNA(cell_line_matches)))
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
    check.names = FALSE,
    row.names = NULL
  )
  return(results_df)
}



ProcessDEMETERDataDf <- function(data_df, check_replacements = TRUE) {
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
  if (check_replacements) {
    stopifnot(!(any(replaced_names %in% depmap_samples_df[["stripped_cell_line_name"]])))
  }

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



# Functions for exporting data --------------------------------------------

ReplaceEssentialNAs <- function(input_df) {

  NA_columns <- c("Achilles_common", "CRISPR_common",
                  "Hart_3_or_more_lines", "Hart_HeLa",
                  "BlomenHart_intersect", "BlomenHart_intersect_DepMap",
                  "Blomen_HAP1_KBM7_intersect"
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




# Functions for creating gene essentiality histograms ---------------------

DrawHistogram <- function(numeric_vec,
                          add = FALSE,
                          hist_color = brewer.pal(9, "Blues")[[8]],
                          use_breaks = 1000
                          ) {
  points_alpha <- 0.5
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
  if (!(is.na(hist_color))) {
    use_color <- paste0(hist_color, alpha_hex)
  } else {
    use_color <- NA
  }
  hist_results <- hist(numeric_vec,
                       breaks = use_breaks,
                       col    = use_color,
                       border = NA,
                       main   = "Depmap \u2013 all cell lines",
                       xlab   = "CRISPR knockout fitness effect",
                       mgp    = c(2.5, 0.5, 0),
                       freq   = TRUE,
                       add    = add,
                       axes   = FALSE,
                       ylab   = ""
                       )
  if (!(add)) {
    box(bty = "l")
    x_axis_pos <- pretty(par("usr")[c(1, 2)], n = 10)
    axis(1, at = x_axis_pos, mgp = c(2.5, 0.55, 0), tcl = -0.45)
  }
  return(invisible(hist_results))
}



DrawEssentialityHistograms <- function(PNG_dir = file_output_directory) {

  for (make_PNG in c(TRUE, FALSE)) {

    if (make_PNG) {
      png(filename = file.path(PNG_dir, "Histograms - gene effects.png"),
          height = 6, width = 8, units = "in", res = 900
          )
    }

    hist_breaks <- DrawHistogram(as.matrix(CRISPR_effects_df[, 4:ncol(CRISPR_effects_df)]),
                                 hist_color = NA
                                 )[["breaks"]]

    abline(v = seq(-0.1, 0.1, by = 0.1), col = c("gray75", "gray50"), lty = "dashed")


    DrawHistogram(as.matrix(CRISPR_effects_df[, 4:ncol(CRISPR_effects_df)]),
                  hist_color = brewer.pal(9, "Greys")[[4]],
                  add = TRUE, use_breaks = hist_breaks
                  )


    DrawHistogram(as.matrix(CRISPR_effects_df[categ_mat[, "Non-essential"], 4:ncol(CRISPR_effects_df)]),
                  add = TRUE, hist_color = brewer.pal(9, "Greens")[[8]], use_breaks = hist_breaks
                  )

    DrawHistogram(as.matrix(CRISPR_effects_df[categ_mat[, "Intermediate"], 4:ncol(CRISPR_effects_df)]),
                  add = TRUE, hist_color = "#aa6c39", use_breaks = hist_breaks
                  )

    DrawHistogram(as.matrix(CRISPR_effects_df[categ_mat[, "Essential"], 4:ncol(CRISPR_effects_df)]),
                  add = TRUE, hist_color = brewer.pal(9, "Reds")[[8]], use_breaks = hist_breaks
                  )


    legend("topleft",
           legend     = c("All genes", "Non-essential", "Intermediate", "Essential"),
           fill       = c(brewer.pal(9, "Greys")[[4]],
                          brewer.pal(9, "Greens")[[8]],
                          "#aa6c39",
                          brewer.pal(9, "Reds")[[8]]
                          ),
           border    = NA,
           bty       = "n",
           y.intersp = 1.1
           )

    if (make_PNG) {
      dev.off()
    }
  }
  return(invisible(NULL))
}

