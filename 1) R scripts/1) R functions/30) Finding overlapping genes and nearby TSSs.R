### 20th July 2020 ###


# Import packages and source code -----------------------------------------

# Packages required for the MergeCommonElements function
library("Matrix")
library("igraph")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R")) # For MakeLocationStrings
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R")) # For GetCutLocations



# Define functions for CRISPRi/CRISPRa ------------------------------------

PrepareTSSDf <- function(TSS_df,
                         distance_before             = 1001L,
                         distance_after              = 1001L,
                         only_consistent_chromosomes = TRUE,
                         only_protein_coding         = FALSE,
                         only_best_TSS               = FALSE,
                         check_entrezs_and_symbols   = TRUE
                         ) {

  all_chromosomes <- paste0("chr", c(1:23, "X", "Y", "M"))
  stopifnot(all(TSS_df[["Chromosome"]] %in% all_chromosomes))

  if (check_entrezs_and_symbols) {
    num_symbols_vec <- lengths(SplitCommas(TSS_df[["Gene_symbol"]]))
    num_entrezs_vec <- lengths(SplitCommas(TSS_df[["Entrez_ID"]]))
    are_discrepant <- num_symbols_vec != num_entrezs_vec
    if (any(are_discrepant)) {
      mapped_df <- MapToEntrezs(entrez_IDs_vec = TSS_df[["Entrez_ID"]][are_discrepant])
      stopifnot(identical(mapped_df[["Gene_symbol"]],
                          TSS_df[["Gene_symbol"]][are_discrepant]
                          )
                )
      stopifnot(identical(lengths(SplitCommas(mapped_df[["Entrez_ID"]])),
                          num_symbols_vec[are_discrepant]
                          )
                )
      TSS_df[["Entrez_ID"]][are_discrepant] <- mapped_df[["Entrez_ID"]]
    }
  }

  are_selected <- rep(TRUE, nrow(TSS_df))

  if (only_consistent_chromosomes) {
    are_selected[are_selected] <- TSS_df[["Has_consistent_chromosome"]][are_selected]
  }
  if (only_best_TSS) {
    are_selected[are_selected] <- TSS_df[["Is_chosen_TSS"]][are_selected]
  }
  if (only_protein_coding) {
    are_selected[are_selected] <- grepl("protein-coding", TSS_df[["Gene_type"]][are_selected], fixed = TRUE)
  }

  names(TSS_df)[names(TSS_df) == "Entrez_ID"]       <- "Entrez_IDs"
  names(TSS_df)[names(TSS_df) == "Gene_symbol"]     <- "Gene_symbols"
  names(TSS_df)[names(TSS_df) == "Original_symbol"] <- "Original_symbols"

  TSS_df[["Number_of_Entrez_IDs"]] <- lengths(SplitCommas(TSS_df[["Entrez_IDs"]]))

  select_columns <- c("Entrez_IDs", "Number_of_Entrez_IDs",
                      "Gene_symbols", "Original_symbols", "Gene_type",
                      "Source", "Score", "Is_main_TSS", "TSS",
                      "Entrez_chromosome", "Chromosome", "Strand"
                      )
  if (only_consistent_chromosomes) {
    select_columns <- setdiff(select_columns, "Entrez_chromosome")
  }

  results_df <- data.frame(
    TSS_df[are_selected, select_columns],
    "Start" = TSS_df[["TSS"]][are_selected] - distance_before,
    "End"   = TSS_df[["TSS"]][are_selected] + distance_after,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



FindNearbyTSSs <- function(CRISPR_df,
                           input_TSS_df,
                           only_protein_coding = FALSE,
                           only_best_TSS = FALSE
                           ) {

  stopifnot(!(anyNA(CRISPR_df[["Entrez_ID"]])))
  stopifnot(!(any(grepl(",", CRISPR_df[["Entrez_ID"]], fixed = TRUE))))

  TSS_df <- PrepareTSSDf(input_TSS_df,
                         only_protein_coding = only_protein_coding,
                         only_best_TSS = only_best_TSS
                         )

  stopifnot(!(anyNA(TSS_df[["Chromosome"]])))

  rename_TSS_columns <- c(
    "TSS"          = "TSS_position",
    "Score"        = "FANTOM5_score",
    "Source"       = "TSS_source"
  )

  for (column_name in names(rename_TSS_columns)) {
    names(TSS_df)[names(TSS_df) == column_name] <- rename_TSS_columns[[column_name]]
  }

  split_results <- Process0MMLoci(CRISPR_df)
  unique_loci_df <- split_results[["unique_df"]]
  expanded_0MM_df <- split_results[["expanded_df"]]
  rm(split_results)

  combined_df <- FindOverlappingHits(unique_loci_df, TSS_df)
  full_combined_df <- AlignHits(expanded_0MM_df, combined_df)

  names(full_combined_df)[names(full_combined_df) == "Strand"] <- "TSS_strand"

  full_combined_df <- AddEntrezIDAvailable(full_combined_df, TSS_df)

  combined_columns <- c(
    "Locus",
    "Index", "Is_primary_location",
    "Intended_Entrez_ID", "Entrez_ID_available", "Affected_Entrez_IDs",
    "Number_of_Entrez_IDs",
    "Intended_gene_symbol", "Affected_gene_symbols", "Num_loci",
    "Guide_locus", "Gene_locus", "Chromosome",
    "Is_main_TSS", "TSS_source", "TSS_strand", "TSS_position",
    "Gene_type"
  )

  full_combined_df <- full_combined_df[, combined_columns]

  stopifnot(identical(unique(full_combined_df[["Index"]]),
                      seq_len(nrow(CRISPR_df))
                      )
            )


  summary_df <- SummarizeFullDf(full_combined_df)

  stopifnot(nrow(summary_df) == nrow(CRISPR_df))


  summary_df <- SummarizeSummaryDf(summary_df)

  results_list <- list(
    "summary_df" = summary_df,
    "full_df"    = full_combined_df
  )
  return(results_list)
}




# Define functions for CRISPRko -------------------------------------------

PrepareGenesDf <- function(genes_df,
                           only_protein_coding       = FALSE,
                           exclude_pseudogenes       = TRUE,
                           only_known_gene_type      = TRUE,
                           require_Entrez_ID         = FALSE,
                           check_entrezs_and_symbols = TRUE
                           ) {

  all_chromosomes <- paste0("chr", c(1:23, "X", "Y", "M"))
  are_inside <- genes_df[["Chromosome"]] %in% all_chromosomes

  if (any(!(are_inside))) {
    message(paste0(sum(!(are_inside)), " entries lay outside the main",
            " reference chromosomes and were excluded."
            ))
  }

  are_selected <- are_inside

  are_NA_entrezs <- is.na(genes_df[["Entrez_IDs"]])

  if (check_entrezs_and_symbols) {
    have_no_symbol <- is.na(genes_df[["Gene_symbols"]]) & !(are_NA_entrezs)
    entrezs_without_symbol <- unique(genes_df[["Entrez_IDs"]][have_no_symbol])
    num_symbols <- sum(have_no_symbol)
    num_entrezs <- length(entrezs_without_symbol)
    message(paste0(num_symbols, " entries were associated with ",
                   num_entrezs, " unique Entrez ID",
                   if (num_entrezs > 1) "s" else "",
                   if (num_entrezs < 10) paste0(" (", paste0(entrezs_without_symbol, collapse = ", "), ")") else "",
                   " that could not be mapped to a gene symbol.",
                   " These were excluded."
                   ))
    are_selected[are_selected] <- !(have_no_symbol[are_selected])
  }
  if (require_Entrez_ID) {
    are_selected[are_selected] <- !(are_NA_entrezs[are_selected])
  }

  stopifnot(!(any(grepl(",", genes_df[["Gene_types"]], fixed = TRUE))))

  if (only_protein_coding) {
    are_selected[are_selected] <- genes_df[["Gene_types"]][are_selected] %in% "protein-coding"
  } else if (exclude_pseudogenes || only_known_gene_type) {
    if (exclude_pseudogenes) {
      exclude_types <- "pseudogene"
    } else {
      exclude_types <- c()
    }
    if (only_known_gene_type) {
      exclude_types <- c(exclude_types, "other", "unknown")
    }
    are_selected[are_selected] <- !(genes_df[["Gene_types"]][are_selected] %in% exclude_types)
  }

  if (!(all(are_selected))) {
    genes_df <- genes_df[are_selected, ]
    row.names(genes_df) <- NULL
  }

  genes_df[["Number_of_Entrez_IDs"]] <- lengths(SplitCommas(genes_df[["Entrez_IDs"]]))
  genes_df[["Number_of_gene_IDs"]] <- lengths(SplitCommas(genes_df[["Gene_IDs"]]))
  return(genes_df)
}




FindOverlappingGenes <- function(CRISPR_df,
                                 input_genes_df,
                                 only_protein_coding       = FALSE,
                                 exclude_pseudogenes       = TRUE,
                                 only_known_gene_type      = TRUE,
                                 require_Entrez_ID         = FALSE
                                 ) {


  stopifnot(!(anyNA(CRISPR_df[["Entrez_ID"]])))
  stopifnot(!(any(grepl(",", CRISPR_df[["Entrez_ID"]], fixed = TRUE))))

  genes_df <- PrepareGenesDf(input_genes_df,
                             only_protein_coding  = only_protein_coding,
                             exclude_pseudogenes  = exclude_pseudogenes,
                             only_known_gene_type = only_known_gene_type,
                             require_Entrez_ID    = require_Entrez_ID
                             )

  rename_gene_columns <- c(
    "Gene_IDs" = "Affected_gene_IDs",
    "Source"   = "Gene_source"
  )

  for (column_name in names(rename_gene_columns)) {
    names(genes_df)[names(genes_df) == column_name] <- rename_gene_columns[[column_name]]
  }

  split_results <- Process0MMLoci(CRISPR_df)
  unique_loci_df <- split_results[["unique_df"]]
  expanded_0MM_df <- split_results[["expanded_df"]]
  rm(split_results)

  combined_df <- FindOverlappingHits(unique_loci_df, genes_df)
  full_combined_df <- AlignHits(expanded_0MM_df, combined_df)


  stopifnot(identical(lengths(SplitCommas(full_combined_df[["Affected_Entrez_IDs"]])),
                      lengths(SplitCommas(full_combined_df[["Affected_gene_symbols"]]))
                      )
            )

  full_combined_df <- AddEntrezIDAvailable(full_combined_df, genes_df)

  combined_columns <- c(
    "Locus",
    "Index", "Is_primary_location",
    "Intended_Entrez_ID", "Entrez_ID_available", "Affected_Entrez_IDs",
    "Number_of_gene_IDs", "Number_of_Entrez_IDs",
    "Intended_gene_symbol", "Affected_gene_symbols", "Num_loci",
    "Guide_locus", "Gene_locus", "Chromosome",
    "Affected_gene_IDs", "Ensembl_gene_IDs",
    "Ensembl_transcript_ID", "Exon_ID",
    "Strand", "Gene_types", "Gene_source"
  )

  full_combined_df <- full_combined_df[, combined_columns]

  stopifnot(identical(unique(full_combined_df[["Index"]]),
                      seq_len(nrow(CRISPR_df))
                      )
            )


  summary_df <- SummarizeFullDf(full_combined_df)

  stopifnot(nrow(summary_df) == nrow(CRISPR_df))

  summary_df <- SummarizeSummaryDf(summary_df)

  results_list <- list(
    "summary_df" = summary_df,
    "full_df"    = full_combined_df
  )
  return(results_list)
}





# Define functions for CRISPRko 4sg combinations --------------------------

Get4sgProjectedDeletions <- function(CRISPR_df,
                                     only_primary_location = FALSE,
                                     split_large_deletions = TRUE,
                                     deletion_size_limit   = 10^6
                                     ) {

  stopifnot(!(anyNA(CRISPR_df[["Entrez_ID"]])))
  stopifnot(all(table(CRISPR_df[["Entrez_ID"]]) == 4))

  locations_splits <- strsplit(CRISPR_df[["Locations_0MM"]], "; ", fixed = TRUE)

  primary_locus_vec <- MakeLocationStrings(CRISPR_df)
  are_primary_location_list <- lapply(seq_along(primary_locus_vec),
                                      function(x) locations_splits[[x]] %in% primary_locus_vec[[x]]
                                      )
  are_primary <- unlist(are_primary_location_list, use.names = FALSE)

  entrez_fac <- factor(CRISPR_df[["Entrez_ID"]], levels = unique(CRISPR_df[["Entrez_ID"]]))
  locations_split_splits <- split(locations_splits, entrez_fac)
  locations_split_splits <- lapply(locations_split_splits, unlist, use.names = FALSE)

  expanded_df <- ExpandList(locations_split_splits)
  names(expanded_df) <- c("Locus_0MM", "Index")

  expanded_df[["Intended_Entrez_ID"]] <- names(locations_split_splits)[expanded_df[["Index"]]]
  matches_vec <- match(expanded_df[["Intended_Entrez_ID"]], CRISPR_df[["Entrez_ID"]])

  expanded_df[["Intended_gene_symbol"]] <- CRISPR_df[["Gene_symbol"]][matches_vec]
  expanded_df[["Num_loci"]] <- lengths(locations_split_splits, use.names = FALSE)[expanded_df[["Index"]]]

  if (only_primary_location) {
    expanded_df <- expanded_df[are_primary, ]
  } else {
    expanded_df[["Is_primary_location"]] <- are_primary
  }

  are_NA <- is.na(expanded_df[["Locus_0MM"]])
  use_indices <- rep(NA, nrow(expanded_df))
  use_indices[!(are_NA)] <- seq_len(sum(!(are_NA)))
  locations_df <- LocationStringToDf(expanded_df[["Locus_0MM"]][!are_NA])
  expanded_df <- data.frame(expanded_df,
                            locations_df[use_indices, ],
                            stringsAsFactors = FALSE,
                            row.names = NULL
                            )
  expanded_df[["Cut_location"]] <- GetCutLocations(expanded_df)

  new_order <- order(expanded_df[["Index"]],
                     expanded_df[["Chromosome"]],
                     expanded_df[["Cut_location"]]
                     )
  expanded_df <- expanded_df[new_order, ]
  row.names(expanded_df) <- NULL

  index_chr_vec <- paste0(expanded_df[["Index"]], "__",
                          expanded_df[["Chromosome"]]
                          )
  index_chr_fac <- factor(index_chr_vec, levels = unique(index_chr_vec))

  if (split_large_deletions) {
    block_ID_vec <- unlist(tapply(expanded_df[["Cut_location"]],
                                  index_chr_fac,
                                  function(x) {
                                    block_ID <- 1L
                                    current_location <- x[[1]]
                                    results_vec <- rep(NA, length(x))
                                    for (i in seq_along(x)) {
                                      if (all(is.na(x))) {
                                        return(rep(1L, length(x)))
                                      }
                                      if ((x[[i]] - current_location) > deletion_size_limit) {
                                        block_ID <- block_ID + 1L
                                      }
                                      results_vec[[i]] <- block_ID
                                      current_location <- x[[i]]
                                    }
                                    return(results_vec)
                                  }, simplify = FALSE), use.names = FALSE)
  } else {
    block_ID_vec <- rep(1L, nrow(expanded_df))
  }

  expanded_df[, "Block"] <- block_ID_vec

  index_chr_block_vec <- paste0(index_chr_vec, "__", block_ID_vec)

  are_first <- !(duplicated(index_chr_block_vec))
  are_last  <- !(duplicated(index_chr_block_vec, fromLast = TRUE))

  retain_columns <- c("Intended_Entrez_ID", "Intended_gene_symbol",
                      "Index", "Chromosome", "Block"
                      )
  reduced_df <- data.frame(expanded_df[are_first, retain_columns],
                           "Num_cuts" = tabulate(factor(index_chr_block_vec,
                                                        levels = unique(index_chr_block_vec)
                                                        )
                                                 ),
                           "Num_loci" = NA,
                           "Locus" = NA,
                           "Start" = expanded_df[["Cut_location"]][are_first] - 1L,
                           "End"   = expanded_df[["Cut_location"]][are_last],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  reduced_df[["Locus"]] <- StandardLocationString(reduced_df)
  reduced_df[["Span"]] <- reduced_df[["End"]] - reduced_df[["Start"]]


  num_loci <- tabulate(reduced_df[["Index"]])
  reduced_df[["Num_loci"]] <- rep(num_loci, num_loci)


  results_list <- list("reduced_df" = reduced_df,
                       "expanded_df" = expanded_df
                       )
  return(results_list)
}



FindOverlapsWithDeletions <- function(CRISPR_df,
                                      input_genes_df,
                                      only_protein_coding       = FALSE,
                                      exclude_pseudogenes       = TRUE,
                                      only_known_gene_type      = TRUE,
                                      require_Entrez_ID         = FALSE
                                      ) {

  stopifnot(!(any(grepl(",", CRISPR_df[["Entrez_ID"]], fixed = TRUE))))

  genes_df <- PrepareGenesDf(input_genes_df,
                             only_protein_coding  = only_protein_coding,
                             exclude_pseudogenes  = exclude_pseudogenes,
                             only_known_gene_type = only_known_gene_type,
                             require_Entrez_ID    = require_Entrez_ID
                             )

  rename_gene_columns <- c(
    "Gene_IDs" = "Affected_gene_IDs",
    "Source"   = "Gene_source"
  )
  for (column_name in names(rename_gene_columns)) {
    names(genes_df)[names(genes_df) == column_name] <- rename_gene_columns[[column_name]]
  }

  deletion_results <- Get4sgProjectedDeletions(CRISPR_df)

  reduced_df <- deletion_results[["reduced_df"]]

  are_NA <- is.na(reduced_df[["Start"]])

  combined_df <- FindOverlappingHits(reduced_df[!(are_NA), ], genes_df)
  full_combined_df <- AlignHits(reduced_df[, !(names(reduced_df) %in% c("Chromosome", "Start", "End"))], combined_df)

  stopifnot(identical(lengths(SplitCommas(full_combined_df[["Affected_Entrez_IDs"]])),
                      lengths(SplitCommas(full_combined_df[["Affected_gene_symbols"]]))
                      )
            )

  full_combined_df <- AddEntrezIDAvailable(full_combined_df, genes_df)

  combined_columns <- c(
    "Index", "Intended_Entrez_ID", "Entrez_ID_available", "Affected_Entrez_IDs",
    "Number_of_gene_IDs", "Number_of_Entrez_IDs",
    "Intended_gene_symbol", "Affected_gene_symbols", "Num_loci",
    "Gene_locus", "Chromosome",
    "Affected_gene_IDs", "Ensembl_gene_IDs", "Ensembl_transcript_ID", "Exon_ID",
    "Strand", "Gene_types", "Gene_source",
    "Block", "Span", "Num_cuts",
    "Locus", "Guide_locus"
  )

  full_combined_df <- full_combined_df[, combined_columns]

  num_entrezs <- length(unique(unique(CRISPR_df[["Entrez_ID"]])))
  stopifnot(identical(unique(full_combined_df[["Index"]]), seq_len(num_entrezs)))

  summary_df <- SummarizeFullDf(full_combined_df, tolerate_num_affected = TRUE)

  stopifnot(nrow(summary_df) == num_entrezs)

  summary_df <- SummarizeSummaryDf(summary_df)

  full_combined_df <- full_combined_df[, names(full_combined_df) != "Guide_locus"]
  names(full_combined_df)[names(full_combined_df) == "Locus"] <- "Deletion_locus"

  results_list <- list(
    "summary_df" = summary_df,
    "full_df"    = full_combined_df
  )
  return(results_list)
}





# Define general functions ------------------------------------------------

AlignSummaryDf <- function(UseFunction, CRISPR_df, genes_df, ...) {

  are_valid <- !(is.na(CRISPR_df[["Entrez_ID"]])) &
               !(grepl(",", CRISPR_df[["Entrez_ID"]], fixed = TRUE))

  indices_vec <- rep(NA, nrow(CRISPR_df))
  indices_vec[are_valid] <- seq_len(sum(are_valid))

  results_list <- UseFunction(CRISPR_df[are_valid, ], genes_df, ...)

  results_list[["summary_df"]] <- results_list[["summary_df"]][indices_vec, ]
  row.names(results_list[["summary_df"]]) <- NULL
  return(results_list)
}





SplitCommas <- function(char_vec) {
  strsplit(char_vec, ", ", fixed = TRUE)
}

MergeCommonElements <- function(input_list) {
  # from https://stackoverflow.com/a/47328161
  i <- rep(seq_along(input_list), lengths(input_list))
  j <- factor(unlist(input_list))
  tab <- sparseMatrix(i = i, j = as.integer(j), x = TRUE, dimnames = list(NULL, levels(j)))
  connects <- tcrossprod(tab, boolArith = TRUE)
  group <- clusters(graph_from_adjacency_matrix(as(connects, "lsCMatrix")))[["membership"]]
  results_list <- tapply(input_list,
                         group,
                         function(x) unique(unlist(x)),
                         simplify = FALSE
                         )
  return(results_list)
}






Process0MMLoci <- function(CRISPR_df) {
  loci_splits <- strsplit(CRISPR_df[, "Locations_0MM"], "; ", fixed = TRUE)

  primary_locus_vec <- MakeLocationStrings(CRISPR_df)
  expanded_df <- ExpandList(loci_splits)

  names(expanded_df) <- c("Locus", "Index")

  index_vec <- expanded_df[["Index"]]
  are_primary <- primary_locus_vec[index_vec] == expanded_df[["Locus"]]
  expanded_df[["Is_primary_location"]] <- are_primary

  expanded_df[["Intended_Entrez_ID"]] <- CRISPR_df[["Entrez_ID"]][index_vec]
  expanded_df[["Intended_gene_symbol"]] <- CRISPR_df[["Gene_symbol"]][index_vec]
  expanded_df[["Num_loci"]] <- CRISPR_df[["Num_0MM"]][index_vec]

  are_not_NA <- !(is.na(expanded_df[["Locus"]]))
  unique_loci_vec <- unique(expanded_df[["Locus"]][are_not_NA])
  unique_loci_df <- data.frame("Locus" = unique_loci_vec,
                               LocationStringToDf(unique_loci_vec),
                               stringsAsFactors = FALSE
                               )
  unique_loci_df[["Cut_location"]] <- GetCutLocations(unique_loci_df)
  results_list <- list(
    "expanded_df" = expanded_df,
    "unique_df"   = unique_loci_df
  )
  return(results_list)
}



AddEntrezIDAvailable <- function(full_df, genes_df) {
  available_entrezs <- unique(unlist(strsplit(genes_df[, "Entrez_IDs"], ", ", fixed = TRUE), use.names = FALSE))
  full_df[, "Entrez_ID_available"] <- full_df[, "Intended_Entrez_ID"] %in% available_entrezs
  return(full_df)
}



StandardLocationString <- function(ranges_df) {
  are_NA <- is.na(ranges_df[, "Chromosome"]) |
            is.na(ranges_df[, "Start"]) |
            is.na(ranges_df[, "End"])

  results_vec <- paste0(ranges_df[, "Chromosome"],
                        ":",
                        ranges_df[, "Start"],
                        "-",
                        ranges_df[, "End"]
                        )
  results_vec <- ifelse(are_NA, NA_character_, results_vec)
  if ("Strand" %in% names(ranges_df)) {
    results_vec <- paste0(results_vec, ":", ranges_df[, "Strand"])
    results_vec <- ifelse(is.na(ranges_df[, "Strand"]), NA, results_vec)
  }
  return(results_vec)
}



FindOverlappingHits <- function(loci_df,
                                genes_df,
                                retain_columns = c("Locus", "Guide_locus")
                                ) {

  if ("Cut_location" %in% names(loci_df)) {
    sgRNA_GRanges_object <- CutLocationsToGRangesObject(loci_df)
  } else {
    sgRNA_GRanges_object <- RangesDfToGRangesObject(loci_df)
  }
  genes_GRanges_object <- RangesDfToGRangesObject(genes_df)

  hits_object <- findOverlaps(sgRNA_GRanges_object,
                              genes_GRanges_object,
                              ignore.strand = TRUE,
                              select = "all"
                              )

  loci_df[["Guide_locus"]] <- StandardLocationString(loci_df)
  genes_df[["Gene_locus"]] <- StandardLocationString(genes_df)

  combined_df <- data.frame(
    loci_df[queryHits(hits_object), retain_columns],
    genes_df[subjectHits(hits_object), ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  rename_columns <- c(
    "Entrez_IDs"   = "Affected_Entrez_IDs",
    "Gene_symbols" = "Affected_gene_symbols"
  )
  for (column_name in names(rename_columns)) {
    names(combined_df)[names(combined_df) == column_name] <- rename_columns[[column_name]]
  }
  return(combined_df)
}





AlignHits <- function(expanded_df, comb_df) {

  were_found <- expanded_df[["Locus"]] %in% comb_df[["Locus"]]
  were_found_indices <- unique(expanded_df[["Index"]][were_found])

  are_to_retain <- (were_found) |
                   (!(duplicated(expanded_df[["Index"]])) &
                    !(expanded_df[["Index"]] %in% were_found_indices)
                    )

  indices_list <- split(
    seq_len(nrow(comb_df)),
    factor(comb_df[["Locus"]],
           levels = unique(expanded_df[["Locus"]][are_to_retain]),
           exclude = c()
           )
  )

  indices_list[lengths(indices_list) == 0] <- list(NA_integer_)

  matches_vec <- match(expanded_df[["Locus"]][are_to_retain], names(indices_list))

  num_reps_vec <- lengths(indices_list)[matches_vec]

  expanded_indices <- rep(which(are_to_retain), times = num_reps_vec)

  combined_indices <- unlist(indices_list[matches_vec], use.names = FALSE)

  full_df <- data.frame(expanded_df[expanded_indices, ],
                        comb_df[combined_indices, names(comb_df) != "Locus"],
                        stringsAsFactors = FALSE,
                        row.names = NULL
                        )
  return(full_df)
}




SplitAndRejoin <- function(long_df, use_column) {

  splits_list <- SplitCommas(long_df[, use_column])

  joined_list <- tapply(splits_list,
                        factor(long_df[, "Index"]),
                        function(x) {
                          results_vec <- unique(unlist(x, use.names = FALSE))
                          results_vec <- results_vec[!(is.na(results_vec))]
                          return(results_vec)
                        },
                        simplify = FALSE
                        )

  joined_vec <- vapply(joined_list,
                       function(x) {
                         if (length(x) == 0) {
                           NA_character_
                         } else {
                           paste0(x, collapse = ", ")
                         }
                       }, "")

  results_list <- list(
    "list"   = joined_list,
    "vector" = joined_vec
  )
  return(results_list)
}




CollapseList <- function(input_list, use_sep = ", ") {
  vapply(input_list,
         function(x) {
           if (all(is.na(x))) {
             NA_character_
           } else {
             paste0(x, collapse = use_sep)
           }
         }, "")
}



SummarizeFullDf <- function(full_df, tolerate_num_affected = FALSE) {

  ## Remove the "ncbi:" prefix before Entrez IDs
  use_gene_IDs <- "Affected_gene_IDs" %in% names(full_df)
  if (use_gene_IDs) {
    full_df[["Affected_gene_IDs"]] <- gsub("ncbi:", "", full_df[["Affected_gene_IDs"]], fixed = TRUE)
    ID_column <- "Affected_gene_IDs"
  } else {
    ID_column <- "Affected_Entrez_IDs"
  }

  ## Re-order the data frame to ensure that the intended gene is listed first,
  ## and that all other genes are ordered by their Entrez ID
  unique_entrezs <- unique(full_df[["Affected_Entrez_IDs"]])
  min_unique_entrezs <- GetMinEntrez(unique_entrezs)
  min_entrezs_vec <- min_unique_entrezs[match(full_df[["Affected_Entrez_IDs"]], unique_entrezs)]
  full_df[["Min_Entrez_ID"]] <- min_entrezs_vec
  new_order <- order(
    full_df[["Index"]],
    full_df[["Intended_Entrez_ID"]] != full_df[["Affected_Entrez_IDs"]],
    min_entrezs_vec,
    match(full_df[["Locus"]], full_df[["Locus"]])
  )
  by_gene_df <- full_df
  by_gene_df[["Original_order"]] <- seq_len(nrow(by_gene_df))
  by_gene_df <- by_gene_df[new_order, ]
  row.names(by_gene_df) <- NULL


  ## Create a summary data frame
  gene_ID_df <- AnalyzeGeneIDs(by_gene_df, ID_column)
  symbol_df <- AnalyzeGeneIDs(by_gene_df, "Affected_gene_symbols")
  shared_columns <- intersect(names(gene_ID_df), names(symbol_df))
  if (tolerate_num_affected) {
    use_col <- "Num_affected_genes"
    compare_columns <- setdiff(shared_columns, use_col)
    if (!(identical(gene_ID_df[[use_col]], symbol_df[[use_col]]))) {
      are_equal <- gene_ID_df[[use_col]] == symbol_df[[use_col]]
      message(paste0(
        "For the following ", sum(!(are_equal)), " genes, ",
        " the number of affected genes differs when determined ",
        " using Entrez IDs or gene symbols: ",
        paste0(symbol_df[["Intended_gene_symbol"]][!(are_equal)], collapse = ", ")
      ))
    }
  } else {
    compare_columns <- shared_columns
  }
  assign("delete_gene_ID_df", gene_ID_df, envir = globalenv())
  assign("delete_symbol_df", symbol_df, envir = globalenv())
  assign("delete_shared_columns", shared_columns, envir = globalenv())
  stopifnot(identical(gene_ID_df[, compare_columns], symbol_df[, compare_columns]))

  first_indices <- match(unique(full_df[["Index"]]),
                         full_df[["Index"]]
                         )
  total_loci <- full_df[["Num_loci"]][first_indices]

  summary_df <- data.frame(
    gene_ID_df[, setdiff(names(gene_ID_df), shared_columns)],
    symbol_df[, setdiff(names(symbol_df), shared_columns)],
    gene_ID_df["Num_affected_genes"],
    "Total_loci" = total_loci,
    gene_ID_df[, setdiff(shared_columns, "Num_affected_genes")],
    stringsAsFactors = FALSE
  )


  ## Identify the strand(s) of affected genes
  if ("TSS_strand" %in% names(by_gene_df)) {
    strand_column <- "TSS_strand"
  } else {
    strand_column <- "Strand"
  }
  strands_vec <- tapply(by_gene_df[[strand_column]],
                        by_gene_df[["Index"]],
                        function(x) {
                          has_pos <- "+" %in% x
                          has_neg <- "-" %in% x
                          if (has_pos && has_neg) {
                            "Both"
                          } else if (has_neg) {
                            "-"
                          } else if (has_pos) {
                            "+"
                          } else {
                            NA_character_
                          }
                        })
  summary_df[["Affected_genes_strand"]] <- strands_vec

  ## List all perfect-match loci that target a gene
  loci_list <- tapply(by_gene_df[["Guide_locus"]], by_gene_df[["Index"]], unique)
  summary_df[["Targeting_loci"]] <- CollapseList(loci_list)

  ## List perfect-match loci that target a gene,
  ## but only list a maximum of 3 loci (and truncate the rest)
  num_loci <- lengths(loci_list, use.names = FALSE)
  are_to_truncate <- num_loci > 3
  loci_list[are_to_truncate] <- lapply(loci_list[are_to_truncate], "[", 1:3)
  loci_vec <- CollapseList(loci_list)
  loci_vec[are_to_truncate] <- paste0(num_loci[are_to_truncate],
                                      " perfect-match hits that may target",
                                      " genes were found. The first 3 are: ",
                                      loci_vec[are_to_truncate]
                                      )
  summary_df[["Targeting_loci_truncated"]] <- loci_vec

  return(summary_df)
}




TruncateListList <- function(list_of_lists, max_items = 20L) {
  results_list <- lapply(list_of_lists, function(x) {
    cum_sum <- cumsum(lengths(x))
    are_to_truncate <- cum_sum > max_items
    if (!(any(are_to_truncate))) {
      return(x)
    }
    last_index <- which(are_to_truncate)[[1]]
    if (last_index == 1) {
      use_length <- max_items
      before_list <- list()
    } else {
      use_length <- max_items - cum_sum[[last_index - 1L]]
      before_list <- x[seq_len(last_index - 1L)]
    }
    trunc_list <- x[[last_index]][seq_len(use_length)]
    return(c(before_list, trunc_list))
  })
}





CollapseListList <- function(list_of_lists) {
  vapply(list_of_lists,
         function(x) paste0(vapply(x, paste0, collapse = ", ", ""), collapse = "; "),
         ""
         )
}




AnalyzeGeneIDs <- function(full_df, use_column) {

  assign("delete_my_full_df", full_df, envir = globalenv())

  if (use_column == "Affected_gene_symbols") {
    intended_column <- "Intended_gene_symbol"
  } else if (use_column %in% c("Affected_Entrez_IDs", "Affected_gene_IDs")) {
    intended_column <- "Intended_Entrez_ID"
  } else {
    stop("Unsupported column name!")
  }

  index_loci_vec <- paste0(full_df[["Index"]], "__", full_df[["Locus"]])
  are_duplicated <- duplicated(index_loci_vec)
  locus_fac <- factor(index_loci_vec, levels = index_loci_vec[!(are_duplicated)])
  index_vec <- full_df[["Index"]][!(are_duplicated)]

  are_NA_IDs <- is.na(full_df[[use_column]])

  ID_splits <- SplitCommas(full_df[[use_column]][!(are_NA_IDs)])

  ID_locus_splits <- tapply(ID_splits,
                            locus_fac[!(are_NA_IDs)],
                            function(x) unlist(x, use.names = FALSE),
                            simplify = FALSE
                            )


  second_splits <- split(ID_locus_splits, index_vec)
  second_splits <- lapply(second_splits, function(x) x[lengths(x) > 0])
  second_splits <- lapply(second_splits,
                          function(x) lapply(x, function(y) unique(unlist(SplitCommas(y))))
                          )


  ## Count the number of loci that affect at least one gene
  loci_with_targets <- lengths(second_splits, use.names = FALSE)



  ## Count the number of loci that target the intended gene
  first_indices <- match(unique(full_df[["Index"]]),
                         full_df[["Index"]]
                         )
  intended_IDs <- full_df[[intended_column]][first_indices]
  stopifnot(!(anyNA(intended_IDs)))
  loci_targeting_intended <- vapply(seq_along(second_splits),
                                    function(x) {
                                      sum(vapply(second_splits[[x]], function(y) intended_IDs[[x]] %in% y, logical(1)))
                                    }, integer(1))


  ## Count the number of distinct loci
  ## (i.e. affecting a distinct, non-overlapping set of genes)
  distinct_splits <- lapply(second_splits, unique)

  have_multiple_loci <- loci_with_targets > 1

  distinct_splits[have_multiple_loci] <- lapply(distinct_splits[have_multiple_loci],
                                                MergeCommonElements
                                                )

  distinct_loci <- lengths(distinct_splits, use.names = FALSE)


  ## Count the number of distinct loci that affect more than one gene
  distinct_multiple <- vapply(distinct_splits,
                              function(x) sum(lengths(x) > 1),
                              integer(1)
                              )


  ## Count the number of affected genes
  unlisted_splits <- lapply(distinct_splits, function(x) unlist(x, use.names = FALSE))
  num_genes <- lengths(unlisted_splits, use.names = FALSE)



  ## Report whether the intended gene is available
  are_available <- full_df[["Entrez_ID_available"]][first_indices]



  ## List the genes targeted by each gRNA or gRNA combination
  joined_IDs <- CollapseListList(distinct_splits)


  ## List the genes targeted, but truncate lists that would become too long
  trunc_max <- 20L
  trunc_list <- TruncateListList(distinct_splits, max_items = trunc_max)
  were_truncated <- num_genes > trunc_max
  trunc_IDs <- CollapseListList(trunc_list)
  trunc_IDs[were_truncated] <- paste0(num_genes[were_truncated],
                                      " genes are potential targets. ",
                                      "The first ", trunc_max, " are: ",
                                      trunc_IDs[were_truncated]
                                      )


  ## List only unintended genes targeted by each gRNA or gRNA combination
  unintended_splits <- lapply(seq_along(distinct_splits),
                              function(x)  {
                                results_list <- lapply(distinct_splits[[x]], function(y) {
                                  y[y != intended_IDs[[x]]]
                                })
                                results_list <- results_list[lengths(results_list) > 0]
                                return(results_list)
                              })
  unintended_IDs <- CollapseListList(unintended_splits)


  ## List only unintended genes, and truncate lists that become too long
  trunc_unintended_list <- TruncateListList(unintended_splits, max_items = trunc_max)
  unlisted_unintended_splits <- lapply(unintended_splits, function(x) unlist(x, use.names = FALSE))
  num_unintended <- lengths(unlisted_unintended_splits, use.names = FALSE)
  were_truncated <- num_unintended > trunc_max
  trunc_unintended_IDs <- CollapseListList(trunc_unintended_list)
  trunc_unintended_IDs[were_truncated] <- paste0(num_genes[were_truncated],
                                                 " genes are unintended targets. ",
                                                 "The first ", trunc_max, " are: ",
                                                 trunc_unintended_IDs[were_truncated]
                                                 )


  ## Add a note when the intended gene is not targeted
  target_intended <- loci_targeting_intended >= 1
  target_unintended <- num_unintended >= 1
  intended_prefix <- paste0(intended_IDs,
                            ifelse(are_available,
                                   " is not targeted!",
                                   " is not in the gene location database!"
                                   )
                            )

  intended_postfix <- paste0(" Instead, the following gene",
                             ifelse(num_unintended == 1,
                                    " is an unintended target",
                                    "s are unintended targets"
                                    ),
                             ": "
                             )
  intended_strings <- ifelse(target_intended,
                             "",
                             paste0(intended_prefix,
                                    ifelse(target_unintended, intended_postfix, "")
                                    )
                             )
  intended_strings_truncated <- ifelse(target_intended,
                                      "",
                                      paste0(intended_prefix,
                                             ifelse(target_unintended,
                                                    ifelse(were_truncated,
                                                           "Instead, ",
                                                           intended_postfix
                                                           ),
                                                    ""
                                                    )
                                             )
                                      )
  unintended_IDs <- paste0(intended_strings, unintended_IDs)
  trunc_unintended_IDs <- paste0(intended_strings_truncated, trunc_unintended_IDs)

  results_df <- data.frame(
    "Intended_gene"                       = intended_IDs,
    "Entrez_ID_available"                 = are_available,
    "Affected_genes"                      = joined_IDs,
    "Affected_genes_truncated"            = trunc_IDs,
    "Unintended_genes"                    = unintended_IDs,
    "Unintended_genes_truncated"          = trunc_unintended_IDs,
    "Num_affected_genes"                  = num_genes,
    "Loci_with_targets"                   = loci_with_targets,
    "Loci_targeting_intended_gene"        = loci_targeting_intended,
    "Distinct_loci"                       = distinct_loci,
    "Distinct_loci_with_multiple_targets" = distinct_multiple,
    stringsAsFactors = FALSE
  )
  names(results_df)[names(results_df) == "Intended_gene"] <- intended_column
  names(results_df)[names(results_df) == "Affected_genes"] <- use_column
  names(results_df)[names(results_df) == "Affected_genes_truncated"] <- paste0(use_column, "_truncated")
  column_suffix <- sub("Affected_", "", use_column, fixed = TRUE)
  names(results_df)[names(results_df) == "Unintended_genes"] <- paste0("Unintended_", column_suffix)
  names(results_df)[names(results_df) == "Unintended_genes_truncated"] <- paste0("Unintended_", column_suffix, "_truncated")
  return(results_df)
}







SummarizeSummaryDf <- function(summary_df) {

  summary_df[["Affects_intended_gene"]]       <- ifelse(summary_df[["Entrez_ID_available"]],
                                                        summary_df[["Loci_targeting_intended_gene"]] >= 1,
                                                        NA
                                                        )
  summary_df[["Affects_genes_at_other_loci"]] <- summary_df[["Distinct_loci"]] >= 2

  summary_df[["Affects_unintended_gene"]]     <- (summary_df[["Num_affected_genes"]] > 1) |
                                                 ((summary_df[["Num_affected_genes"]] == 1) &
                                                  (summary_df[["Loci_targeting_intended_gene"]] == 0))

  category_mat <- cbind(
    "Only intended"                               = (summary_df[["Loci_targeting_intended_gene"]] >= 1) &
                                                    (summary_df[["Distinct_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] == 1),

    "Intended and unintended (in the same locus)" = (summary_df[["Total_loci"]] == 1) &
                                                    (summary_df[["Loci_targeting_intended_gene"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] >= 2),

    "Intended and unintended (in the same loci)"  = (summary_df[["Loci_targeting_intended_gene"]] >= 1) &
                                                    (summary_df[["Distinct_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] >= 2) &
                                                    (summary_df[["Total_loci"]] >= 2),

    "Intended and unintended (in other loci)"     = (summary_df[["Loci_targeting_intended_gene"]] >= 1) &
                                                    (summary_df[["Distinct_loci"]] >= 2),

    "No location data for the intended gene"      = !(summary_df[["Entrez_ID_available"]]),

    "No hits in the reference genome"             = summary_df[["Entrez_ID_available"]] & # Some genes may not be located on the
                                                                                          # main chromosomes of the reference genome,
                                                                                          # so the fact of whether or not location data
                                                                                          # is available should take precedence.
                                                    (summary_df[["Total_loci"]] == 0),

    "No targets (single locus)"                   = summary_df[["Entrez_ID_available"]] &
                                                    (summary_df[["Total_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] == 0),

    "Only unintended (in the same locus)"         = summary_df[["Entrez_ID_available"]] &
                                                    (summary_df[["Total_loci"]] == 1) &
                                                    (summary_df[["Distinct_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] >= 1) &
                                                    (summary_df[["Loci_targeting_intended_gene"]] == 0),

    "Only unintended (in the same loci)"          = summary_df[["Entrez_ID_available"]] &
                                                    (summary_df[["Total_loci"]] >= 2) &
                                                    (summary_df[["Distinct_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] >= 1) &
                                                    (summary_df[["Loci_targeting_intended_gene"]] == 0),

    "No targets (multiple loci)"                  = summary_df[["Entrez_ID_available"]] &
                                                    (summary_df[["Total_loci"]] >= 2) &
                                                    (summary_df[["Num_affected_genes"]] == 0),

    "Only unintended (multiple loci)"             = summary_df[["Entrez_ID_available"]] &
                                                    (summary_df[["Distinct_loci"]] >= 2) &
                                                    (summary_df[["Loci_targeting_intended_gene"]] == 0)
  )

  assign("delete_category_mat", category_mat, envir = globalenv())
  assign("delete_summary_df", summary_df, envir = globalenv())

  stopifnot(rowSums(category_mat) == 1)

  index_vec <- rep(0L, nrow(summary_df))
  for (i in seq_len(ncol(category_mat))) {
    index_vec <- index_vec + (i * category_mat[, i])
  }
  summary_df[["Gene_targets_summary"]] <- factor(colnames(category_mat)[index_vec],
                                                 levels = colnames(category_mat)
                                                 )
  return(summary_df)
}




