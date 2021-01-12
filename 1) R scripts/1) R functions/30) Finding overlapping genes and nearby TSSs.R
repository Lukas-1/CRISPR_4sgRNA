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
    num_symbols_vec <- lengths(strsplit(TSS_df[["Gene_symbol"]], ", ", fixed = TRUE))
    num_entrezs_vec <- lengths(strsplit(TSS_df[["Entrez_ID"]], ", ", fixed = TRUE))
    are_discrepant <- num_symbols_vec != num_entrezs_vec
    if (any(are_discrepant)) {
      mapped_df <- MapToEntrezs(entrez_IDs_vec = TSS_df[["Entrez_ID"]][are_discrepant])
      stopifnot(identical(mapped_df[["Gene_symbol"]],
                          TSS_df[["Gene_symbol"]][are_discrepant]
                          )
                )
      stopifnot(identical(lengths(strsplit(mapped_df[["Entrez_ID"]], ", ", fixed = TRUE)),
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

  TSS_df[["Number_of_Entrez_IDs"]] <- lengths(strsplit(TSS_df[["Entrez_IDs"]], ", ", fixed = TRUE))

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

  combined_columns <- c(
    "Locus_0MM",
    "Index", "Is_primary_location", "Intended_Entrez_ID", "Affected_Entrez_IDs",
    "Number_of_Entrez_IDs",
    "Intended_gene_symbol", "Affected_gene_symbols", "Num_loci",
    "Guide_locus", "Gene_locus", "Chromosome",
    "Is_main_TSS", "TSS_source", "TSS_strand", "TSS_position"
  )

  full_combined_df <- full_combined_df[, combined_columns]

  stopifnot(identical(unique(full_combined_df[["Index"]]),
                      seq_len(nrow(CRISPR_df))
                      )
            )


  summary_df <- SummarizeFullDf(full_combined_df)

  stopifnot(nrow(summary_df) == nrow(CRISPR_df))


  TSS_strands_result <- SplitAndRejoin(full_combined_df, "TSS_strand")

  summary_df[["Affected_genes_strand"]] <- ifelse(TSS_strands_result[["vector"]] == "+",
                                                  "+",
                                                  ifelse(TSS_strands_result[["vector"]] == "-",
                                                         "-",
                                                         "Both"
                                                         )
                                                  )

  summary_df <- SummarizeSummaryDf(summary_df)

  return(summary_df)
}




# Define functions for CRISPRko -------------------------------------------

FindOverlappingGenes <- function(ranges_df,
                                 genes_df,
                                 only_protein_coding = FALSE,
                                 exclude_pseudogenes = FALSE
                                 ) {

}




# Define general functions ------------------------------------------------

MergeCommonElements <- function(input_list) {
  # from https://stackoverflow.com/a/47328161
  i <- rep(seq_along(input_list), lengths(input_list))
  j <- factor(unlist(input_list))
  tab <- sparseMatrix(i = i, j = as.integer(j), x = TRUE, dimnames = list(NULL, levels(j)))
  connects <- tcrossprod(tab, boolArith = TRUE)
  group <- clusters(graph_from_adjacency_matrix(as(connects, "lsCMatrix")))[["membership"]]
  results_list <- tapply(input_list,
                         group,
                         function(x) sort(unique(unlist(x))),
                         simplify = FALSE
                         )
  return(results_list)
}


Process0MMLoci <- function(CRISPR_df) {
  loci_splits <- strsplit(CRISPR_df[, "Locations_0MM"], "; ", fixed = TRUE)

  primary_locus_vec <- MakeLocationStrings(CRISPR_df)
  expanded_df <- ExpandList(loci_splits)

  names(expanded_df) <- c("Locus_0MM", "Index")

  index_vec <- expanded_df[["Index"]]
  are_primary <- primary_locus_vec[index_vec] == expanded_df[["Locus_0MM"]]
  expanded_df[["Is_primary_location"]] <- are_primary

  expanded_df[["Intended_Entrez_ID"]] <- CRISPR_df[["Entrez_ID"]][index_vec]
  expanded_df[["Intended_gene_symbol"]] <- CRISPR_df[["Gene_symbol"]][index_vec]
  expanded_df[["Num_loci"]] <- CRISPR_df[["Num_0MM"]][index_vec]

  are_not_NA <- !(is.na(expanded_df[["Locus_0MM"]]))
  unique_loci_vec <- unique(expanded_df[["Locus_0MM"]][are_not_NA])
  unique_loci_df <- data.frame("Locus_0MM" = unique_loci_vec,
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



StandardLocationString <- function(ranges_df) {
  paste0(ranges_df[, "Chromosome"],
         ":",
         ranges_df[, "Start"],
         "-",
         ranges_df[, "End"],
         ":",
         ranges_df[, "Strand"]
         )
}


FindOverlappingHits <- function(loci_df, genes_df) {

  sgRNA_GPos_object <- LocationsDfToGPosObject(loci_df)
  TSS_GRanges_object <- RangesDfToGRangesObject(genes_df)

  hits_object <- findOverlaps(sgRNA_GPos_object,
                              TSS_GRanges_object,
                              ignore.strand = TRUE,
                              select = "all"
                              )

  loci_df[["Guide_locus"]] <- StandardLocationString(loci_df)
  genes_df[["Gene_locus"]] <- StandardLocationString(genes_df)

  combined_df <- data.frame(
    loci_df[queryHits(hits_object), c("Locus_0MM", "Guide_locus")],
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

  were_found <- expanded_df[["Locus_0MM"]] %in% comb_df[["Locus_0MM"]]
  were_found_indices <- unique(expanded_df[["Index"]][were_found])

  are_to_retain <- (were_found) |
                   (!(duplicated(expanded_df[["Index"]])) &
                    !(expanded_df[["Index"]] %in% were_found_indices)
                    )

  indices_list <- split(
    seq_len(nrow(comb_df)),
    factor(comb_df[["Locus_0MM"]],
           levels = unique(expanded_df[["Locus_0MM"]][are_to_retain]),
           exclude = c()
           )
  )

  indices_list[lengths(indices_list) == 0] <- list(NA_integer_)

  matches_vec <- match(expanded_df[["Locus_0MM"]][are_to_retain], names(indices_list))

  num_reps_vec <- lengths(indices_list)[matches_vec]

  expanded_indices <- rep(which(are_to_retain), times = num_reps_vec)

  combined_indices <- unlist(indices_list[matches_vec], use.names = FALSE)

  full_df <- data.frame(expanded_df[expanded_indices, ],
                        comb_df[combined_indices, names(comb_df) != "Locus_0MM"],
                        stringsAsFactors = FALSE,
                        row.names = NULL
                        )
  return(full_df)
}




SplitAndRejoin <- function(long_df, use_column) {

  splits_list <- strsplit(long_df[, use_column], ", ", fixed = TRUE)

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






SummarizeFullDf <- function(full_df) {

  unique_entrezs <- unique(full_df[["Affected_Entrez_IDs"]])
  min_unique_entrezs <- GetMinEntrez(unique_entrezs)
  min_entrezs_vec <- min_unique_entrezs[match(full_df[["Affected_Entrez_IDs"]], unique_entrezs)]

  new_order <- order(
    full_df[["Index"]],
    full_df[["Intended_Entrez_ID"]] != full_df[["Affected_Entrez_IDs"]],
    min_entrezs_vec,
    match(full_df[["Locus_0MM"]], full_df[["Locus_0MM"]])
  )

  reordered_by_gene_df <- full_df
  reordered_by_gene_df[["Original_order"]] <- seq_len(nrow(reordered_by_gene_df))
  reordered_by_gene_df <- full_df[new_order, ]
  row.names(reordered_by_gene_df) <- NULL

  affected_symbols_results <- SplitAndRejoin(reordered_by_gene_df, "Affected_gene_symbols")
  affected_entrezs_results <- SplitAndRejoin(reordered_by_gene_df, "Affected_Entrez_IDs")

  stopifnot(identical(lengths(affected_symbols_results[["list"]]),
                      lengths(affected_entrezs_results[["list"]])
                      ))

  first_indices <- match(unique(full_df[["Index"]]),
                         full_df[["Index"]]
                         )

  intended_entrezs <- full_df[["Intended_Entrez_ID"]][first_indices]
  intended_symbols <- full_df[["Intended_gene_symbol"]][first_indices]
  total_loci       <- full_df[["Num_loci"]][first_indices]

  summary_df <- data.frame(
    "Intended_Entrez_ID"    = intended_entrezs,
    "Affected_Entrez_IDs"   = affected_entrezs_results[["vector"]],
    "Intended_gene_symbol"  = intended_symbols,
    "Affected_gene_symbols" = affected_symbols_results[["vector"]],
    "Num_affected_genes"    = lengths(affected_entrezs_results[["list"]]),
    "Total_loci"            = total_loci,
    stringsAsFactors        = FALSE
  )

  index_loci_vec <- paste0(full_df[["Index"]], "__", full_df[["Locus_0MM"]])

  are_duplicated <- duplicated(index_loci_vec)

  locus_fac <- factor(index_loci_vec, levels = index_loci_vec[!(are_duplicated)])

  are_NA_entrezs <- is.na(full_df[["Affected_Entrez_IDs"]])

  entrezs_locus_splits <- split(full_df[["Affected_Entrez_IDs"]][!(are_NA_entrezs)],
                                locus_fac[!(are_NA_entrezs)]
                                )

  index_vec <- full_df[["Index"]][!(duplicated(index_loci_vec))]

  second_splits <- split(entrezs_locus_splits, index_vec)
  second_splits <- lapply(second_splits, function(x) x[lengths(x) > 0])

  have_targets <- tapply(!(are_NA_entrezs), locus_fac, any)
  summary_df[["Loci_with_targets"]] <- as.vector(tapply(have_targets, index_vec, sum))
  stopifnot(identical(unname(lengths(second_splits)), summary_df[["Loci_with_targets"]]))


  summary_df[["Loci_targeting_intended_gene"]] <- vapply(seq_len(nrow(summary_df)),
                                                         function(x) {
                                                           intended_entrez <- summary_df[["Intended_Entrez_ID"]][[x]]
                                                           sum(vapply(second_splits[[x]], function(y) intended_entrez %in% y, logical(1)))
                                                         }, integer(1))

  distinct_splits <- lapply(second_splits, unique)

  have_multiple_loci <- summary_df[["Loci_with_targets"]] > 1

  distinct_splits[have_multiple_loci] <- lapply(distinct_splits[have_multiple_loci],
                                                MergeCommonElements
                                                )

  split_distinct_splits <- lapply(distinct_splits,
                                  function(x) lapply(x, function(y) unlist(strsplit(y, ", ", fixed = TRUE)))
                                  )


  summary_df[["Distinct_loci"]] <- lengths(distinct_splits, use.names = FALSE)

  summary_df[["Distinct_loci_with_multiple_targets"]] <- vapply(split_distinct_splits,
                                                                function(x) sum(lengths(x) > 1),
                                                                integer(1)
                                                                )

  return(summary_df)
}





SummarizeSummaryDf <- function(summary_df) {

  summary_df[["Affects_intended_gene"]]       <- summary_df[["Loci_targeting_intended_gene"]] >= 1
  summary_df[["Affects_genes_at_other_loci"]] <- summary_df[["Distinct_loci"]] >= 2

  summary_df[["Affects_unintended_gene"]]     <- (summary_df[["Num_affected_genes"]] >= 1) |
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

    "No hits in the reference genome"             = (summary_df[["Total_loci"]] == 0),

    "No targets (single locus)"                   = (summary_df[["Total_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] == 0),

    "Only unintended (single locus)"              = (summary_df[["Total_loci"]] == 1) &
                                                    (summary_df[["Num_affected_genes"]] >= 1) &
                                                    (summary_df[["Loci_targeting_intended_gene"]] == 0),

    "No targets (multiple loci)"                  = (summary_df[["Total_loci"]] >= 2) &
                                                    (summary_df[["Num_affected_genes"]] == 0),

    "Only unintended (multiple loci)"             = (summary_df[["Total_loci"]] >= 2) &
                                                    (summary_df[["Num_affected_genes"]] >= 1) &
                                                    (summary_df[["Loci_targeting_intended_gene"]] == 0)
  )

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




