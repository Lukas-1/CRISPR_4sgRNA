### 7th August 2019 ####





# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("motifRG")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R")) # For ExpandList
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R"))






# Load global variables ---------------------------------------------------

human_genes_GRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)






# General helper functions ------------------------------------------------

RangesDfToGRangesObject <- function(ranges_df) {
  CheckRangesDf(ranges_df)
  GRanges_object <- GRanges(
    seqnames = ranges_df[, "Chromosome"],
    ranges   = IRanges(start = ranges_df[, "Start"], end = ranges_df[, "End"]),
    strand   = ranges_df[, "Strand"]
  )
  return(GRanges_object)
}





# Functions for retrieving sequences from locations -----------------------

RetrieveSequences <- function(ranges_df) {
  GRanges_object <- RangesDfToGRangesObject(ranges_df)
  results_vec <- as.character(motifRG::getSequence(GRanges_object, BSgenome.Hsapiens.UCSC.hg38))
  return(results_vec)
}


GetNGGPAM <- function(ranges_df) {

  CheckRangesDf(ranges_df)

  start_vec <- ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] + 1L, ranges_df[, "Start"] - 3L)
  end_vec   <- ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] + 3L, ranges_df[, "Start"] - 1L)

  assign("delete_ranges_df", ranges_df, envir = globalenv())
  GRanges_object <- GRanges(
    seqnames = ranges_df[, "Chromosome"],
    ranges   = IRanges(start = start_vec, end = end_vec),
    strand   = ranges_df[, "Strand"]
  )
  results_vec <- as.character(motifRG::getSequence(GRanges_object, BSgenome.Hsapiens.UCSC.hg38))
  return(results_vec)
}






# Functions for mapping from genomic locations to genes -------------------

ProcessHitsObject <- function(Hits_object, num_queries, gene_models_GRanges) {

  entrez_IDs_vec <- mcols(gene_models_GRanges)[, "gene_id"]
  hits_df <- as.data.frame(Hits_object)

  hits_df[, "Entrez_ID"] <- entrez_IDs_vec[hits_df[, "subjectHits"]]
  integer_entrezs <- as.integer(hits_df[, "Entrez_ID"])
  new_order <- order(hits_df[, "queryHits"], integer_entrezs)
  hits_df <- hits_df[new_order, ]

  entrez_matches <- match(hits_df[, "Entrez_ID"], entrez_to_symbol_df[, "Entrez_ID"])
  hits_df[, "Gene_symbol"] <- ifelse(is.na(entrez_to_symbol_df[entrez_matches, "Symbol_Org_Hs_eg_db"]),
                                     entrez_to_symbol_df[entrez_matches, "Symbol_NCBI_Hs_info"],
                                     entrez_to_symbol_df[entrez_matches, "Symbol_Org_Hs_eg_db"]
                                     )

  query_fac <- factor(hits_df[, "queryHits"], unique(hits_df[, "queryHits"]))

  entrezs_list <- split(hits_df[, "Entrez_ID"], query_fac)
  entrezs_vec <- vapply(entrezs_list, function(x) paste0(x, collapse = ", "), "")

  symbols_list <- split(hits_df[, "Gene_symbol"], query_fac)
  symbols_vec <- vapply(symbols_list, function(x) paste0(x, collapse = ", "), "")

  results_df <- data.frame(
    "Entrez_ID"   = entrezs_vec,
    "Gene_symbol" = symbols_vec,
    "Num_genes"   = lengths(entrezs_list),
    stringsAsFactors = FALSE
  )

  if ("distance" %in% colnames(hits_df)) {
    results_df[, "Distance"] <- hits_df[unique(match(hits_df[, "queryHits"], hits_df[, "queryHits"])), "distance"]
  }

  match_matches_vec <- match(seq_len(num_queries), as.integer(names(entrezs_list)))
  results_df <- results_df[match_matches_vec, ]
  rownames(results_df) <- NULL
  return(results_df)
}



FindOverlappingGenes <- function(ranges_df, gene_models_GRanges = human_genes_GRanges, ignore_strand = TRUE) {
  ### This function requires the 'human_genes_GRanges' object in the global environment ###

  message("Finding genes that overlap with the specified genomic ranges...")

  GRanges_object <- RangesDfToGRangesObject(ranges_df)

  hits_object <- findOverlaps(GRanges_object, gene_models_GRanges, ignore.strand = ignore_strand, select = "all")
  hits_df <- ProcessHitsObject(hits_object, nrow(ranges_df), gene_models_GRanges)

  results_df <- data.frame(
    ranges_df,
    "Overlapping_Entrez_IDs" = hits_df[, "Entrez_ID"],
    "Overlapping_symbols"    = hits_df[, "Gene_symbol"],
    "Num_overlapping_genes"  = hits_df[, "Num_genes"],
    stringsAsFactors         = FALSE,
    row.names                = NULL
  )
  return(results_df)
}



FindNearestGenes <- function(ranges_df) {
  ### This function requires the 'human_genes_GRanges' object in the global environment ###

  message("Finding the closest genes to the specified genomic ranges...")

  GRanges_object <- RangesDfToGRangesObject(ranges_df)

  hits_object <- distanceToNearest(GRanges_object, human_genes_GRanges, ignore.strand = TRUE, select = "all")
  hits_df <- ProcessHitsObject(hits_object, nrow(ranges_df), human_genes_GRanges)

  results_df <- data.frame(
    ranges_df,
    "Nearest_Entrez_IDs" = hits_df[, "Entrez_ID"],
    "Nearest_symbols"    = hits_df[, "Gene_symbol"],
    "Distance"           = hits_df[, "Distance"],
    "Num_nearest"        = hits_df[, "Num_genes"],
    stringsAsFactors     = FALSE,
    row.names            = NULL
  )
  return(results_df)
}







# Functions for summarizing mapped sequences ------------------------------

ReturnClearMatches <- function(found_seq_df, have_clear_match, clear_match_indices_vec) {
  columns_list <- lapply(seq_len(ncol(found_seq_df)), function(x) {
    if (typeof(found_seq_df[, x]) == "integer") {
      NA_vec <- rep.int(NA_integer_, nrow(found_seq_df))
    } else if (typeof(found_seq_df[, x]) == "double") {
      NA_vec <- rep.int(NA_real_, nrow(found_seq_df))
    } else if (typeof(found_seq_df[, x]) == "character") {
      NA_vec <- rep.int(NA_character_, nrow(found_seq_df))
    }
    results_vec <- ifelse(have_clear_match, found_seq_df[clear_match_indices_vec, x], NA_vec)
    return(results_vec)
  })
  results_df <- do.call(data.frame, c(columns_list, list(stringsAsFactors = FALSE)))
  colnames(results_df) <- colnames(found_seq_df)
  return(results_df)
}



MakeLocationStrings <- function(ranges_df) {
  CheckRangesDf(ranges_df)
  results_vec <- ifelse(is.na(ranges_df[, "Start"]),
                        NA_character_,
                        paste0(ranges_df[, "Chromosome"], "(", ranges_df[, "Strand"], "):", ranges_df[, "Start"], "-", ranges_df[, "End"])
                        )
  return(results_vec)
}



PasteIndices <- function(found_seq_df, indices_list, column_name, use_separator) {
  vapply(indices_list, function(x) {
    if (length(x) == 0) {
      return(NA_character_)
    } else {
      use_vec <- found_seq_df[x, column_name]
      if (all(is.na(use_vec))) {
        return(NA_character_)
      } else {
        return(paste0(found_seq_df[x, column_name], collapse = use_separator))
      }
    }
  }, "")
}


SummarizeFoundSequencesDf <- function(found_seq_df, all_sequences = NULL, use_separator = "; ") {

  are_0MM       <- found_seq_df[, "Num_MM"] == 0
  # are_old_MM5primeG <- !(are_0MM) & (substr(found_seq_df[, "Sequence"], 2, nchar(found_seq_df[, "Sequence"])) == substr(found_seq_df[, "Reference"], 2, nchar(found_seq_df[, "Reference"])))
  are_MM5primeG <- !(are_0MM) & (substr(found_seq_df[, "Reference"], 1, 1) == "G") &
                   (substr(found_seq_df[, "Sequence"], 2, nchar(found_seq_df[, "Sequence"])) == substr(found_seq_df[, "Reference"], 2, nchar(found_seq_df[, "Reference"])))
  are_1MM       <- !(are_0MM | are_MM5primeG)

  are_potential_PAM <- substr(found_seq_df[, "PAM"], 2, 3) %in% c("GG", "GA", "AG")

  positions_vec <- MakeLocationStrings(found_seq_df)

  references_fac <- factor(found_seq_df[, "Reference"], levels = unique(found_seq_df[, "Reference"]))
  indices_list <- split(seq_len(nrow(found_seq_df)), references_fac)

  which_0MM_list <- tapply(are_0MM & are_potential_PAM, references_fac, which, simplify = FALSE)
  which_0MM_indices_list <- lapply(seq_len(nlevels(references_fac)), function(x) indices_list[[x]][which_0MM_list[[x]]])
  locations_0MM_vec <- vapply(which_0MM_indices_list, function(x) if (length(x) == 0) NA_character_ else paste0(positions_vec[x], collapse = use_separator), "")

  which_MM5primeG_list <- tapply(are_MM5primeG & are_potential_PAM, references_fac, which, simplify = FALSE)
  which_MM5primeG_indices_list <- lapply(seq_len(nlevels(references_fac)), function(x) indices_list[[x]][which_MM5primeG_list[[x]]])

  which_1MM_list <- tapply((are_1MM | are_MM5primeG) & are_potential_PAM, references_fac, which, simplify = FALSE)
  which_1MM_indices_list <- lapply(seq_len(nlevels(references_fac)), function(x) indices_list[[x]][which_1MM_list[[x]]])
  locations_1MM_vec <- vapply(which_1MM_indices_list, function(x) if (length(x) == 0) NA_character_ else paste0(positions_vec[x],            collapse = use_separator), "")
  sequences_1MM_vec <- vapply(which_1MM_indices_list, function(x) if (length(x) == 0) NA_character_ else paste0(found_seq_df[x, "Sequence"], collapse = use_separator), "")

  sum_0MM_and_PAM_vec   <- tapply(are_0MM       & are_potential_PAM, references_fac, sum)
  sum_1MM_and_PAM_vec   <- tapply(are_1MM       & are_potential_PAM, references_fac, sum)
  sum_5G_MM_and_PAM_vec <- tapply(are_MM5primeG & are_potential_PAM, references_fac, sum)

  sum_0MM_vec   <- tapply(are_0MM, references_fac, sum)
  sum_5G_MM_vec <- tapply(are_MM5primeG, references_fac, sum)
  have_clear_match <- (sum_0MM_and_PAM_vec == 1) | ((sum_0MM_and_PAM_vec == 0) & (sum_5G_MM_and_PAM_vec == 1))

  use_overlapping_genes <- all(c("Overlapping_Entrez_IDs", "Overlapping_symbols") %in% colnames(found_seq_df))
  if (use_overlapping_genes) {
    entrez_column <- "Overlapping_Entrez_IDs"
    symbol_column <- "Overlapping_symbols"
  } else {
    entrez_column <- "Nearest_Entrez_IDs"
    symbol_column <- "Nearest_symbols"
  }

  entrezs_0MM_vec <- PasteIndices(found_seq_df, which_0MM_indices_list, entrez_column, use_separator)
  symbols_0MM_vec <- PasteIndices(found_seq_df, which_0MM_indices_list, symbol_column, use_separator)
  entrezs_1MM_vec <- PasteIndices(found_seq_df, which_1MM_indices_list, entrez_column, use_separator)
  symbols_1MM_vec <- PasteIndices(found_seq_df, which_1MM_indices_list, symbol_column, use_separator)

  # The following vector is relevant mainly for duplicated genes, which or may not all have the same PAM
  PAM_vec <- vapply(which_0MM_indices_list, function(x) {
    if (length(x) == 0) {
      return(NA_character_)
    } else {
      unique_PAMs <- unique(found_seq_df[x, "PAM"])
      if (length(unique_PAMs) > 1) {
        return(NA_character_)
      } else {
        return(unique_PAMs)
      }
    }
  }, "")
  PAM_0MM_vec <- PasteIndices(found_seq_df, which_0MM_indices_list, "PAM", use_separator)
  PAM_1MM_vec <- PasteIndices(found_seq_df, which_1MM_indices_list, "PAM", use_separator)

  clear_match_indices_vec <- vapply(seq_len(nlevels(references_fac)), function(x) {
    if (sum_0MM_and_PAM_vec[[x]] == 1) {
      which_0MM_indices_list[[x]]
    } else if (have_clear_match[[x]]) {
      which_MM5primeG_indices_list[[x]]
    } else {
      NA_integer_
    }
  }, integer(1))


  SNP_column_names <- c(
    "SNP_IDs_vcf",     "SNP_AFs_1kGenomes", "SNP_AF_max_1kGenomes", "SNP_AF_sum_1kGenomes",
    "SNP_AFs_TOPMED",  "SNP_AF_sum_TOPMED", "SNP_AF_max_TOPMED",
    "SNP_AFs_Kaviar",  "SNP_AF_sum_Kaviar", "SNP_AF_max_Kaviar",
    "SNP_IDs_1kG_ph1", "SNP_AFs_1kG_ph1",   "SNP_AF_max_1kG_ph1",   "SNP_AF_sum_1kG_ph1",
    "SNP_IDs_1kG_ph3", "SNP_AFs_1kG_ph3",   "SNP_AF_max_1kG_ph3",   "SNP_AF_sum_1kG_ph3",
    "SNP_IDs_gnomAD",  "SNP_AFs_gnomAD",    "SNP_AF_max_gnomAD",    "SNP_AF_sum_gnomAD"
  )
  SNP_column_names <- c(paste0("sgRNA_", SNP_column_names), paste0("PAM_", SNP_column_names), paste0("all23_", SNP_column_names))
  SNP_column_names <- SNP_column_names[SNP_column_names %in% colnames(found_seq_df)]

  all_column_names <- c("Chromosome", "Strand", "Start", "End", SNP_column_names)#, "PAM", )

  assign("delete_all_column_names", all_column_names, envir = globalenv())
  assign("delete_found_seq_df", found_seq_df, envir = globalenv())

  results_df <- data.frame(
    "Sequence"           = levels(references_fac),
    "PAM"                = PAM_vec,
    ReturnClearMatches(found_seq_df[, all_column_names], have_clear_match, clear_match_indices_vec),
    "Num_0MM"            = sum_0MM_and_PAM_vec,
    "Num_5G_MM"          = sum_5G_MM_and_PAM_vec,
    "Num_1MM"            = sum_1MM_and_PAM_vec,
    "Num_0MM_no_PAM"     = sum_0MM_vec - sum_0MM_and_PAM_vec,
    "Num_5G_MM_no_PAM"   = sum_5G_MM_vec - sum_5G_MM_and_PAM_vec,
    "Num_1MM_no_PAM"     = sum_5G_MM_vec - sum_1MM_and_PAM_vec,
    "PAM_0MM"            = PAM_0MM_vec,
    "Locations_0MM"      = locations_0MM_vec,
    "Entrez_nearest_0MM" = entrezs_0MM_vec,
    "Symbol_nearest_0MM" = symbols_0MM_vec,
    "PAM_1MM"            = PAM_1MM_vec,
    "Locations_1MM"      = locations_1MM_vec,
    "Sequences_1MM"      = sequences_1MM_vec,
    "Entrez_nearest_1MM" = entrezs_1MM_vec,
    "Symbol_nearest_1MM" = symbols_1MM_vec,
    stringsAsFactors     = FALSE,
    row.names            = NULL
  )

  if (use_overlapping_genes) {
    old_column_names <- c("Entrez_nearest_0MM", "Symbol_nearest_0MM", "Entrez_nearest_1MM", "Symbol_nearest_1MM")
    new_column_names <- sub("_nearest_", "_overlapping_", old_column_names, fixed = TRUE)
    colnames(results_df)[match(old_column_names, colnames(results_df))] <- new_column_names
  }

  if (!(is.null(all_sequences))) {
    my_matches <- match(toupper(all_sequences), toupper(results_df[, "Sequence"]))
    results_df <- results_df[my_matches, ]
    results_df[is.na(my_matches), "Sequence"] <- all_sequences[is.na(my_matches)]
    rownames(results_df) <- NULL
    for (my_column in c("Num_0MM", "Num_5G_MM", "Num_1MM")) {
      results_df[, my_column] <- ifelse(is.na(my_matches), 0L, results_df[, my_column])
    }
    results_df[, "Found"] <- ifelse(is.na(my_matches), "No", "Yes")
  }

  return(results_df)
}





# Other functions ---------------------------------------------------------

LiftOverAndAnnotate <- function(ranges_df) {
  # Depends on the object 'hg19tohg38_chain' in the global environment

  stopifnot(all((ranges_df[, "End"] - ranges_df[, "Start"]) == 19))

  GRanges_object_sgRNAs <- RangesDfToGRangesObject(ranges_df)

  sequences_vec <- as.character(motifRG::getSequence(GRanges_object_sgRNAs, BSgenome.Hsapiens.UCSC.hg19))

  sgRNAs_lifted_over_GRanges <- liftOver(GRanges_object_sgRNAs, hg19tohg38_chain)

  sgRNAs_lifted_over_df <- as.data.frame(sgRNAs_lifted_over_GRanges)

  location_columns <- c("seqnames" = "Chromosome",
                        "strand"   = "Strand",
                        "start"    = "Start",
                        "end"      = "End"
                        )
  for (column_name in names(location_columns)) {
    colnames(sgRNAs_lifted_over_df)[colnames(sgRNAs_lifted_over_df) == column_name] <- location_columns[[column_name]]
  }


  sgRNAs_lifted_over_df <- sgRNAs_lifted_over_df[sgRNAs_lifted_over_df[, "width"] == 20 , ]

  if (any(duplicated(sgRNAs_lifted_over_df[, "group"]))) {
    stop("Some sgRNAs mapped to multiple locations in the liftOver!")
  }

  sgRNAs_lifted_over_GRanges_20mers <- RangesDfToGRangesObject(sgRNAs_lifted_over_df)

  sgRNAs_lifted_over_df[, "Sequence_liftOver"] <- as.character(motifRG::getSequence(sgRNAs_lifted_over_GRanges_20mers, BSgenome.Hsapiens.UCSC.hg38))

  sgRNAs_lifted_over_df[, "PAM_liftOver"] <- GetNGGPAM(sgRNAs_lifted_over_df)

  lifted_over_matches <- match(seq_len(nrow(ranges_df)), sgRNAs_lifted_over_df[, "group"])

  lifted_over_matched_df <- sgRNAs_lifted_over_df[lifted_over_matches, c("Sequence_liftOver", "PAM_liftOver", location_columns)]
  colnames(lifted_over_matched_df)[3:6] <- paste0(location_columns, "_liftOver")

  results_df <- data.frame(
    "Sequence_hg19"          = sequences_vec,
    lifted_over_matched_df,
    stringsAsFactors         = FALSE,
    row.names                = NULL
  )
  return(results_df)
}



ReassignEntrezsByLocations <- function(confirmed_locations_df) {

  are_ambiguous <- grepl(",", confirmed_locations_df[, "Entrez_ID"], fixed = TRUE)

  ambiguous_df <- confirmed_locations_df[are_ambiguous, ]

  new_entrezs_vec <- rep.int(NA_character_, nrow(ambiguous_df))
  assignment_vec <- rep.int(NA_character_, nrow(ambiguous_df))

  for (combined_ID in unique(ambiguous_df[, "Combined_ID"])) {
    are_this_ID <- ambiguous_df[, "Combined_ID"] == combined_ID
    entrezs_vec <- strsplit(ambiguous_df[which(are_this_ID)[[1]], "Entrez_ID"], ", ", fixed = TRUE)[[1]]
    chromosomes_vec <- entrez_to_symbol_df[match(entrezs_vec, entrez_to_symbol_df[, "Entrez_ID"]), "Chromosome"]
    correct_chromosome <- unique(ambiguous_df[are_this_ID, "Chromosome"])
    stopifnot(length(correct_chromosome) == 1)
    if (!(any(duplicated(chromosomes_vec))) && length(correct_chromosome == 1)) {
      assign("delete_entrezs_vec", entrezs_vec, envir = globalenv())
      assign("delete_correct_chromosome", correct_chromosome, envir = globalenv())
      assign("delete_chromosomes_vec", chromosomes_vec, envir = globalenv())
      new_entrezs_vec[are_this_ID] <- entrezs_vec[[match(correct_chromosome, chromosomes_vec)]]
      assignment_vec[are_this_ID] <- "Unambiguous chromosome"
    } else {
      overlapping_entrezs <- unique(unlist(strsplit(ambiguous_df[, "Overlapping_Entrez_IDs"], ", ", fixed = TRUE), use.names = FALSE))
      if (all(is.na(overlapping_entrezs))) {
        assignment_vec[are_this_ID] <- "Ambiguous chromosome, and no overlaps with genes"
      } else {
        overlapping_entrezs <- overlapping_entrezs[!(is.na(overlapping_entrezs))]
        entrezs_intersect <- intersect(entrezs_vec, overlapping_entrezs)
        if (length(entrezs_intersect) == 1) {
          new_entrezs_vec[are_this_ID] <- entrezs_intersect
          assignment_vec[are_this_ID] <- "Overlaps with gene"
        } else {
          assignment_vec[are_this_ID] <- "Ambiguous chromosome, and ambiguous overlaps"
        }
      }
    }
  }

  results_df <- confirmed_locations_df

  results_df[, "Old_Entrez"] <- results_df[, "Entrez_ID"]
  results_df[, "Old_symbol"] <- results_df[, "Gene_symbol"]

  # Replace the old ambiguous Entrez IDs with the new ones
  were_replaced <- !(is.na(new_entrezs_vec))
  results_df[are_ambiguous, "Entrez_ID"][were_replaced] <- new_entrezs_vec[were_replaced]

  results_df[, "Were_replaced"] <- NA
  results_df[are_ambiguous, "Were_replaced"] <- were_replaced

  results_df[, "Entrez_assignment"] <- NA_character_
  results_df[are_ambiguous, "Entrez_assignment"] <- assignment_vec

  results_df[are_ambiguous, "Combined_ID"] <- ifelse(is.na(results_df[are_ambiguous, "Entrez_ID"]),
                                                     toupper(results_df[are_ambiguous, "Original_symbol"]),
                                                     results_df[are_ambiguous, "Entrez_ID"]
                                                     )

  results_df[are_ambiguous, "Gene_symbol"][were_replaced] <- MapToEntrezs(entrez_IDs_vec = new_entrezs_vec[were_replaced])[, "Gene_symbol"]

  return(results_df)
}


