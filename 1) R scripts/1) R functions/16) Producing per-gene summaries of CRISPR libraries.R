### 17th October 2019 ###





# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "09) Constants and settings.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))






# Define functions --------------------------------------------------------

CollapseOriginal <- function(char_vec) {
  paste0(unique(char_vec[char_vec != ""]), collapse = "")
}


MeetCriteria <- function(CRISPR_df, allow_curated = FALSE) {
  are_to_exclude <- ((CRISPR_df[, preferred_AF_max_column] > SNP_frequency_cutoff) %in% TRUE) |
                    (is.na(CRISPR_df[, "GuideScan_specificity"]) & is.na(CRISPR_df[, "CRISPOR_4MM_specificity"])) |
                    grepl("TTTT", CRISPR_df[, "sgRNA_sequence"], ignore.case = TRUE) |
                    !(((substr(CRISPR_df[, "PAM"], 2, 3) == "GG") %in% TRUE)) |
                    (((CRISPR_df[, "GuideScan_specificity"] < 0.2) %in% TRUE) |
                      (is.na(CRISPR_df[, "GuideScan_specificity"]) & ((CRISPR_df[, "CRISPOR_3MM_specificity"] < 0.2) %in% TRUE))) |
                    (("Exon_number_GPP" %in% colnames(CRISPR_df)) & (CRISPR_df[, "CRISPOR_Graf_status"] %in% c("ggc", "tt")))
  if (!(allow_curated)) {
    are_to_exclude <- are_to_exclude | (CRISPR_df[, "Source"] == "Curated")
  }
  results_vec <- !(are_to_exclude)
  return(results_vec)
}


AggregateSpecificityScores <- function(scores_vec) {
  1 / (1 + sum((1 / scores_vec) - 1))
}


SummarizeCRISPRDf <- function(CRISPR_df) {
  split_indices <- split(seq_len(nrow(CRISPR_df)), factor(CRISPR_df[, "Combined_ID"], levels = unique(CRISPR_df[, "Combined_ID"])))
  include_top4nonoverlapping <- "Best_combination_rank" %in% colnames(CRISPR_df)
  if (include_top4nonoverlapping) {
    CRISPR_df[, "Non_overlapping"] <- !(is.na(CRISPR_df[, "Best_combination_rank"]))
  } else {
    CRISPR_df[, "Non_overlapping"] <- NA
  }
  CRISPR_df[, "Overlap_with_SNP"] <- (CRISPR_df[, preferred_AF_max_column] > SNP_frequency_cutoff) %in% TRUE
  CRISPR_df[, "Meet_criteria"] <- MeetCriteria(CRISPR_df)
  if (include_top4nonoverlapping) {
    CRISPR_df[, "Meet_criteria"] <- CRISPR_df[, "Meet_criteria"] & CRISPR_df[, "Non_overlapping"]
  }
  CRISPR_df[, "Are_unspecific"] <- ((CRISPR_df[, "Num_0MM"] > 1) | (CRISPR_df[, "Num_1MM"] > 0)) %in% TRUE
  include_TSS                   <- "TSS_regions" %in% colnames(CRISPR_df)
  include_searched_by_GuideScan <- "TSS_searched_by_GuideScan" %in% colnames(CRISPR_df)
  include_transcripts           <- "TSS_ID" %in% colnames(CRISPR_df)
  if (include_transcripts && include_top4nonoverlapping) {
    are_incomplete <- (CRISPR_df[, "Spacing"] %in% 0) & (CRISPR_df[, "Rank"] %in% 1:4)
  }
  results_list_list <- lapply(split_indices, function(x) {
    results_list <- list(
      "Combined_ID"                 = unique(CRISPR_df[x, "Combined_ID"]),
      "Entrez_ID"                   = unique(CRISPR_df[x, "Entrez_ID"]),
      "Gene_symbol"                 = unique(CRISPR_df[x, "Gene_symbol"]),
      "Original_entrez"             = CollapseOriginal(CRISPR_df[x, "Original_entrez"]),
      "Original_symbol"             = CollapseOriginal(CRISPR_df[x, "Original_symbol"]),
      "Num_hCRISPRa_v2_transcripts" = NA_integer_,
      "Num_transcripts"             = NA_integer_,
      "Num_overlapping_transcripts" = NA_integer_,
      "Num_unspaced_transcripts"    = NA_integer_,
      "Num_incomplete_transcripts"  = NA_integer_,
      "Spacing"                     = NA_character_,
      "GuideScan_specificity"       = NA_real_,
      "CRISPOR_3MM_specificity"     = NA_real_,
      "CRISPOR_4MM_specificity"     = NA_real_,
      "Longest_subsequence"         = NA_integer_,
      "Num_total"                   = length(x),
      "Num_overlaps"                = NULL,
      "Num_meeting_criteria"        = sum(CRISPR_df[x, "Meet_criteria"]),
      "Num_without_GuideScan"       = sum(is.na(CRISPR_df[x, "GuideScan_specificity"])),
      "Num_unspecific"              = sum(CRISPR_df[x, "Are_unspecific"]),
      "Num_overlapping_with_SNP"    = sum(CRISPR_df[x, "Overlap_with_SNP"]),
      "Submitted_to_GuideScan"      = NA_character_,
      "TSS_regions"                 = NA_character_
      )
    if (include_searched_by_GuideScan) {
      were_searched <- CRISPR_df[x, "TSS_searched_by_GuideScan"]
      if (all(is.na(were_searched))) {
        submitted <- "Not mapped"
      } else if (all(were_searched[!(is.na(were_searched))] %in% c("No", "Not this gene"))) {
        submitted <- "No"
      } else {
        submitted <- "Yes"
      }
    } else {
      submitted <- NULL
    }
    results_list[["Submitted_to_GuideScan"]] <- submitted
    if (include_TSS) {
      results_list[["TSS_regions"]] <- unique(CRISPR_df[x, "TSS_regions"])
    } else {
      results_list[["TSS_regions"]] <- NULL
    }
    if (include_transcripts) {
      transcripts_vec <- CRISPR_df[x, "TSS_ID"]
      if (all(is.na(transcripts_vec))) {
        transcripts_vec <- CRISPR_df[x, "AltTSS_ID"]
      }
      transcripts_fac <- factor(transcripts_vec, levels = unique(transcripts_vec))
      transcripts_indices <- split(seq_along(x), transcripts_fac)
      transcripts_top4_indices <- lapply(transcripts_indices, function(y) y[CRISPR_df[x, "Rank"][y] %in% 1:4])
    }
    if (include_top4nonoverlapping) {
      found_overlap <- FALSE
      for (space in c(50, 45, 40)) {
        are_this_spacing <- CRISPR_df[x, "Spacing"] %in% space
        if (any(are_this_spacing)) {
          num_zero_overlaps <- sum(CRISPR_df[x, "Num_overlaps"][are_this_spacing] %in% 0)
          overlap_string <- paste0(sum(num_zero_overlaps), "*", space)
          found_overlap <- TRUE
          break
        }
      }
      if (!(found_overlap)) {
        if (any(CRISPR_df[x, "Spacing"] %in% 12)) {
          overlap_string <- ">12bp"
        } else {
          overlap_string <- "None"
        }
      }
      results_list[["Spacing"]] <- overlap_string
      are_top_4 <- CRISPR_df[x, "Rank"] %in% 1:4
      results_list[["Num_overlaps"]] <- sum(CRISPR_df[x[are_top_4], "Num_overlaps"])
      if (include_transcripts) {
        guidescan_spec_vec <- vapply(transcripts_top4_indices,
                                     function(y) AggregateSpecificityScores(CRISPR_df[x, "GuideScan_specificity"][y]),
                                     numeric(1)
                                     )
        CRISPOR_3MM_spec_vec <- vapply(transcripts_top4_indices,
                                       function(y) AggregateSpecificityScores(CRISPR_df[x, "CRISPOR_3MM_specificity"][y]),
                                       numeric(1)
                                       )
        CRISPOR_4MM_spec_vec <- vapply(transcripts_top4_indices,
                                       function(y) AggregateSpecificityScores(CRISPR_df[x, "CRISPOR_4MM_specificity"][y]),
                                       numeric(1)
                                       )
        results_list[["GuideScan_specificity"]]   <- if (all(is.na(guidescan_spec_vec)))   NA_real_ else min(guidescan_spec_vec,   na.rm = TRUE)
        results_list[["CRISPOR_3MM_specificity"]] <- if (all(is.na(CRISPOR_3MM_spec_vec))) NA_real_ else min(CRISPOR_3MM_spec_vec, na.rm = TRUE)
        results_list[["CRISPOR_4MM_specificity"]] <- if (all(is.na(CRISPOR_4MM_spec_vec))) NA_real_ else min(CRISPOR_4MM_spec_vec, na.rm = TRUE)
        results_list[["Longest_subsequence"]] <- max(vapply(transcripts_top4_indices,
                                                            function(y) LongestSharedSubsequence(CRISPR_df[x, "sgRNA_sequence"][y]),
                                                            integer(1)
                                                            )
                                                     )
      } else {
        results_list[["GuideScan_specificity"]]   <- AggregateSpecificityScores(CRISPR_df[x[are_top_4], "GuideScan_specificity"])
        results_list[["CRISPOR_3MM_specificity"]] <- AggregateSpecificityScores(CRISPR_df[x[are_top_4], "CRISPOR_3MM_specificity"])
        results_list[["CRISPOR_4MM_specificity"]] <- AggregateSpecificityScores(CRISPR_df[x[are_top_4], "CRISPOR_4MM_specificity"])
        results_list[["Longest_subsequence"]]     <- LongestSharedSubsequence(CRISPR_df[x[are_top_4], "sgRNA_sequence"])
      }
    } else {
      results_list[["Spacing"]]                 <- NULL
      results_list[["Num_overlaps"]]            <- NULL
      results_list[["GuideScan_specificity"]]   <- NULL
      results_list[["CRISPOR_3MM_specificity"]] <- NULL
      results_list[["CRISPOR_4MM_specificity"]] <- NULL
      results_list[["Longest_subsequence"]]     <- NULL
    }
    if (include_transcripts) {
      transcripts_vec <- CRISPR_df[x, "TSS_ID"]
      if (all(is.na(transcripts_vec))) {
        transcripts_vec <- CRISPR_df[x, "AltTSS_ID"]
      }
      unique_transcripts <- unique(transcripts_vec)
      unique_transcripts <- unique_transcripts[!(is.na(unique_transcripts))]
      unique_hC_transcripts <- unique(CRISPR_df[x, "hCRISPRa_v2_transcript"])
      unique_hC_transcripts <- unique_hC_transcripts[!(is.na(unique_hC_transcripts))]
      results_list[["Num_transcripts"]] <- length(unique_transcripts)
      results_list[["Num_hCRISPRa_v2_transcripts"]] <- length(unique_hC_transcripts)
    } else {
      results_list[["Num_transcripts"]] <- NULL
      results_list[["Num_hCRISPRa_v2_transcripts"]] <- NULL
    }
    if (include_transcripts && include_top4nonoverlapping) {
      results_list[["Num_unspaced_transcripts"]]    <- sum(tapply(are_incomplete[x], transcripts_fac, any))
      results_list[["Num_incomplete_transcripts"]]  <- sum(!(lengths(transcripts_top4_indices) != 4))
      results_list[["Num_overlapping_transcripts"]] <- sum(!(tapply(CRISPR_df[x, "Num_overlaps"], transcripts_fac, function(x) sum(x %in% 0) == 4)))
    } else {
      results_list[["Num_unspaced_transcripts"]]    <- NULL
      results_list[["Num_incomplete_transcripts"]]  <- NULL
      results_list[["Num_overlapping_transcripts"]] <- NULL
    }
    return(results_list)
  })

  assign("delete_results_list_list", results_list_list, envir = globalenv())

  summary_df <- do.call(rbind.data.frame, c(results_list_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(summary_df)
}




ReorganizeSummaryDf <- function(summary_df, reference_IDs) {
  summary_mat <- as.matrix(summary_df[, grep("^Num_", colnames(summary_df), value = TRUE)])
  summary_matches <- match(reference_IDs, summary_df[, "Combined_ID"])
  matched_summary_mat <- summary_mat[summary_matches, ]
  # matched_summary_mat[is.na(matched_summary_mat)] <- 0L
  results_df <- data.frame("Combined_ID" = reference_IDs,
                           "Gene_present" = ifelse(is.na(summary_matches), "No", "Yes"),
                           matched_summary_mat,
                           summary_df[summary_matches, !(colnames(summary_df) %in% colnames(summary_mat))],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  column_index <- match("Original_symbol", colnames(summary_df))
  colnames_before <- colnames(summary_df)[seq_len(column_index)]
  colnames_after <- colnames(summary_df)[(column_index + 1):ncol(summary_df)]
  results_df <- results_df[, c(colnames_before, "Gene_present", colnames_after)]

  include_spacing <- "Spacing" %in% colnames(results_df)
  include_transcripts <- "TSS_ID" %in% colnames(results_df)
  if (include_spacing) {
    have_no_space <- (results_df[, "Spacing"] %in% c("None", ">12bp")) | is.na(results_df[, "Spacing"])
    spacing_splits <- strsplit(results_df[, "Spacing"], "*", fixed = TRUE)

    assign("delete_results_df", results_df, envir = globalenv())
    assign("delete_spacing", results_df[, "Spacing"], envir = globalenv())
    assign("delete_spacing_splits", spacing_splits, envir = globalenv())

    spacing_vec <- vapply(seq_along(have_no_space), function(x) if (have_no_space[[x]]) 0L else as.integer(spacing_splits[[x]][[2]]), integer(1))
    num_nonoverlapping_vec <- vapply(seq_along(have_no_space), function(x) if (have_no_space[[x]]) 0L else as.integer(spacing_splits[[x]][[1]]), integer(1))
  }

  NA_vec <- rep.int(NA, nrow(results_df))

  final_order <- order(
    results_df[, "Gene_present"],
    if (include_spacing) !(results_df[, "Spacing"] %in% "None")                                                               else NA_vec,
    if (include_spacing) !(results_df[, "Spacing"] %in% ">12bp")                                                              else NA_vec,
    if (include_spacing && include_transcripts) -(results_df[, "Num_incomplete_transcripts"])                                 else NA_vec,
    if (include_spacing && include_transcripts) (results_df[, "Num_transcripts"] != results_df[, "Num_unspaced_transcripts"]) else NA_vec,
    if (include_spacing && include_transcripts) -(results_df[, "Num_unspaced_transcripts"])                                   else NA_vec,
    if (include_spacing) num_nonoverlapping_vec                                                                               else NA_vec,
    if (include_spacing && include_transcripts) -(results_df[, "Num_overlapping_transcripts"])                                else NA_vec,
    if (include_spacing) spacing_vec                                                                                          else NA_vec,
    if (include_spacing) -(results_df[, "Num_overlaps"])                                                                      else NA_vec,
    ifelse(results_df[, "Num_total"] < 4, results_df[, "Num_total"], NA_integer_),
    ifelse(results_df[, "Num_meeting_criteria"] < 4, results_df[, "Num_meeting_criteria"], NA_integer_),
    if (include_spacing) !(is.na(results_df[, "GuideScan_specificity"]))                                                      else NA_vec,
    if (include_spacing) results_df[, "GuideScan_specificity"]                                                                else NA_vec,
    if (include_spacing) !(is.na(results_df[, "CRISPOR_4MM_specificity"]))                                                    else NA_vec,
    if (include_spacing) results_df[, "CRISPOR_4MM_specificity"]                                                              else NA_vec
  )
  results_df <- results_df[final_order, ]
  rownames(results_df) <- NULL
  return(results_df)
}




FixSymbolsForSummaryDf <- function(reorganized_df) {
  have_no_symbol <- is.na(reorganized_df[, "Gene_symbol"])
  entrez_symbols_df <- MapToEntrezs(reorganized_df[have_no_symbol, "Entrez_ID"])
  could_be_mapped <- !(is.na(entrez_symbols_df[, "Gene_symbol"]))
  results_df <- reorganized_df
  results_df[have_no_symbol, ][could_be_mapped, "Gene_symbol"] <- entrez_symbols_df[could_be_mapped, "Gene_symbol"]
  results_df[is.na(results_df[, "Original_entrez"]), "Original_entrez"] <- ""
  results_df[is.na(results_df[, "Original_symbol"]), "Original_symbol"] <- ""
  return(results_df)
}






FindSharedsgRNAs <- function(CRISPR_df) {
  all_shared_guides_df <- SharedsgRNAsDf(CRISPR_df)
  only_top4_shared_df <- SharedsgRNAsDf(CRISPR_df[CRISPR_df[, "Rank"] %in% 1:4, ])
  results_df <- data.frame(
    all_shared_guides_df[, "Sequence", drop = FALSE],
    "Num_all"  = all_shared_guides_df[, "Num_occurrences"],
    "Num_top4" = only_top4_shared_df[match(all_shared_guides_df[, "Sequence"], only_top4_shared_df[, "Sequence"]), "Num_occurrences"],
    all_shared_guides_df[, c("Gene_symbol", "Entrez_ID", "AltTSS_ID")],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  results_df <- results_df[order(-(results_df[, "Num_top4"])), ]
  rownames(results_df) <- NULL
  return(results_df)
}





SharedsgRNAsDf <- function(CRISPR_df) {

  all_sequences <- toupper(CRISPR_df[, "sgRNA_sequence"])
  num_occurrences <- table(all_sequences)

  multiples_vec <- num_occurrences[num_occurrences >= 2]

  results_df <- data.frame(
    "Sequence"        = names(multiples_vec),
    "Num_occurrences" = as.integer(multiples_vec),
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  filtered_df <- CRISPR_df[toupper(CRISPR_df[, "sgRNA_sequence"]) %in% results_df[, "Sequence"], ]
  capitalized_sequences <- toupper(filtered_df[, "sgRNA_sequence"])

  for (column in c("Gene_symbol", "Entrez_ID", "AltTSS_ID")) {

    results_df[, column] <- vapply(results_df[, "Sequence"],
                                   function(x) {
                                     are_this_sequence <- capitalized_sequences == x
                                     paste0(unique(filtered_df[are_this_sequence, column]), collapse = ", ")
                                   }, "")


  }

  results_df <- results_df[order(-results_df[, "Num_occurrences"],
                                 match(results_df[, "Sequence"], filtered_df[, "sgRNA_sequence"])
                                 ), ]
  rownames(results_df) <- NULL

  return(results_df)
}

























