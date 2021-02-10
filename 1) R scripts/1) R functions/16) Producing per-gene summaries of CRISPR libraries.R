### 17th October 2019 ###





# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "09) Constants and settings.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # For RoundNumericColumns
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R")) # For CRISPRkoAreTop4Mat





# Define constants --------------------------------------------------------

TF_annotation_columns <- c("Gene_symbol", "Ensembl_gene_ID", "Entrez_ID", "Original_symbol", "Original_entrez", "Sources")
all_genes_annotation_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "Original_entrez", "Sources")

TSS_columns <- c("Num_hCRISPR_v2_transcripts", "Num_transcripts", "Num_overlapping_transcripts", "Num_incomplete_transcripts")

selected_metrics <- c("Num_overlaps", "Spacing", "Longest_subsequence", "GuideScan_specificity",
                      "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity",
                      "Num_top4_outside_criteria", "Num_total",
                      "Specific_guides_available", "Num_no_perfect_match", "Gene_annotation_status"
                      )




# General-use helper functions --------------------------------------------

MeetCriteria <- function(CRISPR_df, allow_curated = FALSE) {
  are_to_exclude <- ((CRISPR_df[[preferred_AF_max_column]] > SNP_frequency_cutoff) %in% TRUE) |
                    (is.na(CRISPR_df[["GuideScan_specificity"]]) & is.na(CRISPR_df[["CRISPOR_4MM_specificity"]])) |
                    grepl("TTTT", CRISPR_df[["sgRNA_sequence"]], ignore.case = TRUE) |
                    !(((substr(CRISPR_df[["PAM"]], 2, 3) == "GG") %in% TRUE)) |
                    (((CRISPR_df[["GuideScan_specificity"]] < 0.2) %in% TRUE) |
                      (is.na(CRISPR_df[["GuideScan_specificity"]]) & ((CRISPR_df[["CRISPOR_3MM_specificity"]] < 0.2) %in% TRUE))) |
                    (("Exon_number_GPP" %in% names(CRISPR_df)) & (CRISPR_df[["CRISPOR_Graf_status"]] %in% c("ggc", "tt")))
  stopifnot(length(are_to_exclude) == nrow(CRISPR_df))
  if (!(allow_curated)) {
    are_to_exclude <- are_to_exclude | (CRISPR_df[["Source"]] == "Curated")
  }
  results_vec <- !(are_to_exclude)
  return(results_vec)
}


AggregateSpecificityScores <- function(scores_vec) {
  stopifnot(length(scores_vec) > 0)
  1 / (1 + sum((1 / scores_vec) - 1))
}




# Functions for re-formatting library source information ------------------

libraries_order <- c(
  "Calabrese",
  "hCRISPRa-v2",
  "Dolcetto",
  "hCRISPRi-v2.1",
  "hCRISPRi-v2.0",
  "Brunello",
  "TKOv3",
  "GPP",
  "Caprano",
  "mCRISPRa-v2",
  "Dolomiti"
)

libraries_two_letters <- c(
  "Calabrese"     = "Ca",
  "hCRISPRa-v2"   = "v2",
  "Dolcetto"      = "Do",
  "hCRISPRi-v2.1" = "v2",
  "hCRISPRi-v2.0" = "v2",
  "Brunello"      = "B",
  "TKOv3"         = "T",
  "GPP"           = "G",
  "Caprano"       = "Ca",
  "mCRISPRa-v2"   = "v2",
  "Dolomiti"      = "Do"
)


ReformatSourceToFactor <- function(source_vec, strip_curated = TRUE) {
  if (strip_curated) {
    source_vec <- sub("Curated, ", "", source_vec, fixed = TRUE)
  }
  unique_sources <- unique(source_vec)
  reordered_unique_sources <- vapply(strsplit(unique_sources, ", ", fixed = TRUE),
                                     function(x) paste0(x[order(x == "GPP")], collapse = ", "),
                                     ""
                                     )
  for (i in which(unique_sources != reordered_unique_sources)) {
    source_vec[source_vec == unique_sources[[i]]] <- reordered_unique_sources[[i]]
  }
  levels_present <- libraries_order[libraries_order %in% source_vec]
  stopifnot(length(levels_present) == 3)
  levels_expanded <- c(apply(combn(levels_present, 2), 2, paste0, collapse = ", "),
                       paste0(levels_present, collapse = ", ")
                       )
  source_fac <- factor(source_vec, levels = c(levels_present, levels_expanded))
  return(source_fac)
}


SummarizeSources <- function(source_sub_vec) {
  unique_sources <- unique(unlist(strsplit(source_sub_vec, ", ", fixed = TRUE)))
  unique_sources <- unique_sources[order(match(unique_sources, names(libraries_two_letters)))]
  unique_sources_abbr <- libraries_two_letters[match(unique_sources, names(libraries_two_letters))]
  result_string <- paste0(unique_sources_abbr, collapse = ", ")
  return(result_string)
}





# Functions for generating overview tables --------------------------------

CollapseOriginal <- function(char_vec) {
  paste0(unique(char_vec[char_vec != ""]), collapse = "")
}


SummarizeCRISPRDf <- function(CRISPR_df, sublibraries_entrezs_list) {

  are_curated <- CRISPR_df[["Source"]] == "Curated"
  if (any(are_curated)) {
    CRISPR_df <- CRISPR_df[!(are_curated), ]
  }

  split_indices <- split(seq_len(nrow(CRISPR_df)), factor(CRISPR_df[["Combined_ID"]], levels = unique(CRISPR_df[["Combined_ID"]])))
  include_top4nonoverlapping <- "Best_combination_rank" %in% names(CRISPR_df)
  if (include_top4nonoverlapping) {
    CRISPR_df[["Non_overlapping"]] <- !(is.na(CRISPR_df[["Best_combination_rank"]]))
  } else {
    CRISPR_df[["Non_overlapping"]] <- NA
  }
  CRISPR_df[["Overlap_with_SNP"]] <- (CRISPR_df[[preferred_AF_max_column]] > SNP_frequency_cutoff) %in% TRUE
  CRISPR_df[["Meet_criteria"]] <- MeetCriteria(CRISPR_df)
  if (include_top4nonoverlapping) {
    CRISPR_df[["Meet_criteria"]] <- CRISPR_df[["Meet_criteria"]] & CRISPR_df[["Non_overlapping"]]
  }
  CRISPR_df[["Are_unspecific"]] <- ((CRISPR_df[["Num_0MM"]] > 1) | (CRISPR_df[["Num_1MM"]] > 0)) %in% TRUE
  include_TSS                   <- "TSS_regions" %in% names(CRISPR_df)
  include_searched_by_GuideScan <- "TSS_searched_by_GuideScan" %in% names(CRISPR_df)
  include_transcripts           <- "TSS_ID" %in% names(CRISPR_df)
  if (include_transcripts && include_top4nonoverlapping) {
    are_incomplete <- (CRISPR_df[["Spacing"]] %in% 0) & (CRISPR_df[["Rank"]] %in% 1:4)
  }

  are_valid_4sg <- Are4sg(CRISPR_df, sublibraries_entrezs_list)

  results_list_list <- lapply(split_indices, function(x) {
    results_list <- list(
      "Combined_ID"                 = unique(CRISPR_df[["Combined_ID"]][x]),
      "Entrez_ID"                   = unique(CRISPR_df[["Entrez_ID"]][x]),
      "Gene_symbol"                 = unique(CRISPR_df[["Gene_symbol"]][x]),
      "Original_entrez"             = CollapseOriginal(CRISPR_df[["Original_entrez"]][x]),
      "Original_symbol"             = CollapseOriginal(CRISPR_df[["Original_symbol"]][x]),
      "Sources"                     = SummarizeSources(CRISPR_df[["Source"]][x]),
      "In_4sg_library"              = if (sum(are_valid_4sg[x]) >= 4) "Yes" else "No",
      "Num_hCRISPR_v2_transcripts"  = NA_integer_,
      "Num_transcripts"             = NA_integer_,
      "Num_overlapping_transcripts" = NA_integer_,
      "Num_unspaced_transcripts"    = NA_integer_,
      "Num_incomplete_transcripts"  = NA_integer_,
      "Spacing"                     = NA_character_,
      "Deletion_size"               = NA_integer_,
      "GuideScan_specificity"       = NA_real_,
      "CRISPOR_3MM_specificity"     = NA_real_,
      "CRISPOR_4MM_specificity"     = NA_real_,
      "Longest_subsequence"         = NA_integer_,
      "Num_total"                   = length(x),
      "Num_overlaps"                = NA_integer_,
      "Num_top4_outside_criteria"   = NA_integer_,
      "Num_meeting_criteria"        = sum(CRISPR_df[["Meet_criteria"]][x]),
      "Num_without_GuideScan"       = sum(is.na(CRISPR_df[["GuideScan_specificity"]][x])),
      "Num_unspecific"              = sum(CRISPR_df[["Are_unspecific"]][x]),
      "Num_overlapping_with_SNP"    = sum(CRISPR_df[["Overlap_with_SNP"]][x]),
      "Num_no_perfect_match"        = sum(CRISPR_df[["Num_0MM"]][x] == 0),
      "Num_top4_no_perfect_match"   = NA_integer_,
      "Specific_guides_available"   = if (all(is.na(CRISPR_df[["Start"]][x]) & (CRISPR_df[["Num_0MM"]][x] != 1))) "No" else "Yes",
      "Submitted_to_GuideScan"      = NA_character_,
      "TSS_regions"                 = NA_character_
      )
    if (include_searched_by_GuideScan) {
      were_searched <- CRISPR_df[["TSS_searched_by_GuideScan"]][x]
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
      results_list[["TSS_regions"]] <- unique(CRISPR_df[["TSS_regions"]][x])
    } else {
      results_list[["TSS_regions"]] <- NULL
    }
    if (include_transcripts) {
      results_list[["Deletion_size"]] <- NULL
      transcripts_vec <- CRISPR_df[["TSS_ID"]][x]
      if (all(is.na(transcripts_vec))) {
        transcripts_vec <- CRISPR_df[["AltTSS_ID"]][x]
      }
      transcripts_fac <- factor(transcripts_vec, levels = unique(transcripts_vec))
      transcripts_indices <- split(seq_along(x), transcripts_fac)
      transcripts_top4_indices <- lapply(transcripts_indices, function(y) y[CRISPR_df[["Rank"]][x][y] %in% 1:4])
    }
    if (include_top4nonoverlapping) {
      is_control <- all(CRISPR_df[["Is_control"]][x] == "Yes")
      if (!(is_control)) {
        found_overlap <- FALSE
        for (space in c(50, 45, 40)) {
          are_this_spacing <- CRISPR_df[["Spacing"]][x] %in% space
          if (any(are_this_spacing)) {
            num_zero_overlaps <- sum(CRISPR_df[["Num_overlaps"]][x][are_this_spacing] %in% 0)
            overlap_string <- paste0(sum(num_zero_overlaps), "*", space)
            found_overlap <- TRUE
            break
          }
        }
        if (!(found_overlap)) {
          if (any(CRISPR_df[["Spacing"]][x] %in% 12)) {
            overlap_string <- ">12bp"
          } else {
            overlap_string <- "None"
          }
        }
        results_list[["Spacing"]] <- overlap_string
        are_top_4 <- CRISPR_df[["Rank"]][x] %in% 1:4
        results_list[["Num_overlaps"]] <- sum(CRISPR_df[["Num_overlaps"]][x[are_top_4]])
        results_list[["Num_top4_outside_criteria"]] <- sum(are_top_4) - sum(CRISPR_df[["Meet_criteria"]][x[are_top_4]])
        results_list[["Num_top4_no_perfect_match"]] <- sum(CRISPR_df[["Num_0MM"]][x][are_top_4] == 0)
        if (include_transcripts) {
          guidescan_spec_vec <- vapply(transcripts_top4_indices,
                                       function(y) AggregateSpecificityScores(CRISPR_df[["GuideScan_specificity"]][x][y]),
                                       numeric(1)
                                       )
          CRISPOR_3MM_spec_vec <- vapply(transcripts_top4_indices,
                                         function(y) AggregateSpecificityScores(CRISPR_df[["CRISPOR_3MM_specificity"]][x][y]),
                                         numeric(1)
                                         )
          CRISPOR_4MM_spec_vec <- vapply(transcripts_top4_indices,
                                         function(y) AggregateSpecificityScores(CRISPR_df[["CRISPOR_4MM_specificity"]][x][y]),
                                         numeric(1)
                                         )
          results_list[["GuideScan_specificity"]]   <- if (all(is.na(guidescan_spec_vec)))   NA_real_ else min(guidescan_spec_vec,   na.rm = TRUE)
          results_list[["CRISPOR_3MM_specificity"]] <- if (all(is.na(CRISPOR_3MM_spec_vec))) NA_real_ else min(CRISPOR_3MM_spec_vec, na.rm = TRUE)
          results_list[["CRISPOR_4MM_specificity"]] <- if (all(is.na(CRISPOR_4MM_spec_vec))) NA_real_ else min(CRISPOR_4MM_spec_vec, na.rm = TRUE)
          subsequence_lengths <- vapply(transcripts_top4_indices,
                                        function(y) LongestSharedSubsequence(CRISPR_df[["sgRNA_sequence"]][x][y]),
                                        integer(1)
                                        )
          if (all(is.na(subsequence_lengths))) {
            results_list[["Longest_subsequence"]] <- NA_integer_
          } else {
            results_list[["Longest_subsequence"]] <- max(subsequence_lengths, na.rm = TRUE)
          }
        } else {
          cut_locations_vec <- CRISPR_df[["Cut_location"]][x[are_top_4]]
          chromosomes_vec <- CRISPR_df[["Chromosome"]][x[are_top_4]]
          assign("delete_cut_locations_vec", cut_locations_vec, envir = globalenv())
          assign("delete_are_top_4", are_top_4, envir = globalenv())
          assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
          assign("delete_x", x, envir = globalenv())
          if (anyNA(chromosomes_vec) || (length(unique(chromosomes_vec)) > 1)) {
            results_list[["Deletion_size"]] <- NA
          } else {
            results_list[["Deletion_size"]] <- max(cut_locations_vec) - min(cut_locations_vec)
          }
          results_list[["GuideScan_specificity"]]   <- AggregateSpecificityScores(CRISPR_df[["GuideScan_specificity"]][x[are_top_4]])
          results_list[["CRISPOR_3MM_specificity"]] <- AggregateSpecificityScores(CRISPR_df[["CRISPOR_3MM_specificity"]][x[are_top_4]])
          results_list[["CRISPOR_4MM_specificity"]] <- AggregateSpecificityScores(CRISPR_df[["CRISPOR_4MM_specificity"]][x[are_top_4]])
          results_list[["Longest_subsequence"]]     <- LongestSharedSubsequence(CRISPR_df[["sgRNA_sequence"]][x[are_top_4]])
        }
      }
    } else {
      results_list[["Deletion_size"]]             <- NULL
      results_list[["Spacing"]]                   <- NULL
      results_list[["Num_overlaps"]]              <- NULL
      results_list[["Num_top4_outside_criteria"]] <- NULL
      results_list[["GuideScan_specificity"]]     <- NULL
      results_list[["CRISPOR_3MM_specificity"]]   <- NULL
      results_list[["CRISPOR_4MM_specificity"]]   <- NULL
      results_list[["Longest_subsequence"]]       <- NULL
    }
    if (include_transcripts) {
      transcripts_vec <- CRISPR_df[["TSS_ID"]][x]
      if (all(is.na(transcripts_vec))) {
        transcripts_vec <- CRISPR_df[["AltTSS_ID"]][x]
      }
      unique_transcripts <- unique(transcripts_vec)
      unique_transcripts <- unique_transcripts[!(is.na(unique_transcripts))]
      transcript_column <- grep("_v2_transcript", colnames(CRISPR_df), fixed = TRUE, value = TRUE)
      unique_hC_transcripts <- unique(CRISPR_df[[transcript_column]][x])
      unique_hC_transcripts <- unique_hC_transcripts[!(is.na(unique_hC_transcripts))]
      results_list[["Num_transcripts"]] <- length(unique_transcripts)
      results_list[["Num_hCRISPR_v2_transcripts"]] <- length(unique_hC_transcripts)
    } else {
      results_list[["Num_transcripts"]] <- NULL
      results_list[["Num_hCRISPR_v2_transcripts"]] <- NULL
    }
    if (include_transcripts && include_top4nonoverlapping) {
      results_list[["Num_unspaced_transcripts"]]    <- sum(tapply(are_incomplete[x], transcripts_fac, any))
      results_list[["Num_incomplete_transcripts"]]  <- sum(lengths(transcripts_top4_indices) != 4)
      results_list[["Num_overlapping_transcripts"]] <- sum(!(tapply(CRISPR_df[["Num_overlaps"]][x], transcripts_fac, function(z) sum(z %in% 0) == 4)))
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
  are_duplicated <- duplicated(reference_IDs)
  if (any(are_duplicated)) {
    stop(paste0("The following IDs were duplicated: ", paste0(reference_IDs[are_duplicated], collapse = ", "), "!"))
  }
  summary_mat <- as.matrix(summary_df[, grep("^Num_", names(summary_df), value = TRUE)])
  summary_matches <- match(reference_IDs, summary_df[["Combined_ID"]])
  matched_summary_mat <- summary_mat[summary_matches, ]
  # matched_summary_mat[is.na(matched_summary_mat)] <- 0L
  results_df <- data.frame("Combined_ID" = reference_IDs,
                           "Gene_present" = ifelse(is.na(summary_matches), "No", "Yes"),
                           matched_summary_mat,
                           summary_df[summary_matches, !(names(summary_df) %in% colnames(summary_mat))],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  column_index <- match("Original_symbol", names(summary_df))
  colnames_before <- names(summary_df)[seq_len(column_index)]
  colnames_after <- names(summary_df)[(column_index + 1):ncol(summary_df)]
  results_df <- results_df[, c(colnames_before, "Gene_present", colnames_after)]

  include_spacing <- "Spacing" %in% names(results_df)
  include_transcripts <- "TSS_ID" %in% names(results_df)
  if (include_spacing) {
    have_no_space <- (results_df[["Spacing"]] %in% c("None", ">12bp")) | is.na(results_df[["Spacing"]])
    spacing_splits <- strsplit(results_df[["Spacing"]], "*", fixed = TRUE)
    spacing_vec <- vapply(seq_along(have_no_space), function(x) if (have_no_space[[x]]) 0L else as.integer(spacing_splits[[x]][[2]]), integer(1))
    num_nonoverlapping_vec <- vapply(seq_along(have_no_space), function(x) if (have_no_space[[x]]) 0L else as.integer(spacing_splits[[x]][[1]]), integer(1))
  }

  NA_vec <- rep.int(NA, nrow(results_df))

  final_order <- order(
    results_df[["Gene_present"]],
    ifelse(results_df[["Num_total"]] < 4, results_df[["Num_total"]], 4L),
    if (include_spacing) !(results_df[["Spacing"]] %in% "None")                                                               else NA_vec,
    if (include_spacing) !(results_df[["Spacing"]] %in% ">12bp")                                                              else NA_vec,
    if (include_spacing && include_transcripts) -(results_df[["Num_incomplete_transcripts"]])                                 else NA_vec,
    if (include_spacing) -(results_df[["Num_top4_outside_criteria"]])                                                         else NA_vec,
    if (include_spacing && include_transcripts) (results_df[["Num_transcripts"]] != results_df[["Num_unspaced_transcripts"]]) else NA_vec,
    if (include_spacing && include_transcripts) -(results_df[["Num_unspaced_transcripts"]])                                   else NA_vec,
    if (include_spacing && include_transcripts) -(results_df[["Num_overlapping_transcripts"]])                                else NA_vec,
    if (include_spacing) num_nonoverlapping_vec                                                                               else NA_vec,
    if (include_spacing) spacing_vec                                                                                          else NA_vec,
    if (include_spacing) -(results_df[["Num_overlaps"]])                                                                      else NA_vec,
    ifelse(results_df[["Num_meeting_criteria"]] < 4, results_df[["Num_meeting_criteria"]], NA_integer_),
    if (include_spacing) !(is.na(results_df[["GuideScan_specificity"]]) & is.na(results_df[["CRISPOR_4MM_specificity"]]))     else NA_vec,
    if (include_spacing) !(is.na(results_df[["GuideScan_specificity"]]))                                                      else NA_vec,
    if (include_spacing) results_df[["GuideScan_specificity"]]                                                                else NA_vec,
    if (include_spacing) !(is.na(results_df[["CRISPOR_4MM_specificity"]]))                                                    else NA_vec,
    if (include_spacing) results_df[["CRISPOR_4MM_specificity"]]                                                              else NA_vec
  )
  results_df <- results_df[final_order, ]
  entrezs_vec <- results_df[["Entrez_ID"]]
  are_not_present <- results_df[["Gene_present"]] == "No"
  entrezs_vec[are_not_present] <- results_df[["Combined_ID"]][are_not_present]
  results_df[["Gene_annotation_status"]] <- collected_entrezs_df[["Category"]][match(entrezs_vec, collected_entrezs_df[["Entrez_ID"]])]
  results_df[["Gene_annotation_status"]] <- ifelse(entrezs_vec %in% collected_entrez_IDs,
                                                   results_df[["Gene_annotation_status"]],
                                                   "Not protein-coding"
                                                   )
  row.names(results_df) <- NULL
  return(results_df)
}




ProduceGenomeOverviewDf <- function(strict_CRISPR_df, sublibraries_entrezs_list, lax_CRISPR_df = NULL, use_lax_df = FALSE, is_mouse = FALSE) {
  ## requires 'collected_entrez_IDs' in the global environment

  ## Collect all Entrez IDs from various sources
  CRISPR_df_entrez_IDs <- unique(strict_CRISPR_df[["Entrez_ID"]])
  CRISPR_df_entrez_IDs <- CRISPR_df_entrez_IDs[!(is.na(CRISPR_df_entrez_IDs))]
  CRISPR_df_entrez_IDs <- unique(unlist(strsplit(CRISPR_df_entrez_IDs, ", ", fixed = TRUE)))
  unique_entrez_IDs <- union(collected_entrez_IDs, CRISPR_df_entrez_IDs)

  ## Create an sgRNA overview data frame, with one row per gene
  if (use_lax_df) {
    sgRNAs_summary_df <- SummarizeCRISPRDf(lax_CRISPR_df, sublibraries_entrezs_list)
  } else {
    sgRNAs_summary_df <- SummarizeCRISPRDf(strict_CRISPR_df, sublibraries_entrezs_list)
  }
  sgRNAs_all_genes_df <- ReorganizeSummaryDf(sgRNAs_summary_df, unique_entrez_IDs)
  sgRNAs_all_genes_df[["Entrez_ID"]] <- sgRNAs_all_genes_df[["Combined_ID"]]
  sgRNAs_all_genes_df <- sgRNAs_all_genes_df[, names(sgRNAs_all_genes_df) != "Combined_ID"]
  sgRNAs_overview_df <- FixSymbolsForSummaryDf(sgRNAs_all_genes_df, is_mouse = is_mouse)
  if (!(is.null(lax_CRISPR_df))) {
    are_different <- DifferUsingRelaxedLocations(sgRNAs_overview_df[["Entrez_ID"]], strict_CRISPR_df, lax_CRISPR_df)
    sgRNAs_overview_df[["Lax_locations_differ"]] <- ifelse(are_different, "Yes", "No")
  }
  return(sgRNAs_overview_df)
}




# Functions for exporting overview tables ---------------------------------

FixSymbolsForSummaryDf <- function(reorganized_df, is_mouse = FALSE) {
  have_no_symbol <- is.na(reorganized_df[["Gene_symbol"]])
  entrez_symbols_df <- MapToEntrezs(reorganized_df[["Entrez_ID"]][have_no_symbol], is_mouse = is_mouse)
  could_be_mapped <- !(is.na(entrez_symbols_df[["Gene_symbol"]]))
  results_df <- reorganized_df
  results_df[["Gene_symbol"]][have_no_symbol][could_be_mapped] <- entrez_symbols_df[["Gene_symbol"]][could_be_mapped]
  results_df[["Original_entrez"]][is.na(results_df[["Original_entrez"]])] <- ""
  results_df[["Original_symbol"]][is.na(results_df[["Original_symbol"]])] <- ""
  return(results_df)
}


FormatOverviewDfForExport <- function(overview_df) {
  for (column_name in intersect(c("Num_no_perfect_match", "Num_top4_no_perfect_match"), colnames(overview_df))) {
    overview_df[[column_name]] <- ifelse(overview_df[[column_name]] == 0,
                                         "",
                                         as.character(overview_df[[column_name]])
                                         )
  }
  if ("Lax_locations_differ" %in% colnames(overview_df)) {
    overview_df[["Lax_locations_differ"]] <- ifelse(is.na(overview_df[["Num_total"]]), NA_character_, overview_df[["Lax_locations_differ"]])
  }
  overview_df[["Num_total"]] <- ifelse(is.na(overview_df[["Num_total"]]), 0L, overview_df[["Num_total"]])
  results_df <- RoundNumericColumns(overview_df)
  for (i in seq_along(results_df)) {
    results_df[[i]] <- ifelse(is.na(results_df[[i]]), "", as.character(results_df[[i]]))
  }
  return(results_df)
}


WriteOverviewDfToDisk <- function(overview_df, file_name, use_directory = file_output_directory) {
  write.table(FormatOverviewDfForExport(overview_df),
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )
  return(invisible(NULL))
}







# Functions for identifying sgRNAs that are shared between genes ----------

FindSharedsgRNAs <- function(CRISPR_df) {
  all_shared_guides_df <- SharedsgRNAsDf(CRISPR_df)
  only_top4_shared_df <- SharedsgRNAsDf(CRISPR_df[CRISPR_df[["Rank"]] %in% 1:4, ])
  results_df <- data.frame(
    all_shared_guides_df["Sequence"],
    "Num_all"  = all_shared_guides_df[["Num_occurrences"]],
    "Num_top4" = only_top4_shared_df[["Num_occurrences"]][match(all_shared_guides_df[["Sequence"]], only_top4_shared_df[["Sequence"]])],
    all_shared_guides_df[, c("Gene_symbol", "Entrez_ID", "AltTSS_ID")],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  results_df <- results_df[order(-(results_df[["Num_top4"]])), ]
  row.names(results_df) <- NULL
  return(results_df)
}



SharedsgRNAsDf <- function(CRISPR_df) {

  all_sequences <- toupper(CRISPR_df[["sgRNA_sequence"]])
  num_occurrences <- table(all_sequences)
  multiples_vec <- num_occurrences[num_occurrences >= 2]

  results_df <- data.frame(
    "Sequence"        = names(multiples_vec),
    "Num_occurrences" = as.integer(multiples_vec),
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  filtered_df <- CRISPR_df[toupper(CRISPR_df[["sgRNA_sequence"]]) %in% results_df[["Sequence"]], ]
  capitalized_sequences <- toupper(filtered_df[["sgRNA_sequence"]])

  for (column in c("Gene_symbol", "Entrez_ID", "AltTSS_ID")) {
    results_df[[column]] <- vapply(results_df[["Sequence"]],
                                   function(x) {
                                     are_this_sequence <- capitalized_sequences == x
                                     paste0(unique(filtered_df[[column]][are_this_sequence]), collapse = ", ")
                                   }, "")
  }

  results_df <- results_df[order(-results_df[["Num_occurrences"]],
                                 match(results_df[["Sequence"]], filtered_df[["sgRNA_sequence"]])
                                 ), ]
  row.names(results_df) <- NULL
  return(results_df)
}






# Functions for finding genes affected by relaxed/strict locations --------

DifferUsingRelaxedLocations <- function(combined_IDs, strict_df, lax_df, filter_for_top4 = TRUE) {

  results_vec <- rep(NA, length(combined_IDs))
  are_not_NA <- !(is.na(combined_IDs))
  combined_IDs <- combined_IDs[are_not_NA]

  if (any(duplicated(combined_IDs))) {
    stop("Duplicated entries are not allowed for the 'combined_IDs' parameter!")
  }

  are_included_lax <- lax_df[["Combined_ID"]] %in% combined_IDs
  are_included_strict <- strict_df[["Combined_ID"]] %in% combined_IDs

  if (filter_for_top4) {
    are_included_lax <- are_included_lax & lax_df[["Rank"]] %in% 1:4
    are_included_strict <- are_included_strict & strict_df[["Rank"]] %in% 1:4
  }

  lax_df <- lax_df[are_included_lax, ]
  strict_df <- strict_df[are_included_strict, ]

  lax_seq_list <- split(lax_df[["sgRNA_sequence"]], factor(lax_df[["Combined_ID"]], levels = combined_IDs))
  strict_seq_list <- split(strict_df[["sgRNA_sequence"]], factor(strict_df[["Combined_ID"]], levels = combined_IDs))

  are_identical <- mapply(identical, lax_seq_list, strict_seq_list)
  results_vec[are_not_NA] <- !(are_identical)
  return(results_vec)
}





