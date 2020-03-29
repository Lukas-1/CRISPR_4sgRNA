 ### 2nd December 2019 ###



# Import packages and source code -----------------------------------------

library("parallel")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"

source(file.path(general_functions_directory, "09) Constants and settings.R"))
source(file.path(general_functions_directory, "12) Re-ordering sgRNAs based on their genomic location.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))





# Define functions --------------------------------------------------------

CreateCombinations <- function(sub_df_reordered,
                               were_included,
                               num_overlaps_allowed           = 999L,
                               tolerate_NA_overlaps           = FALSE,
                               min_space                      = 50L,
                               num_sgRNAs                     = 4L,
                               meet_strict_criteria           = NULL
                               ) {

  assign("delete_sub_df_reordered", sub_df_reordered, envir = globalenv())
  assign("delete_meet_strict_criteria", meet_strict_criteria, envir = globalenv())

  if (is.null(meet_strict_criteria)) {
    meet_strict_criteria <- rep.int(TRUE, nrow(sub_df_reordered))
  }

  have_one_0MM <- sub_df_reordered[["Num_0MM"]] == 1L

  are_core_library <- grepl("Calabrese|hCRISPRa-v2|Brunello|TKOv3", sub_df_reordered[["Source"]])

  is_CRISPRa <- "hCRISPRa_v2_rank" %in% names(sub_df_reordered)
  if (is_CRISPRa) {
    are_preferred <- grepl("Calabrese", sub_df_reordered[["Source"]], fixed = TRUE) |
                     (sub_df_reordered[["hCRISPRa_v2_rank"]] %in% 1:5)
  } else {
    are_preferred <- are_core_library | (sub_df_reordered[["GPP_rank"]] %in% 1:10)
  }
  combination_indices_mat <- combn(which(were_included), num_sgRNAs)

  cut_sites_vec <- sub_df_reordered[["Cut_location"]]
  sgRNA_sequences_vec <- toupper(sub_df_reordered[["sgRNA_sequence"]])

  guide_list <- lapply(seq_len(num_sgRNAs), function(x) x == seq_len(num_sgRNAs))

  combinations_list <- lapply(
    seq_len(ncol(combination_indices_mat)),
    function(x) {
      indices_vec <- combination_indices_mat[, x]
      locations_vec <- cut_sites_vec[indices_vec]
      overlap_numbers_vec <- vapply(guide_list,
                                    function(x) sum(abs(locations_vec[x] - locations_vec[!(x)]) < min_space),
                                    integer(1)
                                    )
      total_overlaps <- sum(overlap_numbers_vec) / 2
      too_many_overlaps <- isTRUE(total_overlaps > num_overlaps_allowed)
      results_list <- list(
        "Indices"                  = indices_vec,
        "Overlap_numbers"          = overlap_numbers_vec,
        "Total_overlaps"           = total_overlaps,
        "Too_many_overlaps"        = too_many_overlaps,
        "Num_homologies"           = if (!(too_many_overlaps)) NumHomologousPairs(sgRNA_sequences_vec[indices_vec]) else NA_integer_, # Save computational time by skipping the determination of homologies if the number of overlaps exceeds the limit
        "Mean_rank"                = mean(sub_df_reordered[["Rank"]][indices_vec]),
        "GuideScan_specificity"    = 1 / (1 + sum((1 / sub_df_reordered[["GuideScan_specificity"]][indices_vec]) - 1)),
        "CRISPOR_4MM_specificity"  = 1 / (1 + sum((1 / sub_df_reordered[["CRISPOR_4MM_specificity"]][indices_vec]) - 1)),
        "Num_meet_strict_criteria" = sum(meet_strict_criteria[indices_vec]),
        "Num_exactly_one_match"    = sum(have_one_0MM[indices_vec]),
        "Num_core_library"         = sum(are_core_library[indices_vec]),
        "Num_preferred"            = sum(are_preferred[indices_vec])
        )
      return(results_list)
    }
  )

  combinations_mat <- t(sapply(combinations_list, function(x) unlist(x[3:length(x)])))

  meet_criteria <- !((combinations_mat[, "Num_homologies"] > 0) %in% TRUE)
  if (!(tolerate_NA_overlaps)) {
    meet_criteria <- meet_criteria & !(combinations_mat[, "Too_many_overlaps"])
  }
  combinations_mat <- combinations_mat[meet_criteria, , drop = FALSE]
  combinations_list <- combinations_list[meet_criteria]

  combinations_order <- order(combinations_mat[, "Num_meet_strict_criteria"],
                              combinations_mat[, "Num_exactly_one_match"],
                              -(combinations_mat[, "Total_overlaps"]),
                              combinations_mat[, "Num_core_library"],
                              combinations_mat[, "Num_preferred"],
                              combinations_mat[, "GuideScan_specificity"],
                              combinations_mat[, "CRISPOR_4MM_specificity"],
                              -(combinations_mat[, "Mean_rank"]),
                              decreasing = TRUE
                              )
  combinations_mat <- combinations_mat[combinations_order, , drop = FALSE]
  combinations_list <- combinations_list[combinations_order]

  seq_vec <- seq_len(nrow(combinations_mat))

  indices_mat <- do.call(rbind, lapply(combinations_list, function(x) x[[1]]))

  is_contained_list_list <- lapply(seq_len(nrow(sub_df_reordered)),
                                   function(x) vapply(seq_vec,
                                                      function(y) x %in% indices_mat[y, ],
                                                      logical(1)
                                                      )
                                   )

  rank_vec <- vapply(is_contained_list_list, function(x) {
    if (any(x)) {
      seq_vec[which(x)[[1]]]
    } else {
      NA_integer_
    }
  }, integer(1))

  sub_df_reordered[["Best_combination_rank"]] <- rank_vec
  are_best_combo <- sub_df_reordered[["Best_combination_rank"]] %in% 1

  num_overlaps_vec <- rep.int(NA_integer_, nrow(sub_df_reordered))
  if (any(are_best_combo)) {
    num_overlaps_vec[are_best_combo] <- combinations_list[[1]][["Overlap_numbers"]]
  }

  if (tolerate_NA_overlaps && all(is.na(num_overlaps_vec[are_best_combo]))) {
    spacing <- 20L - 8L
  } else {
    spacing <- min_space
  }

  sub_df_reordered[["Spacing"]] <- ifelse(were_included,
                                          ifelse(are_best_combo, spacing, NA_integer_),
                                          NA_integer_
                                          )
  sub_df_reordered[["Overlaps_tolerance"]] <- ifelse(were_included,
                                                     ifelse(are_best_combo, num_overlaps_allowed, NA_integer_),
                                                     NA_integer_
                                                     )
  sub_df_reordered[["Num_overlaps"]] <- num_overlaps_vec
  sub_df_final <- sub_df_reordered[order(!(are_best_combo), sub_df_reordered[["Rank"]]), ]

  return(sub_df_final)
}




SortCombinations <- function(CRISPR_sub_df,
                             min_spaces      = 50L,
                             num_sgRNAs      = 4L,
                             only_top_24_GPP = FALSE,
                             max_overlaps    = 999L, # warning: If this is changed to, for instance, 0:12, then the number of overlaps is prioritized over the number of unambiguous guides (with exactly one 0MM location!)
                             tolerate_divergent_chromosomes = FALSE
                             ) {

  MessageID(CRISPR_sub_df)

  if (("Num_TSSs" %in% names(CRISPR_sub_df)) &&
       all(CRISPR_sub_df[["Num_TSSs"]] >= 2) &&
       all(is.na(CRISPR_sub_df[["Start"]]))
      ) {
    CRISPR_sub_df[["Best_combination_rank"]] <- NA_integer_
    CRISPR_sub_df[["Spacing"]]               <- NA_integer_
    CRISPR_sub_df[["Overlaps_tolerance"]]    <- NA_integer_
    CRISPR_sub_df[["Num_overlaps"]]          <- NA_integer_
    CRISPR_sub_df[["Original_rank"]]         <- CRISPR_sub_df[["Rank"]]
    CRISPR_sub_df[["Rank"]]                  <- NA_integer_
    return(CRISPR_sub_df)
  }

  reordered_list <- ReorderSubDfByLocation(CRISPR_sub_df, tolerate_divergent_chromosomes = tolerate_divergent_chromosomes)
  sub_df_reordered <- reordered_list[["reordered_df"]]
  were_mapped <- reordered_list[["were_mapped_vec"]]

  are_polyT          <- grepl("TTTT", sub_df_reordered[["sgRNA_sequence"]], ignore.case = TRUE)
  have_canonical_PAM <- (substr(sub_df_reordered[["PAM"]], 2, 3) == "GG") %in% TRUE
  are_curated        <- sub_df_reordered[["Source"]] == "Curated"
  are_specific       <- ((sub_df_reordered[["GuideScan_specificity"]] < 0.2) %in% FALSE) |
                        (is.na(sub_df_reordered[["GuideScan_specificity"]]) &
                         ((sub_df_reordered[["CRISPOR_3MM_specificity"]] < 0.2) %in% FALSE)
                         )

  if ("Exon_number_GPP" %in% names(sub_df_reordered)) { # ==> CRISPRko
    violate_Graf_criteria <- sub_df_reordered[["CRISPOR_Graf_status"]] %in% c("ggc", "tt")
  } else {
    violate_Graf_criteria <- rep.int(FALSE, nrow(sub_df_reordered))
  }

  were_included <- were_mapped &
                   are_specific &
                   !(are_polyT) &
                   have_canonical_PAM &
                   !((sub_df_reordered[[preferred_AF_max_column]] > SNP_frequency_cutoff) %in% TRUE) &
                   !(violate_Graf_criteria) &
                   !(are_curated)
                   #((sub_df_reordered[["Source"]] != "GPP") | (sub_df_reordered[["GPP_rank"]] %in% 1:24))

  if (only_top_24_GPP) {
    are_GPP_top_24 <- (sub_df_reordered[["Source"]] != "GPP") | (sub_df_reordered[["GPP_rank"]] %in% 1:24)
    were_included <- were_included & are_GPP_top_24
  }

  sgRNAs_found <- FALSE
  FoundsgRNA <- function(sub_df) sum(sub_df[["Best_combination_rank"]] %in% 1) >= num_sgRNAs
  for (max_overlap in max_overlaps) {
    for (min_space in min_spaces) {
      if (!(sgRNAs_found) && (sum(were_included) >= num_sgRNAs)) {
        sub_df_final <- CreateCombinations(sub_df_reordered,
                                           were_included,
                                           min_space            = min_space,
                                           num_sgRNAs           = num_sgRNAs,
                                           num_overlaps_allowed = max_overlap
                                           )
        sgRNAs_found <- FoundsgRNA(sub_df_final)
        if (sgRNAs_found) {
          break
        }
      }
    }
  }

  if (!(sgRNAs_found)) {
    if (sum(!(are_polyT | are_curated)) >= num_sgRNAs) {
      sub_df_final <- CreateCombinations(sub_df_reordered,
                                         were_included        = !(are_polyT | are_curated),
                                         min_space            = min(min_spaces),
                                         num_sgRNAs           = num_sgRNAs,
                                         num_overlaps_allowed = max_overlap,
                                         tolerate_NA_overlaps = TRUE,
                                         meet_strict_criteria = were_included
                                         )
      sgRNAs_found <- FoundsgRNA(sub_df_final)
    }
  }

  if (!(sgRNAs_found)) {
    are_core_library <- grepl("Calabrese|hCRISPRa-v2", sub_df_reordered[["Source"]])
    sub_df_final <- sub_df_reordered[order(!(are_core_library), sub_df_reordered[["Rank"]]), ]
    sub_df_final[["Best_combination_rank"]] <- NA_integer_
    sub_df_final[["Spacing"]]               <- ifelse(sub_df_final[["Rank"]] %in% 1:4, 0L, NA_integer_)
    sub_df_final[["Overlaps_tolerance"]]    <- NA_integer_
    sub_df_final[["Num_overlaps"]]          <- NA_integer_
    sgRNAs_found <- FALSE
  }

  sub_df_final[["Original_rank"]] <- sub_df_final[["Rank"]]
  sub_df_final[["Rank"]] <- seq_len(nrow(sub_df_final))

  row.names(sub_df_final) <- NULL
  return(sub_df_final)
}




NonOverlappingDfForControls <- function(CRISPR_df) {
   are_controls <- CRISPR_df[["Is_control"]] == "Yes"
   if (any(are_controls)) {
     controls_df <- CRISPR_df[are_controls, ]
     controls_df[["Best_combination_rank"]] <- NA_integer_
     controls_df[["Original_rank"]]         <- NA_integer_
     controls_df[["Spacing"]]               <- NA_integer_
     controls_df[["Overlaps_tolerance"]]    <- NA_integer_
     controls_df[["Num_overlaps"]]          <- NA_integer_
   } else {
     controls_df <- NULL
   }
  return(controls_df)
}



GetUniqueIDs <- function(CRISPR_df, ID_column = "AltTSS_ID") {
  are_controls <- CRISPR_df[["Is_control"]] == "Yes"
  unique_IDs <- unique(CRISPR_df[[ID_column]][!(are_controls)])
  return(unique_IDs)
}




PrioritizeNonOverlapping <- function(CRISPR_df,
                                     ID_column       = "AltTSS_ID",
                                     max_overlaps    = 999L,
                                     min_spaces      = 50L,
                                     only_top_24_GPP = FALSE,
                                     parallel_mode   = TRUE,
                                     num_cores       = NULL,
                                     tolerate_divergent_chromosomes = FALSE
                                     ) {

  unique_IDs <- GetUniqueIDs(CRISPR_df, ID_column)
  controls_df <- NonOverlappingDfForControls(CRISPR_df)

  if (parallel_mode) {
    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores() - 2
    }
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterExport(cl,
                            varlist = c("SortCombinations", "CreateCombinations", "MessageID",
                                        "ReorderSubDfByLocation", "NumHomologousPairs", "SplitIntoSubstrings",
                                        "CRISPR_df", "tolerate_divergent_chromosomes",
                                        "preferred_AF_max_column", "SNP_frequency_cutoff",
                                        "max_overlaps", "min_spaces", "only_top_24_GPP"
                                        ),
                            envir = environment()
                            )
    reordered_df_list <- parallel::parLapply(cl,
                                             unique_IDs,
                                             function(x) SortCombinations(CRISPR_df[CRISPR_df[[ID_column]] == x, , drop = FALSE],
                                                                          max_overlaps    = max_overlaps,
                                                                          min_spaces      = min_spaces,
                                                                          only_top_24_GPP = only_top_24_GPP,
                                                                          tolerate_divergent_chromosomes = tolerate_divergent_chromosomes
                                                                          )
                                             )
    parallel::stopCluster(cl)

  } else {
    reordered_df_list <- lapply(unique_IDs,
                                function(x) SortCombinations(CRISPR_df[CRISPR_df[[ID_column]] == x, , drop = FALSE],
                                                             min_spaces      = min_spaces,
                                                             only_top_24_GPP = only_top_24_GPP,
                                                             max_overlaps    = max_overlaps,
                                                             tolerate_divergent_chromosomes = tolerate_divergent_chromosomes
                                                             )
                                )
  }

  results_df <- do.call(rbind.data.frame,
                        c(reordered_df_list,
                          list(controls_df),
                          list(stringsAsFactors = FALSE, make.row.names = FALSE)
                          )
                        )
  return(results_df)
}


