### 1st January 2020 ###







# Define functions --------------------------------------------------------

RandomizeAllIndices <- function(n_total = NULL, n_per_plate_vec = NULL, n_per_plate = 384L) {

  if (is.null(n_total) && is.null(n_per_plate_vec)) {
    stop("Either n_total or n_per_plate_vec must be specified!")
  }

  if (is.null(n_per_plate_vec)) {
    num_full_plates <- floor(n_total / n_per_plate)
    last_plate_n <- n_total - (num_full_plates * n_per_plate)
    n_per_plate_vec <- rep(n_per_plate, num_full_plates)
    if (last_plate_n > 0) {
      n_per_plate_vec <- c(n_per_plate_vec, last_plate_n)
    }
  } else {
    n_total <- sum(n_per_plate_vec)
  }

  indices_pool <- seq_len(n_total)
  num_plates <- length(n_per_plate_vec)
  results_list <- vector(mode = "list", length = num_plates)

  for (plate_number in seq_len(num_plates - 1L)) {
    this_sample <- sample(indices_pool, n_per_plate_vec[[plate_number]])
    results_list[[plate_number]] <- this_sample
    indices_pool <- setdiff(indices_pool, this_sample)
  }
  results_list[[num_plates]] <- sample(indices_pool, n_per_plate_vec[[num_plates]])
  return(results_list)
}


RetrieveIndices <- function(CRISPR_df, indices_vec, ID_column) {
  unique_combined_IDs <- unique(CRISPR_df[[ID_column]])
  combined_ID_matches <- match(CRISPR_df[[ID_column]], unique_combined_IDs)
  results_df <- CRISPR_df[combined_ID_matches %in% indices_vec, ]
  return(results_df)
}



ReportProblematic <- function(are_top4_mat, CRISPR_df) {
  if (any(are_top4_mat[, "Are_problematic_top4"])) {
    problematic_IDs <- unique(CRISPR_df[["Combined_ID"]][are_top4_mat[, "Are_problematic_top4"]])
    assign("delete_problematic_IDs", problematic_IDs, envir = globalenv())
    message(paste0(length(problematic_IDs), " genes had an invalid 4sg ",
                   "combination (3 or fewer guides, or 4 guides that shared ",
                   "identical subsequences of 8 bp or more!)"
                   )
            )
  }
}


CRISPRaAreTop4Mat <- function(CRISPRa_df) {

  are_top4 <- CRISPRa_df[["Rank"]] %in% 1:4
  only_one_TSS <- CRISPRa_df[["Num_TSSs"]] == 1
  are_not_homologous <- !(is.na(CRISPRa_df[["Spacing"]])) & !(CRISPRa_df[["Spacing"]] %in% 0)

  altTSS_IDs_fac <- factor(CRISPRa_df[["AltTSS_ID"]], levels = unique(CRISPRa_df[["AltTSS_ID"]]))
  CheckThatFactorIsInOrder(altTSS_IDs_fac)
  rank_splits <- split(CRISPRa_df[["Rank"]], altTSS_IDs_fac)
  are_complete <- vapply(rank_splits, function(x) all(1:4 %in% x), logical(1))
  have_complete_guides <- rep(are_complete, lengths(rank_splits))

  are_valid_top4 <- are_top4 & are_not_homologous & have_complete_guides
  are_valid_or_only_top4 <- are_valid_top4 | (are_top4 & only_one_TSS)

  have_non_homologous_guides <- rep(tapply(are_not_homologous, altTSS_IDs_fac, sum) >= 4, lengths(rank_splits))

  results_mat <- cbind(
    "Are_chosen_4sg"       = are_valid_or_only_top4,
    "Are_top4"             = are_top4,
    "Are_problematic_top4" = are_valid_or_only_top4 & !(are_valid_top4),
    "Have_complete_guides" = have_complete_guides,
    "Have_spaced_guides"   = have_non_homologous_guides,
    "Have_valid_guides"    = have_non_homologous_guides & have_complete_guides
  )
  ReportProblematic(results_mat, CRISPRa_df)
  return(results_mat)
}




CRISPRkoAreTop4Mat <- function(CRISPRko_df) {

  are_top4 <- CRISPRko_df[["Rank"]] %in% 1:4
  are_not_homologous <- !(is.na(CRISPRko_df[["Spacing"]])) & !(CRISPRko_df[["Spacing"]] %in% 0)

  combined_IDs_fac <- factor(CRISPRko_df[["Combined_ID"]], levels = unique(CRISPRko_df[["Combined_ID"]]))
  CheckThatFactorIsInOrder(combined_IDs_fac)

  rank_splits <- split(CRISPRko_df[["Rank"]], combined_IDs_fac)
  are_complete <- vapply(rank_splits, function(x) all(1:4 %in% x), logical(1))
  have_complete_guides <- rep(are_complete, lengths(rank_splits))

  have_non_homologous_guides <- rep(tapply(are_not_homologous, combined_IDs_fac, sum) >= 4, lengths(rank_splits))

  results_mat <- cbind(
    "Are_chosen_4sg"             = are_top4, # sic
    "Are_top4"                   = are_top4,
    "Are_problematic_top4"       = are_top4 & !(are_spaced & have_complete_guides),
    "Have_complete_guides"       = have_complete_guides,
    "Have_non_homologous_guides" = have_non_homologous_guides,
    "Have_valid_guides"          = have_non_homologous_guides & have_complete_guides
  )
  ReportProblematic(results_mat, CRISPRko_df)
  return(results_mat)
}



# AreCompleteTranscripts <- function(CRISPR_df, must_be_spaced = FALSE) {
#   are_complete <- rep.int(NA, nrow(CRISPR_df))
#   unique_TSS_IDs <- unique(CRISPR_df[["AltTSS_ID"]][CRISPR_df[["Is_control"]] == "No"])
#   for (unique_TSS_ID in unique_TSS_IDs) {
#     are_this_TSS <- CRISPR_df[["AltTSS_ID"]] %in% unique_TSS_ID
#     are_unmapped_TSS <- all(CRISPR_df[["Num_TSSs"]][are_this_TSS] >= 2) && all(is.na(CRISPR_df[["Start"]][are_this_TSS]))
#     if (!(are_unmapped_TSS)) { # AltTSS_IDs corresponding to unmapped sgRNAs are left to be 'NA'
#       this_TSS_rank_vec <- CRISPR_df[["Rank"]][are_this_TSS]
#       if (must_be_spaced) {
#         this_TSS_rank_vec <- this_TSS_rank_vec[!(CRISPR_df[["Spacing"]][are_this_TSS] %in% 0)]
#       }
#       is_complete <- all(1:4 %in% this_TSS_rank_vec)
#       are_complete[are_this_TSS] <- is_complete
#     }
#   }
#   return(are_complete)
# }



ReturnProblematicGenes <- function(library_summary_df) {
  are_available_genes        <- library_summary_df[["Gene_present"]] == "Yes"
  are_overlapping_genes      <- library_summary_df[["Num_overlaps"]] != 0
  are_overlappingNoNA_genes  <- (are_overlapping_genes %in% TRUE) | (is.na(are_overlapping_genes) & are_available_genes)
  do_not_meet_criteria_genes <- library_summary_df[["Num_top4_outside_criteria"]] > 0

  are_problematic_genes <- are_overlappingNoNA_genes | do_not_meet_criteria_genes

  results_list <- list(
    "problematic"   = library_summary_df[["Entrez_ID"]][are_problematic_genes %in% TRUE],
    "unproblematic" = library_summary_df[["Entrez_ID"]][are_problematic_genes %in% FALSE],
    "overlap"       = library_summary_df[["Entrez_ID"]][are_overlappingNoNA_genes],
    "fail_criteria" = library_summary_df[["Entrez_ID"]][do_not_meet_criteria_genes %in% TRUE]
  )
  return(results_list)
}




AreGoodControls <- function(CRISPR_df) {
  are_controls <- CRISPR_df[["Is_control"]] == "Yes"
  if (any(duplicated(toupper(CRISPR_df[["sgRNA_sequence"]][are_controls])))) {
    stop("Error: Duplicated control sgRNA sequences found!")
  }
  have_no_issues <- (CRISPR_df[["Num_0MM"]] %in% 0) &
                    (CRISPR_df[["Num_1MM"]] %in% 0) &
                    !(grepl("TTTT", CRISPR_df[["sgRNA_sequence"]], ignore.case = TRUE))
  results_vec <- are_controls & have_no_issues
  return(results_vec)
}



Make4sgControlsList <- function(sequences_vec) {
  indices_list <- RandomizeAllIndices(n_total = length(sequences_vec), n_per_plate = 4)
  indices_list <- indices_list[lengths(indices_list) == 4]
  guides_pool_list <- lapply(indices_list, function(x) sequences_vec[x])
  num_homologies_vec <- vapply(guides_pool_list, NumHomologousPairs, integer(1))
  guides_pool_list <- guides_pool_list[num_homologies_vec == 0]
  return(guides_pool_list)
}



MakeControlGuidesDf <- function(control_sequences_list, controls_CRISPR_df) {

  stopifnot(all(lengths(control_sequences_list) == 4L))

  sequences_vec    <- unlist(control_sequences_list)
  well_numbers_vec <- rep(seq_along(control_sequences_list), each = 4)
  sg_numbers_vec   <- rep(1:4, times = length(control_sequences_list))

  guides_pool_matches <- match(toupper(sequences_vec), toupper(controls_CRISPR_df[["sgRNA_sequence"]]))
  results_df <- controls_CRISPR_df[guides_pool_matches, ]

  results_df[["Combined_ID"]] <- paste0("Control_", FormatFixedWidthInteger(well_numbers_vec))
  results_df[["Rank"]] <- sg_numbers_vec

  row.names(results_df) <- NULL
  return(results_df)
}






RandomlyShufflePlates <- function(plate_df_list) {

  shuffled_indices_list <- lapply(plate_df_list, function(x) sample(seq_len(nrow(x) / 4)))

  plate_df_shuffled_list <- lapply(seq_along(shuffled_indices_list), function(x) {
    assign("delete_x", x, envir = globalenv())
    my_df <- plate_df_list[[x]]
    if ("AltTSS_ID" %in% names(my_df)) {
      IDs_vec <- ifelse(my_df[["Is_control"]] == "Yes", my_df[["Combined_ID"]], my_df[["AltTSS_ID"]])
    } else {
      IDs_vec <- my_df[["Combined_ID"]]
    }
    unique_IDs_vec <- unique(IDs_vec)
    my_indices <- shuffled_indices_list[[x]]
    stopifnot(length(unique_IDs_vec) == length(my_indices))
    IDs_vec_shuffled <- unique_IDs_vec[my_indices]
    my_df <- my_df[order(match(IDs_vec, IDs_vec_shuffled)), ]
    return(my_df)
  })
  return(plate_df_shuffled_list)
}





CombinePlateDfList <- function(plate_df_list) {
  plate_df_list <- lapply(seq_along(plate_df_list), function(x) {
    my_df <- plate_df_list[[x]]
    results_df <- data.frame(
      "Plate_number" = x,
      "Well_number" = rep(seq_len(nrow(my_df) / 4), each = 4),
      my_df,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  plates_df <- do.call(rbind.data.frame, c(plate_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  # Re-number the control wells
  are_controls <- plates_df[["Is_control"]] == "Yes"
  plates_df[["Combined_ID"]][are_controls] <- paste0("Control_", rep(seq_len(sum(are_controls) / 4), each = 4))
  return(plates_df)
}















