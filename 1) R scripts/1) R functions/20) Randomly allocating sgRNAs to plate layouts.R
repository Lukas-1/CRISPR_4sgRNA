### 14th April 2020 ###



# Functions for selecting the 4sg guides to export ------------------------

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



# Functions for splitting genes into sublibraries -------------------------

AddSublibrary <- function(CRISPR_df, sublibrary_list) {
  if (is.null(names(sublibrary_list))) {
    stop("The 'sublibrary_list' parameter must be a named list!")
  }
  sublibrary_vec <- rep(NA_character_, CRISPR_df)
  are_eligible <- (CRISPR_df[["Is_control"]] == "No") &
                  (!(is.na(CRISPR_df[["Entrez_ID"]])))
  for (sublibrary in names(sublibrary_list)) {
    are_this_library <- CRISPR_df[["Entrez_ID"]][are_eligible] %in% sublibrary_list[[sublibrary]]
    sublibrary_vec[are_eligible][are_this_library] <- sublibrary
  }
  CRISPR_df[["Sublibrary_4sg"]] <- factor(sublibrary_vec, levels = names(sublibrary_list))
  return(sublibrary_vec)
}






# Functions for picking control guides ------------------------------------

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


CheckControlsNumber <- function(are_selected, num_controls) {
  if (sum(are_selected) < num_controls) {
    stop(paste0("Not enough controls were available in the supplied data frame",
                "to fill the desired quota of ", num_controls, " control guides!"
                )
         )
  }
}


RandomlyPickControls <- function(CRISPR_df, num_controls = NULL) {
  set.seed(1)
  if (is.null(num_controls)) {
    are_not_controls <- CRISPR_df[["Is_control"]] == "No"
    num_controls <- round(sum(are_not_controls) * 0.01)
    message(paste0(num_controls, " guides will be randomly selected."))
  }
  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
  are_good_controls <- AreGoodControls(CRISPR_df)
  if (is_CRISPRko) {
    are_Brunello <- grepl("Brunello", CRISPR_df[["Source"]], fixed = TRUE)
    are_selected <- are_good_controls & are_Brunello
    CheckControlsNumber(are_selected, num_controls)
    chosen_indices <- sample(which(are_selected), num_controls)
  } else {
    are_Doench <- grepl("Calabrese|Dolcetto", CRISPR_df[["Source"]])
    are_hCRISPR_v2 <- grepl("hCRISPR", CRISPR_df[["Source"]], fixed = TRUE)
    are_Doench_controls <- are_Doench & are_good_controls
    are_hCRISPR_v2_controls <- are_hCRISPR_v2 & are_good_controls
    num_to_select_Doench <- min(sum(are_Doench_controls), round(num_controls / 2))
    num_to_select_hCRISPR_v2 <- num_controls - num_to_select_Doench
    CheckControlsNumber(are_hCRISPR_v2_controls, num_to_select_hCRISPR_v2)
    indices_Doench <- sample(which(are_Doench), num_to_select_Doench)
    indices_hCRISPR_v2 <- sample(which(are_hCRISPR_v2), num_to_select_hCRISPR_v2)
    chosen_indices <- sample(c(indices_Doench, indices_hCRISPR_v2))
  }
  results_df <- CRISPR_df[chosen_indices, ]
  row.names(results_df) <- NULL
  return(results_df)
}




# Functions for assigning guides to 384-well plates -----------------------

AssignToPlates <- function(CRISPR_df, randomize_order = TRUE) {
  if (randomize_wells) {
    set.seed(1)
  }
  use_sublibrary <- "Sublibrary_4sg" %in% CRISPR_df
  if (use_sublibrary) {
    sublibraries <- levels(CRISPR_df[["Sublibrary_4sg"]])
  } else {
    sublibraries <- "dummy"
  }
  plates_vec <- rep(NA_integer_, nrow(CRISPR_df))
  wells_vec <- rep(NA_integer_, nrow(CRISPR_df))
  for (sublibrary in sublibraries) {
    if (use_sublibrary) {
      are_eligible <- CRISPR_df[["Sublibrary_4sg"]] %in% sublibrary
    } else {
      are_eligible <- CRISPR_df[["Is_control"]] == "No"
    }
    num_guides <- sum(are_eligible)
    guides_seq <- seq_len(num_guides)
    num_plates <- ceiling(num_guides / 384)
    plates_sub_vec <- rep(num_plates, num_plates)[guides_seq]
    wells_sub_vec <- rep(seq_len(384), times = num_plates)[guides_seq]
    if (randomize_order) {
      plates_sub_vec <- sample(plates_sub_vec)
      wells_sub_vec <- sample(wells_sub_vec)
    }
    plates_vec[are_eligible] <- plates_sub_vec
    wells_vec[are_eligible] <- wells_sub_vec
  }
  results_df <- CRISPR_df
  results_df[["Plate"]] <- plates_vec
  results_df[["Well"]] <- wells_vec
  return(results_df)
}





