### 14th April 2020 ###






# Define maps -------------------------------------------------------------

sublibraries_short_names <- c(
  "Transcription factors"             = "TF",
  "GPCRs"                             = "GPCR",
  "Secretome"                         = "Secretome",
  "Membrane Proteins"                 = "Membrane",
  "Kinases/Phosphatases/Drug Targets" = "Kinase",
  "Mitochondria/Trafficking/Motility" = "Mito",
  "Stress/Proteostasis"               = "Stress",
  "Cancer/Apoptosis"                  = "Cancer",
  "Gene Expression"                   = "Expression",
  "Unassigned"                        = "Misc",
  "None"                              = "None"
)




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
  have_complete_guides <- rep.int(are_complete, lengths(rank_splits))

  are_valid_top4 <- are_top4 & are_not_homologous & have_complete_guides
  are_valid_or_only_top4 <- are_valid_top4 | (are_top4 & only_one_TSS)

  have_non_homologous_guides <- rep.int(tapply(are_not_homologous, altTSS_IDs_fac, sum) >= 4, lengths(rank_splits))

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
  have_complete_guides <- rep.int(are_complete, lengths(rank_splits))

  have_non_homologous_guides <- rep.int(tapply(are_not_homologous, combined_IDs_fac, sum) >= 4, lengths(rank_splits))

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





# Functions for reporting problematic genes -------------------------------

ShowProblematicGuides <- function(targeting_CRISPR_df, top4_mat) {

  are_invalid_top4 <- top4_mat[, "Are_chosen_4sg"] & !(top4_mat[, "Have_valid_guides"])

  problematic_df <- targeting_CRISPR_df[are_invalid_top4, ]
  problematic_df <- data.frame(
    problematic_df,
    top4_mat[are_invalid_top4, c("Have_complete_guides", "Have_spaced_guides")],
    "Are_protein_coding" = problematic_df[["Entrez_ID"]] %in% collected_entrez_IDs,
    "Annotation_category" = collected_entrezs_df[["Category"]][match(problematic_df[["Entrez_ID"]], collected_entrezs_df[["Entrez_ID"]])],
    stringsAsFactors = FALSE
  )

  show_columns <- c(
    "Combined_ID",  "Entrez_ID", "Gene_symbol", "Original_symbol",
    "Source",
    "sgRNA_sequence", "PAM", "Cut_location", "Chromosome", "Entrez_chromosome",
    "Num_0MM", "Num_1MM",

    "TSS_number", "Allocated_TSS", "Num_TSSs",
    "TSS_ID", "AltTSS_ID",
    "CRISPOR_3MM_specificity",
    "Rank", "Best_combination_rank",
    "Spacing", "Overlaps_tolerance", "Num_overlaps", "Original_rank",
    "Have_complete_guides", "Have_spaced_guides",
    "Are_protein_coding", "Annotation_category"
  )
  message(paste0("The following ", nrow(problematic_df), " guides were not ",
                 "part of valid 4sg combinations and should be excluded!"
                 )
          )
  print(problematic_df[, intersect(show_columns, colnames(problematic_df))])
  return(invisible(problematic_df))
}




# Functions for splitting genes into sublibraries -------------------------

AddSublibrary <- function(CRISPR_df, sublibrary_list, combine_unassigned = TRUE) {
  if (is.null(names(sublibrary_list))) {
    stop("The 'sublibrary_list' parameter must be a named list!")
  }
  sublibrary_vec <- rep.int(NA_character_, nrow(CRISPR_df))
  are_eligible <- (CRISPR_df[["Is_control"]] == "No") &
                  (!(is.na(CRISPR_df[["Entrez_ID"]])))
  for (sublibrary in names(sublibrary_list)) {
    are_this_library <- CRISPR_df[["Entrez_ID"]][are_eligible] %in% sublibrary_list[[sublibrary]]
    stopifnot(all(is.na(sublibrary_vec[are_eligible][are_this_library])))
    sublibrary_vec[are_eligible][are_this_library] <- sublibrary
  }
  use_levels <- names(sublibrary_list)
  if (combine_unassigned) {
    sublibrary_vec[sublibrary_vec == "None"] <- "Unassigned"
    use_levels <- setdiff(use_levels, "None")
  }
  CRISPR_df[["Sublibrary_4sg"]] <- factor(sublibrary_vec, levels = use_levels)
  return(CRISPR_df)
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
    stop(paste0("Not enough controls were available in the supplied data frame ",
                "to fill the desired quota of ", num_controls, " control guides!"
                )
         )
  }
}



ShuffleProblematic <- function(control_combos_list, are_problematic) {
  stopifnot(all(lengths(control_combos_list) == 4))
  num_additional <- min(length(control_combos_list), sum(are_problematic) * 10L)
  additional_indices <- sample(which(!(are_problematic)), num_additional)
  all_indices <- c(which(are_problematic), additional_indices)
  shuffled_vec <- unlist(control_combos_list[all_indices])
  shuffled_vec <- sample(shuffled_vec)
  new_groups_vec <- rep(seq_len(length(all_indices)), each = 4)
  results_list <- control_combos_list
  results_list[all_indices] <- split(shuffled_vec, new_groups_vec)
  return(results_list)
}


AddRandomized4sgControls <- function(CRISPR_df, num_control_wells = NULL) {
  set.seed(1)
  if (is.null(num_control_wells)) {
    num_control_wells <- 384L * 2L
    message(paste0(num_control_wells * 4, " guides for ", num_control_wells,
                   " control wells will be randomly selected."
                   )
            )
  }
  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
  are_good_controls <- AreGoodControls(CRISPR_df)
  are_duplicated <- duplicated(toupper(CRISPR_df[["sgRNA_sequence"]][are_good_controls]))
  if (any(are_duplicated)) {
    stop("Duplicated control sequences found!")
  }
  if (is_CRISPRko) {
    are_Brunello <- grepl("Brunello", CRISPR_df[["Source"]], fixed = TRUE)
    are_selected <- are_good_controls & are_Brunello
    CheckControlsNumber(are_selected, num_control_wells * 4)
    chosen_sequences <- CRISPR_df[["sgRNA_sequence"]][are_selected]
  } else {
    are_Doench <- grepl("Calabrese|Dolcetto", CRISPR_df[["Source"]])
    are_hCRISPR_v2 <- grepl("hCRISPR", CRISPR_df[["Source"]], fixed = TRUE)
    are_Doench_controls <- are_Doench & are_good_controls
    are_hCRISPR_v2_controls <- are_hCRISPR_v2 & are_good_controls
    num_guides_Doench <- min(sum(are_Doench_controls), round((num_control_wells * 4) / 2))
    num_guides_hCRISPR_v2 <- (num_control_wells * 4) - num_guides_Doench
    CheckControlsNumber(are_hCRISPR_v2_controls, num_guides_hCRISPR_v2)
    indices_Doench <- sample(which(are_Doench_controls), num_guides_Doench)
    indices_hCRISPR_v2 <- sample(which(are_hCRISPR_v2_controls), num_guides_hCRISPR_v2)
    chosen_indices <- sample(c(indices_Doench, indices_hCRISPR_v2))
    chosen_sequences <- CRISPR_df[["sgRNA_sequence"]][chosen_indices]
  }

  stopifnot((length(chosen_sequences) / 4) == num_control_wells)

  chosen_sequences <- toupper(chosen_sequences)
  control_combos_vec <- rep(seq_len(num_control_wells), each = 4L)
  control_combos_vec <- sample(control_combos_vec)
  controls_list <- split(chosen_sequences, control_combos_vec)

  are_problematic <- TRUE
  num_tries <- 0L
  while(any(are_problematic)) {
    are_problematic <- vapply(controls_list, NumHomologousPairs, integer(1)) > 0
    if (any(are_problematic)) {
      num_tries <- num_tries + 1L
      message(paste0("Control guide groupings will be re-shuffled to resolve ",
                     "any remaining combinations that share subsequences ",
                     "longer than 8 bp. Round ", num_tries, "..."
                     )
              )
      controls_list <- ShuffleProblematic(controls_list, are_problematic)
    }
  }
  message("Control guide groupings are now unproblematic!")

  sequences_vec <- unlist(controls_list, use.names = FALSE)

  groups_vec <- rep(seq_along(controls_list), each = 4)

  matches_vec <- match(sequences_vec, toupper(CRISPR_df[["sgRNA_sequence"]][are_good_controls]))
  stopifnot(!(anyNA(matches_vec)))

  results_df <- CRISPR_df
  results_df[["Control_group_4sg"]] <- NA_integer_
  results_df[["Control_group_4sg"]][are_good_controls][matches_vec] <- groups_vec

  return(results_df)
}



AssignControlsToPlates <- function(CRISPR_df, num_wells_per_plate = 384L) {

  stopifnot("Control_group_4sg" %in% colnames(CRISPR_df))

  control_IDs_seq <- CRISPR_df[["Control_group_4sg"]]

  are_control_guides <- !(is.na(control_IDs_seq))

  num_control_wells <- sum(are_control_guides) / 4
  wells_seq <- seq_len(num_control_wells)

  stopifnot(identical(wells_seq, sort(unique(control_IDs_seq[are_control_guides]))))

  num_plates <- ceiling(num_control_wells / num_wells_per_plate)
  plates_vec <- rep(seq_len(num_plates), each = num_wells_per_plate)[wells_seq]
  wells_vec <- rep.int(seq_len(num_wells_per_plate), times = num_plates)[wells_seq]

  results_df <- CRISPR_df
  for (column_name in c("Sublibrary_4sg", "Plate_ID", "Plate_number", "Well_number")) {
    if (!(column_name %in% names(results_df))) {
      results_df[[column_name]] <- NA_integer_
    }
  }
  control_matches <- match(control_IDs_seq[are_control_guides], wells_seq)
  results_df[["Plate_number"]][are_control_guides] <- plates_vec[control_matches]
  results_df[["Well_number"]][are_control_guides] <- wells_vec[control_matches]
  results_df[["Plate_ID"]][are_control_guides] <- paste0("Control_", results_df[["Plate_number"]][are_control_guides])
  return(results_df)
}




# Functions for assigning guides to 384-well plates -----------------------

AssignToPlates <- function(CRISPR_df, randomize_order = TRUE, num_wells_per_plate = 384L) {
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  stopifnot(length(unique(table(CRISPR_df[["Rank"]]))) == 1)
  if (randomize_order) {
    set.seed(1)
  }
  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
  if (is_CRISPRko) {
    ID_column <- "Combined_ID"
  } else {
    ID_column <- "AltTSS_ID"
  }
  stopifnot(all(table(CRISPR_df[[ID_column]]) == 4))
  use_sublibraries <- "Sublibrary_4sg" %in% colnames(CRISPR_df)
  if (use_sublibraries) {
    sublibraries <- levels(CRISPR_df[["Sublibrary_4sg"]])
    have_no_sublibrary <- is.na(CRISPR_df[["Sublibrary_4sg"]])
    if (any(have_no_sublibrary)) {
      warning(paste0(sum(have_no_sublibrary), " guides did not belong to a ",
                     "sublibrary and could not be allocated to a plate!"
                     )
              )
    }
  } else {
    sublibraries <- "dummy"
  }
  plates_vec <- rep.int(NA_integer_, nrow(CRISPR_df))
  wells_vec <- rep.int(NA_integer_, nrow(CRISPR_df))
  for (sublibrary in sublibraries) {
    if (use_sublibraries) {
      are_eligible <- as.character(CRISPR_df[["Sublibrary_4sg"]]) %in% sublibrary
    } else {
      are_eligible <- CRISPR_df[["Is_control"]] == "No"
    }
    are_eligible <- are_eligible & (CRISPR_df[["Rank"]] %in% 1:4)
    CheckThatFactorIsInOrder(CRISPR_df[["ID_column"]][are_eligible])
    num_wells <- sum(are_eligible) / 4
    wells_seq <- seq_len(num_wells)
    num_plates <- ceiling(num_wells / num_wells_per_plate)
    plates_sub_vec <- rep(seq_len(num_plates), each = num_wells_per_plate)[wells_seq]
    if (randomize_order) {
      wells_sub_vec <- rep(NA_integer_, num_wells)
      plates_sub_vec <- sample(plates_sub_vec)
      for (plate_number in seq_len(num_plates)) {
        are_this_plate <- plates_sub_vec == plate_number
        wells_sub_vec[are_this_plate] <- sample(seq_len(sum(are_this_plate)))
      }
    } else {
      wells_sub_vec <- rep.int(seq_len(num_wells_per_plate), times = num_plates)[wells_seq]
    }
    plates_vec[are_eligible] <- rep(plates_sub_vec, each = 4)
    wells_vec[are_eligible] <- rep(wells_sub_vec, each = 4)
  }
  results_df <- CRISPR_df
  if (use_sublibraries) {
    assign("delete_CRISPR_df_33", CRISPR_df, envir = globalenv())
    short_sublibraries_fac <- CRISPR_df[["Sublibrary_4sg"]]
    levels(short_sublibraries_fac) <- sublibraries_short_names[levels(short_sublibraries_fac)]
    results_df[["Plate_ID"]] <- paste0(as.character(short_sublibraries_fac), "_", plates_vec)
  } else {
    results_df[["Plate_ID"]] <- plates_vec
  }
  results_df[["Plate_number"]] <- plates_vec
  results_df[["Well_number"]] <- wells_vec
  return(results_df)
}




# Functions that encapsulate the workflow for export ----------------------

AssignAllGuides <- function(CRISPR_df, sublibraries_entrezs_list, num_control_wells = NULL, reorder_df = FALSE) {

  CRISPR_df <- AddRandomized4sgControls(CRISPR_df, num_control_wells = num_control_wells)
  CRISPR_df <- AssignControlsToPlates(CRISPR_df)

  controls_df <- CRISPR_df[!(is.na(CRISPR_df[["Control_group_4sg"]])), ]
  controls_df <- AssignControlsToPlates(controls_df)

  targeting_df <- CRISPR_df[CRISPR_df[["Entrez_ID"]] %in% unlist(sublibraries_entrezs_list, use.names = FALSE), ]

  are_top4_mat <- CRISPRaAreTop4Mat(targeting_df)

  ShowProblematicGuides(targeting_df, are_top4_mat)

  are_valid_chosen <- are_top4_mat[, "Are_chosen_4sg"] & are_top4_mat[, "Have_valid_guides"]

  targeting_df <- targeting_df[are_valid_chosen, ]

  targeting_df <- AddSublibrary(targeting_df, sublibraries_entrezs_list)
  targeting_df <- AssignToPlates(targeting_df)

  results_df <- rbind.data.frame(targeting_df,
                                 controls_df,
                                 make.row.names = FALSE,
                                 stringsAsFactors = FALSE
                                 )

  plate_IDs_vec <- results_df[["Plate_ID"]]
  new_order <- order(results_df[["Sublibrary_4sg"]],
                     results_df[["Plate_number"]],
                     results_df[["Well_number"]],
                     results_df[["Rank"]]
                     )
  plate_IDs_vec <- plate_IDs_vec[new_order]
  results_df[["Plate_ID"]] <- factor(results_df[["Plate_ID"]], levels = unique(plate_IDs_vec))
  if (reorder_df) {
    results_df <- results_df[new_order, ]
    are_controls <- results_df[["Is_control"]] == "Yes"
    results_df[["Rank"]][are_controls] <- rep.int(1:4, sum(are_controls) / 4)
    row.names(results_df) <- NULL
  }
  return(results_df)
}















