### 14th April 2020 ###






# Define maps -------------------------------------------------------------

sublibraries_short_names <- c(
  "Controls & changed TFs"            = "ChangedTF",
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


export_columns <- c(
  "Sublibrary_4sg", "Plate_string", "Well_number",
  "Is_obsolete",
  "Entrez_ID", "Other_target_Entrez_IDs", "Other_Entrez_IDs_4sg",
  "Gene_symbol", "Other_target_symbols", "Other_symbols_4sg", "Original_symbol",
  "Exon_number_Brunello", "Exon_number_TKOv3", "Exon_number_GPP",
  "Transcript_ID", "Genomic_sequence_ID",
  "TSS_ID", "Is_main_TSS",
  "Rank", "Num_overlaps", "Source",
  "sgRNA_sequence", "PAM", "Sequence_with_primers",
  "Calabrese_rank", "Dolcetto_rank",
  "Caprano_rank", "Dolomiti_rank",
  "GPP_rank",
  "hCRISPRa_v2_rank", "hCRISPRi_v2_rank",
  "mCRISPRa_v2_rank", "mCRISPRi_v2_rank",
  "Predicted_score", "Empirical_score",
  "Chromosome", "Strand", "Cut_location", "Distance_from_TSS",
  "GuideScan_efficiency", "CRISPOR_Doench_efficacy",
  "CRISPOR_Graf_status",
  "GuideScan_specificity",
  "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity", "CRISPOR_CFD_specificity",
  "Num_0MM", "Num_1MM", "GuideScan_Num_2MM", "GuideScan_Num_3MM",
  "GuideScan_offtarget_category", "CRISPOR_Num_0MM", "CRISPOR_Num_1MM",
  "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
  "all22_SNP_IDs_vcf", "all22_SNP_AF_max_Kaviar",
  "Locations_0MM", "Locations_1MM", "Sequences_1MM"
)




# Functions for selecting the 4sg guides to export ------------------------

ReportProblematic <- function(are_top4_mat, CRISPR_df) {
  if (any(are_top4_mat[, "Are_problematic_top4"])) {
    problematic_IDs <- unique(CRISPR_df[["Combined_ID"]][are_top4_mat[, "Are_problematic_top4"]])
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
    "Are_chosen_4sg"             = are_valid_or_only_top4,
    "Are_top4"                   = are_top4,
    "Are_problematic_top4"       = are_valid_or_only_top4 & !(are_valid_top4),
    "Have_complete_guides"       = have_complete_guides,
    "Have_non_homologous_guides" = have_non_homologous_guides,
    "Have_valid_guides"          = have_non_homologous_guides & have_complete_guides
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
    "Are_problematic_top4"       = are_top4 & !(are_not_homologous & have_complete_guides),
    "Have_complete_guides"       = have_complete_guides,
    "Have_non_homologous_guides" = have_non_homologous_guides,
    "Have_valid_guides"          = have_non_homologous_guides & have_complete_guides
  )
  ReportProblematic(results_mat, CRISPRko_df)
  return(results_mat)
}




Are4sg <- function(CRISPR_df, sublibraries_entrezs_list, show_messages = TRUE) {
  are_targeting <- CRISPR_df[["Entrez_ID"]] %in% unlist(sublibraries_entrezs_list, use.names = FALSE)
  targeting_df <- CRISPR_df[are_targeting, ]
  if ("Num_TSSs" %in% names(CRISPR_df)) {
    are_top4_mat <- CRISPRaAreTop4Mat(targeting_df)
  } else {
    are_top4_mat <- CRISPRkoAreTop4Mat(targeting_df)
  }
  ShowProblematicGuides(targeting_df, are_top4_mat, show_messages = show_messages)
  are_4sg <- are_targeting
  are_4sg[are_targeting] <- are_top4_mat[, "Are_chosen_4sg"] & are_top4_mat[, "Have_valid_guides"]
  return(are_4sg)
}



# Functions for reporting problematic genes -------------------------------

ShowProblematicGuides <- function(targeting_CRISPR_df, top4_mat, show_messages = TRUE) {

  are_invalid_top4 <- top4_mat[, "Are_chosen_4sg"] & !(top4_mat[, "Have_valid_guides"])

  if (any(are_invalid_top4)) {
    problematic_df <- targeting_CRISPR_df[are_invalid_top4, ]
    problematic_df <- data.frame(
      problematic_df,
      top4_mat[are_invalid_top4, c("Have_complete_guides", "Have_non_homologous_guides")],
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
      "Have_complete_guides", "Have_non_homologous_guides",
      "Are_protein_coding", "Annotation_category"
    )
    if (show_messages) {
      message(paste0("\nThe following ", nrow(problematic_df), " guides were not ",
                     "part of valid 4sg combinations and should be excluded!"
                     )
              )
      print(problematic_df[, intersect(show_columns, names(problematic_df))])
    }
  } else {
    if (show_messages) {
      message("Valid guides were found for all genes in the sgRNA dataframe!")
    }
    problematic_df <- NULL
  }

  are_valid_top4 <- top4_mat[, "Are_chosen_4sg"] & top4_mat[, "Have_valid_guides"]
  num_unique_genes <- length(unique(targeting_CRISPR_df[["Entrez_ID"]][are_valid_top4]))
  num_combos <- sum(are_valid_top4) / 4
  show_message <- paste0("\nThe 4sg library targets ", num_unique_genes,
                         " unique genes"
                         )

  is_CRISPRko <- "Exon_number_GPP" %in% names(targeting_CRISPR_df)
  if (is_CRISPRko) {
    stopifnot(num_unique_genes == num_combos)
  } else {
    num_unique_TSSs <- length(unique(targeting_CRISPR_df[["AltTSS_ID"]][are_valid_top4]))
    stopifnot(num_unique_TSSs == num_combos)
    show_message <- paste0(show_message, " and ", num_unique_TSSs, " unique ",
                           "transcription start sites (TSSs)"
                           )
  }
  show_message <- paste0(show_message, "!\n")
  if (show_messages) {
    message(show_message)
  }
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



AddRandomized4sgControls <- function(CRISPR_df, num_control_wells = NULL, previously_used_controls = NULL, num_wells_per_plate = 384L) {
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  assign("delete_previously_used_controls", previously_used_controls, envir = globalenv())
  set.seed(1)
  if (is.null(num_control_wells)) {
    num_control_wells <- num_wells_per_plate * 2L
    message(paste0(num_control_wells * 4, " guides for ", num_control_wells,
                   " control wells will be randomly selected."
                   )
            )
  }
  is_CRISPRko <- "Exon_number_GPP" %in% names(CRISPR_df)
  are_good_controls <- AreGoodControls(CRISPR_df)
  are_duplicated <- duplicated(toupper(CRISPR_df[["sgRNA_sequence"]][are_good_controls]))
  if (any(are_duplicated)) {
    stop("Duplicated control sequences found!")
  }
  if (!(is.null(previously_used_controls))) {
    are_previously_used <- toupper(CRISPR_df[["sgRNA_sequence"]][are_good_controls]) %in% toupper(previously_used_controls)
    if (any(are_previously_used)) {
      message(paste0(sum(are_previously_used), " previously used control ",
                     "guide sequences were excluded."
                     )
              )
    }
    are_good_controls[are_previously_used] <- FALSE
  }
  if (is_CRISPRko) {
    are_Doench <- grepl("Brunello|Brie", CRISPR_df[["Source"]])
    are_selected <- are_good_controls & are_Doench
    CheckControlsNumber(are_selected, num_control_wells * 4)
    chosen_sequences <- sample(CRISPR_df[["sgRNA_sequence"]][are_selected], num_control_wells * 4)
  } else {
    are_Doench <- grepl("Calabrese|Dolcetto|Caprano|Dolomiti", CRISPR_df[["Source"]])
    are_hCRISPR_v2 <- grepl("[hm]CRISPR", CRISPR_df[["Source"]])
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
                     "8 bp or longer. Round ", num_tries, "..."
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


RenameControls <- function(CRISPR_df) {
  stopifnot("Control_group_4sg" %in% names(CRISPR_df))
  are_controls <- CRISPR_df[["Is_control"]] == "Yes"
  are_chosen <- !(is.na(CRISPR_df[["Control_group_4sg"]]))
  CRISPR_df[["Combined_ID"]] <- ifelse(are_controls,
                                       ifelse(are_chosen, paste0("Control_", CRISPR_df[["Control_group_4sg"]]), "Control"),
                                       CRISPR_df[["Combined_ID"]]
                                       )
  if ("AltTSS_ID" %in% names(CRISPR_df)) {
    CRISPR_df[["AltTSS_ID"]] <- ifelse(are_controls,
                                       CRISPR_df[["Combined_ID"]],
                                       CRISPR_df[["AltTSS_ID"]]
                                       )
  }
  return(CRISPR_df)
}



AssignControlsToPlates <- function(CRISPR_df, num_wells_per_plate = 384L) {
  stopifnot("Control_group_4sg" %in% names(CRISPR_df))

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
  control_indices <- unlist(lapply(wells_seq, function(x) which(control_IDs_seq[are_control_guides] == x)))
  results_df[["Plate_number"]][are_control_guides][control_indices] <- rep(plates_vec, each = 4)
  results_df[["Well_number"]][are_control_guides][control_indices] <- rep(wells_vec, each = 4)
  results_df[["Plate_ID"]][are_control_guides] <- paste0("Control_", results_df[["Plate_number"]][are_control_guides])
  results_df[["Sublibrary_4sg"]][are_control_guides] <- "Controls"
  return(results_df)
}


ShuffleControlRanks <- function(CRISPR_df) {
  set.seed(1)
  control_groups_vec <- CRISPR_df[["Control_group_4sg"]]
  are_controls <- !(is.na(control_groups_vec))
  control_groups_vec <- control_groups_vec[are_controls]
  stopifnot(all(table(control_groups_vec) == 4))
  unique_groups <- unique(control_groups_vec)
  ranks_vec <- rep(NA_integer_, length(control_groups_vec))
  for (group in unique_groups) {
    are_this_group <- control_groups_vec == group
    ranks_vec[are_this_group] <- sample(1:4)
  }
  CRISPR_df[["Rank"]][are_controls] <- ranks_vec
  return(CRISPR_df)
}




# Functions for assigning guides to 384-well plates -----------------------

AssignToPlates <- function(CRISPR_df, randomize_order = TRUE, num_wells_per_plate = 384L) {
  stopifnot(length(unique(table(CRISPR_df[["Rank"]]))) == 1)
  if (randomize_order) {
    set.seed(1)
  }
  is_CRISPRko <- "Exon_number_GPP" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    ID_column <- "Combined_ID"
  } else {
    ID_column <- "AltTSS_ID"
  }
  stopifnot(all(table(CRISPR_df[[ID_column]]) == 4))
  use_sublibraries <- "Sublibrary_4sg" %in% names(CRISPR_df)
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



RenumberPlatesContinuously <- function(CRISPR_df,
                                       num_wells_per_plate = 384L,
                                       start_at_plate_number = 1L
                                       ) {
  gene_or_TSS_IDs <- GetGeneOrTSSIDs(CRISPR_df)
  stopifnot(all(table(gene_or_TSS_IDs) == 4))
  total_num_wells <- nrow(CRISPR_df) / 4
  total_num_plates <- ceiling(total_num_wells / num_wells_per_plate)
  wells_seq <- seq_len(num_wells_per_plate)
  plates_seq <- seq(from = start_at_plate_number,
                    to = start_at_plate_number + total_num_plates - 1L
                    )
  total_seq <- seq_len(total_num_wells)
  plates_vec <- rep(plates_seq, each = num_wells_per_plate)[total_seq]
  wells_vec <- rep(wells_seq, times = total_num_plates)[total_seq]

  results_df <- CRISPR_df
  results_df[["Plate_number"]] <- rep(plates_vec, each = 4)
  results_df[["Well_number"]] <- rep(wells_vec, each = 4)
  results_df[["Plate_ID"]] <- NA_character_
  return(results_df)
}




PlaceCandidateGenesTogether <- function(CRISPR_df, candidate_entrezs) {

  CRISPR_df[["Is_candidate_gene"]] <- CRISPR_df[["Entrez_ID"]] %in% candidate_entrezs
  sublibrary_df_list <- split(CRISPR_df, CRISPR_df[["Sublibrary_4sg"]])

  sublibrary_reordered_df_list <- lapply(sublibrary_df_list, function(sub_df) {
    are_candidates <- sub_df[["Entrez_ID"]] %in% candidate_entrezs
    if (any(are_candidates)) {
      are_well_1 <- sub_df[["Well_number"]] %in% 1
      first_well_index <- which(are_well_1)[[1]]

      plate_layout_columns <- c("Plate_ID", "Plate_number", "Well_number")
      other_columns <- setdiff(names(CRISPR_df), plate_layout_columns)
      plate_layout_df <- sub_df[, plate_layout_columns]
      data_df <- sub_df[, other_columns]

      candidates_df <- data_df[are_candidates, ]
      not_candidates_df <- data_df[!(are_candidates), ]

      ## Notice that this assumes that there are not very many candidate genes,
      ## so that we can start them on a "fresh" plate and still have space!
      if (first_well_index == 1) {
        preceding_indices <- c()
        preceding_df <- NULL
      } else {
        preceding_indices <- seq_len(first_well_index - 1L)
        preceding_df <- not_candidates_df[preceding_indices, ]
      }
      following_indices <- setdiff(seq_len(nrow(not_candidates_df)), preceding_indices)
      following_df <- not_candidates_df[following_indices, ]

      recombined_data_df <- rbind.data.frame(preceding_df,
                                             candidates_df,
                                             following_df,
                                             stringsAsFactors = FALSE,
                                             make.row.names = FALSE
                                             )
      recombined_df <- data.frame(recombined_data_df, plate_layout_df)
      recombined_df <- recombined_df[, names(sub_df)]
      return(recombined_df)
    } else {
      return(sub_df)
    }
  })

  reordered_df <- do.call(rbind.data.frame,
                          c(sublibrary_reordered_df_list,
                            list(stringsAsFactors = FALSE,
                                 make.row.names = FALSE
                                 )
                            )
                          )

  return(reordered_df)
}





# Functions that encapsulate the workflow for export ----------------------

MakeMiscPlate <- function(CRISPR_df, num_wells_per_plate = 384L) {

  sublibrary_vec <- as.character(CRISPR_df[["Sublibrary_4sg"]])
  are_TF_or_controls <- sublibrary_vec == "Controls & changed TFs"
  are_controls <- CRISPR_df[["Is_control"]][are_TF_or_controls] == "Yes"
  sublibrary_vec[are_TF_or_controls][are_controls] <- "Misc / controls"
  sublibrary_vec[are_TF_or_controls][!(are_controls)] <- "Misc / changed TFs"

  first_plate_number <- unique(CRISPR_df[["Plate_number"]][are_TF_or_controls])
  stopifnot(length(first_plate_number) == 1)
  last_plate_number <- max(CRISPR_df[["Plate_number"]])
  are_last_plate <- CRISPR_df[["Plate_number"]] == last_plate_number

  last_plate_is_misc <- ((sum(are_TF_or_controls) + sum(are_last_plate)) / 4) <= num_wells_per_plate

  results_df <- CRISPR_df
  results_df[["Sublibrary_4sg"]] <- sublibrary_vec

  if (last_plate_is_misc) {

    are_other_plates <- !(are_TF_or_controls | are_last_plate)

    stopifnot(all(results_df[["Sublibrary_4sg"]][are_last_plate] == "Unassigned"))
    results_df[["Sublibrary_4sg"]][are_last_plate] <- "Misc / unassigned"

    new_first_plate_df <- rbind.data.frame(results_df[are_TF_or_controls, ],
                                           results_df[are_last_plate, ],
                                           make.row.names = FALSE,
                                           stringsAsFactors = FALSE
                                           )
    new_first_plate_df <- RenumberPlatesContinuously(new_first_plate_df)
    new_first_plate_df[["Plate_ID"]] <- NULL
    results_df <- rbind.data.frame(new_first_plate_df,
                                   results_df[are_other_plates, ],
                                   make.row.names = FALSE,
                                   stringsAsFactors = FALSE
                                   )
  }
  return(results_df)
}



AssignPlateStrings <- function(CRISPR_df, use_prefix = "h") {
  unique_plate_numbers <- sort(unique(CRISPR_df[["Plate_number"]]))
  plates_vec <- as.character(CRISPR_df[["Plate_number"]])
  if ((length(unique_plate_numbers) > 1) &&
      (unique_plate_numbers[[1]] == 1) &&
      (unique_plate_numbers[[2]] == 6)
      ) {
    plates_vec[CRISPR_df[["Plate_number"]] == 1] <- "5+"
  }
  if ("Exon_number_GPP" %in% names(CRISPR_df)) {
    modality_string <- "o"
  } else if (("Entrez_source_Calabrese" %in% names(CRISPR_df)) || ("Entrez_source_Caprano" %in% names(CRISPR_df))) {
    modality_string <- "a"
  } else if (("Entrez_source_Dolcetto" %in% names(CRISPR_df)) || ("Entrez_source_Dolomiti" %in% names(CRISPR_df))) {
    modality_string <- "i"
  }
  plates_vec <- paste0(use_prefix, modality_string, "_", plates_vec,
                       "_sg", CRISPR_df[["Rank"]]
                       )
  CRISPR_df[["Plate_string"]] <- plates_vec
  return(CRISPR_df)
}





AllocateAllGuides_v2 <- function(CRISPR_df,
                                 sublibraries_entrezs_list,
                                 previous_version_CRISPR_df,
                                 candidate_entrezs,
                                 num_control_wells        = 96L,
                                 num_wells_per_plate      = 384L
                                 ) {

  ## Preserve the original order

  CRISPR_df[["Original_index"]] <- seq_len(nrow(CRISPR_df))


  ## Choose control guides

  are_previous_controls <- previous_version_CRISPR_df[["Is_control"]] == "Yes"
  previous_control_sequences <- previous_version_CRISPR_df[["sgRNA_sequence"]][are_previous_controls]

  CRISPR_df <- AddRandomized4sgControls(CRISPR_df,
                                        num_control_wells = 96,
                                        previously_used_controls = previous_control_sequences,
                                        num_wells_per_plate = num_wells_per_plate
                                        )
  controls_df <- CRISPR_df[!(is.na(CRISPR_df[["Control_group_4sg"]])), ]
  controls_df <- AssignControlsToPlates(controls_df, num_wells_per_plate = num_wells_per_plate)
  controls_df <- RenameControls(controls_df)
  controls_df <- ShuffleControlRanks(controls_df)
  controls_df <- ReorderPlates(controls_df)
  controls_df[["Is_candidate_gene"]] <- FALSE


  ## Choose targeting guides

  targeting_df <- CRISPR_df[CRISPR_df[["Entrez_ID"]] %in% unlist(sublibraries_entrezs_list, use.names = FALSE), ]
  message("")
  is_CRISPRko <- "Exon_number_GPP" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    are_top4_mat <- CRISPRkoAreTop4Mat(targeting_df)
  } else {
    are_top4_mat <- CRISPRaAreTop4Mat(targeting_df)
  }
  ShowProblematicGuides(targeting_df, are_top4_mat)

  are_valid_chosen <- are_top4_mat[, "Are_chosen_4sg"] & are_top4_mat[, "Have_valid_guides"]
  targeting_df <- targeting_df[are_valid_chosen, ]
  targeting_df <- AddSublibrary(targeting_df, sublibraries_entrezs_list)
  targeting_df <- AssignToPlates(targeting_df, num_wells_per_plate = num_wells_per_plate)
  targeting_df <- ReorderPlates(targeting_df)


  ## Find genes that change (compared to version 1.0 of the TF library)

  are_TFs <- targeting_df[["Sublibrary_4sg"]] %in% "Transcription factors"
  new_TF_df <- targeting_df[are_TFs, ]
  old_TF_df <- previous_version_CRISPR_df[!(are_previous_controls), ]

  new_TF_IDs <- Get4sgIDs(new_TF_df)
  old_TF_IDs <- Get4sgIDs(old_TF_df)

  changed_new_TF_IDs <- setdiff(new_TF_IDs, old_TF_IDs)
  changed_entrez_IDs <- unique(new_TF_df[["Entrez_ID"]][new_TF_IDs %in% changed_new_TF_IDs])
  changed_TF_df <- new_TF_df[new_TF_df[["Entrez_ID"]] %in% changed_entrez_IDs, ]

  changed_TF_df <- RenumberPlatesContinuously(changed_TF_df)
  reordered_changed_TF_df <- PlaceCandidateGenesTogether(changed_TF_df, candidate_entrezs)


  ## Define the first plate ("5+")

  assign("delete_controls_df", controls_df, envir = globalenv())
  assign("delete_reordered_changed_TF_df", reordered_changed_TF_df, envir = globalenv())

  first_plate_df <- rbind.data.frame(controls_df,
                                     reordered_changed_TF_df,
                                     stringsAsFactors = FALSE,
                                     make.row.names = FALSE
                                     )
  stopifnot((nrow(first_plate_df) / 4) < num_wells_per_plate)


  ## Define the other plates

  other_plates_df <- targeting_df[!(are_TFs), ]


  ## Combine all the plates

  levels(other_plates_df[["Sublibrary_4sg"]])[[1]] <- "Controls & changed TFs"
  first_plate_df[["Sublibrary_4sg"]] <- "Controls & changed TFs"
  first_plate_df[["Sublibrary_4sg"]] <- factor(first_plate_df[["Sublibrary_4sg"]],
                                               levels = levels(other_plates_df[["Sublibrary_4sg"]])
                                               )
  assign("delete_first_plate_df", first_plate_df, envir = globalenv())
  first_plate_df <- RenumberPlatesContinuously(first_plate_df)
  other_plates_df <- RenumberPlatesContinuously(other_plates_df, start_at_plate_number = 6)

  reordered_other_plates_df <- PlaceCandidateGenesTogether(other_plates_df, candidate_entrezs)

  all_plates_df <- rbind.data.frame(first_plate_df,
                                    reordered_other_plates_df,
                                    stringsAsFactors = FALSE,
                                    make.row.names = FALSE
                                    )

  ## Check for shared subsequences 8bp or more in length
  all_plates_df[["Plate_ID"]] <- all_plates_df[["Plate_number"]]
  CheckForSharedSubsequences(all_plates_df)

  all_plates_df[["Plate_ID"]] <- NULL

  all_plates_df <- MakeMiscPlate(all_plates_df)
  all_plates_df <- AssignPlateStrings(all_plates_df)

  return(all_plates_df)
}





AllocateAllGuidesToPlates <- function(CRISPR_df,
                                      sublibraries_entrezs_list,
                                      num_control_wells = NULL,
                                      reorder_df = FALSE
                                      ) {

  if (num_control_wells == 0) {
    CRISPR_df <- CRISPR_df[CRISPR_df[["Is_control"]] == "No", ]
    controls_df <- NULL
  } else {
    CRISPR_df <- AddRandomized4sgControls(CRISPR_df, num_control_wells = num_control_wells)

    controls_df <- CRISPR_df[!(is.na(CRISPR_df[["Control_group_4sg"]])), ]
    controls_df <- AssignControlsToPlates(controls_df)
    controls_df <- RenameControls(controls_df)
    controls_df <- ShuffleControlRanks(controls_df)
  }

  targeting_df <- CRISPR_df[CRISPR_df[["Entrez_ID"]] %in% unlist(sublibraries_entrezs_list, use.names = FALSE), ]

  is_CRISPRko <- "Exon_number_GPP" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    are_top4_mat <- CRISPRkoAreTop4Mat(targeting_df)
  } else {
    are_top4_mat <- CRISPRaAreTop4Mat(targeting_df)
  }

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
                     results_df[["Plate_number"]]
                     )
  plate_IDs_vec <- plate_IDs_vec[new_order]
  results_df[["Plate_ID"]] <- factor(results_df[["Plate_ID"]], levels = unique(plate_IDs_vec))

  CheckForSharedSubsequences(results_df, check_gene_IDs = FALSE)
  if (reorder_df) {
    results_df <- ReorderPlates(results_df)
  }
  return(results_df)
}




CheckForSharedSubsequences <- function(CRISPR_df, check_gene_IDs = TRUE) {
  assign("delete_CRISPR_df_sss", CRISPR_df, envir = globalenv())

  unique_IDs_vec <- paste0(as.character(CRISPR_df[["Plate_ID"]]), "__", CRISPR_df[["Well_number"]])
  unique_IDs_fac <- factor(unique_IDs_vec, levels = unique(unique_IDs_vec))

  sequences_vec <- toupper(CRISPR_df[["sgRNA_sequence"]])
  sequences_list <- split(sequences_vec, unique_IDs_fac)

  if (check_gene_IDs) {
    gene_or_TSS_fac <- GetGeneOrTSSIDs(CRISPR_df)
    stopifnot(identical(unname(sequences_list), unname(split(sequences_vec, gene_or_TSS_fac))))
  }

  shared_subsequence_lengths <- vapply(sequences_list, LongestSharedSubsequence, integer(1))
  if (any(shared_subsequence_lengths >= 8)) {
    stop("Invalid 4sg combination found!")
  } else {
    message(paste0("The library did not contain any 4sg combinations that ",
                   "shared subsequences >7bp in length."
                   )
            )
  }
  return(invisible(NULL))
}



ReorderPlates <- function(CRISPR_df) {
  new_order <- order(CRISPR_df[["Sublibrary_4sg"]],
                     CRISPR_df[["Plate_number"]],
                     CRISPR_df[["Well_number"]],
                     CRISPR_df[["Rank"]]
                     )
  CRISPR_df <- CRISPR_df[new_order, ]
  row.names(CRISPR_df) <- NULL
  return(CRISPR_df)
}




RestoreOriginalOrder <- function(CRISPR_df) {
  are_controls <- CRISPR_df[["Is_control"]] %in% "Yes"
  new_order <- order(ifelse(are_controls, NA, CRISPR_df[["Original_index"]]))
  results_df <- CRISPR_df[new_order, ]
  row.names(results_df) <- NULL
  return(results_df)
}




# Functions for exporting the plate layouts -------------------------------

ExportPlates <- function(export_df,
                         file_name,
                         sub_folder,
                         add_padding_between_plates = FALSE,
                         add_primers = TRUE,
                         add_colors  = TRUE,
                         no_modality = FALSE,
                         add_remove_columns = c()
                         ) {
  export_df[["Source"]] <- sub("Curated, ", "", export_df[["Source"]], fixed = TRUE)
  use_columns <- intersect(export_columns, names(export_df))
  remove_columns <- setdiff(names(export_df), use_columns)
  remove_columns <- union(remove_columns, add_remove_columns)
  is_CRISPRko <- "Exon_number_GPP" %in% names(export_df)
  if (!(is_CRISPRko)) {
    remove_columns <- c(remove_columns, "CRISPOR_Graf_status")
  }
  if (length(unique(export_df[["Rank"]])) == 1) {
    remove_columns <- c(remove_columns, "Rank")
  }
  DfToTSV(export_df[, union(use_columns, names(export_df))],
          file_name                  = file.path(sub_folder, file_name),
          add_primers                = add_primers,
          add_colors                 = add_colors,
          remove_columns             = remove_columns,
          add_padding_between_plates = add_padding_between_plates,
          no_modality                = no_modality
          )
}





# Functions for merging the TF sub-library with the rest ------------------

MergeTFWithRest <- function(sg4_by_well_df, TF_by_well_df) {

  sg4_by_gene_df <- RestoreOriginalOrder(sg4_by_well_df)


  ## Add plate strings
  TF_by_well_df <- AssignPlateStrings(TF_by_well_df)
  TF_by_well_df[["Plate_string"]] <- sub("^ha_", "ha_tf", TF_by_well_df[["Plate_string"]])
  TF_by_well_df[["Plate_string"]] <- sub("^ho_", "ho_tf", TF_by_well_df[["Plate_string"]])


  ## Add data on obsolete wells
  are_obsolete <- TF_by_well_df[["Entrez_ID"]] %in% sg4_by_well_df[["Entrez_ID"]]
  TF_by_well_df[["Is_obsolete"]] <- ifelse(is.na(TF_by_well_df[["Entrez_ID"]]),
                                           NA,
                                           ifelse(are_obsolete, "Yes", "No")
                                           )
  sg4_by_well_df[["Is_obsolete"]] <- NA
  sg4_by_gene_df[["Is_obsolete"]] <- NA


  ## Adjust the names of non-targeting controls

  are_TF_controls <- TF_by_well_df[["Is_control"]] %in% "Yes"
  TF_by_well_df[["Combined_ID"]][are_TF_controls] <- sub("Control_", "Control_TF", TF_by_well_df[["Combined_ID"]][are_TF_controls], fixed = TRUE)

  TF_by_well_df[["Sublibrary_4sg"]] <- ifelse(are_TF_controls, "Controls", "Transcription factors")


  ## Merge the two data frames (ordered by well)

  shared_columns <- intersect(names(sg4_by_well_df), names(TF_by_well_df))

  full_4sg_by_well_df <- rbind.data.frame(
    TF_by_well_df[, shared_columns],
    sg4_by_well_df[, shared_columns],
    stringsAsFactors = FALSE,
    make.row.names = FALSE
  )


  ## Merge the two data frames (ordered by gene)

  full_4sg_by_gene_df <- rbind.data.frame(
    TF_by_well_df[!(are_TF_controls), shared_columns],
    sg4_by_gene_df[sg4_by_gene_df[["Is_control"]] %in% "No", shared_columns],
    TF_by_well_df[are_TF_controls, shared_columns],
    sg4_by_gene_df[sg4_by_gene_df[["Is_control"]] %in% "Yes", shared_columns],
    stringsAsFactors = FALSE,
    make.row.names = FALSE
  )



  are_obsolete <- full_4sg_by_gene_df[["Is_obsolete"]] %in% "Yes"
  assign("delete_full_4sg_by_gene_df", full_4sg_by_gene_df, envir = globalenv())
  if ("TSS_number" %in% names(full_4sg_by_gene_df)) {
    TSS_vec <- full_4sg_by_gene_df[["TSS_number"]]
    are_NA_TSS <- is.na(TSS_vec)
    stopifnot(all(full_4sg_by_gene_df[["Num_TSSs"]][are_NA_TSS] == 1))
    TSS_vec[are_NA_TSS] <- 1L
    TSS_vec <- TSS_vec[!(are_obsolete)]
  } else {
    TSS_vec <- rep(NA, sum(!(are_obsolete)))
  }
  by_gene_order <- order(GetMinEntrez(full_4sg_by_gene_df[["Entrez_ID"]][!(are_obsolete)]),
                         TSS_vec
                         )
  use_indices <- which(!(are_obsolete))[by_gene_order]
  full_4sg_by_gene_df <- full_4sg_by_gene_df[use_indices, ]
  row.names(full_4sg_by_gene_df) <- NULL

  full_4sg_by_gene_df <- full_4sg_by_gene_df[, names(full_4sg_by_gene_df) != "Is_obsolete"]

  full_4sg_by_gene_df[["Sublibrary_4sg"]][full_4sg_by_gene_df[["Sublibrary_4sg"]] == "Misc / controls"]    <- "Controls"
  full_4sg_by_gene_df[["Sublibrary_4sg"]][full_4sg_by_gene_df[["Sublibrary_4sg"]] == "Misc / changed TFs"] <- "Transcription Factors"
  full_4sg_by_gene_df[["Sublibrary_4sg"]][full_4sg_by_gene_df[["Sublibrary_4sg"]] == "Misc / unassigned"]  <- "Unassigned"

  results_list <- list(
    "full_4sg_by_well_df" = full_4sg_by_well_df,
    "full_4sg_by_gene_df" = full_4sg_by_gene_df
  )
  return(results_list)
}





# Functions for examining different versions of the library ---------------

GetGeneOrTSSIDs <- function(CRISPR_df) {
  stopifnot(length(unique(table(CRISPR_df[["Rank"]]))) == 1)
  if ("AltTSS_ID" %in% names(CRISPR_df)) {
    gene_or_TSS_IDs <- CRISPR_df[["AltTSS_ID"]]
  } else {
    gene_or_TSS_IDs <- CRISPR_df[["Combined_ID"]]
  }
  gene_or_TSS_IDs_fac <- factor(gene_or_TSS_IDs,
                                levels = unique(gene_or_TSS_IDs)
                                )
  CheckThatFactorIsInOrder(gene_or_TSS_IDs_fac)
  return(gene_or_TSS_IDs_fac)
}


Get4sgIDs <- function(CRISPR_df) {
  gene_or_TSS_IDs_fac <- GetGeneOrTSSIDs(CRISPR_df)
  sgRNAs_list <- split(toupper(CRISPR_df[["sgRNA_sequence"]]),
                       gene_or_TSS_IDs_fac
                       )
  sgRNAs_vec <- vapply(sgRNAs_list, paste0, collapse = "_", "")
  foursg_IDs <- paste0(as.character(gene_or_TSS_IDs_fac), "__",
                       rep.int(sgRNAs_vec, lengths(sgRNAs_list))
                       )
  return(foursg_IDs)
}









