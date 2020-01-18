


# Import packages and source code -----------------------------------------

library("data.table")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R")) # For GetMinEntrez
source(file.path(general_functions_directory, "09) Constants and settings.R")) # For GetMinEntrez






# Define functions --------------------------------------------------------

OrderControlsVec <- function(CRISPR_df) {
  is_CRISPRko <- "Exon_number_GPP" %in% names(CRISPR_df)
  NA_vec <- rep(NA, nrow(CRISPR_df))
  if (is_CRISPRko) {
    controls_vec <- NA_vec
    are_Brunello_controls <- (CRISPR_df[, "Is_control"] == "Yes") & (CRISPR_df[, "Source"] == "Brunello")
    controls_vec[are_Brunello_controls] <- as.integer(sub("^Control_", "", CRISPR_df[are_Brunello_controls, "Combined_ID"]))
  } else {
    controls_vec <- NA_vec
  }
  return(controls_vec)
}




RankCRISPRDf <- function(CRISPR_df, reorder_by_rank = TRUE, allow_5pG_MM = FALSE, ID_column = "Combined_ID") {
  # Requires the constants 'preferred_AF_max_column' and 'SNP_frequency_cutoff' in the global environment

  MakeFactor <- function(x) factor(x, levels = unique(x))

  pre_order <- order(CRISPR_df[, "Is_control"] == "Yes",
                     GetMinEntrez(CRISPR_df[, "Entrez_ID"]),
                     OrderControlsVec(CRISPR_df),
                     CRISPR_df[, "Combined_ID"],
                     MakeFactor(CRISPR_df[, ID_column])
                     )
  CRISPR_df <- CRISPR_df[pre_order, ]

  CRISPR_df[, "Rank"] <- unlist(tapply(seq_len(nrow(CRISPR_df)),
                                       MakeFactor(CRISPR_df[, ID_column]),
                                       function(x) RankDf(CRISPR_df[x, , drop = FALSE], allow_5pG_MM = allow_5pG_MM),
                                       simplify = FALSE
                                       )
                                )

  if (reorder_by_rank) {
    new_order <- order(CRISPR_df[, "Is_control"] == "Yes",
                       GetMinEntrez(CRISPR_df[, "Entrez_ID"]),
                       OrderControlsVec(CRISPR_df),
                       CRISPR_df[, "Combined_ID"],
                       MakeFactor(CRISPR_df[, ID_column]),
                       CRISPR_df[, "Rank"]
                       )
    CRISPR_df <- CRISPR_df[new_order, ]
    row.names(CRISPR_df) <- NULL
  }
  return(CRISPR_df)
}



RankDf <- function(CRISPR_sub_df, allow_5pG_MM = FALSE) {

  exactly_one_vec <- (CRISPR_sub_df[, "Num_0MM"] == 1)
  if (allow_5pG_MM) {
    are_5pG_MM <- (CRISPR_sub_df[, "Num_0MM"] == 0) & (CRISPR_sub_df[, "Num_5G_MM"] == 1)
    exactly_one_vec <- exactly_one_vec | are_5pG_MM
  }

  SNP_data_present <- (preferred_AF_max_column %in% names(CRISPR_sub_df))

  NA_vec <- rep.int(NA, nrow(CRISPR_sub_df))

  is_CRISPRko <- "Exon_number_GPP" %in% names(CRISPR_sub_df)

  list_for_ranking <- list(CRISPR_sub_df[, "Is_control"] == "No",
                           !(grepl("TTTT", CRISPR_sub_df[, "sgRNA_sequence"], ignore.case = TRUE)),
                           !(is.na(CRISPR_sub_df[, "Start"])),
                           substr(CRISPR_sub_df[, "PAM"], 2, 3) == "GG",
                           if (SNP_data_present) !((CRISPR_sub_df[, preferred_AF_max_column] > SNP_frequency_cutoff) %in% TRUE)                   else NA_vec,
                           !(is.na(CRISPR_sub_df[, "GuideScan_specificity"])),
                           if ("CRISPOR_4MM_specificity" %in% names(CRISPR_sub_df)) !(is.na(CRISPR_sub_df[, "CRISPOR_4MM_specificity"]))          else NA_vec,

                           exactly_one_vec, # prefers exactly one match
                           ifelse(exactly_one_vec, TRUE, CRISPR_sub_df[, "Num_0MM"] != 0), # ensures that Num_0MM >= 2 is preferred over Num_0MM == 0
                           ifelse(CRISPR_sub_df[, "Num_0MM"] >= 2, -(CRISPR_sub_df[, "Num_0MM"]), NA_integer_),
                           if (allow_5pG_MM) ifelse(exactly_one_vec, 1, -(CRISPR_sub_df[, "Num_5G_MM"]))                                          else NA_vec, # penalize multiple 5' G-mismatched locations
                           if (allow_5pG_MM || !("Num_5G_MM" %in% names(CRISPR_sub_df))) -(CRISPR_sub_df[, "Num_1MM"]) else -(rowSums(CRISPR_sub_df[, c("Num_5G_MM", "Num_1MM")])),

                           if (is_CRISPRko && ("CRISPOR_Graf_status" %in% names(CRISPR_sub_df))) (CRISPR_sub_df[, "CRISPOR_Graf_status"] %in% "GrafOK") else NA_vec,

                           CRISPR_sub_df[, "GuideScan_specificity"],
                           if ("CRISPOR_4MM_specificity" %in% names(CRISPR_sub_df)) CRISPR_sub_df[, "CRISPOR_4MM_specificity"]                    else NA_vec,
                           CRISPR_sub_df[, "GuideScan_efficiency"],
                           if ("CRISPOR_Doench_efficacy" %in% names(CRISPR_sub_df)) CRISPR_sub_df[, "CRISPOR_Doench_efficacy"]                    else NA_vec,
                           if ("hCRISPRa_v2_rank" %in% names(CRISPR_sub_df)) -(suppressWarnings(as.integer(CRISPR_sub_df[, "hCRISPRa_v2_rank"]))) else NA_vec,
                           if ("Calabrese_rank" %in% names(CRISPR_sub_df)) -(match(CRISPR_sub_df[, "Calabrese_rank"], c("1/2/3", "4/5/6")))       else NA_vec,
                           -(CRISPR_sub_df[, "GPP_rank"])
                           )
  assign("delete_list_for_ranking", list_for_ranking, envir = globalenv())
  my_rank <- data.table::frank(x = lapply(list_for_ranking, "*", -1L), ties.method = "average")
  return(my_rank)
}




