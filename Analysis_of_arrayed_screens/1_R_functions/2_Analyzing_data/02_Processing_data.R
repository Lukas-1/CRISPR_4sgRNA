# 2022-01-17



# Define functions --------------------------------------------------------

NormalizeWithNTControls <- function(use_df, norm_method = "all NT") {

  use_df[, "CellTiterGlo_foldNT"] <- NormPlates(use_df, "CellTiterGlo_raw", foldNT = TRUE, norm_method = norm_method)

  are_after <- seq_len(ncol(use_df)) > which(names(use_df) == "CellTiterGlo_raw")
  use_columns <- unique(c(names(use_df)[!(are_after)], "CellTiterGlo_foldNT",
                          names(use_df)[are_after]
                          ))
  use_df <- use_df[, use_columns]

  for (ri in paste0("_rep", 1:2)) {
    use_df[, paste0("DeltaNT", ri)]             <- NormPlates(use_df, paste0("Raw", ri), norm_method = norm_method)
    use_df[, paste0("FoldNT", ri)]              <- NormPlates(use_df, paste0("Raw", ri), foldNT = TRUE, norm_method = norm_method)
    use_df[, paste0("PercActivation", ri)]      <- NormPlates(use_df, paste0("Raw", ri), percent_activation = TRUE, norm_method = norm_method)
    use_df[, paste0("Raw_log2", ri)]            <- log2(use_df[, paste0("Raw", ri)])
    use_df[, paste0("Log2FC", ri)]              <- NormPlates(use_df, paste0("Raw", ri), take_log2 = TRUE, norm_method = norm_method)
    use_df[, paste0("PercActivation_log2", ri)] <- NormPlates(use_df, paste0("Raw", ri), percent_activation = TRUE, take_log2 = TRUE, norm_method = norm_method)

    use_df[, paste0("Raw_Glo", ri)]                 <- use_df[, paste0("Raw", ri)] / use_df[, "CellTiterGlo_foldNT"]
    use_df[, paste0("DeltaNT_Glo", ri)]             <- NormPlates(use_df, paste0("Raw_Glo", ri), norm_method = norm_method)
    use_df[, paste0("FoldNT_Glo", ri)]              <- NormPlates(use_df, paste0("Raw_Glo", ri), foldNT = TRUE, norm_method = norm_method)
    use_df[, paste0("PercActivation_Glo", ri)]      <- NormPlates(use_df, paste0("Raw_Glo", ri), percent_activation = TRUE, norm_method = norm_method)
    use_df[, paste0("Raw_log2_Glo", ri)]            <- log2(use_df[, paste0("Raw_Glo", ri)])
    use_df[, paste0("Log2FC_Glo", ri)]              <- NormPlates(use_df, paste0("Raw_Glo", ri), take_log2 = TRUE, norm_method = norm_method)
    use_df[, paste0("PercActivation_log2_Glo", ri)] <- NormPlates(use_df, paste0("Raw_Glo", ri), percent_activation = TRUE, take_log2 = TRUE, norm_method = norm_method)
  }

  stripped_columns <- sub("_rep[12]", "", names(use_df))
  use_df <- use_df[, order(match(stripped_columns, stripped_columns))]

  return(use_df)
}



RunSSMDStats <- function(use_df, norm_method = "all NT") {

  ## Calculate SSMD
  use_df[, "SSMD_deltaNT"]      <- Calculate_SSMD(use_df, "Raw_rep1", norm_method = norm_method)
  use_df[, "SSMD_act"]          <- Calculate_SSMD(use_df, "Raw_rep1", percent_activation = TRUE, norm_method = norm_method)
  use_df[, "SSMD_log2"]         <- Calculate_SSMD(use_df, "Raw_rep1", take_log2 = TRUE, norm_method = norm_method)
  use_df[, "SSMD_act_log2"]     <- Calculate_SSMD(use_df, "Raw_rep1", percent_activation = TRUE, take_log2 = TRUE, norm_method = norm_method)

  use_df[, "SSMD_deltaNT_Glo"]  <- Calculate_SSMD(use_df, "Raw_Glo_rep1", norm_method = norm_method)
  use_df[, "SSMD_act_Glo"]      <- Calculate_SSMD(use_df, "Raw_Glo_rep1", percent_activation = TRUE, norm_method = norm_method)
  use_df[, "SSMD_log2_Glo"]     <- Calculate_SSMD(use_df, "Raw_Glo_rep1", take_log2 = TRUE, norm_method = norm_method)
  use_df[, "SSMD_act_log2_Glo"] <- Calculate_SSMD(use_df, "Raw_Glo_rep1", percent_activation = TRUE, take_log2 = TRUE, norm_method = norm_method)


  ## Calculate p value
  use_df[, "p_value_deltaNT"]      <- Calculate_P(use_df, "Raw_rep1", norm_method = norm_method)
  use_df[, "p_value_act"]          <- Calculate_P(use_df, "Raw_rep1", percent_activation = TRUE, norm_method = norm_method)
  use_df[, "p_value_log2"]         <- Calculate_P(use_df, "Raw_rep1", take_log2 = TRUE, norm_method = norm_method)
  use_df[, "p_value_act_log2"]     <- Calculate_P(use_df, "Raw_rep1", percent_activation = TRUE, take_log2 = TRUE, norm_method = norm_method)

  use_df[, "p_value_deltaNT_Glo"]  <- Calculate_P(use_df, "Raw_Glo_rep1", norm_method = norm_method)
  use_df[, "p_value_act_Glo"]      <- Calculate_P(use_df, "Raw_Glo_rep1", percent_activation = TRUE, norm_method = norm_method)
  use_df[, "p_value_log2_Glo"]     <- Calculate_P(use_df, "Raw_Glo_rep1", take_log2 = TRUE, norm_method = norm_method)
  use_df[, "p_value_act_log2_Glo"] <- Calculate_P(use_df, "Raw_Glo_rep1", percent_activation = TRUE, take_log2 = TRUE, norm_method = norm_method)


  ## Calculate hit strength
  meanFC <- rowMeans(use_df[, c("Log2FC_rep1", "Log2FC_rep2")])
  use_df[, "Hit_strength_deltaNT"]      <- meanFC * -log10(use_df[, "p_value_deltaNT"])
  use_df[, "Hit_strength_act"]          <- meanFC * -log10(use_df[, "p_value_act"])
  use_df[, "Hit_strength_log2"]         <- meanFC * -log10(use_df[, "p_value_log2"])
  use_df[, "Hit_strength_act_log2"]     <- meanFC * -log10(use_df[, "p_value_act_log2"])

  meanFC_Glo <- rowMeans(use_df[, c("Log2FC_Glo_rep1", "Log2FC_Glo_rep2")])
  use_df[, "Hit_strength_deltaNT_Glo"]  <- meanFC_Glo * -log10(use_df[, "p_value_deltaNT_Glo"])
  use_df[, "Hit_strength_act_Glo"]      <- meanFC_Glo * -log10(use_df[, "p_value_act_Glo"])
  use_df[, "Hit_strength_log2_Glo"]     <- meanFC_Glo * -log10(use_df[, "p_value_log2_Glo"])
  use_df[, "Hit_strength_act_log2_Glo"] <- meanFC_Glo * -log10(use_df[, "p_value_act_log2_Glo"])

  return(use_df)
}



CreateHitLists <- function(input_df,
                           log2fc_column,
                           p_value_column,
                           hit_strength_column,
                           p_value_cutoff = 0.05,
                           log2fc_cutoff = log2(2)
                           ) {

  ## Establish criteria for defining hits
  are_gene <- !(is.na(input_df[, "Entrez_ID"]))
  if (grepl("_rep[0-9]+$", log2fc_column)) {
    if (!(grepl("_rep1", log2fc_column, fixed = TRUE))) {
      stop("Unexpected value for 'log2fc_column!")
    }
    rep2_column <- sub("_rep1", "_rep2", log2fc_column, fixed = TRUE)
    log2fc_vec <- rowMeans(input_df[, c(log2fc_column, rep2_column)])
  } else {
    log2fc_vec <- input_df[, log2fc_column]
  }

  meet_p_val_cutoff  <- (input_df[, p_value_column] < p_value_cutoff)
  meet_log2fc_cutoff <- abs(log2fc_vec) > log2fc_cutoff

  meet_criteria <- meet_p_val_cutoff & meet_log2fc_cutoff

  ## Define additional columns that are useful for exporting the hit list
  input_df[, "Passes_cutoffs"]    <- meet_criteria
  input_df[, "Mean_log2FC"]       <- log2fc_vec
  input_df[, "p_value_used"]      <- input_df[, p_value_column]
  input_df[, "Hit_strength_used"] <- input_df[, hit_strength_column]


  ## Re-order columns to emphasize the data that was used for choosing hits,
  ## and re-order genes by their rank.
  are_after <- seq_len(ncol(input_df)) > which(names(input_df) == "Is_pos_ctrl")
  use_columns <- unique(c(names(input_df)[!(are_after)],
                          "Passes_cutoffs", "Mean_log2FC", "p_value_used",
                          "Hit_strength_used",
                          names(input_df)[are_after]
                          ))
  input_df <- input_df[, use_columns]

  new_order <- order(meet_criteria,
                     abs(input_df[, hit_strength_column]),
                     decreasing = TRUE
                     )
  reordered_df <- input_df[new_order, ]

  are_gene <- !(is.na(reordered_df[, "Entrez_ID"]))
  are_selected <- are_gene | reordered_df[, "Is_NT_ctrl"]
  reordered_df <- reordered_df[are_selected, ]
  row.names(reordered_df) <- NULL


  ## Create a data frame containing only hit genes
  are_hit <- (!(is.na(reordered_df[, "Entrez_ID"]))) & reordered_df[, "Passes_cutoffs"]
  hits_df <- reordered_df[are_hit, ]
  row.names(hits_df) <- NULL

  results_list <- list(
    "original_df"  = input_df,
    "reordered_df" = reordered_df,
    "hits_df"      = hits_df
  )
  return(results_list)
}


