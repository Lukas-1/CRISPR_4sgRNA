### 10th September 2019 ###




# Define lookup maps ------------------------------------------------------

CRISPR_library_sources <- c(
  "Curated", "GPP",
  "Calabrese", "hCRISPRa-v2",
  "Brunello", "TKOv3"
)




# Define functions --------------------------------------------------------

FindIdenticalsgRNAs <- function(entrez_ID, symbol, sgRNA, CRISPR_df) {
  are_this_gene <- (CRISPR_df[, "Entrez_ID"] %in% entrez_ID) & (CRISPR_df[, "Gene_symbol"] %in% symbol)
  if (!(any(are_this_gene))) {
    are_this_gene <- is.na(CRISPR_df[, "Entrez_ID"]) & is.na(CRISPR_df[, "Gene_symbol"]) & (CRISPR_df[, "Original_symbol"] %in% symbol)
  }
  if (!(any(are_this_gene))) {
    my_result <- NULL
  } else {
    are_identical_sgRNAs <- toupper(CRISPR_df[are_this_gene, "sgRNA_sequence"]) == toupper(sgRNA)
    if (any(are_identical_sgRNAs)) {
      my_result <- CRISPR_df[are_this_gene, , drop = FALSE][are_identical_sgRNAs, , drop = FALSE]
    } else {
      my_result <- NULL
    }
  }
  return(my_result)
}



ResolveDf <- function(replicates_df, drop_columns, concatenate_columns) {

  assign("delete_replicates_df", replicates_df, envir = globalenv())

  results_df <- replicates_df

  unique_symbols <- replicates_df[, "Original_symbol"]
  unique_symbols <- unique(unique_symbols[unique_symbols != ""])
  if (length(unique_symbols) == 0) {
    results_df[, "Original_symbol"] <- ""
  } else if (length(unique_symbols) == 1) {
    results_df[, "Original_symbol"] <- unique_symbols
  } else {
    unique_symbols <- unique(paste0(replicates_df[, "Original_symbol"], " (", replicates_df[, "Source"], ")")[replicates_df[, "Original_symbol"] %in% unique_symbols])
    results_df[, "Original_symbol"] <- paste0(unique_symbols, collapse = " | ")
  }

  unique_sources <- unique(unlist(strsplit(results_df[, "Source"], ", ", fixed = TRUE)))
  unique_sources <- unique_sources[order(match(unique_sources, CRISPR_library_sources))]
  results_df[, "Source"] <- paste0(unique_sources, collapse = ", ")

  if ("Entrez_ID_assignment" %in% colnames(replicates_df)) {
    assignments_order <- c("Unambiguous",
                           "The gene symbol was ambiguous; mapped using GuideScan",
                           "The gene symbol was ambiguous; mapped by proximity to TSS",
                           "The gene symbol was ambiguous; could not be mapped unambiguously",
                           "The gene symbol was ambiguous; could not be mapped to any gene",
                           "No Entrez ID was found for the gene symbol; mapped using GuideScan",
                           "No Entrez ID was found for the gene symbol; mapped by proximity to TSS",
                           "No Entrez ID was found for the gene symbol; could not be mapped unambiguously",
                           "No Entrez ID was found for the gene symbol; could not be mapped to any gene"
                           )
    unique_assignments <- unique(results_df[, "Entrez_ID_assignment"])
    unique_assignments <- unique_assignments[order(match(unique_assignments, assignments_order))]
    results_df[, "Entrez_ID_assignment"] <- paste0(unique_assignments, collapse = " / ")
  }

  if ("Exchanged_5pG" %in% colnames(replicates_df)) {
    if (any(replicates_df[, "Exchanged_5pG"] %in% "Yes")) {
      results_df[, "Exchanged_5pG"] <- "Yes"
    } else {
      results_df[, "Exchanged_5pG"] <- "No"
    }
  }

  if ("GuideScan_offtarget_category" %in% colnames(replicates_df)) {
    offtarget_vec <- unique(replicates_df[, "GuideScan_offtarget_category"])
    if (length(offtarget_vec) > 1) {
      offtarget_vec <- offtarget_vec[!(offtarget_vec %in% "Unknown")]
    }
    if (length(offtarget_vec) == 1) {
      results_df[, "GuideScan_offtarget_category"] <- offtarget_vec
    }
  }

  for (column_name in grep("^Entrez_source_", colnames(replicates_df), value = TRUE)) {
    entrez_source_vec <- unique(replicates_df[, column_name])
    entrez_source_vec <- entrez_source_vec[!(is.na(entrez_source_vec))]
    if (length(entrez_source_vec) == 0) {
      results_df[, column_name] <- NA_integer_
    } else {
      results_df[, column_name] <- min(entrez_source_vec)
    }
  }
  for (column_name in concatenate_columns) {
    assign("delete_replicates_df", replicates_df, envir = globalenv())
    assign("delete_column_name", column_name, envir = globalenv())
    character_vec <- unique(replicates_df[, column_name])
    character_vec <- character_vec[!(is.na(character_vec))]
    if (length(character_vec) == 0) {
      results_df[, column_name] <- NA_character_
    } else {
      results_df[, column_name] <- paste0(character_vec, collapse = "; ")
    }
  }

  results_df[, "sgRNA_sequence"] <- results_df[, "sgRNA_sequence"][order(results_df[, "sgRNA_sequence"] == toupper(results_df[, "sgRNA_sequence"]))][[1]]
  other_columns <- colnames(replicates_df)[!(colnames(replicates_df) %in% c("Original_symbol", "sgRNA_sequence", "Original_index", "Original_source"))]

  for (column in other_columns) {
    my_values <- unique(replicates_df[, column])
    my_values <- my_values[!(is.na(my_values))]
    my_values <- my_values[my_values != ""]
    if (length(my_values) == 1) {
      results_df[, column] <- my_values
    }
  }

  any_resolved <- FALSE
  are_duplicated <- duplicated(results_df[, !(colnames(results_df) %in% drop_columns)])
  if (any(are_duplicated)) {
    results_df <- results_df[!(are_duplicated), ]
    any_resolved <- TRUE
  }

  ### Deal with triplicates (two transcripts in hCRISPRa-v2, and also found in Calabrese) ###
  if ((nrow(results_df) == 3) && (all(c("hCRISPRa-v2", "Calabrese") %in% results_df[, "Original_source"])) && (all(results_df[, "Original_source"] %in% c("hCRISPRa-v2", "Calabrese")))) {
    results_df <- results_df[results_df[, "Original_source"] == "hCRISPRa-v2", ]
    any_resolved <- TRUE
  }

  if (any_resolved) {
    return(results_df)
  } else {
    return(replicates_df)
  }
}




ResolveDuplicates <- function(CRISPR_df, concatenate_columns = c("Sublibrary", "hCRISPRa_v2_ID")) {

  ### Make columns with helper vectors ###

  CRISPR_df[, "sgID"] <- paste0(CRISPR_df[, "Combined_ID"], "__", toupper(CRISPR_df[, "sgRNA_sequence"]))
  CRISPR_df[, "Original_index"] <- seq_len(nrow(CRISPR_df))
  CRISPR_df[, "Original_source"] <- CRISPR_df[, "Source"]
  drop_columns <- c("sgID", "Original_index", "Original_source")


  ### Split the input data frame ###

  split_df_list <- split(CRISPR_df, factor(CRISPR_df[, "sgID"], levels = unique(CRISPR_df[, "sgID"])))


  ### Find multiplicates ###

  are_multiplicates <- vapply(split_df_list, nrow, integer(1)) > 1
  num_multiples <- sum(are_multiplicates)

  if (num_multiples == 0) {
    message("No duplicates were found!")
    results_df <- CRISPR_df
  } else {

    ### Resolve multiplicates ###

    resolved_df_list <- lapply(seq_len(num_multiples), function(x) {
      if ((((x <= 100) || (num_multiples < 500)) && ((x %% 50) == 0)) || ((x <= 1000) && ((x %% 500) == 0)) || ((x %% 1000) == 0)) {
        message(paste0("#", x, " out of ", num_multiples, " duplicate sgRNAs is being checked..."))
      }
      ResolveDf(split_df_list[[which(are_multiplicates)[[x]]]], drop_columns = drop_columns, concatenate_columns = concatenate_columns)
    })
    message(paste0("All ", num_multiples, " duplicate sgRNAs were checked."))

    num_replicates <- vapply(resolved_df_list, nrow, integer(1))
    if (any(num_replicates >= 3)) {
      message(paste0("Out of ", sum(num_replicates >= 2), " replicates that remained, ",
                     sum(num_replicates >= 3), " had 3 or more identical sgRNAs! Please check these so they can be resolved."
                     )
              )
    }

    assign("delete_num_replicates", num_replicates, envir = globalenv())
    assign("delete_resolved_df_list", resolved_df_list, envir = globalenv())
    assign("delete_split_df_list", split_df_list, envir = globalenv())


    ### Create the results data frame ###

    split_df_list[are_multiplicates] <- resolved_df_list
    results_df <- do.call(rbind.data.frame, c(split_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))


    ### Report results ###
    message(paste0(nrow(CRISPR_df) - nrow(results_df), " duplicate entries were combined!"))
  }

  ### Re-order the results data frame ###

  results_df <- results_df[order(results_df[, "Is_control"] == "Yes",
                                 GetMinEntrez(results_df[, "Entrez_ID"]),
                                 results_df[, "Combined_ID"],
                                 ifelse(results_df[, "Is_control"] == "Yes", NA, results_df[, "Original_symbol"]),
                                 results_df[, "Original_index"]
                                 ), ]
  rownames(results_df) <- NULL


  ### Remove superfluous columns ###

  results_df <- results_df[, !(colnames(results_df) %in% drop_columns)]

  return(results_df)
}







TidyGPPCRISPRaDf <- function(my_GPP_CRISPRa_df) {

  keep_columns <- c("Input", "Target Gene ID", "Target Gene Symbol",
                    "Reference Sequence", "Strand of Target", "TSS Position",
                    "Strand of sgRNA", "sgRNA Cut Position",
                    "sgRNA Sequence", "sgRNA Context Sequence", "PAM Sequence",
                    "sgRNA Cut Site TSS Offset",
                    "# Off-Target Tier I Match Bin I Matches",   "# Off-Target Tier II Match Bin I Matches",   "# Off-Target Tier III Match Bin I Matches",
                    "# Off-Target Tier I Match Bin II Matches",  "# Off-Target Tier II Match Bin II Matches",  "# Off-Target Tier III Match Bin II Matches",
                    "# Off-Target Tier I Match Bin III Matches", "# Off-Target Tier II Match Bin III Matches", "# Off-Target Tier III Match Bin III Matches",
                    "# Off-Target Tier I Match Bin IV Matches",  "# Off-Target Tier II Match Bin IV Matches",  "# Off-Target Tier III Match Bin IV Matches",
                    "On-Target Efficacy Score", "DHS Score",
                    "On-Target Rank", "Off-Target Rank", "Combined Rank",
                    "Pick Order", "Picking Round", "Picking Notes"
                    )

  were_not_found <- grepl("^ERROR: Gene .+ not found", my_GPP_CRISPRa_df[, "Picking Notes"]) & is.na(my_GPP_CRISPRa_df[, "Quota"])

  message(paste0(sum(were_not_found), " genes were not found by the Broad Institute's Genetic Perturbation Platform (GPP) CRISPRa sgRNA picker tool and were omitted from the data frame!"))

  results_df <- my_GPP_CRISPRa_df[!(were_not_found), keep_columns]
  rownames(results_df) <- NULL

  return(results_df)
}





TidyGPPCRISPRkoDf <- function(my_GPP_CRISPRko_df) {

  keep_columns <- c("Input", "Target Gene ID", "Target Gene Symbol",
                    "Target Transcript",
                    "Reference Sequence", "Strand of Target",
                    "Strand of sgRNA", "Orientation", "sgRNA Cut Position (1-based)",
                    "sgRNA Sequence", "sgRNA Context Sequence", "PAM Sequence",
                    "Exon Number",

                    "Target Cut Length", "Target Total Length", "Target Cut %",


                    "# Off-Target Tier I Match Bin I Matches",   "# Off-Target Tier II Match Bin I Matches",   "# Off-Target Tier III Match Bin I Matches",   "# Off-Target Tier IV Match Bin I Matches",
                    "# Off-Target Tier I Match Bin II Matches",  "# Off-Target Tier II Match Bin II Matches",  "# Off-Target Tier III Match Bin II Matches",  "# Off-Target Tier IV Match Bin II Matches",
                    "# Off-Target Tier I Match Bin III Matches", "# Off-Target Tier II Match Bin III Matches", "# Off-Target Tier III Match Bin III Matches", "# Off-Target Tier IV Match Bin III Matches",
                    "# Off-Target Tier I Match Bin IV Matches",  "# Off-Target Tier II Match Bin IV Matches",  "# Off-Target Tier III Match Bin IV Matches",  "# Off-Target Tier IV Match Bin IV Matches",
                    "On-Target Efficacy Score",
                    "On-Target Rank", "Off-Target Rank", "Combined Rank",
                    "Pick Order", "Picking Round", "Picking Notes"
                    )

  were_not_found <- grepl("^ERROR: Gene .+ not found", my_GPP_CRISPRko_df[, "Picking Notes"]) & is.na(my_GPP_CRISPRko_df[, "Quota"])

  message(paste0(sum(were_not_found), " genes were not found by the Broad Institute's Genetic Perturbation Platform (GPP) CRISPRa sgRNA picker tool and were omitted from the data frame!"))

  results_df <- my_GPP_CRISPRko_df[!(were_not_found), keep_columns]
  rownames(results_df) <- NULL

  return(results_df)
}













