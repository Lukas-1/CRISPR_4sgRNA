### 10th September 2019 ###




# Define lookup maps ------------------------------------------------------

CRISPR_library_sources <- c(
  "Curated",
  "GPP",
  "Calabrese",
  "hCRISPRa-v2",
  "Dolcetto",
  "hCRISPRi-v2.1",
  "hCRISPRi-v2.0",
  "Brunello",
  "TKOv3",

  "Caprano",
  "mCRISPRa-v2",
  "Dolomiti",
  "mCRISPRi-v2"
)




# Define functions --------------------------------------------------------

AddHorlbeckTSSSource <- function(sgRNA_df, TSS_df) {
  sgRNA_df_transcript_IDs <- paste0(sgRNA_df[["gene"]], "__", sgRNA_df[["transcript"]])
  TSS_df_transcript_IDs <- paste0(TSS_df[["gene"]], "__", TSS_df[["transcript"]])
  TSS_source_matches <- match(sgRNA_df_transcript_IDs, TSS_df_transcript_IDs)
  sgRNA_df[["TSS_source"]] <- TSS_df[["TSS source"]][TSS_source_matches]
  return(sgRNA_df)
}


Exchange5PrimeG <- function(CRISPR_df) {
  are_5prime_G <- !(is.na(CRISPR_df[["Start"]])) &
                  (CRISPR_df[["Num_0MM"]] == 0) & (CRISPR_df[["Num_5G_MM"]] == 1) &
                  (grepl("[hm]CRISPR", CRISPR_df[["Source"]])) &
                  (CRISPR_df[["Is_control"]] != "Yes") &
                  (CRISPR_df[["Exchanged_5pG"]] %in% "No")
  GRanges_object <- RangesDfToGRangesObject(CRISPR_df[are_5prime_G, ])
  nucleotide_5p_vec <- substr(as.character(motifRG::getSequence(GRanges_object, BSgenome.Hsapiens.UCSC.hg38)), 1, 1)
  CRISPR_df[["Exchanged_5pG"]] <- ifelse(are_5prime_G, "Yes", CRISPR_df[["Exchanged_5pG"]])
  CRISPR_df[["sgRNA_sequence"]][are_5prime_G] <- paste0(nucleotide_5p_vec,
                                                        substr(CRISPR_df[["sgRNA_sequence"]][are_5prime_G],
                                                               2,
                                                               nchar(CRISPR_df[["sgRNA_sequence"]][are_5prime_G])
                                                               )
                                                        )
  return(CRISPR_df)
}


FindIdenticalsgRNAs <- function(entrez_ID, symbol, sgRNA, CRISPR_df) {
  are_this_gene <- (CRISPR_df[["Entrez_ID"]] %in% entrez_ID) & (CRISPR_df[["Gene_symbol"]] %in% symbol)
  if (!(any(are_this_gene))) {
    are_this_gene <- is.na(CRISPR_df[["Entrez_ID"]]) & is.na(CRISPR_df[["Gene_symbol"]]) & (CRISPR_df[["Original_symbol"]] %in% symbol)
  }
  if (!(any(are_this_gene))) {
    my_result <- NULL
  } else {
    are_identical_sgRNAs <- toupper(CRISPR_df[["sgRNA_sequence"]][are_this_gene]) == toupper(sgRNA)
    if (any(are_identical_sgRNAs)) {
      my_result <- CRISPR_df[are_this_gene, , drop = FALSE][are_identical_sgRNAs, , drop = FALSE]
    } else {
      my_result <- NULL
    }
  }
  return(my_result)
}


ResolveDf <- function(replicates_df, drop_columns, concatenate_columns) {
  results_df <- replicates_df
  unique_symbols <- replicates_df[["Original_symbol"]]
  unique_symbols <- unique(unique_symbols[unique_symbols != ""])
  if (length(unique_symbols) == 0) {
    results_df[["Original_symbol"]] <- ""
  } else if (length(unique_symbols) == 1) {
    results_df[["Original_symbol"]] <- unique_symbols
  } else {
    unique_symbols <- unique(paste0(replicates_df[["Original_symbol"]], " (", replicates_df[["Source"]], ")")[replicates_df[["Original_symbol"]] %in% unique_symbols])
    results_df[["Original_symbol"]] <- paste0(unique_symbols, collapse = " | ")
  }

  unique_sources <- unique(unlist(strsplit(results_df[["Source"]], ", ", fixed = TRUE)))
  unique_sources <- unique_sources[order(match(unique_sources, CRISPR_library_sources))]
  results_df[["Source"]] <- paste0(unique_sources, collapse = ", ")

  if ("Entrez_ID_assignment" %in% names(replicates_df)) {
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
    unique_assignments <- unique(results_df[["Entrez_ID_assignment"]])
    unique_assignments <- unique_assignments[order(match(unique_assignments, assignments_order))]
    results_df[["Entrez_ID_assignment"]] <- paste0(unique_assignments, collapse = " / ")
  }

  if ("Exchanged_5pG" %in% names(replicates_df)) {
    if (any(replicates_df[["Exchanged_5pG"]] %in% "Yes")) {
      results_df[["Exchanged_5pG"]] <- "Yes"
    } else {
      results_df[["Exchanged_5pG"]] <- "No"
    }
  }

  if ("GuideScan_offtarget_category" %in% names(replicates_df)) {
    offtarget_vec <- unique(replicates_df[["GuideScan_offtarget_category"]])
    if (length(offtarget_vec) > 1) {
      offtarget_vec <- offtarget_vec[!(offtarget_vec %in% "Unknown")]
    }
    if (length(offtarget_vec) == 1) {
      results_df[["GuideScan_offtarget_category"]] <- offtarget_vec
    }
  }

  for (column_name in grep("^Entrez_source_", names(replicates_df), value = TRUE)) {
    entrez_source_vec <- unique(replicates_df[[column_name]])
    entrez_source_vec <- entrez_source_vec[!(is.na(entrez_source_vec))]
    if (length(entrez_source_vec) == 0) {
      results_df[[column_name]] <- NA_integer_
    } else {
      results_df[[column_name]] <- min(entrez_source_vec)
    }
  }
  if ("Exon_number_GPP" %in% concatenate_columns) {
    results_df <- results_df[order(results_df[["Exon_number_GPP"]]), ]
  }
  for (column_name in intersect(concatenate_columns, colnames(results_df))) {
    character_vec <- unique(replicates_df[[column_name]])
    character_vec <- character_vec[!(is.na(character_vec))]
    if (length(character_vec) == 0) {
      results_df[[column_name]] <- NA_character_
    } else {
      results_df[[column_name]] <- paste0(character_vec, collapse = "; ")
    }
  }

  results_df[["sgRNA_sequence"]] <- results_df[["sgRNA_sequence"]][order(results_df[["sgRNA_sequence"]] == toupper(results_df[["sgRNA_sequence"]]))][[1]]
  other_columns <- setdiff(names(replicates_df), c("Original_symbol", "sgRNA_sequence", "Original_index", "Original_source"))

  for (column in other_columns) {
    my_values <- unique(replicates_df[[column]])
    my_values <- my_values[!(is.na(my_values))]
    my_values <- my_values[my_values != ""]
    if (length(my_values) == 1) {
      results_df[[column]] <- my_values
    }
  }

  any_resolved <- FALSE
  are_duplicated <- duplicated(results_df[, !(names(results_df) %in% drop_columns)])

  if ("sgRNA_context_sequence" %in% colnames(results_df)) {
    are_otherwise_duplicated <- !(are_duplicated) & duplicated(results_df[, !(names(results_df) %in% c("sgRNA_context_sequence", drop_columns))])
    if (any(are_otherwise_duplicated)) {
      results_df[["sgRNA_context_sequence"]] <- NA_character_
      are_duplicated <- are_duplicated | are_otherwise_duplicated
    }
  }
  if (any(are_duplicated)) {
    results_df <- results_df[!(are_duplicated), ]
    any_resolved <- TRUE
  }

  ### Deal with triplicates (two transcripts in the hCRISPR/mCRISPR library, and also found in one of the other libraries) ###
  if ((nrow(results_df) == 3) && any(c("hCRISPRa_v2_rank", "hCRISPRi_v2_rank", "mCRISPRa_v2_rank", "mCRISPRi_v2_rank") %in% colnames(results_df))) {
    are_Horlbeck <- results_df[["Original_source"]] %in% c("hCRISPRa-v2", "hCRISPRi-v2.1", "mCRISPRa-v2", "mCRISPRi-v2")
    are_not_Horlbeck <- results_df[["Original_source"]] %in% c("GPP", "Calabrese", "GPP, Calabrese", "Dolcetto", "GPP, Dolcetto",
                                                               "Caprano", "GPP, Caprano", "Dolomiti", "GPP, Dolomiti"
                                                               )
    if ((sum(are_Horlbeck) == 2) && (sum(are_not_Horlbeck) == 1)) {
      non_Horlbeck_columns <- c("GPP_rank", "Original_PAM",
                                "Entrez_source_Calabrese", "Calabrese_rank",
                                "Entrez_source_Dolcetto", "Dolcetto_rank",
                                "Entrez_source_Caprano", "Caprano_rank",
                                "Entrez_souce_Dolomiti", "Dolomiti_rank"
                                )
      for (column_name in non_Horlbeck_columns) {
        results_df[[column_name]] <- results_df[[column_name]][are_not_Horlbeck]
      }
      results_df <- results_df[are_Horlbeck, ]
      any_resolved <- TRUE
    }
  }

  if (any_resolved) {
    return(results_df)
  } else {
    return(replicates_df)
  }
}




ResolveDuplicates <- function(CRISPR_df, concatenate_columns) {

  if ("Rank" %in% names(CRISPR_df)) {
    warning("Warning: The 'Rank' column will cause trouble, i.e., duplicates being ignored! Do not rank sgRNAs until all duplicates have been eliminated.")
  }

  ### Make columns with helper vectors ###

  CRISPR_df[["sgID"]] <- paste0(CRISPR_df[["Combined_ID"]], "__", toupper(CRISPR_df[["sgRNA_sequence"]]))
  CRISPR_df[["Original_index"]] <- seq_len(nrow(CRISPR_df))
  CRISPR_df[["Original_source"]] <- CRISPR_df[["Source"]]
  drop_columns <- c("sgID", "Original_index", "Original_source")


  ### Split the input data frame ###

  split_df_list <- split(CRISPR_df, factor(CRISPR_df[["sgID"]], levels = unique(CRISPR_df[["sgID"]])))


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
    assign("delete_num_replicates",   num_replicates,   envir = globalenv())
    assign("delete_resolved_df_list", resolved_df_list, envir = globalenv())

    ### Create the results data frame ###

    split_df_list[are_multiplicates] <- resolved_df_list
    results_df <- do.call(rbind.data.frame, c(split_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

    ### Report results ###
    message(paste0(nrow(CRISPR_df) - nrow(results_df), " duplicate entries were combined!"))
  }

  ### Re-order the results data frame ###

  results_df <- results_df[order(results_df[["Is_control"]] == "Yes",
                                 GetMinEntrez(results_df[["Entrez_ID"]]),
                                 results_df[["Combined_ID"]],
                                 ifelse(results_df[["Is_control"]] == "Yes", NA, results_df[["Original_symbol"]]),
                                 results_df[["Original_index"]]
                                 ), ]
  row.names(results_df) <- NULL

  ### Remove superfluous columns ###

  results_df <- results_df[, !(names(results_df) %in% drop_columns)]

  return(results_df)
}







