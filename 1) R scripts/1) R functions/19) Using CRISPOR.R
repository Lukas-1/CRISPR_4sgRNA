### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("data.table") # For data.table::fread (optional)




# Define maps -------------------------------------------------------------

rename_CRISPOR_columns_vec <- c(
  "mitSpecScore"        = "CRISPOR_MIT_specificity",
  "cfdSpecScore"        = "CRISPOR_CFD_specificity",
  "offtargetCount"      = "CRISPOR_off_target_count",
  "Doench '16-Score"    = "CRISPOR_Doench_efficacy",
  "Moreno-Mateos-Score" = "CRISPOR_Moreno_Mateos",
  "Out-of-Frame-Score"  = "CRISPOR_out_of_frame",
  "Lindel-Score"        = "CRISPOR_lindel_score",
  "GrafEtAlStatus"      = "CRISPOR_Graf_status"
)


specificity_demo_columns <- c(
  "Entrez_ID", "Gene_symbol", "Location_ID", "sgRNA_sequence", "PAM", "Original_PAM",
  "Num_0MM", "Num_1MM", "GuideScan_Num_2MM", "GuideScan_Num_3MM",
  "GuideScan_specificity",
  "CRISPOR_CFD_specificity", "CRISPOR_MIT_specificity", "CRISPOR_off_target_count",
  "CRISPOR_Doench_efficacy", "CRISPOR_Graf_status",
  "CRISPOR_3MM_specificity" # This will be NA
)







# Functions for exporting CRISPR databases for input to CRISPOR -----------

MakeBedDf <- function(CRISPR_df, combined_IDs) {
  are_selected <- (CRISPR_df[["Combined_ID"]] %in% combined_IDs) &
                  !(is.na(CRISPR_df[["Start"]]))
  CRISPR_bed_df <- CRISPR_df[are_selected, ]
  export_bed_df <- data.frame(
    CRISPR_bed_df[c("Combined_ID", "Chromosome")],
    "Start"  = ifelse(CRISPR_bed_df[["Strand"]] == "+", CRISPR_bed_df[["Start"]],    CRISPR_bed_df[["Start"]] - 3L) - 1L,
    "End"    = ifelse(CRISPR_bed_df[["Strand"]] == "+", CRISPR_bed_df[["End"]] + 3L, CRISPR_bed_df[["End"]]),
    "Names"  = ".",
    "Scores" = ".",
    CRISPR_bed_df["Strand"],
    stringsAsFactors = FALSE
  )
  are_duplicated <- duplicated(export_bed_df[, !(colnames(export_bed_df) %in% c("sgRNA_sequence", "Combined_ID"))])
  export_bed_df <- export_bed_df[!(are_duplicated), ]
  row.names(export_bed_df) <- NULL
  return(export_bed_df)
}



BreakIntoChunks <- function(UseFunction, CRISPR_df, combined_IDs_list) {
  export_df <- UseFunction(CRISPR_df, combined_IDs = unlist(combined_IDs_list, use.names = FALSE))
  export_df_list <- lapply(combined_IDs_list, function(x) {
    are_this_chunk <- export_df[["Combined_ID"]] %in% x
    return(export_df[are_this_chunk, colnames(export_df) != "Combined_ID"])
  })
  return(export_df_list)
}



PAMorOriginalPAM <- function(CRISPR_df) {
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  are_validated_PAMs <- !(is.na(CRISPR_df[["Original_PAM"]])) &
                        mapply(function(x, y) x %in% strsplit(y, "; ", fixed = TRUE)[[1]],
                               CRISPR_df[["Original_PAM"]],
                               CRISPR_df[["PAM_0MM"]]
                               )
  PAM_vec <- ifelse(is.na(CRISPR_df[["PAM"]]),
                    ifelse(are_validated_PAMs, CRISPR_df[["Original_PAM"]], NA_character_),
                    CRISPR_df[["PAM"]]
                    )
  PAM_vec <- ifelse(substr(PAM_vec, 2, 3) == "GG", PAM_vec, NA_character_)
  return(toupper(PAM_vec))
}



MakeFASTADf <- function(CRISPR_df, combined_IDs) {
  are_selected <- (CRISPR_df[["Combined_ID"]] %in% combined_IDs) &
                  is.na(CRISPR_df[["Start"]]) &
                  (CRISPR_df[["Num_0MM"]] > 0)
  if (!(any(are_selected))) {
    message("No sgRNAs were found that required submission as a FASTA sequence!")
    return(NULL)
  } else {
    results_df <- CRISPR_df[are_selected, ]
    PAM_vec <- PAMorOriginalPAM(results_df)
    results_df <- data.frame(
      results_df[, c("Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "sgRNA_sequence")],
      "PAM" = PAM_vec,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    results_df <- results_df[!(is.na(PAM_vec)), ]
    are_duplicated <- duplicated(results_df[, c("sgRNA_sequence", "PAM")])
    results_df <- results_df[!(are_duplicated), ]
    row.names(results_df) <- NULL
    return(results_df)
  }
}


MakeFASTAvec <- function(use_FASTA_df) {
  full_sequences_vec <- unique(paste0(use_FASTA_df[["sgRNA_sequence"]], use_FASTA_df[["PAM"]]))
  FASTA_list <- lapply(seq_along(full_sequences_vec), function(x) c(paste0(">seq", x), full_sequences_vec[[x]], ""))
  FASTA_vec <- unlist(FASTA_list)
  return(FASTA_vec)
}


WriteCRISPORInputFiles <- function(CRISPOR_chunk_list, file_ending, CRISPOR_input_directory, file_prefix = "Input_for_CRISPOR__") {
  for (i in seq_along(CRISPOR_chunk_list)) {
    file_name <- paste0(file_prefix, names(CRISPOR_chunk_list)[[i]], file_ending)
    write.table(CRISPOR_chunk_list[[i]],
                file = file.path(CRISPOR_input_directory, file_name),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
                )
  }
  return(invisible(NULL))
}






# Functions for processing output from CRISPOR ----------------------------

IdentifyCRISPOROutputFiles <- function() {
  # Requires 'CRISPOR_files_directory' in the global environment

  output_folder_bed   <- "Output_bed"
  output_folder_FASTA <- "Output_FASTA"

  output_files_bed   <- grep("^CRISPOR_output_", list.files(file.path(CRISPOR_files_directory, output_folder_bed)),   value = TRUE)
  output_files_FASTA <- grep("^CRISPOR_output_", list.files(file.path(CRISPOR_files_directory, output_folder_FASTA)), value = TRUE)

  are_offtargets_bed   <- grepl("offs\\.tsv", output_files_bed)
  are_offtargets_FASTA <- grepl("offs\\.tsv", output_files_FASTA)

  results_list <- list(
    "bed"              = file.path(output_folder_bed, output_files_bed[!(are_offtargets_bed)]),
    "FASTA"            = file.path(output_folder_FASTA, output_files_FASTA[!(are_offtargets_FASTA)]),
    "bed_offtargets"   = file.path(output_folder_bed, output_files_bed[are_offtargets_bed]),
    "FASTA_offtargets" = file.path(output_folder_FASTA, output_files_FASTA[are_offtargets_FASTA])
  )
  return(results_list)
}



ReadCRISPOROutputFiles <- function(file_names, is_FASTA, show_messages = FALSE) {
  df_list <- lapply(file_names, function(x) ReadCRISPOROutput(x, use_fread = TRUE, show_messages = show_messages))
  results_df <- rbindlist(df_list)
  rm(df_list) # Save RAM
  is_offtarget <- !("Doench '16-Score" %in% colnames(results_df))
  if (!(is_offtarget)) { # i.e. NOT off-targets output
    for (column_name in c("Doench '16-Score", "Moreno-Mateos-Score", "Out-of-Frame-Score", "Lindel-Score")) {
      if (!(is.integer(results_df[[column_name]]))) {
        are_invalid <- results_df[[column_name]] %in% "NotEnoughFlankSeq"
        message(paste0(sum(are_invalid), " entries in the '", column_name,
                       "' column were invalid (i.e. 'NotEnoughFlankSeq')!"
                       )
                )
        results_df[[column_name]][are_invalid] <- NA_character_
        results_df[[column_name]] <- as.integer(results_df[[column_name]])
      }
    }
    have_no_specificity <- results_df[["mitSpecScore"]] %in% "None"
    stopifnot(identical(have_no_specificity, results_df[["cfdSpecScore"]] %in% "None"))
    if (any(have_no_specificity)) {
      message(paste0(sum(have_no_specificity), " entries had invalid specificity scores (i.e. 'None')!"))
      for (column_name in c("mitSpecScore", "cfdSpecScore")) {
        results_df[[column_name]][have_no_specificity] <- NA_character_
        results_df[[column_name]] <- as.integer(results_df[[column_name]])
      }
      results_df[["offtargetCount"]] <- NA_integer_
    }
  }
  use_columns <- setdiff(colnames(results_df), c("mitOfftargetScore", "cfdOfftargetScore")) # This is to avoid floating-point issues that leads to some duplicates being missed
  if (is_FASTA) {
    use_columns <- setdiff(use_columns, c("seqId", "#seqId"))
  }
  are_duplicated <- duplicated(results_df[, ..use_columns])
  are_forward <- results_df[["guideId"]] %in% "21forw"
  are_to_keep <- are_forward & (!(are_duplicated))
  if (any(are_duplicated)) {
    message(paste0(sum(are_duplicated), " duplicated entries were excluded!"))
    rownames(results_df) <- NULL
  }
  results_df <- results_df[are_to_keep, ]
  keep_columns <- setdiff(colnames(results_df), "guideId")
  results_df <- results_df[, ..keep_columns]
  data.table::setDF(results_df)
  return(results_df)
}


ReadCRISPOROutput <- function(file_name, use_fread = TRUE, show_messages = FALSE) {
  # Requires 'CRISPOR_files_directory' in the global workspace
  if (show_messages) {
    message(paste0("Reading the file: '", file_name, "'..."))
  }
  file_path <- file.path(CRISPOR_files_directory, file_name)
  if (use_fread) {
    results_df <- data.table::fread(file = file_path, sep = "\t",
                                    header = TRUE, quote = "",
                                    data.table = TRUE
                                    )
  } else {
    results_df <- read.table(file_path, sep = "\t", header = TRUE,
                             row.names = NULL, quote = "", comment.char = "",
                             stringsAsFactors = FALSE, check.names = FALSE
                             )
  }
  if (show_messages) {
    message("")
  }
  return(results_df)
}



SummarizeOfftargets <- function(offtargets_df, is_FASTA) {
  assign("delete_offtargets_df", offtargets_df, envir = globalenv())
  offtargets_df[["mismatchCount"]] <- factor(offtargets_df[["mismatchCount"]], levels = 0:4, ordered = TRUE)
  offtargets_df[["cfdOfftargetScore"]] <- as.numeric(ifelse(offtargets_df[["cfdOfftargetScore"]] == "None", NA, offtargets_df[["cfdOfftargetScore"]]))

  if (is_FASTA) {
    ID_column <- "guideSeq"
  } else {
    ID_column <- "seqId"
  }
  assign("delete_ID_column", ID_column, envir = globalenv())
  assign("delete_offtargets_df", offtargets_df, envir = globalenv())

  results_list <- tapply(seq_len(nrow(offtargets_df)),
                         factor(offtargets_df[[ID_column]], levels = unique(offtargets_df[[ID_column]])),
                         function(x) {
                           assign("delete_x", x, envir = globalenv())
                           mismatch_table <- as.integer(table(offtargets_df[["mismatchCount"]][x]))
                           names(mismatch_table) <- paste0("CRISPOR_Num_", 0:4, "MM")
                           specificity_unrounded <- 1 / (1 + sum(offtargets_df[["cfdOfftargetScore"]][x], na.rm = TRUE))
                           specificity_upto3MM <- 1 / (1 + sum(offtargets_df[["cfdOfftargetScore"]][x][offtargets_df[["mismatchCount"]][x] %in% 0:3], na.rm = TRUE))
                           return(c(
                             list("Target_sequence" = offtargets_df[["guideSeq"]][[x[[1]]]]),
                             if (!(is_FASTA)) list("Location_ID" = offtargets_df[["seqId"]][[x[[1]]]]) else NULL,
                             as.list(mismatch_table),
                             list(
                               "CRISPOR_3MM_specificity" = specificity_upto3MM,
                               "CRISPOR_4MM_specificity" = specificity_unrounded
                             )
                           ))
                         })

  results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  results_df[["CRISPOR_Num_2or3MM"]] <- as.integer(rowSums(as.matrix(results_df[, c("CRISPOR_Num_2MM", "CRISPOR_Num_3MM")])))
  return(results_df)
}



ResolveSpecificityOfNone <- function(CRISPOR_df, CRISPR_df = NULL, show_columns = NULL) {
  # Relies on the variable 'specificity_demo_columns' in the global environment

  are_invalid <- (CRISPOR_df[["CRISPOR_CFD_specificity"]] %in% "None") | (CRISPOR_df[["CRISPOR_MIT_specificity"]] %in% "None")
  if (any(are_invalid)) {
    message("\n")
    message("The following sgRNAs had specificity scores of 'None', which were replaced by NA values:")
    if (is.null(CRISPR_df)) {
      demo_df <- CRISPOR_df
    } else {
      demo_df <- cbind.data.frame(CRISPOR_df, CRISPR_df[, !(names(CRISPR_df) %in% names(CRISPOR_df))])
    }
    if (is.null(show_columns)) {
      show_columns <- intersect(specificity_demo_columns, colnames(demo_df))
    }
    print(demo_df[are_invalid, show_columns])
    message("\n")
    for (column_name in c("CRISPOR_CFD_specificity", "CRISPOR_MIT_specificity")) {
      CRISPOR_df[[column_name]] <- as.integer(ifelse(are_invalid, NA_character_, CRISPOR_df[[column_name]]))
    }
    CRISPOR_df[["CRISPOR_off_target_count"]][are_invalid] <- NA_integer_
  }
  return(CRISPOR_df)
}



MakeBedIDsFromRangesDf <- function(ranges_df) {
  results_vec <- paste0(ranges_df[["Chromosome"]],
                        ":",
                        ifelse(ranges_df[["Strand"]] == "+", ranges_df[["Start"]], ranges_df[["Start"]] - 3L) - 1L,
                        "-",
                        ifelse(ranges_df[["Strand"]] == "+", ranges_df[["End"]] + 3L, ranges_df[["End"]]),
                        ":",
                        ranges_df[["Strand"]]
                        )
  return(results_vec)
}




AddCRISPORBedData <- function(CRISPR_df, CRISPOR_output_df, CRISPOR_offtargets_df, resolve_missing_offtargets = TRUE) {

  CRISPR_IDs_vec <- MakeBedIDsFromRangesDf(CRISPR_df)
  CRISPR_IDs_vec[is.na(CRISPR_df[["Start"]])] <- NA_character_

  stopifnot(!(anyNA(CRISPOR_output_df[["#seqId"]])))

  output_matches_vec <- match(CRISPR_IDs_vec, CRISPOR_output_df[["#seqId"]])
  output_matched_df <- CRISPOR_output_df[output_matches_vec, names(rename_CRISPOR_columns_vec)]
  names(output_matched_df) <- rename_CRISPOR_columns_vec

  output_matched_df <- ResolveSpecificityOfNone(output_matched_df, CRISPR_df)

  offtargets_results_df <- SummarizeOfftargets(CRISPOR_offtargets_df, is_FASTA = FALSE)
  offtarget_matches_vec <- match(CRISPR_IDs_vec, offtargets_results_df[["Location_ID"]])

  results_df <- data.frame(CRISPR_df,
                           output_matched_df,
                           offtargets_results_df[offtarget_matches_vec, ],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  results_df[["Location_ID"]] <- CRISPR_IDs_vec

  if (resolve_missing_offtargets) {
    results_df <- ResolveMissingOffTargets(results_df)
  }
  return(results_df)
}



AddCRISPORFASTAData <- function(CRISPR_df, CRISPOR_output_df, CRISPOR_offtargets_df, resolve_missing_offtargets = TRUE) {

  are_mapped <- !(is.na(CRISPR_df[["Start"]]))
  have_scores <- are_mapped #!(is.na(CRISPR_df[["CRISPOR_3MM_specificity"]]))

  ### DELETE THIS!! ###

  if (anyNA(CRISPOR_output_df[["targetSeq"]])) {
    stop("None of the entries in the targetSeq column may be NA, to avoid erroneous matches!")
  }

  if (any(have_scores & !(are_mapped))) {
    stop("A guide with no location, but a CRISPOR score, was unexpectedly found!")
  }
  ### DELETE THIS!! ###

  stopifnot(all(is.na(CRISPR_df[["CRISPOR_off_target_count"]][!(have_scores)])))

  were_not_scored <- !(have_scores) & are_mapped
  if (any(were_not_scored)) {
    message(paste0(sum(were_not_scored), " guides had no CRISPOR scores, ",
                   "even though they had a location! It will be attempted to find ",
                   "CRISPOR scores for the sgRNA + PAM sequence."
                   )
            )
  }

  sequences_vec <- rep(NA_character_, nrow(CRISPR_df))
  PAM_vec <- PAMorOriginalPAM(CRISPR_df[!(have_scores), ])
  sequences_vec[!(have_scores)] <- toupper(ifelse(is.na(PAM_vec),
                                                  NA_character_,
                                                  paste0(CRISPR_df[["sgRNA_sequence"]][!(have_scores)], PAM_vec)
                                                  )
                                           )

  output_matches_vec <- match(sequences_vec, toupper(CRISPOR_output_df[["targetSeq"]]))
  output_matched_df <- CRISPOR_output_df[output_matches_vec, names(rename_CRISPOR_columns_vec)]
  names(output_matched_df) <- rename_CRISPOR_columns_vec

  output_matched_df <- ResolveSpecificityOfNone(output_matched_df, CRISPR_df)

  offtargets_results_df <- SummarizeOfftargets(CRISPOR_offtargets_df, is_FASTA = TRUE)
  stopifnot(!(anyNA(offtargets_results_df[["Target_sequence"]])))
  offtarget_matches_vec <- match(sequences_vec, toupper(offtargets_results_df[["Target_sequence"]]))

  offtargets_df <- data.frame(output_matched_df,
                              offtargets_results_df[offtarget_matches_vec, ],
                              stringsAsFactors = FALSE,
                              row.names = NULL
                              )
  if (resolve_missing_offtargets) {
    offtargets_df <- ResolveMissingOffTargets(offtargets_df)
  }
  for (column_name in setdiff(names(offtargets_df), "Location_ID")) {
    CRISPR_df[[column_name]][!(have_scores)] <- offtargets_df[[column_name]][!(have_scores)]
  }
  return(CRISPR_df)
}



ResolveMissingOffTargets <- function(CRISPR_df, use_for_zero = 0.0001) {

  CFD_columns <- c("CRISPOR_4MM_specificity", "CRISPOR_3MM_specificity")

  lack_detailed_offtargets <- is.na(CRISPR_df[["CRISPOR_4MM_specificity"]]) &
                              !(is.na(CRISPR_df[["CRISPOR_CFD_specificity"]]))

  have_no_offtargets <- lack_detailed_offtargets & (CRISPR_df[["CRISPOR_CFD_specificity"]] %in% 100)
  if (any(have_no_offtargets)) {
    for (Num_MM_column in c(paste0("CRISPOR_Num_", 0:4, "MM"), "CRISPOR_Num_2or3MM")) {
      CRISPR_df[[Num_MM_column]][have_no_offtargets] <- 0L
    }
  }
  for (CFD_column in CFD_columns) {
    CRISPR_df[[CFD_column]][have_no_offtargets] <- 1
  }
  too_many_offtargets <- lack_detailed_offtargets & (CRISPR_df[["CRISPOR_CFD_specificity"]] %in% 0)
  for (CFD_column in CFD_columns) {
    CRISPR_df[[CFD_column]][too_many_offtargets] <- use_for_zero
  }
  return(CRISPR_df)
}





# Functions for filtering out sgRNAs with available CRISPOR data ----------

MakeBedIDsFromInputDf <- function(input_df) {
  paste0(input_df[["Chromosome"]], ":", input_df[["Start"]], "-", input_df[["End"]], ":", input_df[["Strand"]])
}


FilterDfList <- function(use_df_list, are_done_list) {
  postfix_vec <- vapply(are_done_list, function(x) if (all(x)) "all" else if (any(x)) "some" else "none", "")
  results_df_list <- Map(function(x, y) if (all(y)) x else x[!(y), ], use_df_list, are_done_list)
  names(results_df_list) <- paste0("chunk_", sub("chunk_", "", names(results_df_list)), "_", postfix_vec, "_done")
  return(results_df_list)
}


AddCombinedFilteredDf <- function(filtered_df_list) {
  are_to_combine <- !(grepl("_all_done$", names(filtered_df_list)))
  if (any(are_to_combine)) {
    combined_df_list <- CombineDfChunks(filtered_df_list[!(grepl("_all_done$", names(filtered_df_list)))])
    if (length(combined_df_list) == 1) {
      filtered_df_list[["filtered_all_chunks_combined"]] <- combined_df_list[[1]]
    } else  if (length(combined_df_list) != length(filtered_df_list)) {
      combined_names <- sub("chunk_chunk", "chunk", names(combined_df_list), fixed = TRUE)
      combined_names <- gsub("_some_donechunk_", "", combined_names, fixed = TRUE)
      combined_names <- sub("_some_done", "", combined_names, fixed = TRUE)
      combined_names <- paste0("combined_", combined_names)
      names(combined_df_list) <- combined_names
      filtered_df_list <- c(filtered_df_list, combined_df_list)
    }
  } else {
    message("CRISPOR scores for all entries are already available!")
  }
  return(filtered_df_list)
}


FilterBedDfList <- function(use_bed_df_list) {
  output_files_list <- IdentifyCRISPOROutputFiles()
  previous_CRISPOR_bed_df <- ReadCRISPOROutputFiles(output_files_list[["bed"]], is_FASTA = FALSE)
  were_processed_vec_list <- lapply(use_bed_df_list, function(x) MakeBedIDsFromInputDf(x) %in% previous_CRISPOR_bed_df[["#seqId"]])
  output_bed_df_list <- FilterDfList(use_bed_df_list, were_processed_vec_list)
  output_bed_df_list <- AddCombinedFilteredDf(output_bed_df_list)
  return(output_bed_df_list)
}


FilterFASTADfList <- function(use_FASTA_df_list) {
  sequence_vec_list <- lapply(use_FASTA_df_list, function(x) paste0(x[["sgRNA_sequence"]], x[["PAM"]]))
  output_files_list <- IdentifyCRISPOROutputFiles()
  previous_CRISPOR_FASTA_df <- ReadCRISPOROutputFiles(output_files_list[["FASTA"]], is_FASTA = TRUE)
  were_processed_vec_list <- lapply(sequence_vec_list, function(x) toupper(x) %in% toupper(previous_CRISPOR_FASTA_df[["targetSeq"]]))
  output_FASTA_df_list <- FilterDfList(use_FASTA_df_list, were_processed_vec_list)
  output_FASTA_df_list <- AddCombinedFilteredDf(output_FASTA_df_list)
  return(output_FASTA_df_list)
}





# Functions for comparing specificity scores ------------------------------

GetProportion <- function(logical_vec) {
  sum(logical_vec) / length(logical_vec)
}




