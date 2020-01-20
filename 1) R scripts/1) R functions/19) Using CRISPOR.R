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

  are_selected <- (CRISPR_df[, "Combined_ID"] %in% combined_IDs) &
                  !(is.na(CRISPR_df[, "Start"]))

  CRISPR_bed_df <- CRISPR_df[are_selected, ]

  export_bed_df <- data.frame(
    CRISPR_bed_df[, "Chromosome", drop = FALSE],
    "Start"  = ifelse(CRISPR_bed_df[, "Strand"] == "+", CRISPR_bed_df[, "Start"], CRISPR_bed_df[, "Start"] - 3L) - 1L,
    "End"    = ifelse(CRISPR_bed_df[, "Strand"] == "+", CRISPR_bed_df[, "End"] + 3L, CRISPR_bed_df[, "End"]),
    "Names"  = ".",
    "Scores" = ".",
    CRISPR_bed_df[, "Strand", drop = FALSE],
    stringsAsFactors = FALSE
  )
  export_bed_df <- unique(export_bed_df)
  row.names(export_bed_df) <- NULL
  return(export_bed_df)
}




PAMorOriginalPAM <- function(CRISPR_df) {
  are_validated_PAMs <- !(is.na(CRISPR_df[, "Original_PAM"])) &
                        mapply(function(x, y) x %in% strsplit(y, "; ", fixed = TRUE)[[1]],
                               CRISPR_df[, "Original_PAM"],
                               CRISPR_df[, "PAM_0MM"]
                               )
  PAM_vec <- ifelse(is.na(CRISPR_df[, "PAM"]),
                    ifelse(are_validated_PAMs, CRISPR_df[, "Original_PAM"], NA_character_),
                    CRISPR_df[, "PAM"]
                    )
  return(PAM_vec)
}



MakeFASTADf <- function(CRISPR_df, combined_IDs) {
  are_selected <- (CRISPR_df[, "Combined_ID"] %in% combined_IDs) &
                  is.na(CRISPR_df[, "Start"])
  results_df <- CRISPR_df[are_selected, ]
  results_df <- results_df[results_df[, "Num_0MM"] > 0, ]
  PAM_vec <- PAMorOriginalPAM(results_df)
  results_df <- data.frame(
    results_df[, c("Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "sgRNA_sequence")],
    "PAM" = PAM_vec,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  results_df <- results_df[!(is.na(PAM_vec)), ]
  row.names(results_df) <- NULL
  return(results_df)
}



MakeFASTAvec <- function(use_FASTA_df) {
  full_sequences_vec <- unique(paste0(use_FASTA_df[, "sgRNA_sequence"], use_FASTA_df[, "PAM"]))
  FASTA_list <- lapply(seq_along(full_sequences_vec), function(x) c(paste0(">seq", x), full_sequences_vec[[x]], ""))
  FASTA_vec <- unlist(FASTA_list)
  return(FASTA_vec)
}


WriteCRISPORInputFiles <- function(CRISPOR_chunk_list, file_ending, CRISPOR_input_directory) {
  for (i in seq_along(CRISPOR_chunk_list)) {
    file_name <- paste0("Input_for_CRISPOR__chunk_",
                        names(CRISPOR_chunk_list)[[i]], file_ending
                        )
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
  output_files   <- grep("^CRISPOR_output_", list.files(CRISPOR_files_directory), value = TRUE)
  are_FASTA      <- grepl("FASTA", output_files, fixed = TRUE)
  are_offtargets <- grepl("offs\\.tsv", output_files)
  results_list <- list(
    "bed"              = output_files[(!(are_FASTA)) & !(are_offtargets)],
    "FASTA"            = output_files[are_FASTA & !(are_offtargets)],
    "bed_offtargets"   = output_files[(!(are_FASTA)) & are_offtargets],
    "FASTA_offtargets" = output_files[are_FASTA & are_offtargets]
  )
  return(results_list)
}


ReadCRISPOROutputFiles <- function(file_names, use_fread = TRUE, show_messages = FALSE) {
  df_list <- lapply(file_names, function(x) ReadCRISPOROutput(x, use_fread, show_messages))
  results_df <- do.call(rbind.data.frame, c(df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
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
                                    data.table = FALSE
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



SummarizeOfftargets <- function(offtargets_df) {
  offtargets_df <- offtargets_df[offtargets_df[, "guideId"] == "21forw", ]
  offtargets_df[, "mismatchCount"] <- as.ordered(offtargets_df[, "mismatchCount"])
  offtargets_df[, "cfdOfftargetScore"] <- as.numeric(ifelse(offtargets_df[, "cfdOfftargetScore"] == "None", NA, offtargets_df[, "cfdOfftargetScore"]))

  results_list <- tapply(seq_len(nrow(offtargets_df)),
                         factor(offtargets_df[, "seqId"], levels = unique(offtargets_df[, "seqId"])),
                         function(x) {
                           mismatch_table <- as.integer(table(offtargets_df[x, "mismatchCount"]))
                           names(mismatch_table) <- paste0("CRISPOR_Num_", 0:4, "MM")
                           specificity_unrounded <- 1 / (1 + sum(offtargets_df[x, "cfdOfftargetScore"], na.rm = TRUE))
                           specificity_upto3MM <- 1 / (1 + sum(offtargets_df[x, "cfdOfftargetScore"][offtargets_df[x, "mismatchCount"] %in% 0:3], na.rm = TRUE))
                           return(c(
                             list("Location_ID"     = offtargets_df[x[[1]], "seqId"]),
                             list("Target_sequence" = offtargets_df[x[[1]], "guideSeq"]),
                             as.list(mismatch_table),
                             list(
                               "CRISPOR_3MM_specificity" = specificity_upto3MM,
                               "CRISPOR_4MM_specificity" = specificity_unrounded
                             )
                           ))
                         })

  results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  assign("delete_results_df", results_df, envir = globalenv())
  results_df[, "CRISPOR_Num_2or3MM"] <- as.integer(rowSums(as.matrix(results_df[, c("CRISPOR_Num_2MM", "CRISPOR_Num_3MM")])))
  return(results_df)
}



ResolveSpecificityOfNone <- function(CRISPOR_df, CRISPR_df = NULL, show_columns = specificity_demo_columns) {
  are_invalid <- (CRISPOR_df[, "CRISPOR_CFD_specificity"] %in% "None") | (CRISPOR_df[, "CRISPOR_MIT_specificity"] %in% "None")
  if (any(are_invalid)) {
    message("\n")
    message("The following sgRNAs had specificity scores of 'None', which were replaced by NA values:")
    if (is.null(CRISPR_df)) {
      demo_df <- CRISPOR_df
    } else {
      demo_df <- cbind.data.frame(CRISPOR_df, CRISPR_df[, !(names(CRISPR_df) %in% names(CRISPOR_df))])
    }
    print(demo_df[are_invalid, show_columns])
    message("\n")
    for (column_name in c("CRISPOR_CFD_specificity", "CRISPOR_MIT_specificity")) {
      CRISPOR_df[, column_name] <- as.integer(ifelse(are_invalid, NA_character_, CRISPOR_df[, column_name]))
    }
    CRISPOR_df[are_invalid, "CRISPOR_off_target_count"] <- NA_integer_
  }
  return(CRISPOR_df)
}



MakeBedIDsFromRangesDf <- function(ranges_df) {
  results_vec <- paste0(ranges_df[, "Chromosome"],
                        ":",
                        ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "Start"], ranges_df[, "Start"] - 3L) - 1L,
                        "-",
                        ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] + 3L, ranges_df[, "End"]),
                        ":",
                        ranges_df[, "Strand"]
                        )
  return(results_vec)
}




AddCRISPORBedData <- function(CRISPR_df, CRISPOR_output_df, CRISPOR_offtargets_df, resolve_missing_offtargets = TRUE) {

  CRISPOR_output_df <- CRISPOR_output_df[CRISPOR_output_df[, "guideId"] %in% "21forw", ]
  CRISPR_IDs_vec <- MakeBedIDsFromRangesDf(CRISPR_df)
  CRISPR_IDs_vec[is.na(CRISPR_df[, "Start"])] <- NA_character_

  stopifnot(!(anyNA(CRISPOR_output_df[, "#seqId"])))

  output_matches_vec <- match(CRISPR_IDs_vec, CRISPOR_output_df[, "#seqId"])
  output_matched_df <- CRISPOR_output_df[output_matches_vec, names(rename_CRISPOR_columns_vec)]
  names(output_matched_df) <- rename_CRISPOR_columns_vec

  output_matched_df <- ResolveSpecificityOfNone(output_matched_df, CRISPR_df)

  offtargets_results_df <- SummarizeOfftargets(CRISPOR_offtargets_df)
  offtarget_matches_vec <- match(CRISPR_IDs_vec, offtargets_results_df[, "Location_ID"])

  results_df <- data.frame(CRISPR_df,
                           output_matched_df,
                           offtargets_results_df[offtarget_matches_vec, ],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  results_df[, "Location_ID"] <- CRISPR_IDs_vec

  if (resolve_missing_offtargets) {
    results_df <- ResolveMissingOffTargets(results_df)
  }
  return(results_df)
}



AddCRISPORFASTAData <- function(CRISPR_df, CRISPOR_output_df, CRISPOR_offtargets_df, resolve_missing_offtargets = TRUE) {
  CRISPOR_output_df <- CRISPOR_output_df[CRISPOR_output_df[, "guideId"] %in% "21forw", ]

  not_mapped <- is.na(CRISPR_df[, "Start"])
  sequences_vec <- rep(NA_character_, nrow(CRISPR_df))
  sequences_vec[not_mapped] <- toupper(paste0(CRISPR_df[not_mapped, "sgRNA_sequence"], PAMorOriginalPAM(CRISPR_df[not_mapped, ])))

  stopifnot(!(anyNA(CRISPOR_output_df[, "targetSeq"])))
  stopifnot(all(is.na(CRISPR_df[not_mapped, "CRISPOR_off_target_count"])))

  output_matches_vec <- match(sequences_vec, toupper(CRISPOR_output_df[, "targetSeq"]))
  output_matched_df <- CRISPOR_output_df[output_matches_vec, names(rename_CRISPOR_columns_vec)]
  names(output_matched_df) <- rename_CRISPOR_columns_vec

  output_matched_df <- ResolveSpecificityOfNone(output_matched_df, CRISPR_df)

  offtargets_results_df <- SummarizeOfftargets(CRISPOR_offtargets_df)
  stopifnot(!(anyNA(offtargets_results_df[, "Target_sequence"])))
  offtarget_matches_vec <- match(toupper(sequences_vec), toupper(offtargets_results_df[, "Target_sequence"]))

  offtargets_df <- data.frame(output_matched_df,
                              offtargets_results_df[offtarget_matches_vec, ],
                              stringsAsFactors = FALSE,
                              row.names = NULL
                              )
  if (resolve_missing_offtargets) {
    offtargets_df <- ResolveMissingOffTargets(offtargets_df)
  }
  for (column_name in setdiff(names(offtargets_df), "Location_ID")) {
    CRISPR_df[not_mapped, column_name] <- offtargets_df[not_mapped, column_name]
  }
  return(CRISPR_df)
}




ResolveMissingOffTargets <- function(CRISPR_df, use_for_zero = 0.0001) {

  CFD_columns <- c("CRISPOR_4MM_specificity", "CRISPOR_3MM_specificity")

  lack_detailed_offtargets <- is.na(CRISPR_df[, "CRISPOR_4MM_specificity"]) &
                              !(is.na(CRISPR_df[, "CRISPOR_CFD_specificity"]))

  have_no_offtargets <- lack_detailed_offtargets & (CRISPR_df[, "CRISPOR_CFD_specificity"] %in% 100)
  if (any(have_no_offtargets)) {
    for (Num_MM_column in c(paste0("CRISPOR_Num_", 0:4, "MM"), "CRISPOR_Num_2or3MM")) {
      CRISPR_df[have_no_offtargets, Num_MM_column] <- 0L
    }
  }
  for (CFD_column in CFD_columns) {
    CRISPR_df[have_no_offtargets, CFD_column] <- 1
  }
  too_many_offtargets <- lack_detailed_offtargets & (CRISPR_df[, "CRISPOR_CFD_specificity"] %in% 0)
  for (CFD_column in CFD_columns) {
    CRISPR_df[too_many_offtargets, CFD_column] <- use_for_zero
  }
  return(CRISPR_df)
}





# Functions for filtering out sgRNAs with available CRISPOR data ----------

MakeBedIDsFromInputDf <- function(input_df) {
  paste0(input_df[, "Chromosome"], ":", input_df[, "Start"], "-", input_df[, "End"], ":", input_df[, "Strand"])
}

FilterDfList <- function(use_df_list, are_done_list) {
  postfix_vec <- vapply(are_done_list, function(x) if (all(x)) "all" else if (any(x)) "some" else "none", "")
  results_df_list <- Map(function(x, y) if (all(y)) x else x[!(y), ], use_df_list, are_done_list)
  names(results_df_list) <- paste0(names(results_df_list), "_", postfix_vec, "_done")
  return(results_df_list)
}

FilterBedDfList <- function(use_bed_df_list) {
  output_files_list <- IdentifyCRISPOROutputFiles()
  previous_CRISPOR_bed_df <- ReadCRISPOROutputFiles(output_files_list[["bed"]])
  were_processed_vec_list <- lapply(use_bed_df_list, function(x) MakeBedIDsFromInputDf(x) %in% previous_CRISPOR_bed_df[, "#seqId"])
  output_bed_df_list <- FilterDfList(use_bed_df_list, were_processed_vec_list)
  return(output_bed_df_list)
}

FilterFASTADfList <- function(use_FASTA_df_list) {
  sequence_vec_list <- lapply(use_FASTA_df_list, function(x) paste0(x[, "sgRNA_sequence"], x[, "PAM"]))
  output_files_list <- IdentifyCRISPOROutputFiles()
  previous_CRISPOR_FASTA_df <- ReadCRISPOROutputFiles(output_files_list[["FASTA"]])
  were_processed_vec_list <- lapply(sequence_vec_list, function(x) x %in% previous_CRISPOR_FASTA_df[, "targetSeq"])
  output_FASTA_df_list <- FilterDfList(use_FASTA_df_list, were_processed_vec_list)
  return(output_FASTA_df_list)
}




# Functions for comparing specificity scores ------------------------------

GetProportion <- function(logical_vec) {
  sum(logical_vec) / length(logical_vec)
}


ConvertCFDScores <- function(numeric_vec) {
  1 / (1 + (10000 / numeric_vec) - 100)
}


SpecificityScatterPlot <- function(CRISPR_df,
                                   x_column                     = "GuideScan_specificity",
                                   y_column                     = "CRISPOR_CFD_specificity",
                                   x_label                      = "GuideScan specificity score",
                                   y_label                      = "Original CRISPOR CFD specificity score",
                                   identical_axes               = FALSE,
                                   mark_diagonal                = identical_axes,
                                   convert_CRISPOR_to_GuideScan = FALSE,
                                   convert_GuideScan_to_CRISPOR = FALSE,
                                   point_cex                    = 0.5,
                                   point_alpha                  = 0.6,
                                   custom_axis_limits           = NULL,
                                   show_title                   = NULL,
                                   add_jitter                   = FALSE
                                   ) {

  old_par <- par(mar = rep.int(5, 4))
  x_vec <- CRISPR_df[, x_column]
  y_vec <- CRISPR_df[, y_column]

  if (add_jitter) {
    x_vec <- jitter(x_vec, factor = 1.5)
    y_vec <- jitter(y_vec, factor = 1.5)
  }

  if (convert_CRISPOR_to_GuideScan) {
    y_vec <- ConvertCFDScores(y_vec)
  }

  if (convert_GuideScan_to_CRISPOR) {
    x_vec <- (100 / (100 + ((1 / x_vec) - 1))) * 100
  }

  if (identical_axes) {
    if (is.null(custom_axis_limits)) {
      axis_limits <- range(c(x_vec, y_vec), na.rm = TRUE)
    } else {
      axis_limits <- custom_axis_limits
    }
  }

  plot(x_vec,
       y_vec,
       las  = 1,
       mgp  = c(2.7, 0.55, 0),
       tcl  = -0.4,
       type = "n",
       xlab = x_label,
       ylab = y_label,
       xlim = if (identical_axes) axis_limits else NULL,
       ylim = if (identical_axes) axis_limits else NULL
       )

  if (mark_diagonal) {
    abline(a = 0, b = 1, col = "gray88", lwd = 0.5)
  }

  line_color <- "gray80"
  line_type <- "dotted"
  if (x_column %in% c("GuideScan_specificity", "CRISPOR_3MM_specificity")) {
    abline(v = 0.2, col = line_color, lty = line_type)
  }
  if (y_column == "CRISPOR_CFD_specificity") {
    abline(h = 80, col = line_color, lty = line_type)
  } else if (y_column == "CRISPOR_3MM_specificity") {
    abline(h = 0.2, col = line_color, lty = line_type)
  }
  box()

  alpha_hex <- substr(rgb(1, 1, 1, point_alpha), 8, 9)
  points(x_vec,
         y_vec,
         col = paste0(brewer.pal(9, "Blues")[[7]], alpha_hex),
         pch = 16,
         cex = point_cex
         )

  if (!(is.null(show_title))) {
    title(show_title, cex.main = par("cex") * 0.9)
  }
  par(old_par)
  return(invisible(NULL))
}






DrawAllSpecificityScatterPlots <- function(CRISPR_df, append_to_file_name) {

  # Additional plots (not selected for PDF)

  SpecificityScatterPlot(CRISPR_df, point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE)

  SpecificityScatterPlot(CRISPR_df,
                         x_column                     = "CRISPOR_4MM_specificity", # This is just for confirmation
                         y_column                     = "CRISPOR_CFD_specificity",
                         x_label                      = "Unrounded CRISPOR CFD specificity score",
                         convert_GuideScan_to_CRISPOR = TRUE,
                         point_alpha                  = 0.2,
                         point_cex                    = 0.4
                         )

  SpecificityScatterPlot(CRISPR_df,
                         x_column       = "GuideScan_Num_2or3MM",
                         y_column       = "CRISPOR_Num_2or3MM",
                         x_label        = "GuideScan \u2013 number of 2MM or 3MM sites",
                         y_label        = "CRISPOR \u2013 number of 2MM or 3MM sites",
                         point_alpha    = 0.5,
                         point_cex      = 0.2,
                         identical_axes = TRUE
                         )

  SpecificityScatterPlot(CRISPR_df,
                         x_column       = "GuideScan_Num_2MM",
                         y_column       = "CRISPOR_Num_2MM",
                         x_label        = "GuideScan \u2013 number of 2MM sites",
                         y_label        = "CRISPOR \u2013 number of 2MM sites",
                         point_alpha    = 0.5,
                         point_cex      = 0.2,
                         identical_axes = TRUE
                         )



  # Most interesting plots -- used for PDF

  for (make_PDF in c(FALSE, TRUE)) {

    if (make_PDF) {
      plot_dimensions <- 5.75
      pdf(file = file.path(output_plots_directory, paste0("Specificity scores - CRISPOR vs. GuideScan - scatterplots - ", append_to_file_name, ".pdf")),
          width = plot_dimensions, height = plot_dimensions
      )
    }

    SpecificityScatterPlot(CRISPR_df,
                           y_column       = "CRISPOR_3MM_specificity",
                           y_label        = "CRISPOR CFD specificity score: max 3MM",
                           identical_axes = TRUE,
                           point_alpha    = 0.2,
                           point_cex      = 0.4,
                           show_title     = "CRISPOR CFD score (up to 3MM) vs. GuideScan score"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           y_column       = "CRISPOR_4MM_specificity",
                           y_label        = "Modified CRISPOR CFD specificity score",
                           identical_axes = TRUE,
                           mark_diagonal  = FALSE,
                           point_alpha    = 0.2,
                           point_cex      = 0.4,
                           show_title     = "CRISPOR CFD score (unrounded) vs. GuideScan score"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           x_column       = "CRISPOR_3MM_specificity",
                           x_label        = "CRISPOR CFD specificity score: max 3MM",
                           y_column       = "CRISPOR_4MM_specificity",
                           y_label        = "CRISPOR CFD specificity score: max 4MM",
                           identical_axes = TRUE,
                           mark_diagonal  = FALSE,
                           point_alpha    = 0.2,
                           point_cex      = 0.4,
                           show_title     = "CRISPOR CFD scores: 3MM vs. 4MM"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           point_alpha = 0.2,
                           point_cex   = 0.4,
                           show_title  = "Original CRISPOR CFD score vs. GuideScan score"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           x_column           = "GuideScan_Num_2or3MM",
                           y_column           = "CRISPOR_Num_2or3MM",
                           x_label            = "GuideScan \u2013 number of 2MM or 3MM sites",
                           y_label            = "CRISPOR \u2013 number of 2MM or 3MM sites",
                           point_alpha        = 0.4,
                           point_cex          = 0.2,
                           identical_axes     = TRUE,
                           custom_axis_limits = c(0, 1000),
                           show_title         = "Number of mismatches (2-3MM): CRISPOR vs. GuideScan"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           x_column           = "GuideScan_Num_2or3MM",
                           y_column           = "CRISPOR_Num_2or3MM",
                           x_label            = "GuideScan \u2013 number of 2MM or 3MM sites",
                           y_label            = "CRISPOR \u2013 number of 2MM or 3MM sites",
                           point_alpha        = 0.4,
                           point_cex          = 0.2,
                           identical_axes     = TRUE,
                           custom_axis_limits = c(0, 200),
                           show_title         = "Number of mismatches (2-3MM): CRISPOR vs. GuideScan (zoomed in)"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           x_column           = "GuideScan_Num_2MM",
                           y_column           = "CRISPOR_Num_2MM",
                           x_label            = "GuideScan \u2013 number of 2MM sites",
                           y_label            = "CRISPOR \u2013 number of 2MM sites",
                           point_alpha        = 0.4,
                           point_cex          = 0.2,
                           identical_axes     = TRUE,
                           custom_axis_limits = c(0, 300),
                           show_title         = "Number of mismatches (2MM): CRISPOR vs. GuideScan"
                           )

    SpecificityScatterPlot(CRISPR_df,
                           x_column           = "GuideScan_Num_2MM",
                           y_column           = "CRISPOR_Num_2MM",
                           x_label            = "GuideScan \u2013 number of 2MM sites",
                           y_label            = "CRISPOR \u2013 number of 2MM sites",
                           point_alpha        = 0.4,
                           point_cex          = 0.2,
                           identical_axes     = TRUE,
                           custom_axis_limits = c(0, 50),
                           show_title         = "Number of mismatches (2MM): CRISPOR vs. GuideScan (zoomed in)",
                           add_jitter         = TRUE
                           )

    if (make_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}




