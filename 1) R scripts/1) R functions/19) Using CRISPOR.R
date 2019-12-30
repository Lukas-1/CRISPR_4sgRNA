### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")





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
  "Entrez_ID", "Gene_symbol", "Location_ID", "sgRNA_sequence", "PAM",
  "Num_0MM", "Num_1MM", "GuideScan_Num_2MM", "GuideScan_Num_3MM",
  "GuideScan_specificity", "CRISPOR_MIT_specificity", "CRISPOR_CFD_specificity", "CRISPOR_off_target_count",
  "CRISPOR_Doench_efficacy", "CRISPOR_Graf_status",
  "CRISPOR_4MM_specificity"
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
  rownames(export_bed_df) <- NULL

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
  rownames(results_df) <- NULL
  return(results_df)
}



MakeFASTAvec <- function(use_FASTA_df) {
  full_sequences_vec <- unique(paste0(use_FASTA_df[, "sgRNA_sequence"], use_FASTA_df[, "PAM"]))
  FASTA_list <- lapply(seq_along(full_sequences_vec), function(x) c(paste0(">seq", x), full_sequences_vec[[x]], ""))
  FASTA_vec <- unlist(FASTA_list)
  return(FASTA_vec)
}







# Functions for processing output from CRISPOR ----------------------------

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
                             list("Target_sequence" = offtargets_df[x[[1]], "targetSeq"]),
                             as.list(mismatch_table),
                             list(
                               "CRISPOR_3MM_specificity" = specificity_upto3MM,
                               "CRISPOR_4MM_specificity" = specificity_unrounded
                             )
                           ))
                         })

  results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  results_df[, "CRISPOR_Num_2or3MM"] <- rowSums(as.matrix(results_df[, c("CRISPOR_Num_2MM", "CRISPOR_Num_3MM")]))
  return(results_df)
}






AddCRISPORBedData <- function(CRISPR_df, CRISPOR_output_df, CRISPOR_offtargets_df, resolve_missing_offtargets = TRUE) {
  CRISPOR_output_df <- CRISPOR_output_df[CRISPOR_output_df[, "guideId"] %in% "21forw", ]

  CRISPR_IDs_vec <- paste0(CRISPR_df[, "Chromosome"],
                           ":",
                           ifelse(CRISPR_df[, "Strand"] == "+", CRISPR_df[, "Start"], CRISPR_df[, "Start"] - 3L) - 1L,
                           "-",
                           ifelse(CRISPR_df[, "Strand"] == "+", CRISPR_df[, "End"] + 3L, CRISPR_df[, "End"]),
                           ":",
                           CRISPR_df[, "Strand"]
                           )

  output_matches_vec <- match(CRISPR_IDs_vec, CRISPOR_output_df[, "#seqId"])
  output_matched_df <- CRISPOR_output_df[output_matches_vec, names(rename_CRISPOR_columns_vec)]
  colnames(output_matched_df) <- rename_CRISPOR_columns_vec

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
  sequences_vec[not_mapped] <- PAMorOriginalPAM(CRISPR_df[not_mapped, ])

  stopifnot(!(anyNA(CRISPOR_output_df[, "targetSeq"])))
  stopifnot(all(is.na(CRISPR_df[not_mapped, "CRISPOR_off_target_count"])))

  output_matches_vec <- match(sequences_vec, CRISPOR_output_df[, "targetSeq"])
  output_matched_df <- CRISPOR_output_df[output_matches_vec, names(rename_CRISPOR_columns_vec)]
  colnames(output_matched_df) <- rename_CRISPOR_columns_vec

  offtargets_results_df <- SummarizeOfftargets(CRISPOR_offtargets_df)
  stopifnot(!(anyNA(offtargets_results_df[, "Target_sequence"])))
  offtarget_matches_vec <- match(sequences_vec, offtargets_results_df[, "Target_sequence"])

  offtargets_df <- data.frame(output_matched_df,
                              offtargets_results_df[offtarget_matches_vec, ],
                              stringsAsFactors = FALSE,
                              row.names = NULL
                              )

  for (column_name in setdiff(colnames(offtargets_df), "Location_ID")) {
    CRISPR_df[not_mapped, column_name] <- offtargets_df[not_mapped, column_name]
  }

  return(CRISPR_df)
}







ResolveMissingOffTargets <- function(CRISPR_df) {

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
    CRISPR_df[too_many_offtargets, CFD_column] <- 0
  }

  return(CRISPR_df)
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
  # Plot CRISPOR vs. GuideScan specificity scores


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




