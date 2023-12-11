## 2023-10-17


# Functions for preparing data for between-screen comparisons -------------

ScatterInputDf <- function(logfc_1_df, logfc_2_df, choose_rep = NULL) {
  required_objects <- c("essentials_2020Q2_df", "non_essentials_2020Q2_df")
  stopifnot(all(required_objects %in% ls(envir = globalenv())))
  if (is.null(choose_rep)) {
    logfc_column <- "Mean_log2FC"
  } else {
    logfc_column <- paste0("Log2FC_rep", choose_rep)
  }
  common_entrezs <- intersect(
    logfc_1_df[, "Entrez_ID"][(!(is.nan(logfc_1_df[, logfc_column])))],
    logfc_2_df[, "Entrez_ID"][(!(is.nan(logfc_2_df[, logfc_column])))]
  )
  common_entrezs <- sort(common_entrezs)
  logfc_1_df <- logfc_1_df[match(common_entrezs, logfc_1_df[, "Entrez_ID"]), ]
  logfc_2_df <- logfc_2_df[match(common_entrezs, logfc_2_df[, "Entrez_ID"]), ]
  are_essential <- common_entrezs %in% essentials_2020Q2_df[, "Entrez_ID"]
  are_non_essential <- common_entrezs %in% non_essentials_2020Q2_df[, "Entrez_ID"]
  results_df <- data.frame(
    "Entrez_ID"      = common_entrezs,
    "Essentiality"   = ifelse(are_essential,
                              "essential",
                              ifelse(are_non_essential, "non-essential", "neither")
                              ),
    "Is_NT"          = FALSE,
    "Rep1_data"      = logfc_1_df[, logfc_column],
    "Rep2_data"      = logfc_2_df[, logfc_column],
    stringsAsFactors = FALSE
  )
  if ("Min_specificity" %in% names(logfc_1_df)) {
    results_df <- data.frame(
      results_df,
      "Rep1_min_spec"  = logfc_1_df[, "Min_specificity"],
      "Rep2_min_spec"  = logfc_2_df[, "Min_specificity"],
      "Rep1_comb_spec" = logfc_1_df[, "Combined_specificity"],
      "Rep2_comb_spec" = logfc_2_df[, "Combined_specificity"]
    )
  }
  return(results_df)
}



CommonRocDfList <- function(ROC_df_list) {
  common_entrezs <- Reduce(intersect, lapply(ROC_df_list, function(x) x[, "Entrez_ID"]))
  results_df_list <- lapply(ROC_df_list, function(x) {
    are_common <- x[, "Entrez_ID"] %in% common_entrezs
    x <- x[are_common, ]
    row.names(x) <- NULL
    return(x)
  })
  return(results_df_list)
}



CommonPlasmidsRocDfList <- function(logfc_df_list, first_plasmid_only = FALSE) {

  CheckPlasmids <- function(df_list) {
    stopifnot(length(unique(lapply(df_list, function(x) x[, "Plasmid_ID"]))) == 1)
  }

  CheckPlasmids(logfc_df_list)

  if (first_plasmid_only) {
    entrezs_vec <- logfc_df_list[[1]][, "Entrez_ID"]
    are_selected <- !(duplicated(entrezs_vec) | is.na(entrezs_vec))
    logfc_df_list <- lapply(logfc_df_list, function(x) {
      x <- x[are_selected, ]
      row.names(x) <- NULL
      return(x)
    })
  }
  non_NA_plasmids_list <- lapply(logfc_df_list, function(x) {
    are_non_NA <- !(is.na(x[, "Mean_log2FC"]))
    x[, "Plasmid_ID"][are_non_NA]
  })
  common_plasmids <- Reduce(intersect, non_NA_plasmids_list)
  matches_vec <- match(common_plasmids, logfc_df_list[[1]][, "Plasmid_ID"])
  corresponding_entrezs <- logfc_df_list[[1]][, "Entrez_ID"][matches_vec]
  are_selected <- !(duplicated(corresponding_entrezs) | is.na(corresponding_entrezs))
  selected_plasmids <- common_plasmids[are_selected]
  logfc_df_list <- lapply(logfc_df_list, function(x) {
    x <- x[x[, "Plasmid_ID"] %in% selected_plasmids, ]
    row.names(x) <- NULL
    return(x)
  })

  CheckPlasmids(logfc_df_list)

  return(logfc_df_list)
}



LogFcDfListToRocDfList <- function(logfc_df_list, essential_entrezs, non_essential_entrezs) {
  roc_input_df_list <- lapply(logfc_df_list, function(x) ROCInputDf(x, essential_entrezs, non_essential_entrezs))
  roc_df_list <- lapply(roc_input_df_list, function(x) ROCDfForColumn(x, "Mean_log2FC"))
  return(roc_df_list)
}



# Functions for plotting ROC curves ---------------------------------------

IsWholeNumber <- function(x, tolerance = .Machine[["double.eps"]]^0.9) {
  abs(x - round(x)) < tolerance
}


Log10Axis <- function(axis_side = 1, minor_ticks = 1:10, tick_mgp = 0.5, tick_length = 0.4) {

  raw_ticks <- axTicks(axis_side)
  raw_ticks <- raw_ticks[IsWholeNumber(raw_ticks)]

  extended_log10_ticks <- 10^c(raw_ticks[[1]] - 1, raw_ticks, raw_ticks[[length(raw_ticks)]] + 1)

  minor_ticks_list <- lapply(extended_log10_ticks, function(x) x * minor_ticks)
  minor_ticks_vec <- unique(unlist(minor_ticks_list))
  minor_ticks_vec <- minor_ticks_vec[!(minor_ticks_vec %in% 10^raw_ticks)]
  if (axis_side == 1) {
    axis_limits <- par("usr")[1:2]
  } else if (axis_side == 2) {
    axis_limits <- par("usr")[3:4]
  }
  axis_range <- axis_limits[[2]] - axis_limits[[1]]
  side_space <- axis_range * 0.01
  final_limits <- c(axis_limits[[1]] + side_space, axis_limits[[2]] - side_space)
  minor_ticks_vec <- minor_ticks_vec[(log10(minor_ticks_vec) > final_limits[[1]]) & (log10(minor_ticks_vec) < final_limits[[2]])]

  axis_labels <- formatC(10^raw_ticks, format = "e", digits = 0)
  axis_labels <- sub("e+0", "0^", axis_labels, fixed = TRUE)
  axis_labels <- sub("e-0", "0^-", axis_labels, fixed = TRUE)
  axis_labels <- parse(text = axis_labels)

  axis(axis_side, at = raw_ticks, labels = axis_labels,
       mgp = c(3, tick_mgp, 0), tcl = -(tick_length), lwd = NA, lwd.ticks = par("lwd")
       )

  axis(axis_side, at = log10(minor_ticks_vec), labels = rep(NA, length(minor_ticks_vec)),
       lwd = NA, lwd.ticks = par("lwd"), tcl = -(tick_length) * 0.4
       )
  return(invisible(NULL))
}



DrawBottomLegend <- function(labels_list,
                             use_colors,
                             border_colors        = NULL,
                             use_pch              = 16,
                             use_point_size       = 1.2,
                             lines_x_start        = 0.7,
                             lines_y_start        = 1.3,
                             y_mid                = 0.5,
                             small_gap_size       = 1.25,
                             large_gap_multiplier = 1.4,
                             use_lwd              = 2,
                             line_x_distance      = -0.2,
                             text_cex             = 1
                             ) {

  lines_x_start <- lines_x_start - diff(grconvertX(c(0, strwidth(expression(""^"gh"))), from = "user", to = "lines"))

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing legend
  small_gap <- diff(grconvertY(c(0, small_gap_size), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_multiplier

  text_list <- labels_list
  if (all(lengths(text_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(text_list))
    are_first <- rep(TRUE, length(text_list))
  } else {
    are_first <- unlist(lapply(text_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- total_span + diff(grconvertY(c(0, lines_y_start), from = "lines", to = "npc"))
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  all_expressions <- sapply(unlist(labels_list), VerticalAdjust)

  x_text  <- par("usr")[[2]] -
             diff(grconvertX(c(0, lines_x_start), from = "lines", to = "user")) -
             max(strwidth(all_expressions))

  ## Draw legend
  text(x      = x_text,
       y      = y_pos,
       labels = all_expressions,
       adj    = c(0, 0.5),
       cex    = text_cex,
       xpd    = NA
       )
  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  length_in_lines <- 0.5
  line_x_start <- x_text + diff(grconvertX(c(0, line_x_distance), from = "lines", to = "user"))
  segments(x0  = line_x_start,
           x1  = line_x_start + diff(grconvertX(c(0, length_in_lines), from = "lines", to = "user")),
           y0  = tapply(y_pos, groups_vec, mean),
           col = use_colors,
           lwd = par("lwd") * use_lwd,
           xpd = NA
           )

  return(invisible(NULL))
}



MultiLinesROC <- function(ROC_df_list,
                          flip             = TRUE,
                          embed_PNG        = FALSE,
                          only_annotation  = NULL,
                          transparency     = TRUE,
                          use_lwd          = 2,
                          legend_lwd       = use_lwd * 1.25,
                          legend_order     = seq_along(ROC_df_list),
                          x_axis_limits    = c(0, 1),
                          y_axis_limits    = NULL,
                          legend_inside    = TRUE,
                          middle_line      = !(legend_inside),
                          use_colors       = c(brewer.pal(9, "Greys")[[7]],
                                               brewer.pal(9, "Blues")[[6]],
                                               "#B8363F"
                                               ),
                          black_alpha      = 0.55,
                          colors_alpha     = 0.8,
                          long_labels      = FALSE,
                          GapFunction      = NULL,
                          x_label_line     = 1.75,
                          y_label_line     = 2.3,
                          legend_vec       = c("Re-analysis", "CRISPRoff", "T.gonfio"),
                          lines_y_start    = 1.3,
                          draw_legend      = TRUE,
                          AUC_num_digits   = 2,
                          log_FPR          = FALSE,
                          precision_recall = FALSE,
                          ...
                          ) {

  if (log_FPR && precision_recall) {
    stop("The 'logFPR' and 'precision_recall' options are mutually exclusive!")
  }

  if (is.null(y_axis_limits)) {
    if (log_FPR) {
      y_axis_limits <- c(0, 1)
    } else {
      y_axis_limits <- x_axis_limits
    }
  }

  if (transparency) {
    legend_colors <- Palify(use_colors, fraction_pale = 1 - colors_alpha)
    line_colors <- adjustcolor(use_colors, alpha.f = colors_alpha)
    if (length(unique(col2rgb(use_colors[[1]]))) == 1) {
      legend_colors[[1]] <- Palify("black", fraction_pale = 1 - black_alpha)
      line_colors[[1]] <- adjustcolor("black", alpha.f = black_alpha)
    }
  } else {
    legend_colors <- use_colors
    line_colors <- use_colors
  }

  if (precision_recall) {
    xy_mat_list <- lapply(ROC_df_list, function(x) {
      if (flip) {
        cbind("Recall"    = x[, "Specificity"],
              "Precision" = x[, "Precision_flipped"]
              )
      } else {
        cbind("Recall"    = x[, "Sensitivity"],
              "Precision" = x[, "Precision"]
              )
      }
    })
    AUC_vec <- vapply(ROC_df_list, function(x) {
      if ("Log2FC" %in% names(x)) {
        use_numeric_vec <- x[, "Log2FC"]
      } else {
        use_numeric_vec <- x[, "Mean_log2FC"]
      }
      AUCForROCdf(numeric_vec = use_numeric_vec,
                  logical_vec = x[, "Is_essential"],
                  precision_recall = TRUE,
                  flip = flip
                  )
    }, numeric(1))
  } else {
    ROC_mat_list <- lapply(ROC_df_list, GetROCMat)
    AUC_vec <- vapply(ROC_mat_list, GetAUC, numeric(1))

    if (flip) {
      for (i in seq_along(ROC_mat_list)) {
        sens_vec <- ROC_mat_list[[i]][, "Sensitivity"]
        spec_vec <- ROC_mat_list[[i]][, "Specificity"]
        ROC_mat_list[[i]][, "Sensitivity"] <- spec_vec
        ROC_mat_list[[i]][, "Specificity"] <- sens_vec
      }
    }

    xy_mat_list <- lapply(ROC_mat_list, function(x) {
      cbind("FPR"         = 1 - x[, "Specificity"],
            "Sensitivity" = x[, "Sensitivity"]
            )
    })

    if (log_FPR) {
      min_log10_FPR_vec <- vapply(xy_mat_list, function(x) {
        log10_FPR_vec <- log10(x[, "FPR"])
        min(log10_FPR_vec[is.finite(log10_FPR_vec)])
      }, numeric(1))
      x_axis_lower_limit <- min(min_log10_FPR_vec)
      x_axis_lower_limit <- x_axis_lower_limit - (abs(x_axis_lower_limit) * 0.025)
      x_axis_limits <- c(x_axis_lower_limit, 0)

      for (i in seq_along(xy_mat_list)) {
        if (log_FPR) {
          x_vec <- log10(xy_mat_list[[i]][, 1])
          x_vec[x_vec == "-Inf"] <- x_axis_limits[[1]]
          xy_mat_list[[i]][, 1] <- x_vec
        }
      }
    }
  }

  DrawAxes <- function() {
    use_tcl <- -0.36
    x_axis_mgp <- 0.4
    if (precision_recall) {
      x_axis_label <- "Recall"
      y_axis_label <- "Precision"
    } else {
      x_axis_label <- "False positive rate"
      y_axis_label <- "True positive rate"
    }
    if (log_FPR) {
      Log10Axis(1, tick_mgp = x_axis_mgp, tick_length = -(use_tcl))
      x_axis_label <- "False positive rate (log scale)"
    } else {
      axis(1, mgp = c(3, x_axis_mgp, 0), tcl = use_tcl, lwd = par("lwd"))
    }
    mtext(x_axis_label, side = 1, line = x_label_line, cex = par("cex"))
    axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1, lwd = par("lwd"))
    mtext(y_axis_label, side = 2, line = y_label_line, cex = par("cex"))
    box()
  }

  if (embed_PNG) {
    current_device <- StartEmbedPNG(figures_dir)
  }

  MakeEmptyPlot(x_axis_limits, y_axis_limits)

  if (!(isFALSE(only_annotation))) {
    if (middle_line && (!(log_FPR)) && (!(precision_recall))) {
      abline(a = 0, b = 1, col = "gray78", lty = "dashed")
    }
    if (!(is.null(GapFunction))) {
      GapFunction()
    }
    if (!(embed_PNG)) {
      DrawAxes()
    }
  }

  clip_lines <- (!(log_FPR)) && (length(unique(list(c(0, 1), x_axis_limits, y_axis_limits))) != 1)

  if (!(isTRUE(only_annotation))) {
    for (i in seq_along(xy_mat_list)) {
      lines(x   = xy_mat_list[[i]][, 1],
            y   = xy_mat_list[[i]][, 2],
            lwd = use_lwd * par("lwd"),
            col = line_colors[[i]],
            xpd = if (clip_lines) FALSE else NA
            )
    }
  }

  if (embed_PNG) {
    StopEmbedPNG(current_device, figures_dir)
    DrawAxes()
  }

  if (draw_legend && (!(isFALSE(only_annotation)))) {
    ## Draw legend
    if (length(AUC_vec) == 2) {
      legend_vec <- legend_vec[-1]
      if (long_labels) {
        legend_vec <- list(expression("CRISPRoff", "library"), expression("T.gonfio", "library"))
      }
    }
    AUC_legend_vec <- format(round(AUC_vec, digits = AUC_num_digits), nsmall = AUC_num_digits)
    AUC_legend_vec <- sapply(AUC_legend_vec, function(x) {
      as.expression(bquote("(AUC" * scriptscriptstyle(" ") * "=" * scriptscriptstyle(" ") * .(x) * ")"))
    })
    labels_list <- lapply(seq_along(legend_vec), function(x) c(as.expression(legend_vec[[x]]), AUC_legend_vec[[x]]))

    if (legend_inside) {
      DrawBottomLegend(labels_list = labels_list[legend_order],
                       use_colors = legend_colors[legend_order],
                       use_lwd = legend_lwd,
                       lines_y_start = lines_y_start,
                       ...
                       )
    } else {
      DrawSideLegend(labels_list = labels_list[legend_order],
                     use_colors = legend_colors[legend_order],
                     use_lwd = legend_lwd,
                     draw_lines = TRUE,
                     ...
                     )
    }
  }
  return(invisible(NULL))
}



# Plotting differences in discriminatory power between screens ------------

MeanSwarms <- function(rep_list, group_labels = c("CRISPRoff", "T.gonfio"), ...) {

  x_positions <- BeeViolinPlot(lapply(rep_list[c(1, 3, 2, 4)], function(x) x / 10),
                               violin_colors = rep(c(brewer.pal(9, "Purples")[[3]], "#c7e7c0"), each = 2),
                               point_colors  = rep(c("#7c7198", "#5b8669"), each = 2),
                               point_cex     = 0.25,
                               use_spacing   = 0.3,
                               adjust        = 1,
                               lower_bound   = -0.6,
                               upper_bound   = 0.2,
                               y_limits      = c(-0.608, 0.2),
                               groups_vec    = c(rep(1, 2), rep(2, 2)),
                               gap_ratio     = 1,
                               draw_groups_n = FALSE,
                               side_gap      = 0.45,
                               wex           = 0.95,
                               ...
                               )

  mtext(expression("Mean phenotype (" * gamma * ")"), side = 2, line = 2.1,
        cex = par("cex")
        )

  segments(x0  = x_positions[c(1, 3, 5, 7)] - 0.25,
           x1  = x_positions[c(2, 4, 6, 8)] + 0.25,
           y0  = par("usr")[[4]] + diff(grconvertY(c(0, 0.25), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(text = gsub("-", "\uad", c("essential\ngenes", "non-essential\ngenes"), fixed = TRUE),
        at = tapply(x_positions, rep(1:2, each = 2), mean),
        line = 0.4, padj = 0, cex = par("cex")
        )
  text(x      = x_positions + diff(grconvertY(c(0, 0.75), from = "lines", to = "user")),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.5), from = "lines", to = "user")),
       labels = rep(group_labels, 2),
       adj    = c(1, 0.4),
       srt    = 45,
       xpd    = NA
       )
  return(invisible(x_positions))
}



ComparePoints <- function(points_vec,
                          side_gap          = 0.5,
                          left_gap          = side_gap,
                          right_gap         = side_gap,
                          y_upper_limit     = 1.5,
                          y_axis_label      = "SSMD*",
                          use_tcl           = 0.375,
                          y_axis_label_line = 1.8,
                          y_axis_mgp        = 0.55,
                          y_axis_n          = 5,
                          group_labels_y    = 1.25,
                          points_color      = brewer.pal(9, "Blues")[[9]],
                          group_labels      = c("Re-analysis", "CRISPRoff", "T.gonfio")
                          ) {


  ## Determine point positions
  num_groups <- length(points_vec) / 2
  group_positions <- seq_len(num_groups)
  groups_vec <- rep(group_positions, each = 2)
  group_limits <- c((min(group_positions) - left_gap) - (num_groups * 0.04),
                    (max(group_positions) + right_gap) + (num_groups * 0.04)
                     )

  ## Prepare the data axis
  if (is.null(y_upper_limit)) {
    y_upper_limit <- max(pretty(c(0, max(points_vec))))
  }
  numeric_limits <- c(0, y_upper_limit)

  ## Draw lines
  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)
  use_width <- 0.15
  final_width <- use_width * ((max(group_positions) - min(group_positions)) / (num_groups - 1))
  means_vec <- tapply(points_vec, groups_vec, mean)
  segments(x0  = group_positions - final_width,
           x1  = group_positions + final_width,
           y0  = means_vec,
           lwd = 1
           )
  # SEMs <- tapply(points_vec, groups_vec, function(x) sd(x) / length(x))
  # segments(x0  = group_positions,
  #          y0  = means_vec - SEMs,
  #          y1  = means_vec + SEMs
  #          )
  # segments(x0  = group_positions - final_width / 3,
  #          x1  = group_positions + final_width / 3,
  #          y0  = means_vec - SEMs
  #          )
  # segments(x0  = group_positions - final_width / 3,
  #          x1  = group_positions + final_width / 3,
  #          y0  = means_vec + SEMs
  #          )

  ## Draw points
  points(x   = groups_vec,
         y   = points_vec,
         pch = 16,
         col = points_color,
         cex = 0.9
         )

  ## Draw the y axis
  tick_locations <- pretty(numeric_limits, n = y_axis_n)
  axis(2,
       at     = tick_locations,
       labels = format(tick_locations),
       las    = 2,
       mgp    = c(3, y_axis_mgp, 0),
       tcl    = -(use_tcl),
       lwd    = par("lwd")
       )
  if (!(is.null(y_axis_label))) {
    mtext(VerticalAdjust(y_axis_label),
          side = 2,
          line = y_axis_label_line,
          cex  = par("cex")
          )
  }

  if (length(group_positions) == 2) {
    group_labels <- group_labels[-1]
  }
  text(x      = group_positions + diff(grconvertX(c(0, 0.25), from = "lines", to = "user")),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.5), from = "lines", to = "user")),
       labels = group_labels,
       adj    = c(1, 0.5),
       srt    = 45,
       xpd    = NA
       )

  box(bty = "l")
}





# Functions for drawing custom scatter plots ------------------------------

RegressionScatter <- function(x_vec,
                              y_vec,
                              identity_line   = TRUE,
                              same_limits     = FALSE,
                              grid_lines      = FALSE,
                              zero_lines      = TRUE,
                              points_color    = "black",
                              points_alpha    = 0.5,
                              regression_line = TRUE,
                              confint_level   = 0.95,
                              band_color      = "#c0dcfc",
                              line_color      = "#0a5dbd",
                              show_axes       = TRUE,
                              GridFunction    = NULL
                              ) {

  ## Define axis limits
  if (same_limits) {
    x_limits <- range(x_vec, na.rm = TRUE)
    y_limits <- range(y_vec, na.rm = TRUE)
  } else {
    x_limits <- range(c(x_vec, y_vec), na.rm = TRUE)
    y_limits <- x_limits
  }

  if (regression_line) {
    ## Perform linear regression (and compute 95% confidence interval)
    model_df <- data.frame("x_var" = x_vec, "y_var" = y_vec)
    lm_model <- lm(y_var ~ x_var, data = model_df)
    lm_summary <- summary(lm_model)
    assign("delete_lm_summary", lm_summary, envir = globalenv())

    new_seq <- seq(min(x_vec), max(x_vec), length.out = 200)
    new_df <- data.frame("x_var" = new_seq)
    conf_int_mat <- predict(lm_model,
                            newdata = new_df,
                            interval = "confidence",
                            level = confint_level
                            )

    ## Prepare R-squared text
    r_squared <- format(round(lm_summary[["r.squared"]], digits = 2), nsmall = 2)
    p_value <- lm_summary[["coefficients"]][1, 4]
    scientific_split <- strsplit(formatC(p_value, format = "e", digits = 0),
                                 "e", fixed = TRUE
                                 )[[1]]
    power_of_10 <- as.integer(scientific_split[[2]])
    if (abs(power_of_10) <= 4) {
      corr_text <- bquote(italic("R") * ""^2  ~ "=" ~ .(r_squared) *
                          " (" * italic("p") * " = " * .(formatC(p_value, digits = 1, format = "fg")) * ")"
                          )
    } else {
      corr_text <- bquote(italic("R") * ""^2  ~ "=" ~ .(r_squared) *
                            " (" * italic("p") * " = " * .(scientific_split[[1]]) %*% 10^.(power_of_10) * ")"
                          )
    }
  }

  ## Set up plot region
  plot(NA,
       xlim = x_limits,
       ylim = y_limits,
       ann  = FALSE,
       axes = FALSE
       )

  if (!(is.null(GridFunction))) {
    GridFunction()
  }
  if (identity_line) {
    abline(a = 0, b = 1, col = "gray50")
  }
  if (zero_lines) {
    abline(h = 0, v = 0, col = "gray60")
  }

  if (regression_line) {
    ## Draw linear regression line and CI
    polygon(c(new_df[, 1], rev(new_df[, 1])),
            c(conf_int_mat[, 2], rev(conf_int_mat[, 3])),
            col = band_color, border = NA
            )
    lines(new_df[, 1], conf_int_mat[, 1], col = line_color)
    mtext(corr_text, line = 0.5, cex = par("cex"))
  }

  ## Draw points
  points(x_vec, y_vec, pch = 16, col = adjustcolor(points_color, points_alpha), cex = 0.5)

  ## Annotate plot
  if (show_axes) {
    axis(1, tcl = -0.375, mgp = c(3, 0.55, 0), lwd = par("lwd"), gap.axis = 0.25)
    axis(2, tcl = -0.375, mgp = c(3, 0.55, 0), las = 1, lwd = par("lwd"))
  }
  box()

  return(invisible(NULL))
}



GetDensityMat <- function(numeric_vec) {
  density_output <- density(numeric_vec)
  height_vec <- density_output[["y"]]
  values_vec <- density_output[["x"]]

  data_limits <- range(numeric_vec)
  are_within_bounds <- (values_vec >= data_limits[[1]]) &
    (values_vec <= data_limits[[2]])
  values_vec <- values_vec[are_within_bounds]
  height_vec <- height_vec[are_within_bounds]
  results_mat <- cbind(
    "value"  = values_vec,
    "height" = height_vec
  )
  return(results_mat)
}



GetQuantilesMat <- function(numeric_vec, show_quantiles = c(0.25, 0.5, 0.75)) {
  quantile_values <- quantile(numeric_vec, probs = show_quantiles)
  density_mat <- GetDensityMat(numeric_vec)
  quantile_heights <- stats::approxfun(density_mat[, "value"], density_mat[, "height"])(quantile_values)
  results_mat <- cbind(
    "quantile" = show_quantiles,
    "value"    = quantile_values,
    "height"   = quantile_heights
  )
  return(results_mat)
}



TwoCustomTrapezoids <- function(quantiles_list,
                                line_colors = c("#2870af", "#a51d7e"),
                                fill_colors = c(RColorBrewer::brewer.pal(9, "Blues")[[3]],
                                                RColorBrewer::brewer.pal(9, "RdPu")[[4]]
                                                ),
                                overlay     = FALSE
                                ) {

  stopifnot(all(c(length(quantiles_list), length(line_colors), length(fill_colors)) == 2))

  use_lwd <- 0.85

  x_vec <- c(quantiles_list[[1]][["25%"]],
             quantiles_list[[1]][["75%"]],
             quantiles_list[[1]][["75%"]],
             quantiles_list[[1]][["25%"]]
             )
  y_vec <- c(par("usr")[[4]],
             par("usr")[[4]],
             quantiles_list[[1]][["75%"]],
             quantiles_list[[1]][["25%"]]
             )
  polygon(x      = x_vec,
          y      = y_vec,
          border = NA,
          col    = adjustcolor(fill_colors[[1]], alpha.f = if (overlay) 0.1 else 0.5)
          )

  if (!(overlay)) {
    for (use_quantile in c("25%", "50%", "75%")) {
      segments(x0  = quantiles_list[[1]][[use_quantile]],
               y0  = par("usr")[[4]],
               y1  = quantiles_list[[1]][[use_quantile]],
               col = adjustcolor(line_colors[[1]], alpha.f = if (use_quantile == "50%") 0.9 else 0.6),
               lwd = par("lwd") * use_lwd,
               lty = "21"
               )
    }
  }

  x_vec <- c(par("usr")[[2]],
             par("usr")[[2]],
             quantiles_list[[2]][["25%"]],
             quantiles_list[[2]][["75%"]]
             )
  y_vec <- c(quantiles_list[[2]][["75%"]],
             quantiles_list[[2]][["25%"]],
             quantiles_list[[2]][["25%"]],
             quantiles_list[[2]][["75%"]]
             )
  polygon(x      = x_vec,
          y      = y_vec,
          border = NA,
          col    = adjustcolor(fill_colors[[2]], alpha.f = if (overlay) 0.1 else 0.5)
          )

  if (!(overlay)) {
    for (use_quantile in c("25%", "50%", "75%")) {
      segments(x0  = par("usr")[[2]],
               x1  = quantiles_list[[2]][[use_quantile]],
               y0  = quantiles_list[[2]][[use_quantile]],
               col = adjustcolor(line_colors[[2]], alpha.f = if (use_quantile == "50%") 0.9 else 0.6),
               lwd = par("lwd") * use_lwd,
               lty = "21"
               )
    }
  }
  return(invisible(NULL))
}



DensityTrapezoidScatter <- function(vec_list,
                                    trapezoid_fill_colors = c(RColorBrewer::brewer.pal(9, "Blues")[[3]],
                                                              RColorBrewer::brewer.pal(9, "RdPu")[[4]]
                                                              ),
                                    trapezoid_line_colors = c("#2870af", "#a51d7e"),
                                    swarm_fill_colors     = c(trapezoid_fill_colors[[1]],
                                                              Palify(trapezoid_fill_colors[[2]], fraction_pale = 0.45)
                                                              ),
                                    swarm_border_colors   = c(Palify(trapezoid_line_colors[[1]], fraction_pale = 0.1),
                                                              "#cc66af"
                                                              ),
                                    quantile_line_colors  = c(trapezoid_line_colors[[1]], "#a80079")
                                    ) {


  stopifnot(all(lengths(list(vec_list, trapezoid_fill_colors, trapezoid_line_colors,
                             swarm_fill_colors, swarm_border_colors, quantile_line_colors
                             )) == 2))


  tick_locations <- seq(-0.6, 0, by = 0.2)
  three_quantiles_list <- lapply(vec_list, function(x) quantile(x, probs = c(0.25, 0.5, 0.75)))

  UseGridFunction <- function() {
    abline(v = tick_locations, h = tick_locations, col = "gray86")
    TwoCustomTrapezoids(three_quantiles_list)
  }

  RegressionScatter(vec_list[[1]], vec_list[[2]], same_limits = TRUE,
                    points_color = "#4b357e",
                    regression_line = FALSE, show_axes = FALSE, zero_lines = FALSE,
                    GridFunction = UseGridFunction,
                    identity_line = FALSE
                    )

  TwoCustomTrapezoids(three_quantiles_list, overlay = TRUE)
  abline(a = 0, b = 1, col = adjustcolor("black", alpha.f = 0.4))
  axis(1, at = tick_locations, tcl = -0.35, mgp = c(3, 0.35, 0), lwd = par("lwd"))
  axis(2, at = tick_locations, tcl = -0.35, mgp = c(3, 0.5, 0), las = 1, lwd = par("lwd"))
  mtext(VerticalAdjust(expression("Prepool" ~ "(" * gamma * ")")),
        side = 1, line = 1.8, cex = par("cex")
        )
  mtext(VerticalAdjust(expression("Postpool" ~ "(" * gamma * ")")),
        side = 2, line = 2.1, cex = par("cex")
        )


  x_density_mat <- GetDensityMat(vec_list[[1]])
  y_density_mat <- GetDensityMat(vec_list[[2]])

  height_factor <- 1 / max(x_density_mat[, "height"], y_density_mat[, "height"])
  density_gap <- 0.25
  density_height <- 2.3
  beeswarm_spacing <- 0.35
  quantiles_lty <- c("21", "21", "21")
  beeswarm_cex <- 0.2

  start_pos <- par("usr")[[4]] + diff(grconvertY(c(0, density_gap), from = "lines", to = "user"))
  end_pos <- start_pos + diff(grconvertY(c(0, density_gap + density_height), from = "lines", to = "user"))
  height_range <- end_pos - start_pos
  height_scale <- height_factor * height_range
  values_vec <- x_density_mat[, "value"]
  polygon(x      = c(values_vec[[1]], values_vec, values_vec[[length(values_vec)]]),
          y      = start_pos + c(0, x_density_mat[, "height"] * height_scale, 0),
          col    = Palify(trapezoid_fill_colors[[1]], fraction_pale = 0.45),
          border = NA,
          xpd    = NA
          )

  numeric_vec <- vec_list[[1]]
  quantiles_mat <- GetQuantilesMat(numeric_vec)
  line_heights_vec <- quantiles_mat[, "height"] * height_scale - GetHalfLineWidth()

  for (i in seq_len(ncol(quantiles_mat))) {
    segments(x0   = quantiles_mat[, "value"][[i]],
             y0   = start_pos,
             y1   = start_pos + line_heights_vec[[i]],
             lty  = quantiles_lty[[i]],
             col  = quantile_line_colors[[1]],
             lend = "butt",
             xpd  = NA
             )
  }

  set.seed(1)
  swarm_df <- beeswarm(numeric_vec,
                       cex      = beeswarm_cex,
                       spacing  = beeswarm_spacing,
                       side     = 1,
                       priority = "random",
                       do.plot  = FALSE
                       )
  displacement_vec <- swarm_df[, "x"] - 1
  point_radius <- (par("cxy")[2] / pi) * par("cex") * beeswarm_cex
  points(x   = swarm_df[, "y"],
         y   = start_pos + point_radius + displacement_vec,
         cex = beeswarm_cex,
         pch = 21,
         col = swarm_border_colors[[1]],
         bg  = swarm_fill_colors[[1]],
         lwd = par("lwd") * 0.4,
         xpd = NA
         )


  start_pos <- par("usr")[[2]] + diff(grconvertX(c(0, density_gap), from = "lines", to = "user"))
  end_pos <- start_pos + diff(grconvertX(c(0, density_gap + density_height), from = "lines", to = "user"))
  height_range <- end_pos - start_pos
  height_scale <- height_factor * height_range * 0.8
  values_vec <- y_density_mat[, "value"]
  polygon(x      = start_pos + c(0, y_density_mat[, "height"] * height_scale, 0),
          y      = c(values_vec[[1]], values_vec, values_vec[[length(values_vec)]]),
          col    = Palify(trapezoid_fill_colors[[2]], fraction_pale = 0.45),
          border = NA,
          xpd    = NA
          )


  numeric_vec <- vec_list[[2]]
  quantiles_mat <- GetQuantilesMat(numeric_vec)
  line_heights_vec <- quantiles_mat[, "height"] * height_scale - GetHalfLineWidth()
  for (i in seq_len(ncol(quantiles_mat))) {
    segments(x0   = start_pos,
             x1   = start_pos + line_heights_vec[[i]],
             y0   = quantiles_mat[, "value"][[i]],
             lty  = quantiles_lty[[i]],
             col  = quantile_line_colors[[2]],
             lend = "butt",
             xpd  = NA
             )
  }

  set.seed(1)
  swarm_df <- beeswarm(numeric_vec,
                       cex      = beeswarm_cex,
                       spacing  = beeswarm_spacing,
                       side     = 1,
                       priority = "random",
                       do.plot  = FALSE
                       )
  displacement_vec <- swarm_df[, "x"] - 1
  points(x   = start_pos + point_radius + displacement_vec,
         y   = swarm_df[, "y"],
         cex = beeswarm_cex,
         pch = 21,
         col = swarm_border_colors[[2]],
         bg  = swarm_fill_colors[[2]],
         lwd = par("lwd") * 0.4,
         xpd = NA
         )

  return(invisible(NULL))
}


