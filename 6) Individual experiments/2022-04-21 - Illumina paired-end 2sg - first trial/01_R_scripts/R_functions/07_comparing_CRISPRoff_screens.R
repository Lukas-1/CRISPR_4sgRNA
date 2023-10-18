## 2023-10-17



# Define functions --------------------------------------------------------

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



ThreeLinesROC <- function(ROC_df_list,
                          flip            = TRUE,
                          embed_PNG       = FALSE,
                          only_annotation = NULL,
                          transparency    = TRUE,
                          use_lwd         = 2,
                          legend_lwd      = use_lwd * 1.25,
                          legend_order    = seq_along(ROC_df_list),
                          axis_limits     = c(0, 1),
                          legend_inside   = TRUE,
                          middle_line     = !(legend_inside),
                          use_colors      = c(brewer.pal(9, "Greys")[[7]],
                                              brewer.pal(9, "Blues")[[6]],
                                              "#B8363F"
                                              ),
                          black_alpha     = 0.55,
                          colors_alpha    = 0.8,
                          long_labels     = FALSE,
                          GapFunction     = NULL,
                          x_label_line    = 1.75,
                          y_label_line    = 2.3,
                          legend_vec      = c("Re-analysis", "CRISPRoff", "T.gonfio"),
                          ...
                          ) {

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

  if (embed_PNG) {
    current_device <- StartEmbedPNG(figures_dir)
  }

  MakeEmptyPlot(axis_limits, axis_limits)
  if (!(isFALSE(only_annotation))) {
    if (middle_line) {
      abline(a = 0, b = 1, col = "gray78", lty = "dashed")
    }
    if (!(is.null(GapFunction))) {
      GapFunction()
    }
    if (!(embed_PNG)) {
      use_tcl <- -0.36
      axis(1, mgp = c(3, 0.4, 0), tcl = use_tcl, lwd = par("lwd"))
      mtext("False positive rate", side = 1, line = x_label_line, cex = par("cex"))
      axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1, lwd = par("lwd"))
      mtext("True positive rate", side = 2, line = y_label_line, cex = par("cex"))
      box()
    }
  }

  if (!(isTRUE(only_annotation))) {
    for (i in seq_along(ROC_mat_list)) {
      lines(x   = 1 - ROC_mat_list[[i]][, "Specificity"],
            y   = ROC_mat_list[[i]][, "Sensitivity"],
            lwd = use_lwd * par("lwd"),
            col = line_colors[[i]],
            xpd = NA
            )
    }
  }

  if (embed_PNG) {
    StopEmbedPNG(current_device, figures_dir)
  }

  if (embed_PNG) {
    use_tcl <- -0.36
    axis(1, mgp = c(3, 0.4, 0), tcl = use_tcl, lwd = par("lwd"))
    mtext("False positive rate", side = 1, line = x_label_line, cex = par("cex"))
    axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1, lwd = par("lwd"))
    mtext("True positive rate", side = 2, line = y_label_line, cex = par("cex"))
    box()
  }


  if (!(isFALSE(only_annotation))) {
    ## Draw legend
    if (length(AUC_vec) == 2) {
      legend_vec <- legend_vec[-1]
      if (long_labels) {
        legend_vec <- list(expression("CRISPRoff", "library"), expression("T.gonfio", "library"))
      }
    }
    AUC_legend_vec <- format(round(AUC_vec, digits = 2), nsmall = 2)
    AUC_legend_vec <- sapply(AUC_legend_vec, function(x) {
      as.expression(bquote("(AUC" * scriptscriptstyle(" ") * "=" * scriptscriptstyle(" ") * .(x) * ")"))
    })
    labels_list <- lapply(seq_along(legend_vec), function(x) c(as.expression(legend_vec[[x]]), AUC_legend_vec[[x]]))
    assign("delete_labels_list", labels_list, envir = globalenv())

    if (legend_inside) {
      DrawBottomLegend(labels_list = labels_list[legend_order],
                       use_colors = legend_colors[legend_order],
                       use_lwd = legend_lwd,
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
  x_user <- grconvertX(x = x_text, from = "npc", to = "user")
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       labels = all_expressions,
       adj    = c(0, 0.5),
       cex    = text_cex,
       xpd    = NA
       )
  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  length_in_lines <- 0.5
  line_x_start <- x_user + diff(grconvertX(c(0, line_x_distance), from = "lines", to = "user"))
  segments(x0  = line_x_start,
           x1  = line_x_start + diff(grconvertX(c(0, length_in_lines), from = "lines", to = "user")),
           y0  = tapply(y_pos, groups_vec, mean),
           col = use_colors,
           lwd = par("lwd") * use_lwd,
           xpd = NA
           )

  return(invisible(NULL))
}




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
  mtext(text = c("essential\ngenes", "non-essential\ngenes"),
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
  return(invisible(NULL))
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


