### 21st October 2020 ###



# Define plot dimensions --------------------------------------------------

use_height <- 6.5
use_width <- 7





# Define column selections and labels -------------------------------------

source(file.path(R_functions_directory, "01) Define titles and labels.R"))

titles_list <- c(
  titles_list,
  list(
    "All4_sg_cr_pass"    = "All 4 sgRNAs (+ tracrRNA) pass percentage threshold",
    "All4_pr_sg_cr_pass" = "All 4 sgRNAs (+ promoter) pass percentage threshold"
  )
)

titles_list[["Count_all_4_promoters"]] <- "Percentage of reads for which all 4 gRNAs (+ full promoters) are 100% correct"


eCDF_combos_list <- list(
  "4_guides" = list(
    "Count_all_4"       = expression(scriptscriptstyle(" ") * "All 4 gRNAs",
                                     scriptscriptstyle(" ") * "correct in",
                                     scriptscriptstyle(" ") * "same read"
                                     ),
    "All4_sg_cr_pass"   = expression(scriptscriptstyle(" ") * "Correct %",
                                     "" >=  "cutoff for",
                                     scriptscriptstyle(" ") * "all 4 gRNAs"
                                     )#,
    # "Count_mean_sg1to4" = c("Mean correct", "percentage", "(sg1-4)")
  ),
  "Contaminations" = list(
    "Num_contaminated_reads" = c("Contains", "gRNA from", "other wells")
  ),
  "Deletions" = list(
    "Num_reads_with_deletions_exceeding_20bp"     = expression("All deletions", "(" >= "20 bp)"),
    "Num_reads_with_deletions_spanning_tracrRNAs" = c("Deletions", "spanning", "tracrRNAs"),
    "Num_reads_with_deletions_spanning_promoters" = c("Deletions", "spanning", "promoters")
  )

)

eCDF_combos_list[["4_guides_pr"]] <- eCDF_combos_list[["4_guides"]]
names(eCDF_combos_list[["4_guides_pr"]]) <- c("Count_pr_all_4", "All4_pr_sg_cr_pass")#, "Count_mean_pr_sg1to4")

count_metrics <- c(
  "Count_sg1_cr1", "Count_sg2_cr2", "Count_sg3_cr3", "Count_sg4_cr4",
  "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
  "All4_sg_cr_pass", "Count_mean_sg1to4",

  "Count_pr1_sg1_cr1", "Count_pr2_sg2_cr2", "Count_pr3_sg3_cr3", "Count_pr4_sg4_cr4",
  "Count_pr_at_least_1", "Count_pr_at_least_2", "Count_pr_at_least_3", "Count_pr_all_4",
  "All4_pr_sg_cr_pass", "Count_mean_pr_sg1to4",

  "Count_all_4_promoters", "Count_whole_plasmid"
)

percentages_metrics <- c(
  "Num_contaminated_reads",
  "Num_contaminated_reads_aligned",
  "Num_cross_plate_contaminated",
  "Num_under_2kb",
  "Num_reads_with_deletions_exceeding_20bp",
  "Num_reads_with_sgRNA_deletion",
  "Num_reads_with_deletions_spanning_tracrRNAs",
  "Num_reads_with_deletions_spanning_promoters",
  count_metrics
)




# Helper functions --------------------------------------------------------

SideTextAndAxes <- function(side_text,
                            use_lwd               = 0.75,
                            use_cex_axis          = 0.8,
                            sg_label_cex          = 1.3,
                            horizontal_y_label    = TRUE,
                            horizontal_y_lab_pos  = -0.06,
                            vertical_y_label_line = 3,
                            draw_outer_box        = FALSE,
                            many_ticks            = FALSE,
                            sparse_ticks          = TRUE,
                            label_x_axis          = FALSE,
                            x_axis_label          = "Percentage of sequenced genes",
                            x_label_line          = 1.55,
                            x_ticks_line          = 0,
                            x_ticks_length        = 0.2,
                            both_axes_the_same    = FALSE,
                            use_mtext             = FALSE
                            ) {
  if (many_ticks) {
    y_tick_locations <- seq(0, 1, 0.2)
  } else if (sparse_ticks) {
    y_tick_locations <- seq(0, 1, 0.5)
  } else {
    y_tick_locations <- seq(0, 1, 0.25)
  }
  y_tick_labels <- paste0(y_tick_locations * 100, "%")
  axis(2,
       labels   = y_tick_labels,
       at       = y_tick_locations,
       las      = 1,
       mgp      = c(3, 0.38, 0),
       tcl      = -0.3,
       lwd      = use_lwd * par("lwd"),
       cex.axis = use_cex_axis,
       pos      = par("usr")[[1]] - (use_lwd * GetHalfLineWidth())
       )
  if (label_x_axis) {
    if (sparse_ticks && !(many_ticks)) {
      x_tick_locations <- seq(0, 1, 0.2)
    } else if (both_axes_the_same) {
      x_tick_locations <- y_tick_locations
    } else {
      x_tick_locations <- seq(0, 1, 0.1)
    }
    if (both_axes_the_same) {
      x_tick_labels <- y_tick_labels
    } else {
      x_tick_labels <- paste0(x_tick_locations * 100)
    }
    axis(1,
         labels   = x_tick_labels,
         at       = x_tick_locations,
         las      = 1,
         mgp      = c(3, x_ticks_line, 0),
         tcl      = -x_ticks_length,
         lwd      = use_lwd * par("lwd"),
         cex.axis = use_cex_axis,
         pos      = par("usr")[[3]] - (use_lwd * GetHalfLineWidth(y_axis = TRUE))
         )
    text(x      = 0.5,
         y      = par("usr")[[3]] - diff(grconvertY(c(0, x_label_line), from = "lines", to = "user")),
         adj    = if (horizontal_y_label) c(1, 0.5) else 0.5,
         labels = x_axis_label,
         xpd    = NA,
         cex    = sg_label_cex
         )
  }
  if (use_mtext && !(horizontal_y_label)) {
    mtext(side_text, side = 2, line = vertical_y_label_line, cex = par("cex"))
  } else {
  text(x      = if (horizontal_y_label) horizontal_y_lab_pos else par("usr")[[1]] - diff(grconvertX(c(0, vertical_y_label_line), from = "lines", to = "user")),
       y      = 0.5,
       adj    = if (horizontal_y_label) c(1, 0.5) else 0.5,
       labels = side_text,
       xpd    = NA,
       cex    = sg_label_cex,
       srt    = if (horizontal_y_label) 0 else 90
       )
  }

  if (draw_outer_box) {
    DrawOuterBox()
  }
  return(invisible(NULL))
}





# Functions for creating alteration bar plots -----------------------------

DrawAlterationBarplot <- function(summary_df,
                                  main_title           = NULL,
                                  title_color          = "black",
                                  reorder_wells        = FALSE,
                                  gap_weight           = 2L,
                                  space_height         = 1.2,
                                  top_space            = 2.5,
                                  bottom_space         = 1,
                                  show_color_text      = TRUE,
                                  show_color_legend    = FALSE,
                                  maintain_cex         = FALSE,
                                  use_lwd              = 0.75,
                                  use_cex_axis         = 0.8,
                                  sg_label_cex         = 1.3,
                                  horizontal_y_label   = TRUE,
                                  left_space           = 0.13,
                                  right_space          = 0.07,
                                  draw_outer_box       = FALSE,
                                  exclude_blocks       = NULL,
                                  sparse_ticks         = FALSE,
                                  horizontal_y_lab_pos = -0.06,
                                  color_legend_x_pos   = 0.036,
                                  title_y_pos          = 0.65,
                                  color_box_y_pos      = 0.2
                                  ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  percent_columns <- paste0("Perc_sg", 1:4, "_cr", 1:4)
  if (reorder_wells) {
    mean_accuracies <- rowMeans(as.matrix(summary_df[, percent_columns]))
    new_order <- order(summary_df[["Perc_all_4"]],
                       summary_df[["Perc_at_least_3"]],
                       summary_df[["Perc_at_least_2"]],
                       summary_df[["Perc_at_least_1"]],
                       mean_accuracies
                       )
  } else {
    new_order <- seq_len(nrow(summary_df))
  }

  if ("Empty_well" %in% names(sg_sequences_df)) {
    are_to_include <- !(sg_sequences_df[["Empty_well"]])
  } else {
    are_to_include <- rep(TRUE, nrow(sg_sequences_df))
  }
  have_blocks <- "Block" %in% names(sg_sequences_df)
  if (have_blocks && !(is.null(exclude_blocks))) {
    are_to_include[are_to_include] <- !(sg_sequences_df[["Block"]][are_to_include] %in% exclude_blocks)
  }

  use_indices <- new_order[are_to_include[new_order]]
  summary_df <- summary_df[use_indices, ]
  row.names(summary_df) <- NULL

  num_wells <- length(use_indices)


  ## Set up the plot layout

  barplot_height <- 2.8

  layout(cbind(rep(1, 9),
               3:(9 + 3 - 1),
               rep(2, 9)
               ),
         widths  = c(left_space, (1 - left_space - right_space), right_space),
         heights = c(top_space,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     bottom_space
                     )
         )

  use_cex <- par("cex")
  if (maintain_cex) {
    use_cex <- use_cex / 0.66
  }

  old_par <- par(mar = rep(0, 4), cex = use_cex)

  for (i in 1:3) {
    MakeEmptyPlot()
  }
  if (!(is.null(main_title))) {
    text(x      = 0.5,
         y      = title_y_pos,
         labels = main_title,
         col    = title_color,
         cex    = 1.1,
         xpd    = NA
         )

  }

  four_colors <- c("#F9F4EC", "#EE442F", "#63ACBE", "#601A4A")
  four_alterations <- c("Correct", "Mutation", "Deletion", "Contamination")

  four_colors <- rev(four_colors)
  four_alterations <- rev(four_alterations)


  add_gap <- (!(reorder_wells)) && have_blocks
  if (add_gap) {
    positions_vec <- GappedPositionsVec(sg_sequences_df[["Block"]][are_to_include],
                                        gap_weight = gap_weight
                                        )
  } else {
    positions_vec <- seq_len(num_wells)
  }


  ## Draw a vertical barplot for the mutation categories

  column_name_list <- lapply(1:4,
                             function(x) {
                               paste0(four_alterations, "_sg", x, "_cr", x)
                             })

  column_mat_list <- lapply(column_name_list, function(x) {
    use_mat <- as.matrix(summary_df[, x])
    for (i in seq_len(ncol(use_mat))) {
      use_mat[, i] <- use_mat[, i] / summary_df[["Count_total"]]
    }
    return(t(use_mat))
  })

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[1]], four_colors, positions_vec)
  SideTextAndAxes("sg1",
                  use_lwd              = use_lwd,
                  use_cex_axis         = use_cex_axis,
                  sg_label_cex         = sg_label_cex,
                  horizontal_y_label   = horizontal_y_label,
                  draw_outer_box       = draw_outer_box,
                  sparse_ticks         = sparse_ticks,
                  horizontal_y_lab_pos = horizontal_y_lab_pos
                  )


  if (show_color_text) {
    color_text_vec <- c('bold(color1("% plasmids with") * ',
                        'color1(" ") * color2("mutations,") * ',
                        'color1(" ") * color3("deletions,") * ',
                        'color1(" or ") * color4("contaminations"))'
                        )

    use_colors <- c("#000000",
                    rev(four_colors)[2:4]
                    )

    color_indices <- seq_along(use_colors)
    text_indices <- seq_along(color_text_vec)
    color_text <- paste0(color_text_vec[text_indices], collapse = "")

    for (j in color_indices) {
      text(x      = 0.5,
           y      = 1.25,
           labels = VerticalAdjust(parse(text = MakeInvisible(color_text, j))),
           adj    = c(0.5, 0),
           xpd    = NA,
           col    = use_colors[[j]],
           cex    = 1.2
           )
    }
  }

  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[2]], four_colors, positions_vec)
  SideTextAndAxes("sg2",
                  use_lwd              = use_lwd,
                  use_cex_axis         = use_cex_axis,
                  sg_label_cex         = sg_label_cex,
                  horizontal_y_label   = horizontal_y_label,
                  draw_outer_box       = draw_outer_box,
                  sparse_ticks         = sparse_ticks,
                  horizontal_y_lab_pos = horizontal_y_lab_pos
                  )
  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[3]], four_colors, positions_vec)
  SideTextAndAxes("sg3",
                  use_lwd              = use_lwd,
                  use_cex_axis         = use_cex_axis,
                  sg_label_cex         = sg_label_cex,
                  horizontal_y_label   = horizontal_y_label,
                  draw_outer_box       = draw_outer_box,
                  sparse_ticks         = sparse_ticks,
                  horizontal_y_lab_pos = horizontal_y_lab_pos
                  )
  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[4]], four_colors, positions_vec)
  SideTextAndAxes("sg4",
                  use_lwd              = use_lwd,
                  use_cex_axis         = use_cex_axis,
                  sg_label_cex         = sg_label_cex,
                  horizontal_y_label   = horizontal_y_label,
                  draw_outer_box       = draw_outer_box,
                  sparse_ticks         = sparse_ticks,
                  horizontal_y_lab_pos = horizontal_y_lab_pos
                  )
  MakeEmptyPlot()

  if (show_color_legend) {
    plot.window(c(0, 1), c(0, 1), "", asp = 1)
    MakeColorBoxLegend(labels_vec         = rev(four_alterations),
                       colors_vec         = rev(four_colors),
                       x_pos              = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * color_legend_x_pos),
                       y_pos              = color_box_y_pos,
                       use_constant_space = FALSE,
                       vertical_adjust    = 0,
                       x_space_adjust     = 4
                       )
  }
  par(old_par)
  layout(1)
  return(invisible(NULL))
}



ExportAlterationsForManuscript <- function(summary_df, use_prefix) {
  use_width <- 3.5
  use_height <- 4
  use_cex <- 0.7
  use_lwd <- 0.8
  for (use_PDF in TRUE) {
    file_name <- paste0(use_prefix,
                        " - alterations barplot - CCS7 (filtered) - SmrtLink 7",
                        ".pdf"
                        )
    if (use_PDF) {
      pdf(file = file.path(manuscript_directory, file_name),
          width = use_width, height = use_height
          )
    }
    par(lwd = use_lwd, cex = use_cex)
    DrawAlterationBarplot(summary_df,
                          space_height       = 0.5,
                          top_space          = 0.5,
                          bottom_space       = 1,
                          show_color_text    = FALSE,
                          maintain_cex       = FALSE,
                          use_lwd            = 1,
                          use_cex_axis       = 1,
                          sg_label_cex       = 1,
                          horizontal_y_label = FALSE,
                          show_color_legend  = TRUE,
                          left_space         = 0.15,
                          right_space        = 0.02,
                          draw_outer_box     = FALSE,
                          exclude_blocks     = 2,
                          sparse_ticks       = FALSE
                          )
    if (use_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}





# Basic functions for creating both eCDF plots and sand charts ------------

MakeSteps <- function(data_vec) {

  stopifnot(!(anyNA(data_vec)))

  num_obs <- length(data_vec)
  indices_vec <- seq_along(data_vec)

  step_size <- 1 / num_obs
  end_pos   <- step_size * indices_vec
  start_pos <- end_pos - step_size

  data_vec <- rep(data_vec, each = 2)
  fraction_vec <- c(rbind(start_pos, end_pos))

  ## Remove duplicates
  data_rle <- rle(data_vec)
  repeat_groups <- rep(seq_along(data_rle[["lengths"]]), data_rle[["lengths"]])
  are_intervening <- duplicated(repeat_groups) & duplicated(repeat_groups, fromLast = TRUE)
  data_vec <- data_vec[!(are_intervening)]
  fraction_vec <- fraction_vec[!(are_intervening)]

  results_mat <- cbind(
    "data_axis"     = data_vec,
    "fraction_axis" = fraction_vec
  )
  return(results_mat)
}





SingleSandAxes <- function(rotate_axes,
                           data_axis_label,
                           fraction_axis_label,
                           vertical_y_label_line = 3.3,
                           x_label_line = 2.5,
                           x_ticks_line = 0.35,
                           use_mtext = FALSE
                           ) {

  if (rotate_axes) {
    x_axis_label <- fraction_axis_label
    y_axis_label <- data_axis_label
  } else {
    x_axis_label <- data_axis_label
    y_axis_label <- fraction_axis_label
  }

  SideTextAndAxes(y_axis_label,
                  use_lwd               = 1,
                  use_cex_axis          = 1,
                  many_ticks            = TRUE,
                  horizontal_y_label    = FALSE,
                  sg_label_cex          = 1,
                  vertical_y_label_line = vertical_y_label_line,
                  label_x_axis          = TRUE,
                  x_axis_label          = x_axis_label,
                  x_ticks_length        = 0.3,
                  x_ticks_line          = x_ticks_line,
                  x_label_line          = x_label_line,
                  both_axes_the_same    = TRUE,
                  use_mtext             = use_mtext
                  )

  return(invisible(NULL))
}








# Functions for creating eCDF plots ---------------------------------------

ColumnToCDFVec <- function(summary_df, column_name) {

  ## Compute new metrics (not yet present in summary_df)
  new_columns <- c("Count_mean_sg1to4", "Count_mean_pr_sg1to4",
                   "All4_sg_cr_pass", "All4_pr_sg_cr_pass"
                   )
  if (column_name %in% new_columns) {
    if (column_name %in% new_columns[c(1, 3)]) {
      count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
    } else if (column_name %in% new_columns[c(2, 4)]) {
      count_columns <- paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4)
    }
    new_mat <- as.matrix(summary_df[, count_columns])
    if (column_name %in% new_columns[1:2]) {
      summary_df[[column_name]] <- rowMeans(new_mat)
    } else if (column_name %in% new_columns[3:4]) {
      summary_df[[column_name]] <- rowMins(new_mat)
    }
  }

  ## Prepare data for plotting. Exclude missing values (corresponding to wells
  ## with no reads) where appropriate, or replace them with 0.
  use_vec <- summary_df[, column_name] / summary_df[, "Count_total"]
  use_vec <- sort(use_vec, na.last = FALSE)
  if (grepl("mut|del|contam|num_under_2kb", column_name, ignore.case = TRUE)) {
    include_NA <- FALSE
  } else if (grepl("correct|count|all", column_name, ignore.case = TRUE)) {
    include_NA <- TRUE
  } else {
    stop("Unexpected value for 'column_name': ", column_name, "!")
  }
  are_NA <- is.na(use_vec)
  if (include_NA) {
    use_vec[are_NA] <- 0
  } else {
    use_vec <- use_vec[!(are_NA)]
  }

  return(use_vec)
}



FractionsForCutoffs <- function(summary_df, column_name, cutoffs_vec) {
  ## This is a convenience function for interactive use at the console
  sorted_vec <- ColumnToCDFVec(summary_df, column_name)
  results_vec <- vapply(cutoffs_vec, function(x) sum(sorted_vec >= x), integer(1))
  results_vec <- results_vec / length(sorted_vec)
  return(results_vec)
}


ValuesForQuantiles <- function(summary_df, column_name, quantiles_vec) {
  sorted_vec <- ColumnToCDFVec(summary_df, column_name)
  results_vec <- quantile(sorted_vec, probs = quantiles_vec)
  return(results_vec)
}




StepsAUC <- function(xy_mat) {
  xy_mat <- xy_mat[order(xy_mat[, "x"]), ]
  x_starts <- xy_mat[seq_len(nrow(xy_mat) - 1), "x"]
  x_ends <- xy_mat[2:nrow(xy_mat), "x"]
  rect_areas <- (x_ends - x_starts) * xy_mat[seq_len(nrow(xy_mat) - 1), "y"]
  auc_result <- sum(rect_areas)
  return(auc_result)
}



Plot_eCDF <- function(summary_df,
                      show_columns,
                      show_grid            = TRUE,
                      rotate_axes          = FALSE,
                      flip_axis            = FALSE,
                      data_axis_label      = "Percentage of reads from each well",
                      fraction_axis_label  = "Percentile of wells",
                      top_title            = NULL,
                      always_side_legend   = FALSE,
                      reverse_legend_order = FALSE,
                      reverse_colors       = FALSE,
                      use_brewer_pal       = "Blues",
                      line_colors          = NULL,
                      set_mar              = TRUE,
                      legend_y_gap         = 1.25,
                      legend_gap_ratio     = 1.75,
                      legend_x_start       = 1,
                      legend_segment_left  = legend_x_start - 0.3,
                      legend_segment_right = legend_x_start + 0.4,
                      legend_pch           = NULL,
                      point_x_start        = legend_x_start + 0.2,
                      lwd_multiplier       = 2,
                      ...
                      ) {

  if ((length(show_columns) == 1) && (show_columns %in% names(eCDF_combos_list))) {
    legend_labels_list <- eCDF_combos_list[[show_columns]]
    show_columns <- names(eCDF_combos_list[[show_columns]])
    predefined_combo <- TRUE
  } else {
    predefined_combo <- FALSE
  }

  ## Prepare the eCDF curves
  vec_list <- sapply(show_columns, function(x) ColumnToCDFVec(summary_df, x), simplify = FALSE)
  eCDF_mat_list <- lapply(vec_list, function(x) {
    steps_mat <- MakeSteps(x)
    ReorientSteps(steps_mat, rotate_axes = rotate_axes, flip_axis = flip_axis)
  })
  auc_vec <- vapply(eCDF_mat_list, StepsAUC, numeric(1))
  legend_order <- order(auc_vec, decreasing = !(reverse_legend_order))
  ## The curves that deviate most from the identity line should be plotted first
  ## (i.e. "bottom-most layer"):
  plot_order <- order(abs(auc_vec - 0.5), decreasing = TRUE)


  ## Set up the graphical parameters

  if (is.null(line_colors)) {
    if (length(show_columns) == 1) {
      line_colors <- brewer.pal(9, use_brewer_pal)[[7]]
    } else if (length(show_columns) == 2) {
      line_colors <- brewer.pal(9, use_brewer_pal)[c(9, 5)]
    } else {
      line_colors <- colorRampPalette(brewer.pal(9, use_brewer_pal)[9:4])(length(show_columns))
    }
  }

  space_fraction <- 0.0
  axis_limits <- c(-(space_fraction), 1 + space_fraction)

  show_side_legend <- always_side_legend || (length(show_columns) > 1)
  if (set_mar) {
    if (show_side_legend) {
      old_mar <- par("mar" = c(4.7, 5, 4.1, 10))
    } else {
      old_mar <- par("mar" = c(4.7, 5, 4.1, 3))
    }
  }

  ## Set up the plot region
  plot(1, xlim = axis_limits, ylim = axis_limits, xaxs = "i", yaxs = "i",
       axes = FALSE, ann = FALSE, typ = "n"
       )

  if (!(is.null(top_title))) {
    mtext(top_title, line = 0.3, cex = par("cex"))
  } else if (!(show_side_legend)) {
    title(titles_list[[show_columns]], cex.main = 1)
  }

  ## Draw the grid lines
  light_grey <- "gray95"
  darker_grey <- "gray88"
  abline(v = seq(0.05, 0.95, by = 0.1), col = light_grey)
  abline(h = seq(0.05, 0.95, by = 0.1), col = light_grey)
  abline(v = seq(0, 1, by = 0.1), col = darker_grey)
  abline(h = seq(0, 1, by = 0.1), col = darker_grey)


  ## Plot the x and y axes
  SingleSandAxes(rotate_axes, data_axis_label, fraction_axis_label, ...)
  DrawOuterBox()


  ## Plot the eCDF curve
  for (i in plot_order) {
    lines(eCDF_mat_list[[i]],
          lend  = "butt",
          ljoin = "mitre",
          lwd   = par("lwd") * lwd_multiplier,
          col   = line_colors[[i]],
          xpd   = NA
          )
  }


  if (show_side_legend) {
    labels_list <- lapply(show_columns, function(x) {
      if (predefined_combo) {
        legend_labels_list[[x]]
      } else {
        strwrap(titles_list[[x]], 21)
      }
    })
    DrawSideLegend(labels_list     = labels_list[legend_order],
                   use_colors      = line_colors[legend_order],
                   small_y_gap     = legend_y_gap,
                   large_gap_ratio = legend_gap_ratio,
                   lines_x_start   = legend_x_start,
                   point_x_start   = point_x_start,
                   use_pch         = legend_pch
                   )
  }


  if (set_mar) {
    par(old_mar)
  }

  return(invisible(NULL))
}




DrawSideLegend <- function(labels_list,
                           use_colors,
                           top_labels      = NULL,
                           use_pch         = NULL,
                           point_cex       = 1.5,
                           use_line_width  = 3,
                           lines_x_start   = 1,
                           lines_x_title   = lines_x_start - 1,
                           small_y_gap     = 1.25,
                           large_gap_ratio = 1.75,
                           segment_left    = lines_x_start - 0.3,
                           segment_right   = lines_x_start + 0.4,
                           point_x_start   = lines_x_start + 0.2
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing the legend
  y_mid <- 0.5
  small_gap <- diff(grconvertY(c(0, small_y_gap), from = "line", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_ratio
  have_multiline <- any(lengths(labels_list) > 1)

  if (have_multiline) {
    are_first <- unlist(lapply(labels_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  } else {
    gaps_vec <- rep(medium_gap, length(labels_list))
  }

  are_top_labels <- rep(FALSE, length(gaps_vec))
  if (!(is.null(top_labels))) {
    if (have_multiline) {
      are_first <- c(TRUE, rep(FALSE, length(top_labels) - 1))
      top_gaps_vec <- ifelse(are_first, large_gap, small_gap)
    } else {
      top_gaps_vec <- rep(medium_gap, length(top_labels))
    }
    are_top_labels <- c(rep(TRUE, length(top_labels)), are_top_labels)
    gaps_vec[[1]] <- (large_gap + 2 * medium_gap) / 3
    gaps_vec <- c(top_gaps_vec, gaps_vec)
    text_list <- c(as.list(top_labels), labels_list)
  } else {
    text_list <- labels_list
    are_top_labels <- rep(FALSE, length(gaps_vec))
  }
  gaps_vec[[1]] <- 0

  total_span <- sum(gaps_vec)

  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  LinesFromEdge <- function(num_lines) {
    par("usr")[[2]] + diff(grconvertX(c(0, num_lines), from = "lines", to = "user"))
  }

  ## Draw the legend
  text(x      = ifelse(are_top_labels,
                       LinesFromEdge(lines_x_title),
                       LinesFromEdge(lines_x_start)
                       ),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(text_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )

  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  assign("delete_groups_vec", groups_vec, envir = globalenv())
    assign("delete_y_pos", y_pos, envir = globalenv())
  assign("delete_are_top_labels", are_top_labels, envir = globalenv())

  groups_y_pos <- tapply(y_pos[!(are_top_labels)], groups_vec, mean)


  if (!(is.null(use_pch))) {
    points(x   = rep(LinesFromEdge(point_x_start), length(use_colors)),
           y   = groups_y_pos,
           col = "gray30",
           bg  = use_colors,
           pch = use_pch,
           cex = point_cex,
           xpd = NA
           )
  } else {
    segments(x0  = LinesFromEdge(segment_left),
             x1  = LinesFromEdge(segment_right),
             y0  = groups_y_pos,
             lwd = use_line_width * par("lwd"),
             col = use_colors,
             xpd = NA
             )
  }

  return(invisible(NULL))
}









# Functions for creating sand charts --------------------------------------

AddCorner <- function(input_mat) {
  if (input_mat[1, "data_axis"] == 0) {
    first_corner <- NULL
  } else {
    first_corner <- 0
  }
  cbind("data_axis"     = c(first_corner, input_mat[, "data_axis"],     0),
        "fraction_axis" = c(first_corner, input_mat[, "fraction_axis"], 1)
        )
}


DrawReorderedSandPlots <- function(summary_df,
                                   consider_promoters = FALSE,
                                   main_title         = NULL,
                                   title_color        = "black",
                                   space_height       = 1.2,
                                   top_space          = 1.8,
                                   bottom_space       = 2,
                                   maintain_cex       = FALSE,
                                   use_lwd            = 0.5,
                                   use_cex_axis       = 0.8,
                                   sg_label_cex       = 0.8,
                                   left_space         = 0.13 * 1.25,
                                   right_space        = 0.07 * 1.25,
                                   draw_outer_box     = FALSE,
                                   exclude_blocks     = NULL,
                                   sparse_ticks       = FALSE,
                                   use_cex            = par("cex")
                                   ) {


  ## Filter summary_df, if required
  if ("Empty_well" %in% names(sg_sequences_df)) {
    are_to_include <- !(sg_sequences_df[["Empty_well"]])
  } else {
    are_to_include <- rep(TRUE, nrow(sg_sequences_df))
  }
  have_blocks <- "Block" %in% names(sg_sequences_df)
  if (have_blocks && !(is.null(exclude_blocks))) {
    are_to_include[are_to_include] <- !(sg_sequences_df[["Block"]][are_to_include] %in% exclude_blocks)
  }
  summary_df <- summary_df[are_to_include, ]
  row.names(summary_df) <- NULL


  ## Prepare data for plotting

  two_colors <- c("#D2D1CF", "#2A3B6B")
  six_colors <- c("#FFFFFF", # NA
                  "#FCFDCB",  # 0 correct
                  "#FDAC87", "#DB678B", "#934E9B", "#4E4073" # 1-4 correct
                  )

  four_columns <- c("Perc_at_least_1", "Perc_at_least_2", "Perc_at_least_3", "Perc_all_4")
  if (consider_promoters) {
    four_columns <- sub("Perc_", "Perc_pr_", four_columns, fixed = TRUE)
  }
  column_mat <- as.matrix(summary_df[, four_columns]) / 100
  column_mat <- apply(column_mat, 2, sort, na.last = FALSE)
  column_mat_noNA <- column_mat
  column_mat_noNA[is.na(column_mat_noNA)] <- 0

  num_NA <- sum(is.na(column_mat[, 1]))
  num_wells <- nrow(summary_df)
  NA_fraction <- num_NA / num_wells


  ## Prepare graphical parameters
  two_colors <- c("#D2D1CF", "#2A3B6B")
  six_colors <- c("#FFFFFF", # NA
                  "#FCFDCB",  # 0 correct
                  "#FDAC87", "#DB678B", "#934E9B", "#4E4073" # 1-4 correct
                  )

  y_labels_list <- list(
    expression("" >= "1 correct sgRNA"),
    expression("" >= "2 correct sgRNAs"),
    expression("" >= "3 correct sgRNAs"),
    expression("4 correct sgRNAs")
  )


  ## Set up the plot layout

  barplot_height <- 2.8
  layout_mat <- cbind(rep(1, 11),
                      3:(11 + 3 - 1),
                      rep(2, 11)
                      )
  layout_mat[1, ] <- 3

  layout(layout_mat,
         widths  = c(left_space, (1 - left_space - right_space), right_space),
         heights = c(top_space,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     bottom_space
                     )
         )

  if (maintain_cex) {
    use_cex <- (use_cex * par("cex")) / 0.66
  }
  old_par <- par(mar = rep(0, 4), cex = use_cex)

  for (i in 1:3) {
    MakeEmptyPlot()
  }
  if (!(is.null(main_title))) {
    text(x      = 0.5,
         y      = 0.5,
         labels = main_title,
         col    = title_color,
         cex    = 1,
         xpd    = NA
         )

  }

  ## Draw the 4 two-color plots
  for (i in 1:4) {
    MakeEmptyPlot()
    basic_mat <- cbind("data_axis" = c(0, 1, 1), "fraction_axis" = c(0, 0, 1))
    polygon_mat <- SliceFraction(basic_mat, "fraction_axis", 0, NA_fraction)
    DrawPolygon(polygon_mat, "white",
                rotate_axes = TRUE, flip_axis = FALSE
                )

    polygon_mat <- SliceFraction(basic_mat, "fraction_axis", NA_fraction, 1)
    DrawPolygon(polygon_mat, two_colors[[1]],
                rotate_axes = TRUE, flip_axis = FALSE
                )
    steps_mat <- MakeSteps(column_mat_noNA[, i])
    DrawPolygon(AddCorner(steps_mat), two_colors[[2]],
                rotate_axes = TRUE, flip_axis = FALSE
                )
    SideTextAndAxes(y_labels_list[[i]],
                    use_lwd               = use_lwd,
                    use_cex_axis          = use_cex_axis,
                    sg_label_cex          = sg_label_cex,
                    horizontal_y_label    = FALSE,
                    draw_outer_box        = draw_outer_box,
                    sparse_ticks          = sparse_ticks,
                    vertical_y_label_line = 2.8,
                    label_x_axis          = TRUE
                    )
    MakeEmptyPlot()
  }


  ## Draw the four-color plot

  MakeEmptyPlot()
  polygon_mat <- SliceFraction(basic_mat, "fraction_axis", 0, NA_fraction)
  DrawPolygon(polygon_mat, six_colors[[1]],
              rotate_axes = TRUE, flip_axis = FALSE
              )
  polygon_mat <- SliceFraction(basic_mat, "fraction_axis", NA_fraction, 1)
  DrawPolygon(polygon_mat, six_colors[[2]],
              rotate_axes = TRUE, flip_axis = FALSE
              )
  for (i in 1:4) {
    steps_mat <- MakeSteps(column_mat_noNA[, i])
    DrawPolygon(AddCorner(steps_mat), six_colors[[i + 2]],
                rotate_axes = TRUE, flip_axis = FALSE
                )
  }
  SideTextAndAxes("Superimposed",
                  use_lwd               = use_lwd,
                  use_cex_axis          = use_cex_axis,
                  sg_label_cex          = sg_label_cex,
                  horizontal_y_label    = FALSE,
                  draw_outer_box        = draw_outer_box,
                  sparse_ticks          = sparse_ticks,
                  vertical_y_label_line = 2.8,
                  label_x_axis          = TRUE
                  )
  MakeEmptyPlot()

  ## Draw the legend
  plot.window(c(0, 1), c(0, 1), "", asp = 1)
  par("cex" = par("cex") * 0.9)
  MakeColorBoxLegend(labels_vec         = as.character(0:4),
                     colors_vec         = six_colors[2:6],
                     x_pos              = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.16),
                     y_pos              = 0.28,
                     use_constant_space = FALSE,
                     vertical_adjust    = 0.2,
                     x_space_adjust     = 4,
                     after_text         = " correct gRNAs",
                     use_lwd            = use_lwd * 1.5
                     )

  par(old_par)
  layout(1)
  return(invisible(NULL))
}


ReorientSteps <- function(input_mat, rotate_axes = FALSE, flip_axis = FALSE) {
  data_vec <- input_mat[, "data_axis"]
  fraction_vec <- input_mat[, "fraction_axis"]
  if (flip_axis) {
    fraction_vec <- 1 - fraction_vec
  }
  if (rotate_axes) {
    x_vec <- fraction_vec
    y_vec <- data_vec
  } else {
    x_vec <- data_vec
    y_vec <- fraction_vec
  }
  results_mat <- cbind("x" = x_vec, "y" = y_vec)
  return(results_mat)
}



DrawPolygon <- function(input_mat, fill_color = "black", rotate_axes = FALSE, flip_axis = FALSE) {
  xy_mat <- ReorientSteps(input_mat, rotate_axes = rotate_axes, flip_axis = flip_axis)
  polygon(xy_mat, col = fill_color, border = NA, xpd = NA, lend = "butt")
}



PolygonGrid <- function(polygon_mat, grid_color = "gray50", rotate_axes = FALSE, flip_axis = FALSE) {

  lines_seq <- seq(0.05, 0.95, by = 0.05)

  fraction_half_lwd  <- GetHalfLineWidth(y_axis = (!(rotate_axes))) * 0.6
  fraction_line_starts <- lines_seq - fraction_half_lwd
  fraction_line_ends   <- lines_seq + fraction_half_lwd

  data_half_lwd <- GetHalfLineWidth(y_axis = rotate_axes) * 0.6
  data_line_starts <- lines_seq - data_half_lwd
  data_line_ends   <- lines_seq + data_half_lwd

  for (i in seq_along(lines_seq)) {
    sliced_mat <- SliceFraction(polygon_mat, "fraction_axis", fraction_line_starts[[i]], fraction_line_ends[[i]])
    DrawPolygon(sliced_mat, rotate_axes = rotate_axes, flip_axis = flip_axis, fill_color = grid_color)
  }
  for (i in seq_along(lines_seq)) {
    sliced_mat <- SliceFraction(polygon_mat, "data_axis", data_line_starts[[i]], data_line_ends[[i]])
    DrawPolygon(sliced_mat, rotate_axes = rotate_axes, flip_axis = flip_axis, fill_color = grid_color)
  }
  return(invisible(NULL))
}



SliceFraction <- function(input_mat, slice_column, start_value, end_value, add_corners = TRUE) {

  slice_vec <- input_mat[, slice_column]
  is_sorted <- !(is.unsorted(round(slice_vec, digits = 12)))
  stopifnot(is_sorted)
  are_before <- slice_vec < start_value
  are_after <- slice_vec > end_value
  are_between <- !(are_before | are_after)
  sliced_mat <- input_mat[are_between, ]

  if (any(are_before)) {
    before_indices <- which(are_before)
    before_mat <- input_mat[before_indices[[length(before_indices)]], , drop = FALSE]
    before_mat[1, slice_column] <- start_value
    sliced_mat <- rbind(before_mat, sliced_mat)
  }
  if (any(are_after)) {
    after_mat <- input_mat[which(are_after)[[1]], , drop = FALSE]
    after_mat[1, slice_column] <- end_value
    sliced_mat <- rbind(sliced_mat, after_mat)
  }
  if (nrow(sliced_mat) == 0) {
    return(NULL)
  }

  if (add_corners) {
    if (slice_column == "fraction_axis") {
      corner_pos <- 0
    } else if (slice_column == "data_axis") {
      corner_pos <- 1
    } else {
      stop("Unexpected value for 'slice_column'!")
    }
    before_mat <- sliced_mat[1, , drop = FALSE]
    before_mat[, colnames(before_mat) != slice_column] <- corner_pos
    after_mat <- sliced_mat[nrow(sliced_mat), , drop = FALSE]
    after_mat[, colnames(after_mat) != slice_column] <- corner_pos
    sliced_mat <- rbind(before_mat, sliced_mat, after_mat)
  }
  return(sliced_mat)
}




SingleSandPlot <- function(summary_df,
                           consider_promoters  = FALSE,
                           show_grid           = TRUE,
                           invert_x_axis       = TRUE,
                           rotate_axes         = TRUE,
                           data_axis_label     = "Percentage of reads from each well",
                           fraction_axis_label = "Percentile of wells",
                           side_legend         = TRUE,
                           top_title           = NULL,
                           set_mar             = TRUE,
                           NA_in_legend        = TRUE, # Only implemented for the side legend
                           legend_y_gap        = 1,    # For the side legend
                           legend_x_start      = 1.5,  # For the side legend
                           legend_x_title      = 0,    # For the side legend
                           ...
                           ) {

  ## Prepare data for plotting

  four_columns <- c("Perc_at_least_1", "Perc_at_least_2", "Perc_at_least_3", "Perc_all_4")
  if (consider_promoters) {
    four_columns <- sub("Perc_", "Perc_pr_", four_columns, fixed = TRUE)
  }
  column_mat <- as.matrix(summary_df[, four_columns]) / 100
  column_mat <- apply(column_mat, 2, sort, na.last = FALSE)
  column_mat_noNA <- column_mat
  column_mat_noNA[is.na(column_mat_noNA)] <- 0

  num_NA <- sum(is.na(column_mat[, 1]))
  NA_fraction <- num_NA / nrow(column_mat)


  assign("delete_column_mat", column_mat_noNA, envir = globalenv())


  ## Prepare graphical parameters

  polygon_colors <- c("#FFFFFF", # NA
                      "#FCFDCB",  # 0 correct
                      "#FDAC87", "#DB678B", "#934E9B", "#4E4073" # 1-4 correct
                      )

  if (show_grid) {
    grid_colors <- rep(NA, 6)
    grid_colors[1:4] <- vapply(polygon_colors[1:4], Darken, factor = 1.05, "")
    grid_colors[5:6] <- vapply(polygon_colors[5:6], Palify, fraction_pale = 0.1, "")
  }


  if (side_legend) {
    if (set_mar) {
      old_par <- par("mar" = c(5, 4.1, 4, 6.7))
    }
    MakeEmptyPlot()
  } else {
    ## Set up the plot layout (this is necessary to make a legend with
    ## an aspect ratio of 1, see below.)
    layout_mat <- cbind(rep(1, 3),
                        3:(3 + 3 - 1),
                        rep(2, 3)
                        )
    layout(layout_mat,
           widths  = c(0.2, 0.8, 0.2),
           heights = c(0.1, 0.8, 0.3)
           )
    old_par <- par(mar = rep(0, 4), cex = par("cex") * (1 / 0.66))

    for (i in 1:4) {
      MakeEmptyPlot()
    }
    DrawOuterBox(fill = TRUE)
  }

  if (!(is.null(top_title))) {
    mtext(top_title, line = 0.3, cex = par("cex"))
  }


  ## Draw the "NA" rectangle
  basic_mat <- cbind("data_axis" = c(0, 1, 1), "fraction_axis" = c(0, 0, 1))
  polygon_mat <- SliceFraction(basic_mat, "fraction_axis", 0, NA_fraction)
  DrawPolygon(polygon_mat, polygon_colors[[1]],
              rotate_axes = rotate_axes, flip_axis = invert_x_axis
              )
  if (show_grid) {
    PolygonGrid(basic_mat, grid_colors[[1]],
                rotate_axes = rotate_axes, flip_axis = invert_x_axis
                )
  }

  ## Draw the "0 correct" rectangle
  polygon_mat <- SliceFraction(basic_mat, "fraction_axis", NA_fraction, 1)
  DrawPolygon(polygon_mat, polygon_colors[[2]],
              rotate_axes = rotate_axes, flip_axis = invert_x_axis
              )
  if (show_grid) {
    PolygonGrid(basic_mat, grid_colors[[2]],
                rotate_axes = rotate_axes, flip_axis = invert_x_axis
                )
  }

  ## Draw the "1-4 correct" polygon
  for (i in 1:4) {
    steps_mat <- MakeSteps(column_mat_noNA[, i])
    DrawPolygon(AddCorner(steps_mat), polygon_colors[[i + 2]],
                rotate_axes = rotate_axes, flip_axis = invert_x_axis
                )
    if (show_grid) {
      PolygonGrid(steps_mat, grid_colors[[i + 2]],
                  rotate_axes = rotate_axes, flip_axis = invert_x_axis
                  )
    }
  }

  ## Draw the axes
  SingleSandAxes(rotate_axes, data_axis_label, fraction_axis_label, ...)
  DrawOuterBox()

  if (side_legend) {
    legend_labels <- as.list(as.character(0:4))
    legend_colors <- polygon_colors[2:6]
    if (NA_in_legend) {
      legend_labels <- c(legend_labels, "No data")
      legend_colors <- c(legend_colors, polygon_colors[[1]])
    }
    DrawSideLegend(legend_labels, #as.list(as.character(0:4)),
                   use_colors = legend_colors,
                   use_pch = 22, small_y_gap = legend_y_gap,
                   lines_x_start = legend_x_start,
                   lines_x_title = legend_x_title,
                   top_labels = c("Correct", "gRNAs:")
                   )
  } else {
    MakeEmptyPlot()
    plot.window(c(0, 1), c(0, 1), "", asp = 1)
    MakeColorBoxLegend(labels_vec         = as.character(0:4),
                       colors_vec         = polygon_colors[2:6],
                       x_pos              = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.12),
                       y_pos              = 0.2,
                       use_constant_space = FALSE,
                       vertical_adjust    = 0.2,
                       x_space_adjust     = 4,
                       after_text         = " correct gRNAs",
                       use_lwd            = 1
                       )
    layout(1)
  }

  if (set_mar || (!(side_legend))) {
    par(old_par)
  }
  return(invisible(NULL))
}



