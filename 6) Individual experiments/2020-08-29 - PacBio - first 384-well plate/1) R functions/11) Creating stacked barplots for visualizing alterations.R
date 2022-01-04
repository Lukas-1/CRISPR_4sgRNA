### 21st October 2020 ###




# Define plot dimensions --------------------------------------------------

use_height <- 6.5
use_width <- 7




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
                            x_ticks_length        = -0.2
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
    } else {
      x_tick_locations <- seq(0, 1, 0.1)
    }
    x_tick_labels <- paste0(x_tick_locations * 100)
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
  text(x      = if (horizontal_y_label) horizontal_y_lab_pos else par("usr")[[1]] - diff(grconvertX(c(0, vertical_y_label_line), from = "lines", to = "user")),
       y      = 0.5,
       adj    = if (horizontal_y_label) c(1, 0.5) else 0.5,
       labels = side_text,
       xpd    = NA,
       cex    = sg_label_cex,
       srt    = if (horizontal_y_label) 0 else 90
       )
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


  ## Set up the layout

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




# Functions for creating "sand charts" ------------------------------------

SandPlotMatList <- function(summary_df, consider_promoters = FALSE) {
  four_columns <- c("Perc_at_least_1", "Perc_at_least_2", "Perc_at_least_3", "Perc_all_4")
  if (consider_promoters) {
    four_columns <- sub("Perc_", "Perc_pr_", four_columns, fixed = TRUE)
  }
  column_mat_list <- lapply(four_columns, function(x) {
    results_mat <- cbind(summary_df[[x]],
                         100 - summary_df[[x]]
                         )
    results_mat <- results_mat / 100
    colnames(results_mat) <- c(x, "Rest")
    results_mat <- results_mat[order(results_mat[, 1], na.last = FALSE), ]
    return(t(results_mat))
  })
  for (i in 2:4) {
    stopifnot(all(column_mat_list[[i - 1]][1, ] >= column_mat_list[[i]][1, ], na.rm = TRUE))
  }
  return(column_mat_list)
}




DrawReorderedSandPlots <- function(summary_df,
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


  ## Set up the layout

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
    use_cex <- use_cex / 0.66
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

  num_wells <- nrow(summary_df)
  positions_vec <- seq_len(num_wells)

  use_colors <- c("#2A3B6B", "#D2D1CF")
  barplot_colors <- c("#FCFDCB", "#FDAC87", "#DB678B", "#934E9B", "#4E4073")

  column_mat_list <- SandPlotMatList(summary_df)

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[1]], use_colors, positions_vec)
  SideTextAndAxes(expression("" >= "1 correct sgRNA"),
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = FALSE,
                  draw_outer_box     = draw_outer_box,
                  sparse_ticks       = sparse_ticks,
                  vertical_y_label_line = 2.8,
                  label_x_axis       = TRUE
                  )
  MakeEmptyPlot()


  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[2]], use_colors, positions_vec)
  SideTextAndAxes(expression("" >= "2 correct sgRNAs"),
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = FALSE,
                  draw_outer_box     = draw_outer_box,
                  sparse_ticks       = sparse_ticks,
                  vertical_y_label_line = 2.75,
                  label_x_axis       = TRUE
                  )
  MakeEmptyPlot()


  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[3]], use_colors, positions_vec)
  SideTextAndAxes(expression("" >= "3 correct sgRNAs"),
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = FALSE,
                  draw_outer_box     = draw_outer_box,
                  sparse_ticks       = sparse_ticks,
                  vertical_y_label_line = 2.8,
                  label_x_axis       = TRUE
                  )
  MakeEmptyPlot()



  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[4]], use_colors, positions_vec)
  SideTextAndAxes(expression("4 correct sgRNAs"),
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = FALSE,
                  draw_outer_box     = draw_outer_box,
                  sparse_ticks       = sparse_ticks,
                  vertical_y_label_line = 2.8,
                  label_x_axis       = TRUE
                  )
  MakeEmptyPlot()


  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[1]], c(barplot_colors[[2]], barplot_colors[[1]]), positions_vec)
  PlotBarplotMat(column_mat_list[[2]], c(barplot_colors[[3]], NA), positions_vec)
  PlotBarplotMat(column_mat_list[[3]], c(barplot_colors[[4]], NA), positions_vec)
  PlotBarplotMat(column_mat_list[[4]], c(barplot_colors[[5]], NA), positions_vec)
  SideTextAndAxes("Superimposed",
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  draw_outer_box     = draw_outer_box,
                  sparse_ticks       = sparse_ticks,
                  horizontal_y_label = FALSE,
                  sg_label_cex       = sg_label_cex,
                  vertical_y_label_line = 2.8,
                  label_x_axis       = TRUE
                  )

  MakeEmptyPlot()

  plot.window(c(0, 1), c(0, 1), "", asp = 1)
  par("cex" = par("cex") * 0.9)
  MakeColorBoxLegend(labels_vec         = as.character(0:4),
                     colors_vec         = barplot_colors,
                     x_pos              = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.16),
                     y_pos              = 0.28,
                     use_constant_space = FALSE,
                     vertical_adjust    = 0.2,
                     x_space_adjust     = 4,
                     after_text         = " correct gRNAs",
                     use_lwd            = use_lwd * 1.5
                     )


  par(old_par)
  return(invisible(NULL))
}



SingleSandPlot <- function(summary_df,
                           consider_promoters = TRUE,
                           show_grid          = TRUE,
                           invert_x_axis      = TRUE,
                           x_axis_label       = "% wells"
                           ) {


  ## Prepare the plot

  num_wells <- nrow(summary_df)
  positions_vec <- seq_len(num_wells)

  barplot_colors <- c("#FCFDCB", "#FDAC87", "#DB678B", "#934E9B", "#4E4073")

  column_mat_list <- SandPlotMatList(summary_df)

  num_NA <- sum(is.na(summary_df[, "Perc_at_least_1"]))
  NA_vec <- c(rep(1, num_NA), rep(0, num_wells - num_NA))
  NA_mat <- rbind("Fraction_NA" = NA_vec, "Fraction_present" = abs(1 - NA_vec))

  if (invert_x_axis) {
    rev_indices <- rev(seq_len(num_wells))
    column_mat_list <- lapply(column_mat_list, function(x) x[, rev_indices])
    NA_mat <- NA_mat[, rev_indices]
  }


  ## Prepare the grid lines

  four_columns <- c("Perc_at_least_1", "Perc_at_least_2", "Perc_at_least_3", "Perc_all_4")
  if (consider_promoters) {
    four_columns <- sub("Perc_", "Perc_pr_", four_columns, fixed = TRUE)
  }
  barplot_mat <- do.call(rbind, lapply(four_columns, function(x) summary_df[[x]]))
  barplot_mat <- barplot_mat / 100
  barplot_mat <- rbind(
    "none"      = 1 - barplot_mat[1, ],
    "exactly_1" = barplot_mat[1, ] - barplot_mat[2, ],
    "exactly_2" = barplot_mat[2, ] - barplot_mat[3, ],
    "exactly_3" = barplot_mat[3, ] - barplot_mat[4, ],
    "all_4"     = barplot_mat[4, ]
  )


  ## Set up the layout

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

  PlotBarplotMat(column_mat_list[[1]], c(barplot_colors[[2]], barplot_colors[[1]]), positions_vec)
  PlotBarplotMat(column_mat_list[[2]], c(barplot_colors[[3]], NA), positions_vec)
  PlotBarplotMat(column_mat_list[[3]], c(barplot_colors[[4]], NA), positions_vec)
  PlotBarplotMat(column_mat_list[[4]], c(barplot_colors[[5]], NA), positions_vec)
  PlotBarplotMat(NA_mat,               c("white", NA),             positions_vec)

  if (show_grid) {
    grid_colors <- c("#FFFFFF", barplot_colors)
    grid_colors[1:4] <- vapply(grid_colors[1:4], Darken, factor = 1.05, "")
    grid_colors[5:6] <- vapply(grid_colors[5:6], Palify, fraction_pale = 0.1, "")

    pixels_per_line <- 101 # ceilGetHalfLineWidthing(line_width / (1 / num_wells)) + 1
    vert_line_width  <- GetHalfLineWidth(y_axis = FALSE) * 2
    horiz_line_width <- GetHalfLineWidth(y_axis = TRUE) * 2

    vert_line_width  <- vert_line_width * 0.85
    horiz_line_width <- horiz_line_width * 0.85

    vert_rect_width <- vert_line_width / pixels_per_line
    horiz_rect_width <- horiz_line_width / pixels_per_line


    for (use_fraction in seq(0.05, 0.95, by = 0.05)) {
      fractions_seq <- seq(from = use_fraction - (vert_line_width / 2) + (vert_rect_width / 2),
                           to   = use_fraction + (vert_line_width / 2) - (vert_rect_width / 2),
                           by   = vert_rect_width
                           )
      for (sub_fraction in fractions_seq) {
        segments_mat <- StartStopMat(VerticalLineLengths(column_mat_list, use_fraction = sub_fraction))
        for (i in seq_len(nrow(segments_mat))) {
          rect(xleft   = sub_fraction - (vert_rect_width / 2),
               xright  = sub_fraction + (vert_rect_width / 2),
               ybottom = segments_mat[i, 1],
               ytop    = segments_mat[i, 2],
               col     = rev(grid_colors)[[i]],
               border  = NA
               )
        }
      }
    }
    for (use_fraction in seq(0.05, 0.95, by = 0.05)) {
      fractions_seq <- seq(from = use_fraction - (horiz_line_width / 2) + (horiz_rect_width / 2),
                           to   = use_fraction + (horiz_line_width / 2) - (horiz_rect_width / 2),
                           by   = vert_rect_width
                           )
      for (sub_fraction in fractions_seq) {
        segments_mat <- StartStopMat(HorizontalLineLengths(column_mat_list, use_fraction = sub_fraction))
        if (invert_x_axis) {
          segments_mat <- 1 - segments_mat
        }
        for (i in seq_len(nrow(segments_mat))) {
          rect(xleft   = segments_mat[i, 1],
               xright  = segments_mat[i, 2],
               ybottom = sub_fraction - (horiz_rect_width / 2),
               ytop    = sub_fraction + (horiz_rect_width / 2),
               col     = grid_colors[[i]],
               border  = NA
               )
        }
      }
    }

  }

  SideTextAndAxes("% reads from each well",
                  use_lwd               = 1,
                  use_cex_axis          = 1,
                  many_ticks            = TRUE,
                  horizontal_y_label    = FALSE,
                  sg_label_cex          = 1,
                  vertical_y_label_line = 3.3,
                  label_x_axis          = TRUE,
                  x_axis_label          = x_axis_label,
                  x_ticks_length        = 0.3,
                  x_ticks_line          = 0.35,
                  x_label_line          = 2.5
                  )

  DrawOuterBox()
  MakeEmptyPlot()

  plot.window(c(0, 1), c(0, 1), "", asp = 1)
  MakeColorBoxLegend(labels_vec         = as.character(0:4),
                     colors_vec         = barplot_colors,
                     x_pos              = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.12),
                     y_pos              = 0.2,
                     use_constant_space = FALSE,
                     vertical_adjust    = 0.2,
                     x_space_adjust     = 4,
                     after_text         = " correct gRNAs",
                     use_lwd            = 1
                     )

  par(old_par)
  return(invisible(NULL))
}



VerticalLineLengths <- function(column_mat_list, use_fraction = 0.5) {

  lengths_mat <- cbind(
    "NA"        = ifelse(is.na(column_mat_list[[1]][2, ]), 1, 0),
    "None"      = column_mat_list[[1]][2, ],
    "Exactly_1" = column_mat_list[[1]][1, ] - column_mat_list[[2]][1, ],
    "Exactly_2" = column_mat_list[[2]][1, ] - column_mat_list[[3]][1, ],
    "Exactly_3" = column_mat_list[[3]][1, ] - column_mat_list[[4]][1, ],
    "Exactly_4" = column_mat_list[[4]][1, ]
  )
  lengths_mat[is.na(lengths_mat)] <- 0
  stopifnot(all(abs(rowSums(lengths_mat) - 1) < sqrt(.Machine$double.eps)))

  num_wells <- ncol(column_mat_list[[1]])
  use_well <- round(num_wells * use_fraction)
  lengths_vec <- lengths_mat[use_well, ]
  lengths_vec <- rev(lengths_vec)
  return(lengths_vec)
}



HorizontalLineLengths <- function(column_mat_list, use_fraction = 0.5) {

  perc_mat <- t(do.call(rbind, lapply(column_mat_list, function(x) x[1, , drop = FALSE])))
  total_vec <- apply(perc_mat, 2, function(x) sum(x > use_fraction, na.rm = TRUE))
  total_vec <- total_vec / nrow(perc_mat)

  NA_fraction <- sum(is.na(perc_mat[, 1])) / nrow(perc_mat)

  lengths_vec <- c("NA"        = NA_fraction,
                   "None"      = 1 - total_vec[[1]] - NA_fraction,
                   "Exactly_1" = total_vec[[1]] - total_vec[[2]],
                   "Exactly_2" = total_vec[[2]] - total_vec[[3]],
                   "Exactly_3" = total_vec[[3]] - total_vec[[4]],
                   "Exactly_4" = total_vec[[4]]
                   )

  assign("delete_lengths_vec", lengths_vec, envir = globalenv())
  stopifnot((sum(lengths_vec) - 1) < sqrt(.Machine$double.eps))

  return(lengths_vec)
}



StartStopMat <- function(lengths_vec) {
  cumsum_vec <- cumsum(lengths_vec)
  results_mat <- cbind(
    "from" = c(0, cumsum_vec[seq_len(length(lengths_vec) - 1)]),
    "to"   = cumsum_vec
  )
  rownames(results_mat) <- names(lengths_vec)
  return(results_mat)
}



