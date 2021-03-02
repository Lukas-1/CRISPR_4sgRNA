### 21st October 2020 ###




# Define plot dimensions --------------------------------------------------

use_height <- 6.5
use_width <- 7




# Define functions --------------------------------------------------------

SideTextAndAxes <- function(side_text,
                            use_lwd            = 0.75,
                            use_cex_axis       = 0.8,
                            sg_label_cex       = 1.3,
                            horizontal_y_label = TRUE,
                            draw_outer_box     = FALSE
                            ) {

  # tick_locations <- axTicks(2)
  tick_locations <- seq(0, 1, 0.25)
  tick_labels <- paste0(tick_locations * 100, "%")
  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.38, 0),
       tcl      = -0.3,
       lwd      = use_lwd * par("lwd"),
       cex.axis = use_cex_axis,
       pos      = par("usr")[[1]] - (use_lwd * GetHalfLineWidth())
       )

  text(x      = if (horizontal_y_label) -0.06 else par("usr")[[1]] - diff(grconvertX(c(0, 3), from = "lines", to = "user")),
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




DrawAlterationBarplot <- function(summary_df,
                                  main_title         = NULL,
                                  reorder_wells      = FALSE,
                                  gap_weight         = 2L,
                                  space_height       = 1,
                                  top_space          = 3,
                                  bottom_space       = 1,
                                  show_color_text    = TRUE,
                                  show_color_legend  = FALSE,
                                  maintain_cex       = FALSE,
                                  use_lwd            = 0.75,
                                  use_cex_axis       = 0.8,
                                  sg_label_cex       = 1.3,
                                  horizontal_y_label = TRUE,
                                  left_space         = 0.13,
                                  right_space        = 0.07,
                                  draw_outer_box     = FALSE,
                                  exclude_blocks     = NULL
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
    text(x = 0.5,
         y = 0.65,
         labels = main_title,
         cex = 1.1,
         xpd = NA
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
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = horizontal_y_label,
                  draw_outer_box     = draw_outer_box
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
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = horizontal_y_label,
                  draw_outer_box     = draw_outer_box
                  )
  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[3]], four_colors, positions_vec)
  SideTextAndAxes("sg3",
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = horizontal_y_label,
                  draw_outer_box     = draw_outer_box
                  )
  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[4]], four_colors, positions_vec)
  SideTextAndAxes("sg4",
                  use_lwd            = use_lwd,
                  use_cex_axis       = use_cex_axis,
                  sg_label_cex       = sg_label_cex,
                  horizontal_y_label = horizontal_y_label,
                  draw_outer_box     = draw_outer_box
                  )
  MakeEmptyPlot()


  if (show_color_legend) {
    plot.window(c(0, 1), c(0, 1), "", asp = 1)
    MakeColorBoxLegend(labels_vec         = rev(four_alterations),
                       colors_vec         = rev(four_colors),
                       x_pos              = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.036),
                       y_pos              = 0.2,
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
                          exclude_blocks     = 2
                          )
    if (use_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}


