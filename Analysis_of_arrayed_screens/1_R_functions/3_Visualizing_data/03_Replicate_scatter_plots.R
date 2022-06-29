# 2021-12-27


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define functions --------------------------------------------------------

ScatterPlot <- function(x_vec,
                        y_vec,
                        use_limits      = NULL,
                        point_size      = 0.6,
                        label_y_axis    = TRUE,
                        use_color       = NULL,
                        top_label       = NULL,
                        draw_regression = TRUE,
                        use_mgp         = c(2.8, 0.55, 0),
                        x_axis_mgp      = c(2.2, use_mgp[2:3]),
                        xlab_line       = 2.2
                        ) {

  stopifnot(length(x_vec) == length(y_vec))

  ## Determine axis limits
  if (is.null(use_limits)) {
    use_limits = range(c(x_vec, y_vec))
  }
  use_limits <- use_limits + (diff(use_limits) * 0.04 * c(-1, 1))

  if (draw_regression) {
    ## Perform a linear regression (and compute a 95% confidence interval)
    model_df <- data.frame("x_var" = x_vec, "y_var" = y_vec)
    lm_model <- lm(y_var ~ x_var, data = model_df)
    lm_summary <- summary(lm_model)
    new_seq <- seq(use_limits[[1]], use_limits[[2]], length.out = 200)
    new_df <- data.frame("x_var" = new_seq)
    conf_int_mat <- predict(lm_model,
                            newdata = new_df,
                            interval = "confidence",
                            level = 0.95
                            )
    corr_text <- bquote(italic("R") * ""^2  ~ "=" ~
                        .(format(round(lm_summary[["r.squared"]], digits = 2), nsmall = 2))
                        )
  } else {
    pearsons_r <- cor.test(x_vec, y_vec)[["estimate"]][[1]]
    corr_text <- bquote(italic("r")  ~ "=" ~
                        .(format(round(pearsons_r, digits = 2), nsmall = 2))
                        )
  }

  ## Define graphical parameters
  use_tcl <- -0.35
  if (is.null(use_color)) {
    use_color <- "#000000"
  }

  ## Set up the plot region
  plot(1,
       xlim = use_limits,
       ylim = use_limits,
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       ann  = FALSE,
       type = "n"
       )

  ## Draw and label axes
  axis(1, mgp = x_axis_mgp, tcl = use_tcl, gap.axis = 0.5, lwd = par("lwd"))
  mtext(FormatPlotMath("Replicate 1"), side = 1, line = x_axis_mgp[[1]], cex = par("cex"))
  axis(2, las = 1, mgp = use_mgp, tcl = use_tcl, lwd = par("lwd"))
  if (label_y_axis) {
    mtext(FormatPlotMath("Replicate 2"), side = 2, line = use_mgp[[1]], cex = par("cex"))
  }

  ## Draw the plot title
  if (!(is.null(top_label))) {
    mtext(Embolden(VerticalAdjust(top_label)),
          line = 1.6, cex = par("cex"), font = 2
          )
  }
  mtext(VerticalAdjust(as.expression(corr_text)),
        line = 0.05, cex = par("cex"), font = 2
        )

  ## Draw indicator lines
  abline(a = 0, b = 1, col = "grey80", lty = "dashed")
  abline(h = 0, col = "gray90")
  abline(v = 0, col = "gray90")

  ## Draw the regression line and 95% CI region
  if (draw_regression) {
    polygon(c(new_df[, 1], rev(new_df[, 1])),
            c(conf_int_mat[, 2], rev(conf_int_mat[, 3])),
            col = Palify(use_color, fraction_pale = 0.8), border = NA
            )
    lines(new_df[, 1], conf_int_mat[, 1], col = use_color, lwd = 1.5)
  }
  box()

  ## Draw the points of the scatter plot
  points(x_vec,
         y_vec,
         pch = 16,
         col = adjustcolor(use_color, alpha.f = 0.5),
         cex = point_size * par("cex")
         )

  return(invisible(NULL))
}


ReplicateScatter <- function(input_df,
                             rep1_column,
                             rep2_column = NULL,
                             show_title = "Replicate scatter plot",
                             same_scale = TRUE,
                             ...
                             ) {

  if (is.null(rep2_column)) {
    rep2_column <- sub("rep1", "rep2", rep1_column, fixed = TRUE)
  }

  are_gene <- !(is.na(input_df[, "Entrez_ID"]))
  corr_gene <- cor.test(input_df[are_gene, rep1_column],
                        input_df[are_gene, rep2_column]
                        )[["estimate"]][[1]]
  are_NT      <- input_df[, "Is_NT_ctrl"]
  are_posctrl <- input_df[, "Is_pos_ctrl"]
  are_valid <- are_NT | are_posctrl | are_gene

  if (same_scale) {
    axis_limits <- range(input_df[are_valid, c(rep1_column, rep2_column)])
  } else {
    axis_limits <- NULL
  }

  ## Set up the multi-plot layout
  layout_mat <- rbind(c(1, 2, 2, 2, 2, 2, 3),
                      5:11,
                      rep(4, 7)
                      )
  use_heights <- c(0.35, 0.45, 0.2)
  use_widths <- c(0.11, 0.21, 0.1, 0.21, 0.1, 0.21, 0.06)
  layout(layout_mat,
         heights = use_heights,
         widths = use_widths
         )
  old_par <- par(mar = rep(0, 4), cex = par("cex") / 0.66)

  for (i in 1:2) {
    MakeEmptyPlot()
  }
  text(x = 0.5, y = 0.7, labels = show_title, cex = par("cex") * 1.1)
  for (i in 1:3) {
    MakeEmptyPlot()
  }

  ## Draw the 3 scatter plots
  pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
  NT_ctrl_color <- brewer.pal(5, "Blues")[[4]]

  ScatterPlot(input_df[are_gene, rep1_column],
              input_df[are_gene, rep2_column],
              top_label = "Transcription factors",
              use_limits = axis_limits,
              ...
              )
  MakeEmptyPlot()
  ScatterPlot(input_df[are_NT, rep1_column],
              input_df[are_NT, rep2_column],
              top_label = "NT controls",
              use_color = NT_ctrl_color,
              use_limits = axis_limits,
              label_y_axis = FALSE,
              ...
              )
  MakeEmptyPlot()
  ScatterPlot(input_df[are_posctrl, rep1_column],
              input_df[are_posctrl, rep2_column],
              top_label = "Positive controls",
              use_color = pos_ctrl_color,
              use_limits = axis_limits,
              label_y_axis = FALSE,
              ...
              )
  MakeEmptyPlot()

  par(old_par)
  layout(1)

  return(invisible(NULL))
}




ExportAllReplicateScatterPlots <- function(input_df) {

  use_dir <- file.path(output_dir, "Figures", "Replicate scatter plots")

  rep_columns <- grep("_rep", names(column_file_names), value = TRUE, fixed = TRUE)

  plot_height <- 4.5
  plot_ratio <- 0.45 / 0.21

  message("Exporting PDF plots...")

  pdf(file = file.path(use_dir, "Replicate scatter plots - flexible axes.pdf"),
      width = plot_height * plot_ratio, height = plot_height
      )
  for (use_column in rep_columns) {
    ReplicateScatter(input_df, rep1_column = use_column,
                     show_title = FormatPlotMath(long_column_labels[[use_column]]),
                     same_scale = FALSE
                     )
  }
  dev.off()


  pdf(file = file.path(use_dir, "Replicate scatter plots - fixed axes.pdf"),
      width = plot_height * plot_ratio, height = plot_height
      )
  for (use_column in rep_columns) {
    ReplicateScatter(input_df, rep1_column = use_column,
                     show_title = FormatPlotMath(long_column_labels[[use_column]]),
                     same_scale = TRUE
                     )
  }
  dev.off()

  message("Exporting PNG plots...")

  for (fixed_axes in c(FALSE, TRUE)) {
    for (i in seq_along(rep_columns)) {
      use_column <- rep_columns[[i]]
      file_name <- paste0("Replicate scatter plot - ", i,  ") ",
                          column_file_names[[use_column]], " - ",
                          if (fixed_axes) "fixed axes" else "flexible axes",
                          ".png"
                          )
      png(filename = file.path(use_dir, "Replicate scatter plots - PNGs", file_name),
          width = plot_height * plot_ratio, height = plot_height,
          units = "in", res = 600
          )
      ReplicateScatter(input_df, rep1_column = use_column,
                       show_title = FormatPlotMath(long_column_labels[[use_column]]),
                       same_scale = fixed_axes
                       )
      dev.off()
    }
  }
  return(invisible(NULL))
}


