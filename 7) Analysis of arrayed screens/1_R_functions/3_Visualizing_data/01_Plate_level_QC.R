# 2022-01-17


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define functions --------------------------------------------------------

PlotPlateQualities <- function(rep1_vec,
                               rep2_vec,
                               y_limits_include  = NULL,
                               y_axis_label      = "",
                               quality_ranges    = list(c(0.5, 1),
                                                        c(0, 0.5),
                                                        c(-Inf, 0)
                                                        ),
                               use_mai           = c(0.7, 0.82, 0.5, 0.42),
                               y_label_line      = 2.2,
                               reorder_plates    = FALSE,
                               plates_in_order   = NULL,
                               roman_plates      = TRUE,
                               label_plates      = TRUE,
                               plate_labels_line = 0.3,
                               x_label_line      = 1.85,
                               y_axis_ticks      = NULL,
                               y_axis_labels     = NULL,
                               point_cex         = 0.7
                               ) {

  stopifnot(length(rep1_vec) == length(rep2_vec))

  all_plates <- seq_along(rep1_vec)

  if (reorder_plates) {
    if (is.null(plates_in_order)) {
      average_qualities <- rowMeans(cbind(rep1_vec, rep2_vec))
      plates_order <- order(average_qualities)
    } else {
      plates_order <- order(match(all_plates, plates_in_order))
    }
    rep1_vec <- rep1_vec[plates_order]
    rep2_vec <- rep2_vec[plates_order]
    all_plates <- all_plates[plates_order]
  }

  if (roman_plates) {
    plate_names <- as.character(as.roman(all_plates))
  } else {
    plate_names <- as.character(all_plates)
  }
  data_vec <- c(rep1_vec, rep2_vec)

  ## Prepare x axis positions
  x_mids <- seq_along(rep1_vec)
  x_space <- 0.5
  x_positions <- c(x_mids - (x_space / 2), x_mids + (x_space / 2))
  x_space <- 0.5 + length(x_mids) * 0.015
  x_limits <- c(1 - x_space, length(x_mids) + x_space)
  print(x_positions)

  ## Prepare y axis positions
  y_limits <- range(c(y_limits_include, data_vec))
  y_space <- diff(y_limits) * 0.05
  if (y_limits[[1]] > (min(data_vec) - y_space)) {
    y_limits[[1]] <- y_limits[[1]] - y_space
  }
  if (y_limits[[2]] < (max(data_vec) + y_space)) {
    y_limits[[2]] <- y_limits[[2]] + y_space
  }

  ## Set up the plot region
  old_mai <- par("mai" = use_mai)
  plot(1,
       xlim = x_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  if (is.null(y_axis_ticks)) {
    y_axis_ticks <- axTicks(2)
  }
  if (is.null(y_axis_labels)) {
    y_axis_labels = as.character(y_axis_ticks)
  }
  axis(2, at = y_axis_ticks, labels = y_axis_labels,
       las = 1, mgp = c(3, 0.55, 0), tcl = -0.35, lwd = par("lwd")
       )
  mtext(y_axis_label, side = 2, line = y_label_line, cex = par("cex"))

  if (label_plates) {
    mtext(plate_names,
          at   = x_mids,
          side = 1,
          line = plate_labels_line,
          cex  = 0.9 * par("cex")
          )
    mtext(FormatPlotMath("Plate number"), side = 1, line = x_label_line, cex = par("cex"))
  }

  ## Indicate specific y axis ranges with colors
  color_scheme <- c(colorRampPalette(brewer.pal(9, "Greens"))(100)[[26]],
                    brewer.pal(9, "YlOrRd")[[2]],
                    colorRampPalette(brewer.pal(9, "Reds"))(100)[[21]]
                    )
  rect(xleft   = x_limits[[1]],
       xright  = x_limits[[2]],
       ybottom = quality_ranges[[1]][[1]],
       ytop    = y_limits[[2]],
       col     = color_scheme[[1]],
       border  = NA
       )
  rect(xleft   = x_limits[[1]],
       xright  = x_limits[[2]],
       ybottom = quality_ranges[[2]][[1]],
       ytop    = quality_ranges[[2]][[2]],
       col     = color_scheme[[2]],
       border  = NA
       )
  rect(xleft   = x_limits[[1]],
       xright  = x_limits[[2]],
       ybottom = y_limits[[1]],
       ytop    = quality_ranges[[3]][[2]],
       col     = color_scheme[[3]],
       border  = NA
       )

  ## Draw horizontal and vertical indicator lines
  if (y_limits[[2]] != quality_ranges[[1]][[2]]) {
    abline(h = quality_ranges[[1]][[2]], col = "gray50", lty = "dotted")
  }
  abline(v = seq_len(length(x_mids) + 1) - 0.5,
         col = "gray70", lwd = 0.5
         )
  box(lwd = par("lwd"))

  ## Plot the quality control metrics
  points(x   = x_positions,
         y   = data_vec,
         pch = 21,
         bg  = "black",
         cex = point_cex
         )

  par(old_mai)

  return(invisible(all_plates))
}


GetQualityMetric <- function(input_df, UseFunction, filter_NT = FALSE, ...) {
  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  df_list <- split(input_df, plate_numbers_vec)
  rep1_vec <- vapply(df_list, UseFunction, use_column = "Raw_rep1", filter_NT = filter_NT, numeric(1))
  rep2_vec <- vapply(df_list, UseFunction, use_column = "Raw_rep2", filter_NT = filter_NT, numeric(1))
  results_mat <- cbind("rep1" = rep1_vec, "rep2" = rep2_vec)
  return(results_mat)
}


PlotZPrimes <- function(input_df, filter_NT = FALSE, ...) {
  z_primes_mat <- GetQualityMetric(input_df, Calculate_Z_Prime, filter_NT = filter_NT)
  PlotPlateQualities(z_primes_mat[, 1], z_primes_mat[, 2],
                     y_limits_include = c(0, 1, -0.2), y_axis_label = "Z' factor",
                     ...
                     )
}


PlotSSMDControls <- function(input_df,
                             filter_NT = FALSE,
                             y_limits_include = c(0, 7),
                             y_axis_label = "SSMD (pos./neg. controls)",
                             ...
                             ) {
  z_primes_mat <- GetQualityMetric(input_df, Calculate_SSMD_ctrls, filter_NT = filter_NT)
  PlotPlateQualities(z_primes_mat[, 1], z_primes_mat[, 2],
                     y_limits_include = y_limits_include,
                     y_axis_label = y_axis_label,
                     quality_ranges = list(c(5, 7), c(3, 5), c(-Inf, 3)),
                     ...
                     )
}



