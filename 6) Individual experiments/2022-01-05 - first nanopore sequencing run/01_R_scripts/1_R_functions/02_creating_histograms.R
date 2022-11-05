## 2022-04-08


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define functions --------------------------------------------------------

HistogramPolygon <- function(input_df, half_width) {
  x_vec <- rbind(input_df[, "Mids"] - half_width, input_df[, "Mids"] + half_width)
  attributes(x_vec) <- NULL
  y_vec <- rep(input_df[, "Counts"], each = 2)
  results_mat <- cbind(
    x = c(x_vec[[1]], x_vec, x_vec[[length(x_vec)]]),
    y = c(0, y_vec, 0)
  )
  return(results_mat)
}


MakeHistogramPolygons <- function(hist_results) {
  half_width <- diff(hist_results[["mids"]][1:2]) / 2
  are_empty <- hist_results[["counts"]] == 0
  group_lengths <- rle(are_empty)[["lengths"]]
  rle_vec <- rep(seq_along(group_lengths), group_lengths)
  hist_df <- data.frame(
    "Mids"     = hist_results[["mids"]],
    "Counts"   = hist_results[["counts"]],
    "Is_empty" = are_empty,
    "Group"    = rle_vec,
    stringsAsFactors = FALSE
  )
  hist_df <- hist_df[!(are_empty), ]
  row.names(hist_df) <- NULL
  hist_df[, "Group"] <- match(hist_df[, "Group"], unique(hist_df[, "Group"]))
  split_df_list <- split(hist_df, hist_df[, "Group"])
  mat_list <- lapply(split_df_list, HistogramPolygon, half_width = half_width)
  return(mat_list)
}



DrawHistogram <- function(numeric_input,
                          truncation_limit    = NULL,
                          num_breaks          = 200L,
                          title_text          = "",
                          x_axis_label        = "",
                          y_axis_label        = "Count",
                          x_axis_upper_limit  = NULL,
                          y_axis_upper_limit  = NULL,
                          fixed_y_upper_limit = FALSE,
                          y_axis_limits       = NULL,
                          x_axis_limits       = NULL,
                          x_axis_space        = 0.03,
                          y_axis_space        = 0.03,
                          hist_colors         = c(brewer.pal(9, "Blues")[[7]],
                                                  brewer.pal(9, "Purples")[[7]]
                                                  ),
                          title_font          = 2,
                          draw_box            = TRUE,
                          show_outline        = NULL,
                          use_lwd             = 1.5,
                          use_tcl             = 0.35,
                          x_axis_label_line   = 2.2,
                          y_axis_label_line   = 2.9,
                          x_axis_mgp          = 0.5,
                          y_axis_mgp          = 0.5
                          ) {

  if (!(is.list(numeric_input))) {
    numeric_list <- list(numeric_input)
    if (is.null(show_outline)) {
      show_outline <- FALSE
    }
  } else {
    numeric_list <- numeric_input
    if (is.null(show_outline)) {
      show_outline <- length(numeric_list) >= 2
    }
  }

  if (!(is.null(truncation_limit))) {
    numeric_list <- lapply(numeric_list, function(x) {
      ifelse(x > truncation_limit,
             truncation_limit,
             x
             )})
  }

  ## Calculate histogram
  hist_list <- lapply(numeric_list, function(x) {
    hist(x, breaks = num_breaks, plot = FALSE)
  })
  assign("delete_hist_list", hist_list, envir = globalenv())

  ## Prepare axes
  counts_max <- max(vapply(hist_list, function(x) max(x[["counts"]]), integer(1)))
  half_width <- diff(hist_list[[1]][["mids"]][c(1, 2)]) / 2
  numeric_min <- min(vapply(hist_list, function(x) min(x[["mids"]]), numeric(1))) - half_width
  numeric_max <- max(vapply(hist_list, function(x) max(x[["mids"]]), numeric(1))) + half_width
  if (is.null(y_axis_upper_limit)) {
    y_axis_upper_limit <- max(counts_max, y_axis_limits[[2]])
  }
  if (is.null(y_axis_limits)) {
    y_axis_lower_limit <- y_axis_upper_limit * (-y_axis_space)
    if (fixed_y_upper_limit) {
      y_axis_limits <- c(y_axis_lower_limit, y_axis_upper_limit)
    } else {
      y_axis_limits <- c(y_axis_upper_limit * (-y_axis_space), y_axis_upper_limit * (1 + y_axis_space))
    }
  }
  y_axis_ticks <- pretty(c(0, y_axis_upper_limit), n = 6)
  if (all(y_axis_ticks[-1] >= 10^5)) {
    y_axis_labels <- paste0(y_axis_ticks / 1000, "k")
  } else {
    y_axis_labels <- format(y_axis_ticks)
  }
  if (is.null(x_axis_upper_limit)) {
    x_axis_upper_limit <- numeric_max
  }
  if (is.null(x_axis_limits)) {
    x_axis_gap <- (x_axis_upper_limit - numeric_min) * x_axis_space
    x_axis_lower_limit <- min(0, numeric_min - x_axis_gap)
    x_axis_limits <- c(x_axis_lower_limit, x_axis_upper_limit)
  }

  ## Draw histogram
  plot(NA, ann = FALSE, axes = FALSE,
       xlim = x_axis_limits, ylim = y_axis_limits,
       xaxs = "i", yaxs = "i"
       )

  for (i in seq_along(hist_list)) {
    polygon_mat_list <- MakeHistogramPolygons(hist_list[[i]])
    if (show_outline) {
      for (polygon_mat in polygon_mat_list) {
        lines(polygon_mat[, "x"],
              polygon_mat[, "y"],
              col  = hist_colors[[i]],
              lwd  = par("lwd") * use_lwd,
              lend = "butt",
              xpd  = NA
              )
      }
    } else {
      for (polygon_mat in polygon_mat_list) {
        polygon(polygon_mat[, "x"],
                polygon_mat[, "y"],
                col    = hist_colors[[i]],
                border = NA,
                xpd    = NA
                )
      }
    }
  }

  ## Draw x axis
  x_axis_ticks <- axTicks(1)
  x_axis_labels <- format(x_axis_ticks)
  if (is.null(truncation_limit)) {
    x_axis_labels <- format(x_axis_ticks)
  } else {
    x_axis_labels <- ifelse(x_axis_ticks >= truncation_limit,
                            as.expression(bquote("" >= .(as.character(truncation_limit)))),
                            format(x_axis_ticks)
                            )
  }

  axis(1, mgp = c(3, x_axis_mgp, 0), tcl = -(use_tcl),
       at = x_axis_ticks, labels = x_axis_labels, lwd = par("lwd")
       )
  mtext(x_axis_label, side = 1, line = x_axis_label_line, cex = par("cex"))

  ## Draw y axis
  axis(2, las = 1, mgp = c(3, y_axis_mgp, 0), tcl = -(use_tcl),
       at = y_axis_ticks, labels = y_axis_labels, lwd = par("lwd")
       )
  mtext(y_axis_label, side = 2, line = y_axis_label_line, cex = par("cex"))

  ## Final steps
  title(title_text, cex.main = 1, font.main = title_font)
  if (draw_box) {
    box(bty = "l")
  }

  return(invisible(NULL))
}







ReadLengthsHistogram <- function(read_lengths,
                                 read_length_limit = 4000L,
                                 title_text = "Nanopore sequencing of the CRISPRa library",
                                 x_axis_upper_limit = NULL
                                 ) {
  DrawHistogram(read_lengths,
                num_breaks         = 500L,
                truncation_limit   = read_length_limit,
                title_text         = title_text,
                x_axis_label       = "Read length (base pairs)",
                y_axis_label       = "Read count",
                x_axis_upper_limit = x_axis_upper_limit
                )
  return(invisible(NULL))
}



ExportCountHistograms <- function(use_counts_df, title_postfix = "") {

  use_width <- 5
  use_height <- 4
  use_res <- 600

  for (use_device in c("none", "pdf", "png")) {

    if (use_device == "pdf") {
      pdf(file.path(figures_dir, "PDFs", "Histogram - reads per plasmid.pdf"),
          width = use_width, height = use_height
          )
    }


    if (use_device == "png") {
      png(file.path(figures_dir, "PNGs", "Histogram - reads per plasmid - unfiltered.png"),
          width = use_width, height = use_height, units = "in", res = use_res
          )
    }
    DrawHistogram(use_counts_df[, "Count_unfiltered"],
                  truncation_limit = 800,
                  x_axis_upper_limit = 800,
                  num_breaks = 200,
                  x_axis_label = "Number of reads per plasmid",
                  y_axis_label = "Plasmid count",
                  title_text = paste0("Unfiltered reads", title_postfix)
                  )
    if (use_device == "png") {
      dev.off()
    }


    if (use_device == "png") {
      png(file.path(figures_dir, "PNGs", "Histogram - reads per plasmid - 2) sg2 and sg3.png"),
          width = use_width, height = use_height, units = "in", res = use_res
          )
    }
    DrawHistogram(use_counts_df[, "Count_sg2_match_sg3"],
                  truncation_limit = 800,
                  x_axis_upper_limit = 800,
                  num_breaks = 200,
                  x_axis_label = "Number of reads per plasmid",
                  y_axis_label = "Plasmid count",
                  title_text = paste0("Reads using sg2 & sg3", title_postfix)
                  )
    if (use_device == "png") {
      dev.off()
    }

    if (use_device == "png") {
      png(file.path(figures_dir, "PNGs", "Histogram - reads per plasmid - 3) sg3 & sg4.png"),
          width = use_width, height = use_height, units = "in", res = use_res
          )
    }
    DrawHistogram(use_counts_df[, "Count_sg3_match_sg4"],
                  truncation_limit = 800,
                  x_axis_upper_limit = 800,
                  num_breaks = 200,
                  x_axis_label = "Number of reads per plasmid",
                  y_axis_label = "Plasmid count",
                  title_text = paste0("Reads using sg3 & sg4", title_postfix)
                  )
    if (use_device == "png") {
      dev.off()
    }

    if (use_device == "pdf") {
      dev.off()
    }
  }
  return(invisible(NULL))
}



RawCountsHistogram <- function(raw_mat,
                               show_replicates = FALSE,
                               semitransparent_lines = FALSE,
                               title_text = "Plasmid count at baseline",
                               x_axis_label = expression("Log"[10] ~ "(raw count + 1)"),
                               x_axis_space = 0.01,
                               y_axis_space = 0.01,
                               x_start_lines = 0.75,
                               y_start_lines = 1,
                               text_x_lines = 1.1,
                               segment_length_lines = 0.45,
                               short_labels = FALSE,
                               ...
                               ) {
  LogTransform <- function(x) log10(x + 1)
  if (show_replicates) {
    input_list <- rev(list(LogTransform(raw_mat[, 1]), LogTransform(raw_mat[, 2])))
    line_colors <- c(brewer.pal(9, "Blues")[[8]], brewer.pal(9, "Blues")[[5]])
    hist_colors <- line_colors
    if (semitransparent_lines) {
      hist_colors <- adjustcolor(hist_colors, alpha.f = 0.8)
    }
  } else {
    input_list <- LogTransform(raw_mat[, 1:2])
    hist_colors <- brewer.pal(9, "Blues")[[8]]
  }
  DrawHistogram(input_list,
                num_breaks          = 100,
                x_axis_label        = x_axis_label,
                y_axis_label        = "Frequency",
                title_text          = title_text,
                title_font          = 1,
                hist_colors         = rev(hist_colors),
                draw_box            = FALSE,
                x_axis_space        = x_axis_space,
                y_axis_space        = y_axis_space,
                ...
                )
  if (show_replicates) {
    ## Draw the legend
    x_start <- par("usr")[[1]] + diff(grconvertX(c(0, x_start_lines), from = "lines", to = "user"))
    y_start <- par("usr")[[4]] - diff(grconvertY(c(0, y_start_lines), from = "lines", to = "user"))
    timepoints_seq <- seq_along(hist_colors) - 1L
    y_vec <- y_start - diff(grconvertY(c(0, 1.2), from = "lines", to = "user")) * timepoints_seq
    segments(x0  = x_start,
             x1  = x_start + diff(grconvertX(c(0, segment_length_lines), from = "lines", to = "user")),
             y0  = y_vec,
             col = line_colors,
             lwd = par("lwd") * 2,
             xpd = NA
             )
    # points(x   = rep(x_start + diff(grconvertX(c(0, 0.1))), length(y_vec)),
    #        y   = y_vec,
    #        pch = 22,
    #        col = hist_colors,
    #        lwd = par("lwd") * 2
    #        )
    text(x      = x_start + diff(grconvertX(c(0, text_x_lines), from = "lines", to = "user")),
         y      = y_vec,
         labels = if (short_labels) c("R1", "R2") else c("Replicate 1", "Replicate 2"),
         adj    = c(0, 0.5),
         xpd    = NA
         )
  } else {
    ## Show the mean
    mean_count <- mean(raw_mat[, 1:2])
    mean_transf <- LogTransform(mean_count)
    abline(v = mean_transf, col = brewer.pal(9, "Reds")[[6]], lwd = 1.5)
    text(x      = mean_transf + diff(grconvertX(c(0, 0.75), from = "lines", to = "user")),
         y      = par("usr")[[4]] - diff(grconvertY(c(0, 0.5), from = "lines", to = "user")),
         labels = paste0("Mean: ", round(mean_count / 10) * 10),
         adj    = c(0, 0.5),
         xpd    = NA
         )
  }
  box(bty = "l")

  return(invisible(NULL))
}



ManuscriptRawCountsHistogram <- function(raw_counts_mat, title_text) {
  old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
  RawCountsHistogram(raw_counts_mat,
                     y_axis_upper_limit    = 1600,
                     fixed_y_upper_limit   = TRUE,
                     show_replicates       = TRUE,
                     semitransparent_lines = TRUE,
                     x_axis_space          = 0.02,
                     y_axis_space          = 0.02,
                     title_text            = title_text,
                     x_axis_label          = expression("Log"[10] ~ "count at baseline"),
                     use_lwd               = 1.5,
                     segment_length_lines  = 0.4,
                     x_start_lines         = 1.25,
                     y_start_lines         = 2,
                     text_x_lines          = 0.75,
                     short_labels          = TRUE,
                     use_tcl               = 0.3,
                     x_axis_label_line     = 1.7,
                     y_axis_label_line     = 2.68,
                     x_axis_mgp            = 0.35,
                     y_axis_mgp            = 0.45
                     )
  par(old_par)
  return(invisible(NULL))
}



