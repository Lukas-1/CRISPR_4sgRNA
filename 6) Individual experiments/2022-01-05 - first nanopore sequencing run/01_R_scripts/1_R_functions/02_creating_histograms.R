## 2022-04-08


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define functions --------------------------------------------------------

DrawHistogram <- function(numeric_vec,
                          truncation_limit   = NULL,
                          num_breaks         = 200L,
                          title_text         = "",
                          x_axis_label       = "",
                          y_axis_label       = "Count",
                          x_axis_upper_limit = NULL,
                          y_axis_upper_limit = NULL,
                          y_axis_limits      = NULL
                          ) {

  if (!(is.null(truncation_limit))) {
    numeric_vec <- ifelse(numeric_vec > truncation_limit,
                          truncation_limit,
                          numeric_vec
                          )
  }

  ## Calculate histogram
  hist_results <- hist(numeric_vec, breaks = num_breaks, plot = FALSE)

  ## Prepare axes
  if (is.null(y_axis_upper_limit)) {
    y_axis_upper_limit <- max(c(hist_results[["counts"]], y_axis_limits[[2]]))
  }
  if (is.null(y_axis_limits)) {
    y_axis_limits <- c(y_axis_upper_limit * (-0.03), y_axis_upper_limit * 1.03)
  }
  y_axis_ticks <- pretty(c(0, y_axis_upper_limit), n = 6)
  if (all(y_axis_ticks[-1] >= 10^5)) {
    y_axis_labels <- paste0(y_axis_ticks / 1000, "k")
  } else {
    y_axis_labels <- format(y_axis_ticks)
  }
  if (is.null(x_axis_upper_limit)) {
    x_axis_upper_limit <- max(numeric_vec)
  }

  ## Draw histogram
  plot(1, type = "n", ann = FALSE, axes = FALSE,
       xlim = c(0, x_axis_upper_limit), ylim = y_axis_limits,
       xaxs = "i", yaxs = "i"
       )
  half_width <- (hist_results[["mids"]][[2]] - hist_results[["mids"]][[1]]) / 2
  are_not_zero <- hist_results[["counts"]] != 0
  mids_vec <- hist_results[["mids"]][are_not_zero]
  rect(xleft   = mids_vec - half_width,
       xright  = mids_vec + half_width,
       ybottom = 0,
       ytop    = hist_results[["counts"]][are_not_zero],
       col     = brewer.pal(9, "Blues")[[7]],
       border  = NA,
       lwd     = 0.5,
       xpd     = NA
       )

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

  axis(1, mgp = c(2.6, 0.5, 0), tcl = -0.35,
       at = x_axis_ticks, labels = x_axis_labels
       )
  mtext(x_axis_label, side = 1, line = 2.2)

  ## Draw y axis
  axis(2, las = 1, mgp = c(2.6, 0.5, 0), tcl = -0.35,
       at = y_axis_ticks, labels = y_axis_labels
       )
  mtext(y_axis_label, side = 2, line = 2.9)

  ## Final steps
  title(title_text, cex.main = 1)
  box(bty = "l")

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


