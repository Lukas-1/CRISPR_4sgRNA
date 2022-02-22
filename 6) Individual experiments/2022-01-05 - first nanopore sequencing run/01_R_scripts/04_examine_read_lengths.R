## 2022-02-03


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define paths ------------------------------------------------------------

project_dir <- file.path("~", "NP1_4sg")
rdata_dir   <- file.path(project_dir, "03_R_objects")
figures_dir <- file.path(project_dir, "04_output_data", "Figures")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_combine_aligned_reads.RData"))



# Define functions --------------------------------------------------------

ReadLengthsHistogram <- function(read_lengths, read_length_limit = 4000L) {

  ## Truncate very long read lengths
  read_lengths <- ifelse(read_lengths > read_length_limit,
                         read_length_limit,
                         read_lengths
                         )

  ## Draw the histogram
  hist_results <- hist(read_lengths, breaks = 500, plot = FALSE)
  y_max <- max(hist_results[["counts"]])
  y_axis_limits <- c(y_max * (-0.03), y_max * 1.03)
  plot(1, type = "n", ann = FALSE, axes = FALSE,
       xlim = c(0, max(read_lengths)), ylim = y_axis_limits,
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

  ## Draw the x axis
  x_axis_ticks <- axTicks(1)
  x_axis_labels <- ifelse(x_axis_ticks >= read_length_limit,
                          as.expression(bquote("" >= .(as.character(read_length_limit)))),
                          format(x_axis_ticks)
                          )
  axis(1, mgp = c(2.6, 0.5, 0), tcl = -0.35,
       at = x_axis_ticks, labels = x_axis_labels
       )
  mtext("Read length (base pairs)", side = 1, line = 2.2)

  ## Draw the y axis
  y_axis_ticks <- axTicks(2)
  if (all(y_axis_ticks[-1] >= 10^5)) {
    y_axis_labels <- paste0(y_axis_ticks / 1000, "k")
  } else {
    y_axis_labels <- format(y_axis_ticks)
  }
  axis(2, las = 1, mgp = c(2.6, 0.5, 0), tcl = -0.35,
       at = y_axis_ticks, labels = y_axis_labels
       )
  mtext("Read count", side = 2, line = 2.9)

  ## Final steps
  title("Nanopore sequencing of the CRISPRa library", cex.main = 1)
  box(bty = "l")

  return(invisible(NULL))
}



# Prepare data ------------------------------------------------------------

alignments_df[, "Read_length"] <- nchar(alignments_df[, "Read_sequence"])



# Draw a histogram of read lengths ----------------------------------------

ReadLengthsHistogram(alignments_df[, "Read_length"])


pdf(file = file.path(figures_dir, "PDFs", "Read length histogram.pdf"),
    width = 5, height = 4
    )
ReadLengthsHistogram(alignments_df[, "Read_length"])
dev.off()


png(filename = file.path(figures_dir, "PNGs", "Read length histogram.png"),
    width = 5, height = 4, units = "in", res = 600
    )
ReadLengthsHistogram(alignments_df[, "Read_length"])
dev.off()





