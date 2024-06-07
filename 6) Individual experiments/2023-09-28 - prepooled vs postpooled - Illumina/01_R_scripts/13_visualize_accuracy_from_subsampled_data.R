### 2024-05-05


# Load packages and source code -------------------------------------------

library("beeswarm")
library("RColorBrewer")



# Define paths ------------------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")
rdata_dir             <- file.path(project_dir, "03_R_objects")
figures_dir           <- file.path(project_dir, "04_output", "Figures")
PDFs_dir              <- file.path(figures_dir, "Subsampling")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "12_compute_accuracy_from_subsampled_data.RData"))



# Define functions --------------------------------------------------------

GetMinorTicks <- function(ticks_vec) {
  setdiff(seq(min(ticks_vec), max(ticks_vec), (ticks_vec[[2]] - ticks_vec[[1]]) / 2),
          ticks_vec
          )
}


PlotColumn <- function(show_column) {
  stopifnot("subsampling_df" %in% ls(envir = globalenv()))

  if (grepl("_AUC_", show_column, fixed = TRUE)) {
    y_limits <- c(0.5, 1)
    y_label <- "AUC"
  } else if (grepl("_SSMD_", show_column, fixed = TRUE)) {
    y_limits <- c(0, 2.5)
    y_label <- "SSMD*"
  } else {
    stop("Unexpected value for 'show_column'!")
  }

  plot.new()
  plot.window(xlim = c(0, 100), ylim = y_limits, xaxs = "i", yaxs = "i")

  x_ticks <- axTicks(1)
  y_ticks <- axTicks(2)

  grid_color <- "gray78"
  light_color <- "gray92"

  segments(x0  = par("usr")[[1]],
           x1  = par("usr")[[2]],
           y0  = y_ticks,
           col = grid_color,
           xpd = NA
           )
  segments(x0  = par("usr")[[1]],
           x1  = par("usr")[[2]],
           y0  = GetMinorTicks(y_ticks),
           col = light_color,
           xpd = NA
           )
  segments(x0  = x_ticks,
           y0  = par("usr")[[3]],
           y1  = par("usr")[[4]],
           col = grid_color,
           xpd = NA
           )
  segments(x0  = GetMinorTicks(x_ticks),
           y0  = par("usr")[[3]],
           y1  = par("usr")[[4]],
           col = light_color,
           xpd = NA
           )

  axis(1, mgp = c(3, 0.55, 0), tcl = -0.4, at = x_ticks, labels = paste0(x_ticks, "%"), col = grid_color)
  axis(2, mgp = c(3, 0.8, 0), tcl = -0.6, las = 1, col = grid_color)

  groups_vec <- subsampling_df[, "Fraction_sampled"]
  numeric_vec <- abs(subsampling_df[, show_column])

  use_cex <- 0.4

  group_means <- tapply(numeric_vec, groups_vec, mean)
  segments(x0  = as.numeric(names(group_means)) - 1,
           x1  = as.numeric(names(group_means)) + 1,
           y0  = group_means,
           xpd = NA,
           lwd = 1.5
           )

  beeswarm_df <- beeswarm(split(numeric_vec, groups_vec), do.plot = FALSE, cex = use_cex)

  x_vec <- beeswarm_df[, "x"] -
    match(beeswarm_df[, "x.orig"], unique(beeswarm_df[, "x.orig"])) +
    as.numeric(beeswarm_df[, "x.orig"])

  points(x   = x_vec,
         y   = beeswarm_df[, "y"],
         pch = 16,
         cex = use_cex,
         col = adjustcolor(brewer.pal(9, "Blues")[[9]], alpha.f = 0.5)
         )

  title(ylab = y_label, xlab = "Fraction of subsampled reads", mgp = c(2.5, 1, 0))
  return(invisible(NULL))

}


pdf(file = file.path(PDFs_dir, "Accuracy for subsampled reads.pdf"),
    width = 5.5, height = 4.5
    )
par(mar = c(4, 4, 4, 2))

PlotColumn("Postpool_AUC_both")
title("Prepool \u2013 both replicates", cex.main = 1, font.main = 1)

PlotColumn("Prepool_AUC_rep1")
title("Prepool \u2013 replicate 1", cex.main = 1, font.main = 1)

PlotColumn("Prepool_AUC_rep2")
title("Prepool \u2013 replicate 2", cex.main = 1, font.main = 1)


PlotColumn("Postpool_AUC_both")
title("Postpool \u2013 both replicates", cex.main = 1, font.main = 1)

PlotColumn("Postpool_AUC_rep1")
title("Postpool \u2013 replicate 1", cex.main = 1, font.main = 1)

PlotColumn("Postpool_AUC_rep2")
title("Postpool \u2013 replicate 2", cex.main = 1, font.main = 1)


PlotColumn("Prepool_SSMD_rep1")
title("Prepool \u2013 replicate 1", cex.main = 1, font.main = 1)
PlotColumn("Prepool_SSMD_rep2")
title("Prepool \u2013 replicate 2", cex.main = 1, font.main = 1)

PlotColumn("Postpool_SSMD_rep1")
title("Postpool \u2013 replicate 1", cex.main = 1, font.main = 1)
PlotColumn("Postpool_SSMD_rep2")
title("Postpool \u2013 replicate 2", cex.main = 1, font.main = 1)

dev.off()


