#2022-01-04


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define labels -----------------------------------------------------------

controls_labels <- list(
  "Gene" = c("Genes in ", "CRISPRa", "library"),
  "NT"   = c("Non-targeting", "controls"),
  "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)"))
)



# Define functions --------------------------------------------------------

PlotHistResults <- function(hist_results, fill_color, border_color) {
  half_width <- (hist_results[["mids"]][[2]] - hist_results[["mids"]][[1]]) / 2
  are_not_zero <- hist_results[["counts"]] != 0
  mids_vec <- hist_results[["mids"]][are_not_zero]
  rect(xleft   = mids_vec - half_width,
       xright  = mids_vec + half_width,
       ybottom = 0,
       ytop    = hist_results[["counts"]][are_not_zero],
       col     = fill_color,
       border  = border_color,
       lwd     = 0.5,
       xpd     = NA
       )
}


ThreeHistograms <- function(input_df,
                            use_column,
                            x_axis_label      = NULL,
                            gene_border_color = "gray55",
                            gene_fill_color   = "black",
                            gene_fill_alpha   = 0.35,
                            use_mai           = c(1.0, 0.92, 0.76, 1.5),
                            legend_x_start    = 0.75,
                            legend_y_mid      = 0.5,
                            use_mgp           = c(2.7, 0.55, 0),
                            x_axis_mgp        = use_mgp,
                            ...
                            ) {

  if (is.null(x_axis_label)) {
    x_axis_label <- FormatPlotMath(long_column_labels[[use_column]])
  }

  ## Prepare data
  are_gene  <- !(is.na(input_df[, "Entrez_ID"]))
  are_NT    <- input_df[, "Is_NT_ctrl"]
  are_pos   <- input_df[, "Is_pos_ctrl"]
  are_valid <- are_NT | are_pos | are_gene

  has_replicates <- grepl("_rep", use_column, fixed = TRUE)
  if (has_replicates) {
    rep2_column <- sub("_rep1", "_rep2", use_column, fixed = TRUE)
    numeric_vec <- rowMeans(input_df[, c(use_column, rep2_column)])
  } else {
    numeric_vec <- input_df[, use_column]
  }

  ## Determine x axis limits
  side_space_fraction <- 0.03
  data_range <- range(numeric_vec[are_valid])
  x_limits   <- DataAxisLimits(numeric_vec[are_valid], space_fraction = side_space_fraction)

  ## Compute the 3 histograms
  use_breaks <- seq(from = data_range[[1]], to = data_range[[2]], length.out = 70)
  gene_hist <- hist(numeric_vec[are_gene], breaks = use_breaks, plot = FALSE)
  NT_hist   <- hist(numeric_vec[are_NT],   breaks = use_breaks, plot = FALSE)
  pos_hist  <- hist(numeric_vec[are_pos],  breaks = use_breaks, plot = FALSE)

  ## Determine y axis limits
  count_max <- max(c(gene_hist[["counts"]], NT_hist[["counts"]], pos_hist[["counts"]]))
  y_limits <- c(count_max * -(side_space_fraction), count_max * (1 + side_space_fraction))

  ## Prepare graphical parameters
  fill_colors <- c(
    gene_fill_color,                                     # genes
    colorRampPalette(brewer.pal(9, "Blues"))(100)[[80]], # NT control
    brewer.pal(5, "Reds")[[4]]                           # positive control
  )
  for (i in 1:3) {
    if (i == 1) {
      use_alpha <- gene_fill_alpha
    } else {
      use_alpha <- 0.75
    }
    fill_colors[[i]] <- adjustcolor(fill_colors[[i]], alpha.f = use_alpha)
  }
  border_colors <- c(
    gene_border_color,           # genes
    brewer.pal(5, "Blues")[[5]], # NT control
    brewer.pal(5, "Reds")[[5]]   # positive control
  )

  use_tcl <- -0.35

  ## Prepare the plot region
  old_mai <- par(mai = use_mai)
  plot(1,
       xlim = x_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       xlab = "",
       ylab = "Count",
       mgp  = use_mgp,
       type = "n"
       )
  axis(1, mgp = x_axis_mgp, tcl = use_tcl, lwd = par("lwd"))
  mtext(x_axis_label, side = 1, line = x_axis_mgp[[1]], cex = par("cex"))
  axis(2, las = 2, mgp = use_mgp, tcl = use_tcl, lwd = par("lwd"))
  if (use_column == "CellTiterGlo_foldNT") {
    abline(v = 1, lty = "dashed", col = "gray40")
  }

  ## Add the 3 superimposed histograms
  if (grepl("SSMD|CellTiterGlo", use_column)) {
    PlotHistResults(gene_hist, fill_colors[[1]], border_colors[[1]])
    PlotHistResults(NT_hist,   fill_colors[[2]], border_colors[[2]])
    PlotHistResults(pos_hist,  fill_colors[[3]], border_colors[[3]])
  } else {
    # Prevent the positive controls from obscuring genes,
    ## by placing them in the background.
    PlotHistResults(pos_hist,  fill_colors[[3]], border_colors[[3]])
    PlotHistResults(gene_hist, fill_colors[[1]], border_colors[[1]])
    PlotHistResults(NT_hist,   fill_colors[[2]], border_colors[[2]])

  }
  box(bty = "l")

  # Add legend
  DrawSideLegend(labels_list    = controls_labels,
                 use_pch        = 22,
                 use_colors     = fill_colors,
                 border_colors  = border_colors,
                 use_point_size = 1.4,
                 lines_x_start  = legend_x_start,
                 y_mid          = legend_y_mid,
                 ...
                 )

  par(old_mai)
  return(invisible(NULL))
}




