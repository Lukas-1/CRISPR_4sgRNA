## 2023-11-28


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")
illumina_pooling_dir     <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R")) # For RepositionByGroups
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(illumina_pooling_dir, "01_R_scripts", "R_functions", "02_annotating_plots.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2023-10-05 - prepooled vs postpooled - Nanopore")
rdata_dir   <- file.path(project_dir, "03_R_objects")
output_dir  <- file.path(project_dir, "05_output")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "07_produce_read_counts_and_metrics.RData"))



# Define functions --------------------------------------------------------

StackedBars <- function(bars_mat,
                        groups_vec         = NULL,
                        bar_positions      = NULL,
                        y_axis_label       = "",
                        y_label_line       = 2.3,
                        y_axis_upper_limit = NULL,
                        legend_title_vec   = NULL,
                        brewer_pal         = "Blues",
                        invert_colors      = FALSE,
                        convert_to_percent = FALSE,
                        show_legend        = NULL,
                        draw_box           = TRUE,
                        bar_width          = 0.55,
                        bar_labels_line    = -0.25,
                        set_mar            = TRUE,
                        point_x_start      = 0.0,
                        lines_x_start      = 1.4,
                        title_x_start      = 0.0,
                        side_gap           = 0.5,
                        gap_ratio          = 2,
                        grid_lwd           = 1,
                        legend_lwd         = 1
                        ) {

  ## Prepare data
  if (!(is.matrix(bars_mat))) {
    if (convert_to_percent) {
      bars_mat <- rbind(
        "1MM" = bars_mat,
        "0MM" = 1 - bars_mat
      )
      if (is.null(show_legend)) {
        show_legend <- FALSE
      }
    } else {
      stop(paste0("StackedBars is only intended for displaying multiple ",
                  "categories or a vector of percentages!"
                  ))
    }
  } else {
    if (is.null(show_legend)) {
      show_legend <- TRUE
    }
  }
  if (convert_to_percent) {
    bars_mat <- prop.table(bars_mat, margin = 2)
  }

  ## Define color scheme
  use_colors <- rev(colorRampPalette(brewer.pal(9, brewer_pal)[3:8])(nrow(bars_mat)))
  if (invert_colors) {
    use_colors <- rev(use_colors)
  }

  ## Determine bar positions
  if (is.null(bar_positions)) {
    if (is.null(groups_vec)) {
      bar_positions <- seq_len(ncol(bars_mat))
    } else {
      stopifnot(length(groups_vec) == ncol(bars_mat))
      bar_positions <- RepositionByGroups(groups_vec, gap_ratio = gap_ratio)
    }
  } else {
    stopifnot(length(bar_positions) == ncol(bars_mat))
  }

  num_bars <- length(bar_positions)
  # final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - side_gap) - (num_bars * 0.04),
                     max(bar_positions) + side_gap  + (num_bars * 0.04)
                    )


  ## Prepare the data axis
  use_numeric_limits <- c(0, max(colSums(bars_mat)) * 1.00)
  numeric_axis_pos <- pretty(use_numeric_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  if (!(is.null(y_axis_upper_limit))) {
    numeric_limits[[2]] <- y_axis_upper_limit
  }

  ## Draw the bar plot
  if (show_legend && set_mar) {
    old_mar <- par(mar = c(4, 3.75, 3.75, 6.5))
  }
  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)

  ## Draw the grid
  grid_pos <- pretty(numeric_limits, n = 30)
  if ((use_numeric_limits[[1]] == 0) && (use_numeric_limits[[2]] == 1)) {
    are_major <- vapply(grid_pos, function(x) any(abs(seq(0, 1, by = 0.1) - x) < (10^-12)), logical(1))
  } else {
    are_major <- grid_pos %in% axTicks(2)
  }
  segments(x0  = par("usr")[[1]],
           x1  = par("usr")[[2]],
           y0  = grid_pos,
           col = ifelse(are_major, "gray88", "gray95"),
           lwd = par("lwd") * grid_lwd,
           xpd = NA
           )
  PlotBarplotMat(bars_mat,
                 colors_vec    = use_colors,
                 positions_vec = bar_positions,
                 bar_width     = bar_width
                 )
  if (draw_box) {
    box(bty = "l")
  }

  ## Draw the y axis
  tick_pos <- axTicks(2)
  if (convert_to_percent) {
    tick_labels <- paste0(tick_pos * 100, "%")
  } else if (all(tick_pos[-c(1, 2)] >= 10^6)) {
    tick_labels <- paste0(format(tick_pos / 10^6), "M")
  } else {
    tick_labels <- tick_pos
  }
  axis(2,
       at     = tick_pos,
       labels = tick_labels,
       las    = 2,
       mgp    = c(3, 0.45, 0),
       tcl    = -0.35,
       lwd    = par("lwd")
       )
  mtext(VerticalAdjust(y_axis_label),
        side = 2,
        line = y_label_line,
        cex  = par("cex")
        )

  ## Draw the x axis labels
  x_labels <- sub(" ", "\n", colnames(bars_mat))
  mtext(x_labels, side = 1, padj = 1,
        at = bar_positions, line = bar_labels_line, cex = par("cex")
        )

  ## Draw the legend
  if (show_legend) {
    DrawSideLegend(title_vec     = legend_title_vec,
                   labels_list   = rownames(bars_mat),
                   use_colors    = rev(use_colors),
                   border_colors = "gray50",
                   use_pch       = 22,
                   point_x_start = point_x_start,
                   lines_x_start = lines_x_start,
                   title_x_start = title_x_start,
                   large_gap_multiplier = 1.25,
                   small_gap_size = 1.15,
                   use_point_size = 1.4,
                   border_lwd     = legend_lwd
                   )
  }

  if (show_legend && set_mar) {
    par(old_mar)
  }
  return(invisible(bar_positions))
}




# Export plots ------------------------------------------------------------

CreateEMF <- function(file_name) {
  devEMF::emf(file.path(output_dir, "Figures", "Read-level metrics", paste0(file_name, ".emf")),
              width = 2.8, height = 2.2,
              emfPlus = FALSE, coordDPI = 7000
              )

  par(mai = c(0.52, 0.5, 0.38, 0.6), cex = 0.6, lwd = 0.7)

}




y_label_line <- 1.9

num_reads_mat <- t(cbind("Num_0_sgRNAs" = num_any_mapped_mat[, "Total_num_reads"] - num_any_mapped_mat[, "Num_mapped_reads"],
                         num_any_mapped_mat[, c("Num_1_sgRNA", "Num_2_sgRNAs", "Num_3_sgRNAs", "Num_4_sgRNAs")]
                         ))
colnames(num_reads_mat) <- character(8)
rownames(num_reads_mat) <- 0:4



CreateEMF("01) Number of reads")
par(mai = c(0.52, 0.5, 0.38, 0.6) / 0.6 * par("cex"))
StackedBars(num_reads_mat, legend_title_vec = c("Mapped", "sgRNAs"),
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "Number of reads", y_label_line = y_label_line
            )
AnnotatePrepoolPostpoolPlot()
dev.off()





num_1MM_mat <- rbind(rowMeans(num_any_mapped_mat[, paste0("Num_1MM_sg", 1:4)]),
                     num_any_mapped_mat[, "Num_mapped_reads"]
                     )
colnames(num_1MM_mat) <- character(8)

CreateEMF("02) Percent 1MM mismatch")
StackedBars(prop.table(num_1MM_mat, margin = 2) * 100, invert_colors = TRUE, show_legend = FALSE,
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "% mapped sgRNAs", y_label_line = y_label_line
            )
mtext("Single-base mismatch", line = 0.5, cex = par("cex"))
AnnotatePrepoolPostpoolPlot()
dev.off()




CreateEMF("03) Mismatched sgRNAs")
num_mismatched_mat <- num_all4_mapped_mat[, c("Num_1_mismatched_sgRNA", "Num_2_mismatched_sgRNAs", "Num_3_mismatched_sgRNAs", "Num_4_mismatched_sgRNAs")]
num_mismatched_mat <- rbind("Num_0_mismatched_sgRNAs" = num_all4_mapped_mat[, "Num_perfect_reads"] - rowSums(num_mismatched_mat),
                            t(num_mismatched_mat)
                            )
colnames(num_mismatched_mat) <- character(8)
rownames(num_mismatched_mat) <- 0:4

StackedBars(prop.table(num_mismatched_mat, margin = 2) * 100, legend_title_vec = c("1MM", "sgRNAs"),
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "% full reads", y_label_line = y_label_line
            )
AnnotatePrepoolPostpoolPlot()
dev.off()






CreateEMF("04) Template switches")
use_columns <- c("Num_0_template_switches", "Num_1_template_switch", "Num_2_template_switches", "Num_3_template_switches")
num_switches_mat <- t(num_all4_mapped_mat[, use_columns])
colnames(num_switches_mat) <- character(8)
rownames(num_switches_mat) <- 0:3

StackedBars(prop.table(num_switches_mat, margin = 2) * 100,
            legend_title_vec = c("Template", "switches"),
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "% full reads", y_label_line = y_label_line
            )
AnnotatePrepoolPostpoolPlot()
dev.off()





use_columns <- c("Num_1_targeted_plasmid",  "Num_2_targeted_plasmids", "Num_3_targeted_plasmids", "Num_4_targeted_plasmids")
targeted_plasmids_mat <- t(num_all4_mapped_mat[, use_columns])
colnames(targeted_plasmids_mat) <- character(8)
rownames(targeted_plasmids_mat) <- 1:4

CreateEMF("05) Targeted plasmids - percentage")
StackedBars(prop.table(targeted_plasmids_mat, margin = 2) * 100,
            legend_title_vec = c("Targeted", "plasmids"),
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "% full reads", y_label_line = y_label_line
            )
AnnotatePrepoolPostpoolPlot()
dev.off()


CreateEMF("06) Targeted plasmids - absolute numbers")
StackedBars(targeted_plasmids_mat[4:1,],
            legend_title_vec = c("Targeted", "plasmids"),
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "Number of full reads", y_label_line = 2.4,
            y_axis_upper_limit = 2.5 * 10^6
            )
AnnotatePrepoolPostpoolPlot()
dev.off()





CreateEMF("07) Switch-backs")
use_columns <- c("Num_0_switch_backs",  "Num_1_switch_back", "Num_2_switch_backs")
switch_backs_mat <- t(num_all4_mapped_mat[, use_columns])
colnames(switch_backs_mat) <- character(8)
rownames(switch_backs_mat) <- 0:2
StackedBars(prop.table(switch_backs_mat, margin = 2) * 100,
            legend_title_vec = c("Switch-", "backs"),
            bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
            y_axis_label = "% full reads", y_label_line = y_label_line
            )
AnnotatePrepoolPostpoolPlot()
dev.off()




PlotSwitchRate <- function(sg1, sg2) {

  switch_rates <- num_all4_mapped_mat[, paste0("Num_switch_sg", sg1, "_to_sg", sg2)]
  switch_rat_mat <- rbind(num_all4_mapped_mat[, "Num_perfect_reads"] - switch_rates, switch_rates)
  colnames(switch_rat_mat) <- character(8)
  StackedBars(prop.table(switch_rat_mat, margin = 2) * 100,
              show_legend = FALSE,
              bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
              y_axis_label = "% full reads", y_label_line = y_label_line
              )
  AnnotatePrepoolPostpoolPlot()
  mtext(paste0("sg", sg1, "\u2013", sg2, " switch"), line = 0.5, cex = par("cex"))
}

CreateEMF("08) Switch rate sg1-2")
PlotSwitchRate(1, 2)
dev.off()

CreateEMF("09) Switch rate sg2-3")
PlotSwitchRate(2, 3)
dev.off()

CreateEMF("10) Switch rate sg3-4")
PlotSwitchRate(3, 4)
dev.off()



PlotNumMatching <- function(sg1, sg2) {

  use_columns <- paste0("Num_", 2:4, "_matching_sg", sg1, "andsg", sg2)
  matching_plasmids_mat <- t(num_all4_mapped_mat[, use_columns])
  colnames(matching_plasmids_mat) <- character(8)
  rownames(matching_plasmids_mat) <- 2:4
  StackedBars(prop.table(matching_plasmids_mat, margin = 2) * 100,
              legend_title_vec = c("Matching", "sgRNAs"),
              bar_positions = PrepoolPostpoolBarPositions(), set_mar = FALSE, draw_box = FALSE,
              y_axis_label = "% full reads", y_label_line = y_label_line
              )
  AnnotatePrepoolPostpoolPlot()
  mtext(paste0("sg", sg1, "\u2009=\u2009", sg2), line = 0.5, cex = par("cex"))
}



CreateEMF("11) Switch rate sg1 matches sg2")
PlotNumMatching(1, 2)
dev.off()

CreateEMF("12) Switch rate sg2 matches sg3")
PlotNumMatching(2, 3)
dev.off()

CreateEMF("13) Switch rate sg3 matches sg4")
PlotNumMatching(3, 4)
dev.off()




