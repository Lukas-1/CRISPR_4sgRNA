### 2022-11-16


# Load packages and source code -------------------------------------------

library("devEMF")

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R")) # For RepositionByGroups
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

first_nanopore_dir <- file.path(first_nanopore_dir,
                                "03_R_objects"
                                )
nanopore_Q20_dir   <- file.path(experiments_directory,
                                "2022-04-06 - first Q20 nanopore run",
                                "03_R_objects"
                                )
pacbio_dir         <- file.path(experiments_directory,
                                "2022-04-06 - PacBio pooled 4sg - first trial",
                                "03_R_objects"
                                )
project_dir        <- file.path(experiments_directory, "2022-11-16 - compare long-read sequencing methods")
output_dir         <- file.path(project_dir, "02_output_data")
figures_dir        <- output_dir



# Load data ---------------------------------------------------------------

load(file.path(nanopore_Q20_dir, "06_assign_sgRNAs_to_plasmids.RData"))
Q20_counts_df <- counts_df
Q20_df <- nano_df
Q20_num_reads <- total_num_reads

load(file.path(first_nanopore_dir, "08_assign_sgRNAs_to_plasmids.RData"))
nano_counts_df <- counts_df
nano_num_reads <- total_num_reads

load(file.path(pacbio_dir, "07_assign_sgRNAs_to_plasmids.RData"))
pb_counts_df <- counts_df
pb_num_reads <- total_num_reads

rm(counts_df)
rm(total_num_reads)



# Define functions --------------------------------------------------------

Percent1MM <- function(input_df) {
  num_MM_mat <- as.matrix(input_df[, paste0("Num_MM_sg", 1:4)])
  sum(num_MM_mat, na.rm = TRUE) / sum(!(is.na(num_MM_mat)))
}

NumMissing <- function(input_counts_df) {
  sum(input_counts_df[, "Count_2_or_more"][!(input_counts_df[, "Is_obsolete"] %in% "Yes")] == 0)
}

FractionSwitched <- function(input_df, are_eligible, first_sg, second_sg) {
  sum(input_df[, paste0("Plasmid_sg", first_sg)][are_eligible] !=
      input_df[, paste0("Plasmid_sg", second_sg)][are_eligible]
      ) / sum(are_eligible)
}


NumMatchingGuides <- function(input_df, are_eligible, first_sg, second_sg) {
  are_pair <- input_df[, paste0("Plasmid_sg", first_sg)][are_eligible] ==
              input_df[, paste0("Plasmid_sg", second_sg)][are_eligible]
  other_sgs <- setdiff(1:4, c(first_sg, second_sg))
  match_pair_A <- input_df[, paste0("Plasmid_sg", other_sgs[[1]])][are_eligible][are_pair] ==
                  input_df[, paste0("Plasmid_sg", first_sg)][are_eligible][are_pair]
  match_pair_B <- input_df[, paste0("Plasmid_sg", other_sgs[[2]])][are_eligible][are_pair] ==
                  input_df[, paste0("Plasmid_sg", first_sg)][are_eligible][are_pair]
  num_matching_vec <- ifelse(match_pair_A & match_pair_B,
                             4L,
                             ifelse(match_pair_A | match_pair_B, 3L, 2L)
                             )
  table(num_matching_vec) / sum(are_eligible)
}



StackedBars <- function(bars_mat,
                        groups_vec         = NULL,
                        y_axis_label       = "",
                        y_label_line       = 2.3,
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
  if (is.null(groups_vec)) {
    bar_positions <- seq_len(ncol(bars_mat))
  } else {
    stopifnot(length(groups_vec) == ncol(bars_mat))
    bar_positions <- RepositionByGroups(groups_vec, gap_ratio = gap_ratio)
  }
  num_bars <- length(bar_positions)
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - side_gap) - (num_bars * 0.04),
                     max(bar_positions) + side_gap  + (num_bars * 0.04)
                    )


  ## Prepare the data axis
  use_numeric_limits <- c(0, max(colSums(bars_mat)) * 1.00)
  numeric_axis_pos <- pretty(use_numeric_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])


  ## Draw the bar plot
  if (show_legend && set_mar) {
    old_mar <- par(mar = c(4, 3.75, 3.75, 6.5))
  }
  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)

  ## Draw the grid
  grid_pos <- pretty(use_numeric_limits, n = 30)
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
  } else if (all(tick_pos[-1] >= 10^6)) {
    tick_labels <- paste0(tick_pos / 10^6, "M")
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



Percent1MMBars <- function(bars_vec, ...) {

  bar_positions <- StackedBars(bars_vec,
                               convert_to_percent = TRUE,
                               invert_colors      = TRUE,
                               y_axis_label       = "Percentage of mapped sgRNAs",
                               show_legend        = FALSE,
                               ...
                               )

  ## Label the percentages
  mtext(paste0(format(round(bars_vec * 100, digits = 1), nsmall = 1), "%"),
        at = bar_positions, side = 3, padj = 1, line = 0.7, cex = par("cex")
        )

  ## Draw the title
  mtext("Single-base mismatch", line = 1.6, cex = par("cex"))

  return(invisible(NULL))
}



MissingBarPlot <- function(bars_vec, ...) {
  BarsOrLollipop(bars_vec,
                 y_axis_label = "Number of missed plasmids",
                 ...
                 )
}



BarsOrLollipop <- function(bars_vec,
                           lollipop      = FALSE,
                           bar_color     = brewer.pal(9, "Blues")[[8]],
                           stem_color    = brewer.pal(9, "Blues")[[4]],
                           y_axis_label  = "",
                           y_label_line  = 2.3,
                           use_title     = NULL,
                           show_title    = TRUE,
                           bar_width     = 0.55,
                           side_space    = 0.5,
                           y_upper_limit = NULL
                           ) {

  ## Determine bar positions
  bar_positions <- seq_along(bars_vec)
  num_bars <- length(bar_positions)
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - side_space) - (num_bars * 0.04),
                     max(bar_positions) + side_space  + (num_bars * 0.04)
                    )

  ## Prepare the data axis
  if (is.null(y_upper_limit)) {
    y_upper_limit <- max(bars_vec)
  }
  use_numeric_limits <- c(0, y_upper_limit)
  numeric_axis_pos <- pretty(use_numeric_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])

  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)

  if (lollipop) {
    segments(x0  = bar_positions,
             y0  = 0,
             y1  = bars_vec,
             col = stem_color,
             xpd = NA,
             lwd = par("lwd") * 2
             )
    points(x   = bar_positions,
           y   = bars_vec,
           col = bar_color,
           cex = par("cex") * 1.5,
           pch = 16,
           xpd = NA
           )
  } else {
    PlotBarplotMat(t(bars_vec),
                   colors_vec    = bar_color,
                   positions_vec = bar_positions,
                   bar_width     = bar_width
                   )
  }

  ## Draw the y axis
  axis(2,
       las = 2,
       mgp = c(3, 0.45, 0),
       tcl = -0.35,
       lwd = par("lwd")
       )
  mtext(VerticalAdjust(y_axis_label),
        side = 2,
        line = y_label_line,
        cex  = par("cex")
        )

  ## Draw the bar labels
  mtext(sub(" ", "\n", names(bars_vec)),
        at = bar_positions, side = 1, padj = 1, line = -0.25, cex = par("cex")
        )

  ## Final steps
  box(bty = "l")
  return(invisible(bar_positions))
}



ThreeHistograms <- function(numeric_list, x_axis_label_line = 2, only_annotation = NULL, ...) {

  new_order <- order(match(c("PacBio", "Nanopore Q20+", "Nanopore"), names(numeric_list)))
  numeric_list <- numeric_list[new_order]

  ## Define color scheme
  hist_colors <- c(brewer.pal(9, "Blues")[[7]],
                   brewer.pal(9, "Reds")[[7]],
                   brewer.pal(9, "Greys")[[7]]
                   )

  ## Draw histograms
  DrawHistogram(numeric_list,
                hist_colors         = hist_colors,
                num_breaks          = 250,
                use_exact_breaks    = TRUE,
                x_axis_space        = 0.02,
                x_axis_label        = "Plasmid-level read count",
                y_axis_label        = VerticalAdjust("Frequency"),
                title_font          = 1,
                truncation_limit    = 250,
                # x_axis_upper_limit  = 250,
                y_axis_upper_limit  = 1600,
                fixed_y_upper_limit = TRUE,
                x_axis_mgp          = 0.4,
                y_axis_mgp          = 0.45,
                x_axis_label_line   = x_axis_label_line,
                y_axis_label_line   = 2.3,
                show_both           = TRUE,
                only_annotation     = only_annotation,
                ...
                )

  ## Draw legend
  if (!(isFALSE(only_annotation))) {
    x_start <- par("usr")[[2]] - strwidth("Nanopore Q20") - diff(grconvertX(c(0, 1.8), from = "lines", to = "user"))
    y_start <- par("usr")[[4]] - diff(grconvertY(c(0, 1), from = "lines", to = "user"))
    lines_seq <- seq_along(hist_colors) - 1L
    y_vec <- y_start - diff(grconvertY(c(0, 1.2), from = "lines", to = "user")) * lines_seq
    segments(x0  = x_start,
             x1  = x_start + diff(grconvertX(c(0, 0.5), from = "lines", to = "user")),
             y0  = y_vec,
             col = hist_colors,
             lwd = par("lwd") * 2,
             xpd = NA
             )
    text(x      = x_start + diff(grconvertX(c(0, 0.9), from = "lines", to = "user")),
         y      = y_vec,
         labels = names(numeric_list),
         adj    = c(0, 0.5),
         xpd    = NA
         )
  }

  return(invisible(NULL))
}



TemplateSwitchPerRegion <- function(region_vec_list,
                                    use_title = "Template switches between sgRNAs",
                                    y_axis_label = "Percentage of reads",
                                    title_line = 0.7,
                                    bar_width = 0.65,
                                    NQP_in_gray = FALSE,
                                    ...
                                    ) {

  numeric_mat <- do.call(cbind,
                         lapply(region_vec_list,
                                function(x) rbind("No" = 1 - x,  "Yes" = x)
                                )
                         )
  colnames(numeric_mat) <- rep(NA, ncol(numeric_mat))
  groups_vec <- rep(1:3, each = 3)
  bar_positions <- StackedBars(numeric_mat, groups_vec = groups_vec,
                               y_axis_label = y_axis_label,
                               convert_to_percent = TRUE, bar_width = bar_width,
                               show_legend = FALSE,
                               ...
                               )

  mtext(use_title, line = title_line, cex = par("cex"))

  mtext(rep(c("N", "Q", "P"), times = 3), side = 1, padj = 1,
        at = bar_positions, line = -0.5, cex = par("cex"),
        col = if (NQP_in_gray) "gray55" else "black"
        )
  group_labels <- expression(
    "sg1" * scriptscriptstyle(" ") * "\u20132",
    "sg2" * scriptscriptstyle(" ") * "\u20133",
    "sg3" * scriptscriptstyle(" ") * "\u20134",
  )
  # group_labels <- c("sg1&2", "sg2&3", "sg3&4")
  segments(x0  = tapply(bar_positions, groups_vec, min) - 0.25,
           x1  = tapply(bar_positions, groups_vec, max) + 0.25,
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.15), from = "lines", to = "user")),
           col = "gray60",
           xpd = NA
           )
  mtext(group_labels, at = tapply(bar_positions, groups_vec, mean),
        side = 1, line = 1.5, cex = par("cex")
        )

  return(invisible(NULL))
}



NumSgRNAsPerPair <- function(pairs_vec_list,
                             bar_width = 0.65,
                             NQP_in_gray = FALSE,
                             y_axis_label = "Percentage of reads",
                             ...
                             ) {

  numeric_mat <- do.call(cbind, pairs_vec_list)
  colnames(numeric_mat) <- rep(NA, ncol(numeric_mat))
  groups_vec <- rep(1:3, each = 3)
  bar_positions <- StackedBars(numeric_mat,
                               convert_to_percent = TRUE,
                               legend_title_vec   = c("Matching", "sgRNAs"),
                               y_axis_label       = y_axis_label,
                               groups_vec         = rep(1:3, each = 3),
                               bar_width          = bar_width,
                               ...
                               )
  x_bit <- 0.6
  y_bit <- diff(grconvertY(c(0, 0.1), from = "npc", to = "user"))
  rect(xleft   = tapply(bar_positions, groups_vec, max)[c(1, 2)] + x_bit,
       xright  = tapply(bar_positions, groups_vec, min)[c(2, 3)] - x_bit,
       ybottom = par("usr")[[3]] - y_bit,
       ytop    = par("usr")[[4]] + y_bit,
       col     = "white",
       border  = NA,
       xpd     = NA
       )

  mtext(rep(c("N", "Q", "P"), times = 3), side = 1, padj = 1,
        at = bar_positions, line = -0.5, cex = par("cex"),
        col = if (NQP_in_gray) "gray55" else "black"
        )

  group_labels <- expression(
    "sg1" * scriptscriptstyle(" ") * "=" * scriptscriptstyle(" ") * "2",
    "sg2" * scriptscriptstyle(" ") * "=" * scriptscriptstyle(" ") * "3",
    "sg3" * scriptscriptstyle(" ") * "=" * scriptscriptstyle(" ") * "4",
  )
  segments(x0  = tapply(bar_positions, groups_vec, min) - 0.25,
           x1  = tapply(bar_positions, groups_vec, max) + 0.25,
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.15), from = "lines", to = "user")),
           col = "gray60",
           xpd = NA
           )
  mtext(group_labels, at = tapply(bar_positions, groups_vec, mean),
        side = 1, line = 1.5, cex = par("cex")
        )
  return(invisible(NULL))
}



DistancesBarPlot <- function(distances_vec,
                             add_bp        = TRUE,
                             use_title     = "Distance between sgRNAs",
                             title_line    = 0.7,
                             label_beneath = TRUE,
                             bar_color     = brewer.pal(9, "Blues")[[7]],
                             ...
                             ) {

  names(distances_vec) <- NA
  bar_positions <- BarsOrLollipop(distances_vec,
                                  y_axis_label  = "Base pairs",
                                  y_upper_limit = 1000,
                                  bar_color     = bar_color,
                                  side_space    = 0.6,
                                  ...
                                  )

  mtext(use_title, line = title_line, cex = par("cex"))
  group_labels <- expression(
    "sg1" * scriptscriptstyle(" ") * "\u20132",
    "sg2" * scriptscriptstyle(" ") * "\u20133",
    "sg3" * scriptscriptstyle(" ") * "\u20134"
  )
  mtext(group_labels, at = bar_positions,
        side = 1, line = 0.75, cex = par("cex")
        )
  distance_labels <- paste0(distances_vec, if (add_bp) " bp" else "")
  if (label_beneath) {
    text(x      = bar_positions,
         y      = distances_vec - diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
         labels = distance_labels,
         col    = Palify(bar_color, fraction_pale = 0.85),
         xpd    = NA
         )
  } else {
    text(x      = bar_positions,
         y      = distances_vec + diff(grconvertY(c(0, 0.6), from = "lines", to = "user")),
         labels = distance_labels,
         xpd    = NA
         )
  }

  return(invisible(NULL))
}



# Assemble data -----------------------------------------------------------

read_counts_mat <- cbind(
  "Nanopore"      = c("0" = nano_num_reads - nrow(nano_df),
                      table(nano_df[, "Num_matched_sgRNAs"])
                      ),
  "Nanopore Q20+" = c("0" = Q20_num_reads - nrow(Q20_df),
                      table(Q20_df[, "Num_matched_sgRNAs"])
                      ),
  "PacBio"        = c("0" = pb_num_reads - nrow(pb_df),
                      table(pb_df[, "Num_matched_sgRNAs"])
                      )
)

fraction_1MM_vec <- c(
  "Nanopore"      = Percent1MM(nano_df),
  "Nanopore Q20+" = Percent1MM(Q20_df),
  "PacBio"        = Percent1MM(pb_df)
)

num_missing_vec <- c(
  "Nanopore"      = NumMissing(nano_counts_df),
  "Nanopore Q20+" = NumMissing(Q20_counts_df),
  "PacBio"        = NumMissing(pb_counts_df)
)

counts_vec_list <- list(
  "Nanopore"      = nano_counts_df[!(nano_counts_df[, "Is_obsolete"] %in% "Yes"), "Count_2_or_more"],
  "Nanopore Q20+" = Q20_counts_df[!(Q20_counts_df[, "Is_obsolete"] %in% "Yes"), "Count_2_or_more"],
  "PacBio"        = pb_counts_df[!(pb_counts_df[, "Is_obsolete"] %in% "Yes"), "Count_2_or_more"]
)

are_nano_4 <- nano_df[, "Num_matched_sgRNAs"] == 4
are_Q20_4  <- Q20_df[, "Num_matched_sgRNAs"] == 4
are_pb_4   <- pb_df[, "Num_matched_sgRNAs"] == 4


num_switches_mat <- cbind(
  "Nanopore"      = table(nano_df[, "Num_template_switches"][are_nano_4]),
  "Nanopore Q20+" = table(Q20_df[, "Num_template_switches"][are_Q20_4]),
  "PacBio"        = table(pb_df[, "Num_template_switches"][are_pb_4])
)

num_targeted_plasmids_mat <- cbind(
  "Nanopore"      = table(nano_df[, "Num_targeted_plasmids"][are_nano_4]),
  "Nanopore Q20+" = table(Q20_df[, "Num_targeted_plasmids"][are_Q20_4]),
  "PacBio"        = table(pb_df[, "Num_targeted_plasmids"][are_pb_4])
)

num_switch_backs_mat <- cbind(
  "Nanopore"      = table(nano_df[, "Num_switch_backs"][are_nano_4]),
  "Nanopore Q20+" = table(Q20_df[, "Num_switch_backs"][are_Q20_4]),
  "PacBio"        = table(pb_df[, "Num_switch_backs"][are_pb_4])
)


fraction_switched_vec_list <- list(
  "sg1 != sg2" = c(
    "Nanopore"      = FractionSwitched(nano_df, are_nano_4, 1, 2),
    "Nanopore Q20+" = FractionSwitched(Q20_df, are_Q20_4, 1, 2),
    "PacBio"        = FractionSwitched(pb_df, are_pb_4, 1, 2)
  ),
  "sg2 != sg3" = c(
    "Nanopore"      = FractionSwitched(nano_df, are_nano_4, 2, 3),
    "Nanopore Q20+" = FractionSwitched(Q20_df, are_Q20_4, 2, 3),
    "PacBio"        = FractionSwitched(pb_df, are_pb_4, 2, 3)
  ),
  "sg3 != sg4" = c(
    "Nanopore"      = FractionSwitched(nano_df, are_nano_4, 3, 4),
    "Nanopore Q20+" = FractionSwitched(Q20_df, are_Q20_4, 3, 4),
    "PacBio"        = FractionSwitched(pb_df, are_pb_4, 3, 4)
  )
)


num_matching_mat_list <- list(
  "sg1 == sg2" = cbind(
    "Nanopore"      = NumMatchingGuides(nano_df, are_nano_4, 1, 2),
    "Nanopore Q20+" = NumMatchingGuides(Q20_df, are_Q20_4, 1, 2),
    "PacBio"        = NumMatchingGuides(pb_df, are_pb_4, 1, 2)
  ),
  "sg2 == sg3" = cbind(
    "Nanopore"      = NumMatchingGuides(nano_df, are_nano_4, 2, 3),
    "Nanopore Q20+" = NumMatchingGuides(Q20_df, are_Q20_4, 2, 3),
    "PacBio"        = NumMatchingGuides(pb_df, are_pb_4, 2, 3)
  ),
  "sg3 == sg4" = cbind(
    "Nanopore"      = NumMatchingGuides(nano_df, are_nano_4, 3, 4),
    "Nanopore Q20+" = NumMatchingGuides(Q20_df, are_Q20_4, 3, 4),
    "PacBio"        = NumMatchingGuides(pb_df, are_pb_4, 3, 4)
  )
)


sg_distances_vec <- c(
  "sg1-sg2" = 719L,
  "sg2-sg3" = 339L,
  "sg3-sg4" = 359L
)



# Draw plots --------------------------------------------------------------

par(mar = c(4, 3.75, 3.75, 2))
StackedBars(read_counts_mat,
            legend_title_vec = c("Mapped", "sgRNAs"),
            y_axis_label = "Number of reads"
            )
Percent1MMBars(fraction_1MM_vec)
ThreeHistograms(counts_vec_list)
MissingBarPlot(num_missing_vec)
MissingBarPlot(num_missing_vec, lollipop = TRUE)

StackedBars(num_switches_mat,
            legend_title_vec = c("Template", "switches"),
            convert_to_percent = TRUE,
            y_axis_label = "Percentage of reads"
            )

StackedBars(num_targeted_plasmids_mat,
            legend_title_vec = c("Targeted", "plasmids"),
            convert_to_percent = TRUE,
            y_axis_label = "Percentage of reads",
            draw_box = FALSE
            )

StackedBars(num_switch_backs_mat,
            legend_title_vec = c("Switch-", "backs"),
            convert_to_percent = TRUE,
            y_axis_label = "Percentage of reads"
            )

TemplateSwitchPerRegion(fraction_switched_vec_list)
DistancesBarPlot(sg_distances_vec)

NumSgRNAsPerPair(num_matching_mat_list)



# Export plots ------------------------------------------------------------

scaling_factor <- 20L
correction_factor <- 0.95

legend_par_list <- list(
  "mai" = c(0.45, 0.5, 0.35, 0.6),
  "lwd" = 0.8,
  "cex" = 0.7
)
no_legend_par_list <- c(
  list("mai" = c(0.45, 0.5, 0.35, 0.2)),
  legend_par_list[c("lwd", "cex")]
)

scaled_legend_par_list <- list(
  "mar" = legend_par_list[["mai"]] * 5 / (legend_par_list[["cex"]] * correction_factor),
  "lwd" = legend_par_list[["lwd"]] * scaling_factor,
  "cex" = legend_par_list[["cex"]] * scaling_factor * correction_factor
)
scaled_no_legend_par_list <- c(
  list("mar" = c(0.45, 0.5, 0.35, 0.2) * 5 / (legend_par_list[["cex"]] * correction_factor)),
  scaled_legend_par_list[c("lwd", "cex")]
)


devEMF::emf(file.path(output_dir, "1A) Number of mapped sgRNAs.emf"),
            width = 3 * scaling_factor, height = 2.2 * scaling_factor,
            emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, scaled_legend_par_list)
StackedBars(read_counts_mat,
            legend_title_vec = c("Mapped", "sgRNAs"),
            y_axis_label = "Number of reads", set_mar = FALSE, grid_lwd = 0.7
            )
dev.off()



devEMF::emf(file.path(output_dir, "1B) Percentage 1MM.emf"),
            width = 2.6 * scaling_factor, height = 2.2 * scaling_factor,
            emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, scaled_no_legend_par_list)
Percent1MMBars(fraction_1MM_vec, grid_lwd = 0.7)
dev.off()





devEMF::emf(file.path(output_dir, "1C) Plasmid count histograms - only annotation.emf"),
            width = 3 * scaling_factor, height = 2.2 * scaling_factor,
            emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, scaled_no_legend_par_list)
ThreeHistograms(counts_vec_list, use_lwd = 1, x_axis_label_line = 1.7,
                only_annotation = TRUE
                )
dev.off()



svglite::svglite(file.path(output_dir, "1C) Plasmid count histograms.svg"),
                 width = 3 * scaling_factor, height = 2.2 * scaling_factor,
                 bg = "transparent"
                 )
do.call(par, scaled_no_legend_par_list)
ThreeHistograms(counts_vec_list, use_lwd = 0.9, x_axis_label_line = 1.7,
                only_annotation = FALSE
                )
dev.off()



devEMF::emf(file.path(output_dir, "1D) Number of missing plasmids.emf"),
            width = 2.6 * scaling_factor, height = 2.2 * scaling_factor,
            emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, scaled_no_legend_par_list)
MissingBarPlot(num_missing_vec)
dev.off()



devEMF::emf(file.path(output_dir, "2A) Number of template switches.emf"),
            width = 2.4, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, legend_par_list)
abbr_num_switches_mat <- num_switches_mat
colnames(abbr_num_switches_mat)[[2]] <- "Q20+"
StackedBars(abbr_num_switches_mat,
            legend_title_vec = c("Template", "switches"),
            convert_to_percent = TRUE,
            y_axis_label = "Reads",
            set_mar = FALSE,
            point_x_start = -0.07, legend_lwd = 0.8,
            lines_x_start = 1.2,
            title_x_start = -0.2, y_label_line = 2.2, grid_lwd = 0.7
            )
dev.off()



devEMF::emf(file.path(output_dir, "2B) Template switches between sgRNAs.emf"),
            width = 2.4, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, legend_par_list)
TemplateSwitchPerRegion(fraction_switched_vec_list,
                        use_title = "Switch rate",
                        y_axis_label = "Reads",
                        title_line = 0.5, side_gap = 0.6, gap_ratio = 1.5,
                        bar_width = 0.725, NQP_in_gray = TRUE, y_label_line = 2.2,
                        grid_lwd = 0.7
                        )
dev.off()



devEMF::emf(file.path(output_dir, "2C) Distances between sgRNAs.emf"),
            width = 2, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, no_legend_par_list)
DistancesBarPlot(sg_distances_vec, add_bp = FALSE, use_title = "Distance",
                 title_line = 0.5, y_label_line = 2.2, bar_width = 0.65
                 )
dev.off()



devEMF::emf(file.path(output_dir, "2D) Number of targeted plasmids.emf"),
            width = 2.4, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, legend_par_list)
abbr_num_targeted_plasmids_mat <- num_targeted_plasmids_mat
colnames(abbr_num_targeted_plasmids_mat)[[2]] <- "Q20+"
StackedBars(abbr_num_targeted_plasmids_mat,
            legend_title_vec = c("Targeted", "plasmids"),
            convert_to_percent = TRUE,
            y_axis_label = "Reads",
            draw_box = FALSE, set_mar = FALSE,
            point_x_start = -0.07, legend_lwd = 0.8,
            lines_x_start = 1.2,
            title_x_start = -0.2, y_label_line = 2.2, grid_lwd = 0.7
            )
dev.off()



devEMF::emf(file.path(output_dir, "2E) Number of switch-backs.emf"),
            width = 2.4, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, legend_par_list)
abbr_num_switch_backs_mat <- num_switch_backs_mat
colnames(abbr_num_switch_backs_mat)[[2]] <- "Q20+"
StackedBars(abbr_num_switch_backs_mat,
            legend_title_vec = c("Switch-", "backs"),
            convert_to_percent = TRUE,
            y_axis_label = "Reads", set_mar = FALSE,
            point_x_start = -0.07, legend_lwd = 0.8,
            lines_x_start = 1.2,
            title_x_start = -0.2, y_label_line = 2.2, grid_lwd = 0.7
            )
dev.off()



devEMF::emf(file.path(output_dir, "2F) Number of matching sgRNAs.emf"),
            width = 2.4, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
do.call(par, legend_par_list)
NumSgRNAsPerPair(num_matching_mat_list, set_mar = FALSE,
                 point_x_start = -0.07, legend_lwd = 0.8,
                 lines_x_start = 1.2,
                 title_x_start = -0.2, side_gap = 0.6, gap_ratio = 1.5,
                 bar_width = 0.725, NQP_in_gray = TRUE,
                 y_axis_label = "Reads", y_label_line = 2.2, grid_lwd = 0.7
                 )
dev.off()



# Calculate the hypothetical per-sample costs of sequencing ---------------

sum(are_nano_4)
sum(are_Q20_4)
sum(are_pb_4)

signif(sum(are_nano_4), 2)
signif(sum(are_Q20_4), 2)
signif(sum(are_pb_4), 2)

300 * (22000 / 260000) * 1000
300 / (230000 / 22000) * 500
300 / (170000 / 22000) * 2000

signif(300 * (22000 / 260000) * 1000, 2)
signif(300 / (230000 / 22000) * 500, 2)
signif(300 / (170000 / 22000) * 2000, 2)

sum(nano_counts_df[, "Count_perfect"], na.rm = TRUE)
sum(Q20_counts_df[, "Count_perfect"], na.rm = TRUE)
sum(pb_counts_df[, "Count_perfect"], na.rm = TRUE)




# Calculate additional metrics of interest --------------------------------

colSums(num_targeted_plasmids_mat[3, , drop = FALSE]) / colSums(num_switches_mat[3:4, , drop = FALSE]) * 100





