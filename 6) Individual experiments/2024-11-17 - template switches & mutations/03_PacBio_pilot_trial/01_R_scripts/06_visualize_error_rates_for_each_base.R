## 2024-12-31


# Load packages and source code -------------------------------------------

root_dir      <- "~/CRISPR_4sgRNA"
exper_dir     <- file.path(root_dir, "6) Individual experiments")
illumina_dir  <- file.path(exper_dir, "2022-04-21 - Illumina paired-end 2sg - first trial")
lumi_func_dir <- file.path(illumina_dir, "01_R_scripts", "R_functions")
project_dir   <- file.path(exper_dir, "2024-11-17 - template switches & mutations")

source(file.path(lumi_func_dir, "01_violin_swarm_plots.R")) # For RepositionByGroups
source(file.path(lumi_func_dir, "05_creating_figures_from_count_data.R")) # For DrawSideLegend and VerticalAdjust
source(file.path(project_dir, "01_R_functions", "03_visualizing_error_rates_for_subsequences.R"))


# Define paths ------------------------------------------------------------

plasmids_dir <- file.path(project_dir, "02_PacBio_plasmids", "02_R_objects")
sub_dir      <- file.path(project_dir, "03_PacBio_pilot_trial")
rdata_dir    <- file.path(sub_dir, "02_R_objects")
output_dir   <- file.path(sub_dir, "03_output")



# Load data ---------------------------------------------------------------

load(file.path(plasmids_dir, "03_compute_error_rates_for_each_base.RData"))
load(file.path(rdata_dir, "01_extract_and_categorize_subsequences__features_df.RData"))
df_list_names <- load(file.path(rdata_dir, "05_compute_error_rates_for_each_base.RData"))



# Define labels -----------------------------------------------------------

short_feature_labels <- c(
  "promoter1_hU6"  = "hU6", # (hU6)
  "sg1"            = "sg1",
  "tracrRNA1"      = "cr1",
  "polyT_1"        = "poly(T)",

  "EM7_promoter"   = "EM7 promoter",
  "pre_TpR"        = "unannotated",
  "TpR_DHFR"       = "TpR",
  "polyT_TpR"      = "poly(T)",

  "promoter2_mU6"  = "mU6", # (mU6)
  "sg2"            = "sg2",
  "tracrRNA2"      = "cr2",
  "polyT_2"        = "poly(T)",

  "promoter3_hH1"  = "hH1", # (hH1)
  "sg3"            = "sg3",
  "tracrRNA3"      = "cr3",
  "polyT_3"        = "poly(T)",

  "promoter4_h7SK" = "h7SK", # (h7SK)
  "sg4"            = "sg4",
  "tracrRNA4"      = "cr4",
  "polyT_4"        = "poly(T)"
)



# Define functions --------------------------------------------------------

MakeEmptyPlot <- function(x_limits = c(0, 1), y_limits = c(0, 1)) {
  plot.new()
  plot.window(xlim = x_limits,
              ylim = y_limits,
              xaxs = "i",
              yaxs = "i"
              )
}


RollingOutliers <- function(numeric_vec, use_window = 20, outlier_dist = Inf) {
  fill_vec <- rep(NA, use_window - 1L)
  padded_vec <- c(fill_vec, numeric_vec, fill_vec)
  start_vec <- seq(from = 1, to = length(padded_vec) - use_window + 1L)
  stop_vec <- seq(from = use_window, to = length(padded_vec))
  indices_list <- mapply(function(x, y) seq(x, y), start_vec, stop_vec, SIMPLIFY = FALSE)
  data_list <- lapply(indices_list, function(x) padded_vec[x])
  use_indices <- seq(from = use_window / 2, to = (use_window / 2) + length(numeric_vec) - 1)
  results_df <- data.frame(
    "original"      = numeric_vec,
    "rolling_mean"  = vapply(data_list, mean, na.rm = TRUE, numeric(1))[use_indices],
    "rolling_sd"    = vapply(data_list, sd, na.rm = TRUE, numeric(1))[use_indices]
  )
  results_df[, "rolling_upper"] <- results_df[, "rolling_mean"] + 2 * results_df[, "rolling_sd"]
  results_df[, "is_outlier"] <- numeric_vec > results_df[, "rolling_upper"]
  results_df[, "is_outlier"] <- numeric_vec > (results_df[, "rolling_mean"] + outlier_dist)
  # results_df[, "is_outlier"] <- FALSE
  return(results_df)
}



SixErrorRates <- function(display_df_list, show_title = NULL, y_upper = NULL) {

  stopifnot("features_df" %in% ls(envir = globalenv()))

  ## Define the x axis limits
  range_padding <- 30L
  x_padding_fraction <- 0.02
  show_range <- seq(from = min(features_df[, "Start"]) - range_padding,
                    to   = max(features_df[, "End"]) + range_padding
                    )
  x_spacing <- (max(show_range) - min(show_range)) * x_padding_fraction
  x_limits <- c(min(show_range) - x_spacing, max(show_range) + x_spacing)


  ## Define the y axis limits
  fractions_columns <- c(
    "Fraction_incorrect_nonswitched", "Fraction_incorrect_switched",
    "Fraction_reference"
  )
  y_max <- max(unlist(lapply(display_df_list, function(y) {
    are_within <- y[, "Base_number"] %in% show_range
    max(as.matrix(y[are_within, fractions_columns]), na.rm = TRUE)
  })))
  if (is.null(y_upper)) {
    y_padding_fraction <- x_padding_fraction
    y_upper <- y_max + y_max * y_padding_fraction
  }
  y_limits <- c(0, y_upper)


  ## Prepare for drawing x axis gridlines
  are_sg <- grepl("^sg[1-4]", features_df[, "Feature"])
  are_polyT <- grepl("polyT", features_df[, "Feature"], fixed = TRUE)
  are_to_mark <- are_sg | are_polyT
  x_markers <- c(features_df[, "Start"][are_to_mark] - 0.5,
                 features_df[, "End"][are_to_mark] + 0.5
                 )


  ## Prepare the plot layout
  layout_mat <- cbind(rep(1L, 13), 5:17, c(2L, rep(4L, 11), 3L))
  gap_ratio <- 4
  gap_total <- 0.3
  gap_height <- gap_total / (3 + 2 * gap_ratio)
  plot_height <- (1 - gap_total) / 4
  top_bottom_correction <- 0.01
  layout(layout_mat,
         widths = c(0.125, 0.8, 0.125),
         heights = c((gap_height * gap_ratio) - top_bottom_correction,
                     rep(c(plot_height, gap_height), 5),
                     plot_height,
                     (gap_height * gap_ratio) + top_bottom_correction
                     )
         )
  old_par <- par(mar = rep(0, 4))

  for (i in 1:4) {
    MakeEmptyPlot()
  }

  ## Draw the legend
  legend_labels_list <- list(
    c("No", "switch"),
    c("With", "switch"),
    c("Original", "plasmids")
  )
  reference_color <- "gray25"
  legend_colors <- c(hcl.colors(9, "Blues")[[2]], hcl.colors(9, "Reds")[[2]], reference_color)#, "Purples")
  DrawSideLegend(labels_list          = legend_labels_list,
                 use_colors           = vapply(legend_colors, function(x) adjustcolor(x, alpha.f = 0.5), ""),
                 border_colors        = legend_colors,
                 lines_x_start        = 0.75,
                 point_x_start        = 0.15,
                 use_pch              = 22,
                 border_lwd           = 1.15,
                 use_point_size       = 1.4,
                 small_gap_size       = 1.1,
                 large_gap_multiplier = 1.7,
                 x_starting_point     = 0
                 )


  for (i in 1:6) {
    display_df <- display_df_list[[i]]

    MakeEmptyPlot()
    if ((i == 1) && (!(is.null(show_title)))) {
      text(x      = 0.5,
           y      = 0.5,
           labels = show_title,
           cex    = 1 / 0.66 * 0.8
           )
    }
    MakeEmptyPlot(x_limits = x_limits, y_limits = y_limits)


    ## Draw the x axis markers
    axis_color <- "gray50"
    segments(x0  = x_markers,
             y0  = par("usr")[[3]],
             y1  = par("usr")[[4]],
             col = "gray90"
             )
    box(col = axis_color)


    ## Draw and label the y axis
    tick_pos <- axTicks(2)
    axis(2,
         at     = tick_pos,
         labels = paste0(format(tick_pos * 100), "%"),
         las    = 2,
         mgp    = c(3, 0.45, 0),
         tcl    = -0.375,
         col    = axis_color
         )
    ylab_text <- gsub("sg", "", names(display_df_list)[[i]], fixed = TRUE)
    ylab_text <- gsub("_", "\u2009\u2013\u2009", ylab_text, fixed = TRUE)
    ylab_text <- paste0("sg", ylab_text, " pairing")
    text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 3.65), from = "lines", to = "user")),
         y      = mean(par("usr")[3:4]),
         srt    = 90,
         labels = ylab_text,
         xpd    = NA
         )


    ## Prepare for plotting data
    are_within_range  <- display_df[, "Base_number"] %in% show_range
    x_vec             <- display_df[are_within_range, "Base_number"]
    reference_y_vec   <- display_df[are_within_range, "Fraction_reference"]
    nonswitched_y_vec <- display_df[are_within_range, "Fraction_incorrect_nonswitched"]
    switched_y_vec    <- display_df[are_within_range, "Fraction_incorrect_switched"]
    if (y_max > 0.1) {
      max_distance <- 0.01
    } else {
      max_distance <- 0.001
    }
    reference_outliers_df   <- RollingOutliers(reference_y_vec,   outlier_dist = max_distance)
    nonswitched_outliers_df <- RollingOutliers(nonswitched_y_vec, outlier_dist = max_distance)
    switched_outliers_df    <- RollingOutliers(switched_y_vec,    outlier_dist = max_distance)

    sg_indices <- Map(function(x, y) x:y, features_df[are_sg, "Start"], features_df[are_sg, "End"])
    reference_outliers_df[, "is_outlier"] <- ifelse(x_vec %in% unlist(sg_indices),
                                                    FALSE,
                                                    reference_outliers_df[, "is_outlier"]
                                                    )

    ## Prepare colors
    reference_line_color    <- adjustcolor(reference_color, alpha.f = 0.65)
    nonswitched_line_color  <- adjustcolor(hcl.colors(9, "Blues")[[2]], alpha.f = 0.65)
    switched_line_color     <- adjustcolor(hcl.colors(9, "Reds")[[2]],  alpha.f = 0.65)

    reference_point_color   <- adjustcolor(reference_color, alpha.f = 0.65)
    nonswitched_point_color <- adjustcolor(hcl.colors(9, "Blues")[[2]], alpha.f = 0.65)
    switched_point_color    <- adjustcolor(hcl.colors(9, "Reds")[[2]], alpha.f = 0.65)


    ## Draw lines
    use_lwd <- 1.25

    lines(x   = x_vec[!(reference_outliers_df[, "is_outlier"])],
          y   = reference_y_vec[!(reference_outliers_df[, "is_outlier"])],
          col = reference_line_color,
          lwd = use_lwd
          )
    lines(x   = x_vec[!(nonswitched_outliers_df[, "is_outlier"])],
          y   = nonswitched_y_vec[!(nonswitched_outliers_df[, "is_outlier"])],
          col = nonswitched_line_color,
          lwd = use_lwd
          )
    lines(x   = x_vec[!(switched_outliers_df[, "is_outlier"])],
          y   = switched_y_vec[!(switched_outliers_df[, "is_outlier"])],
          col = switched_line_color,
          lwd = use_lwd
          )

    ## Draw points (for outliers)$
    use_pt_cex <- 0.4
    points(x   = x_vec[reference_outliers_df[, "is_outlier"]],
           y   = reference_y_vec[reference_outliers_df[, "is_outlier"]],
           cex = use_pt_cex,
           pch = 16,
           col = reference_point_color
           )
    points(x   = x_vec[nonswitched_outliers_df[, "is_outlier"]],
           y   = nonswitched_y_vec[nonswitched_outliers_df[, "is_outlier"]],
           cex = use_pt_cex,
           pch = 16,
           col = nonswitched_point_color
           )
    points(x   = x_vec[switched_outliers_df[, "is_outlier"]],
           y   = switched_y_vec[switched_outliers_df[, "is_outlier"]],
           cex = use_pt_cex,
           pch = 16,
           col = switched_point_color
           )
  }

  ###############################
  ## Draw the bottom schematic ##
  ###############################

  ## Prepare colors for the schematic
  colors_vec <- rep("gray90", nrow(features_df))
  colors_vec[grepl("promoter", features_df[, "Feature"], fixed = TRUE)] <- "white"
  colors_vec[are_sg] <- "#ff00ff"
  colors_vec[grepl("^tracrRNA[1-4]", features_df[, "Feature"])] <- "#57ddff" #"#00ccff"
  colors_vec[features_df[, "Feature"] == "TpR_DHFR"] <- "#ccffcc"
  colors_vec[are_polyT] <- "gray60"

  ## Prepare the y location
  y_line <- diff(grconvertY(c(0, 1), from = "lines", to = "user"))
  y_top <- par("usr")[[3]] - y_line * 1.25
  y_bottom <- y_top - y_line


  ## Draw the line segments for the labels that don't fit
  text_cex <- 0.9
  are_to_label <- features_df[, "Feature"] %in% c("sg1", "polyT_1")
  segment_length <- y_line * 0.3
  segments(x0  = rowMeans(features_df[are_to_label, c("Start", "End")]),
           y0  = y_bottom,
           y1  = y_bottom - segment_length,
           xpd = NA
           )

  ## Draw the colored elements and their labels
  for (row_index in seq_len(nrow(features_df))) {
    current_feature <- features_df[, "Feature"][[row_index]]
    if (current_feature %in% "pre_TpR") {
      next
    }
    x_start <- features_df[, "Start"][[row_index]]
    x_end <- features_df[, "End"][[row_index]]
    rect(xleft   = x_start - 0.5,
         xright  = x_end + 0.5,
         ybottom = y_bottom,
         ytop    = y_top,
         border  = NA,
         xpd     = NA,
         col     = colors_vec[[row_index]]
         )
    if ((x_end - x_start) > 50) {
      text(x      = mean(c(x_start, x_end)),
           y      = mean(c(y_bottom, y_top)) + y_line * 0.01,
           labels = VerticalAdjust(short_feature_labels[[current_feature]]),
           cex    = text_cex,
           xpd    = NA
           )
    } else if (current_feature == "sg1") {
      text(x      = mean(c(x_start, x_end)),
           y      = y_bottom - (y_line * 0.25),
           adj    = c(0.5, 1),
           labels = VerticalAdjust(short_feature_labels[[current_feature]]),
           cex    = text_cex,
           col    = "#d100d1",
           xpd    = NA
           )
    } else if (current_feature == "polyT_1") {
      text(x      = mean(c(x_start, x_end)) - strwidth("sg1", units = "user", cex = text_cex) / 2,
           y      = y_bottom - (y_line * 0.25),
           adj    = c(0, 1),
           labels = VerticalAdjust(short_feature_labels[[current_feature]]),
           cex    = text_cex,
           col    = "gray25",
           xpd    = NA
           )
    }
  }

  ## Draw the schematic element borders
  schematic_lwd <- 0.75
  segments(x0  = min(features_df[, "Start"]),
           x1  = max(features_df[, "End"]),
           y0  = c(y_top, y_bottom),
           lwd = schematic_lwd,
           xpd = NA
           )
  all_transitions <- unique(c(features_df[, "Start"] - 0.5, features_df[, "End"] + 0.5))
  segments(x0  = all_transitions,
           y0  = y_top,
           y1  = y_bottom,
           lwd = schematic_lwd,
           xpd = NA
           )

  par(old_par)
  layout(1)
  return(invisible(NULL))
}



# Add comparison data from the original plasmids --------------------------

sg_pairs <- names(error_mat_list)[2:7]

full_reads_deletions_df_list <- lapply(full_reads_deletions_df_list, function(x) {
  x[, "Fraction_reference"] <- error_mat_list[["Full"]][, "Fraction_deleted"]
  x
})

full_reads_errors_df_list <- lapply(full_reads_errors_df_list, function(x) {
  x[, "Fraction_reference"] <- error_mat_list[["Full"]][, "Fraction_incorrect"]
  x
})

all_reads_deletions_df_list <- lapply(1:6, function(x) {
  use_df <- all_reads_deletions_df_list[[x]]
  pair_name <- sg_pairs[[x]]
  print(pair_name)
  use_df[, "Fraction_reference"] <- error_mat_list[[pair_name]][, "Fraction_deleted"]
  use_df
})
names(all_reads_deletions_df_list) <- sg_pairs

all_reads_errors_df_list <- lapply(1:6, function(x) {
  use_df <- all_reads_errors_df_list[[x]]
  pair_name <- sg_pairs[[x]]
  use_df[, "Fraction_reference"] <- error_mat_list[[pair_name]][, "Fraction_incorrect"]
  use_df
})
names(all_reads_errors_df_list) <- sg_pairs



# Display error rates -----------------------------------------------------

for (create_PDF in c(FALSE, TRUE)) {

  if (create_PDF) {
    pdf(file.path(output_dir, "Base-level error rate.pdf"),
        width = 8, height = 10
        )
  }

  SixErrorRates(full_reads_errors_df_list,    show_title = "Error rates \u2013 fully mapped reads only")
  SixErrorRates(full_reads_deletions_df_list, show_title = "Deletion rates \u2013 fully mapped reads only")
  SixErrorRates(all_reads_errors_df_list,     show_title = "Mutation rates \u2013 all reads", y_upper = 0.12)
  SixErrorRates(all_reads_deletions_df_list,  show_title = "Deletion rates \u2013 all reads", y_upper = 0.12)

  if (create_PDF) {
    dev.off()
  }

}


