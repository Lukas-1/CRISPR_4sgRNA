### 19th September 2020 ###



# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))



# Define folder paths -----------------------------------------------------

file_input_directory   <- file.path(file_directory, "2) Input")
file_output_directory  <- file.path(file_directory, "5) Output")
R_objects_directory    <- file.path(file_directory, "3) R objects")

plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "10) Process demultiplexed PacBio reads.RData"))






# Define functions --------------------------------------------------------

MakeEmptyPlot <- function(x_limits = c(0, 1), y_limits = c(0, 1)) {
  plot(1, type = "n", axes = FALSE, ann = FALSE,
       xlim = x_limits, ylim = y_limits, xaxs = "i", yaxs = "i"
       )
}


# BarSubPlot <- function(numeric_vec, show_column, y_limits) {
#
#   is_percentage <- grepl("^(Count|Num)_", show_column)
#
#   MakeEmptyPlot(y_limits = y_limits)
#
#   num_bars <- length(numeric_vec)
#
#   tick_locations <- axTicks(2)
#   if (is_percentage) {
#     tick_labels <- paste0(tick_locations * 100, "%")
#   } else {
#     tick_labels <- TRUE
#   }
#
#   abline(h = tick_locations, col = "gray94")
#   if (is_percentage) {
#     abline(h = seq(0.1, 0.9, 0.2), col = "gray98")
#   }
#
#   for (i in seq_len(num_bars)) {
#     rect(xleft   = (i - 1) / num_bars,
#          xright  = (i / num_bars),
#          ybottom = 0,
#          ytop    = numeric_vec[[i]],
#          col     = brewer.pal(9, "Blues")[[9]],
#          border  = NA
#          )
#   }
#
#
#   axis(2,
#        labels   = tick_labels,
#        at       = tick_locations,
#        las      = 1,
#        mgp      = c(3, 0.45, 0),
#        tcl      = -0.35,
#        lwd      = 0.75,
#        cex.axis = 0.9
#        )
#   box(lwd = 0.75)
# }




BeeBarSubPlot <- function(summary_df,
                          show_column,
                          groups_column,
                          color_scheme = "Blues"
                          ) {

  numeric_vec <- summary_df[[show_column]]

  is_percentage <- grepl("^(Count|Num)_", show_column)

  if (show_column == "Longest_subsequence") {
    y_limits <- c(0, 20)
  } else if (is_percentage) {
    numeric_vec <- numeric_vec / summary_df[["Count_total"]]
    y_limits <- c(0, 1)
  } else {
    y_limits <- c(0, max(numeric_vec) * 1.02)
  }

  split_list <- split(numeric_vec, summary_df[[groups_column]])

  num_groups <- length(split_list)
  side_margin <- 0.3
  margin_increment <- 0.03

  x_limits <- c((1 - side_margin) - (margin_increment * num_groups),
                num_groups + side_margin + (margin_increment * num_groups)
                )

  x_positions <- seq_len(num_groups)

  MakeEmptyPlot(x_limits = x_limits, y_limits = y_limits)

  tick_locations <- axTicks(2)
  if (is_percentage) {
    tick_labels <- paste0(tick_locations * 100, "%")
  } else {
    tick_labels <- TRUE
  }

  if (is_percentage) {
    abline(h = seq(0.1, 0.9, 0.2), col = "gray98", lwd = 0.75)
  }
  abline(v = c(x_positions, x_limits),
         col = "gray98", lwd = 0.75
         )
  abline(h = tick_locations, col = "gray94", lwd = 0.75)
  # box(lwd = 0.75, col = "gray94")

  bar_heights <- vapply(split_list, mean, numeric(1))

  for (i in seq_len(num_groups)) {
    rect(xleft   = x_positions - 0.35,
         xright  = x_positions + 0.35,
         ybottom = 0,
         ytop    = bar_heights,
         col     = brewer.pal(9, color_scheme)[[5]],
         border  = NA
         )
  }

  beeswarm_df <- beeswarm(split_list,
                          priority = "density",
                          cex      = 0.3,
                          spacing  = 0.8,
                          do.plot  = FALSE
                          )
  points(beeswarm_df[["x"]],
         beeswarm_df[["y"]],
         pch = 16,
         col = brewer.pal(9, color_scheme)[[8]],
         xpd = NA,
         cex = 0.3
         )

  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.35,
       lwd      = 0.75,
       cex.axis = 0.9
       )

  text(x      = x_positions,
       y      = par("usr")[[3]] - ((par("usr")[[4]] - par("usr")[[3]]) * 0.065),
       labels = seq_len(num_groups),
       xpd    = NA,
       cex    = 0.8,
       col    = "gray30"
       )

  return(invisible(NULL))
}





WellLayoutBarplot <- function(summary_df,
                              show_column,
                              indicate_homologies = FALSE,
                              number_wells = FALSE,
                              gap_fraction = 0.25
                              ) {

  MakeEmptyPlot()

  numeric_vec <- summary_df[[show_column]]
  if (show_column == "Longest_subsequence") {
    numeric_vec <- numeric_vec / 20
  } else if (grepl("^(Count|Num)_", show_column)) {
    numeric_vec <- numeric_vec / summary_df[["Count_total"]]
  } else {
    numeric_vec <- numeric_vec / max(numeric_vec)
  }

  row_coords_mat <- matrix(rep(rev(seq_len(16)), each = 24),
                           nrow  = 16,
                           ncol  = 24,
                           byrow = TRUE
                           )
  col_coords_mat <- matrix(rep(seq_len(24), times = 16),
                           nrow  = 16,
                           ncol  = 24,
                           byrow = TRUE
                           )

  well_number_mat <- matrix(seq_len(384),
                            nrow  = 16,
                            ncol  = 24,
                            byrow = TRUE
                            )

  well_width  <- 1 / ((25 * gap_fraction) + 24)
  well_height <- 1 / ((17 * gap_fraction) + 16)

  x_gap <- well_width * gap_fraction
  y_gap <- well_height * gap_fraction

  start_x <- x_gap
  end_x   <- 1
  start_y <- y_gap
  end_y   <- 1

  x_starts <- seq(start_x, end_x, by = well_width + x_gap)
  y_starts <- seq(start_y, end_y, by = well_height + y_gap)

  rect(xleft   = 0,
       xright  = 1,
       ybottom = 0,
       ytop    = 1,
       col     = "gray80",
       border  = NA,
       xpd     = NA
       )

  text(x      = -(x_gap),
       y      = y_starts[seq_len(16)] + (well_height / 2),
       labels = rev(seq_len(16)),
       xpd    = NA,
       cex    = 0.5,
       font   = 2,
       col    = "gray30",
       adj    = c(1, 0.5)
       )
  text(x      = x_starts[seq_len(24)] + (well_width / 2),
       y      = -(y_gap),
       labels = seq_len(24),
       xpd    = NA,
       cex    = 0.5,
       font   = 2,
       col    = "gray30",
       adj    = c(0.5, 1)
       )

  empty_colors <- ifelse(sg_sequences_df[["Empty_well"]],
                         brewer.pal(9, "Pastel1")[[6]],
                         "#FFFFFF"
                         )


  homology_color <- brewer.pal(9, "Reds")[[7]]
  no_homology_color <- brewer.pal(9, "Greys")[[7]]
  homology_colors_vec <- ifelse(sg_sequences_df[["Longest_subsequence"]] >= 8,
                                homology_color,
                                no_homology_color
                                )


  if (indicate_homologies) {

    x_coord <- 0.3
    y_coord <- -0.09
    text_y <- y_coord + (well_height * 0.25)

    text(x      = x_coord - 0.22,
         y      = text_y,
         adj    = c(0, 0.5),
         labels = "Longest subsequence:",
         cex    = 0.9,
         xpd    = NA
         )

    rect(xleft   = x_coord,
         xright  = x_coord + (well_width / 2),
         ybottom = y_coord,
         ytop    = y_coord + (well_height / 2),
         border  = no_homology_color,
         lwd     = 1,
         col     = NA,
         xpd     = NA
         )

    text(x      = x_coord + 0.018,
         y      = text_y,
         adj    = c(0, 0.5),
         labels = expression("" <= "7 bp"),
         col    = no_homology_color,
         cex    = 0.9,
         xpd    = NA
         )

    rect(xleft   = x_coord + 0.1,
         xright  = x_coord + 0.1 + (well_width / 2),
         ybottom = y_coord,
         ytop    = y_coord + (well_height / 2),
         border  = homology_color,
         lwd     = 1,
         col     = NA,
         xpd     = NA
         )

    text(x      = x_coord + 0.118,
         y      = text_y,
         adj    = c(0, 0.5),
         labels = expression("" > "8 bp"),
         col    = homology_color,
         cex    = 0.9,
         xpd    = NA
         )

  }

  for (well_number in seq_len(384)) {

    are_this_well_mat <- well_number_mat == well_number
    x_coord <- x_starts[col_coords_mat[are_this_well_mat]]
    y_coord <- y_starts[row_coords_mat[are_this_well_mat]]

    bar_height <- numeric_vec[[well_number]] * well_height

    border_fraction <- 0.1
    if (indicate_homologies) {
      rect(xleft   = x_coord - (well_width * border_fraction),
           xright  = x_coord + well_width + (well_width * border_fraction),
           ybottom = y_coord - (well_height  * border_fraction),
           ytop    = y_coord + well_height + (well_height * border_fraction),
           border  = NA,
           col     = homology_colors_vec[[well_number]],
           )
    }

    rect(xleft   = x_coord,
         xright  = x_coord + well_width,
         ybottom = y_coord,
         ytop    = y_coord + bar_height,
         border  = NA,
         col     = brewer.pal(9, "Blues")[[9]]
         )

    rect(xleft   = x_coord,
         xright  = x_coord + well_width,
         ybottom = y_coord + bar_height,
         ytop    = y_coord + well_height,
         border  = NA,
         col     = empty_colors[[well_number]]
         )

    if (number_wells) {
      text(x      = x_coord + (well_width / 2),
           y      = y_coord + (well_height / 2),
           labels = as.character(well_number),
           xpd    = NA,
           col    = "gray75",
           cex    = 0.4,
           font   = 2
           )
    }


  }

  return(invisible(NULL))
}




CenterText <- function(show_text, y_position = 0.3, use_cex = 1) {
  text(x      = 0.5,
       y      = y_position,
       labels = show_text,
       adj    = c(0.5, 0.5),
       cex    = use_cex,
       xpd    = NA
       )
}



BarPlotPanel <- function(summary_df,
                         show_column,
                         width_to_height_ratio = 2,
                         indicate_homologies = FALSE,
                         well_gap = if (indicate_homologies) 0.4 else 0.25,
                         number_wells = FALSE
                         ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  if (show_column == "Longest_subsequence") {
    summary_df[["Longest_subsequence"]] <- sg_sequences_df[["Longest_subsequence"]]
  } else if (show_column == "Count_mean_sg1to4") {
    count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
    summary_df[["Count_mean_sg1to4"]] <- rowMeans(as.matrix(summary_df[, count_columns]))
  }

  ## Set up the layout

  old_mar <- par("mar" = rep(0, 4))

  inner_space <- 1
  space_height <- 1
  golden_ratio <- (1 + sqrt(5)) / 2
  layout_mat <- cbind(rep(3, 6),
                      c(6, 7, 7, 7, 7, 7),
                      rep(4, 6),
                      8:13,
                      rep(5, 6)
                      )
  layout_mat <- rbind(rep(1, 5),
                      layout_mat,
                      rep(2, 5)
                      )

  widths_vec <- c(1, 12, 1.5, 8, 1)


  space_1 <- 1.5
  space_2 <- 0.5
  space_3 <- 1.5

  heights_vec <- c(1,
                   space_1,
                   space_2,
                   3,
                   space_3,
                   3,
                   0.5,
                   0.8
                   )

  plate_prop_width <- (widths_vec[[2]] / sum(widths_vec)) * width_to_height_ratio

  plate_width_in_wells <- 24 + (well_gap * 25)
  plate_height_in_wells <- 16 + (well_gap * 17)

  plate_prop_height <- plate_prop_width * (plate_height_in_wells / plate_width_in_wells)
  plate_height <- sum(heights_vec[3:7])

  total_height <- plate_height / plate_prop_height

  inflexible_height <- sum(heights_vec[2:7])
  flexible_height <- total_height - inflexible_height

  remaining_props <- heights_vec[c(1, 8)] / sum(heights_vec[c(1, 8)])
  heights_vec[c(1, 8)] <- flexible_height * remaining_props

  layout(layout_mat, widths = widths_vec, heights = heights_vec)

  MakeEmptyPlot()
  CenterText(as.expression(bquote(bold(.(titles_list[[show_column]])))),
             y_position = 0.35,
             use_cex = 1.1
             )

  ## Draw the well layout
  for (i in 1:5) {
    MakeEmptyPlot()
  }

  text_y <- 0.4

  CenterText("Original 384-well plate", y_position = text_y * (1 / space_1))
  WellLayoutBarplot(summary_df,
                    show_column,
                    gap_fraction = well_gap,
                    indicate_homologies = indicate_homologies,
                    number_wells = number_wells
                    )

  ## Draw the barplots (with superimposed beeeswarms)

  summary_df[["Row_groups"]] <- barcode_combos_df[["Row_number"]]
  summary_df[["Column_groups"]] <- barcode_combos_df[["Column_number"]]

  summary_df <- summary_df[!(sg_sequences_df[["Empty_well"]]), ]

  MakeEmptyPlot()
  MakeEmptyPlot()
  CenterText("Summarized by row", y_position = text_y * (1 / space_2))
  BeeBarSubPlot(summary_df,
                show_column,
                groups_column = "Row_groups",
                color_scheme = "Greens"
                )

  MakeEmptyPlot()
  CenterText("Summarized by column", y_position = text_y * (1 / space_3))
  BeeBarSubPlot(summary_df,
                show_column,
                groups_column = "Column_groups",
                color_scheme = "Purples"
                )
  MakeEmptyPlot()

  ## Tidy up

  par(mar = old_mar)
  return(invisible(NULL))
}







# Start loop --------------------------------------------------------------

for (use_filtered in c(FALSE, TRUE)) {
  for (smrtlink_version in c(7, 9)) {
    for (highlight_homologies in c(FALSE, TRUE)) {
      for (show_well_numbers in c(FALSE, TRUE)) {

        if (show_well_numbers && !(highlight_homologies)) {
          next
        }

        file_prefix <- paste0("Quality control - SmrtLink ", smrtlink_version, " - ")

        if (use_filtered) {
          file_postfix <- " - filtered"
          df_name <- "filtered_summary_df"
        } else {
          file_postfix <- " - original"
          df_name <- "original_summary_df"
        }
        if (highlight_homologies) {
          file_postfix <- paste0(file_postfix, " - homologies indicated")
        }
        if (highlight_homologies) {
          file_postfix <- paste0(file_postfix, " - wells numbered")
        }
        file_postfix <- paste0(file_postfix, ".pdf")

        if (smrtlink_version == 7) {
          ccs3_df_list <- sl7_ccs3_df_list
          ccs5_df_list <- sl7_ccs5_df_list
        } else if (smrtlink_version == 9) {
          ccs3_df_list <- sl9_ccs3_df_list
          ccs5_df_list <- sl9_ccs5_df_list
        }


        # Export plots ------------------------------------------------------------

        expansion_factor <- 5

        pdf(file = file.path(plots_output_directory,
                             paste0("SmrtLink ", smrtlink_version),
                             "Quality control",
                             paste0(file_prefix, "CCS3", file_postfix)
                             ),
            width  = 2 * expansion_factor,
            height = 1 * expansion_factor
            )
        for (i in seq_along(titles_list)) {
          BarPlotPanel(ccs3_df_list[[df_name]],
                       names(titles_list)[[i]],
                       indicate_homologies = highlight_homologies,
                       number_wells = show_well_numbers
                       )
        }
        dev.off()


        pdf(file = file.path(plots_output_directory,
                             paste0("SmrtLink ", smrtlink_version),
                             "Quality control",
                             paste0(file_prefix, "CCS5", file_postfix)
                             ),
            width  = 2 * expansion_factor,
            height = 1 * expansion_factor
            )
        for (i in seq_along(titles_list)) {
          BarPlotPanel(ccs5_df_list[[df_name]],
                       names(titles_list)[[i]],
                       indicate_homologies = highlight_homologies,
                       number_wells = show_well_numbers
                       )
        }
        dev.off()




        # End loop ----------------------------------------------------------------

      }
    }
  }
}





BarPlotPanel(ccs3_df_list[["original_summary_df"]],
             "Longest_subsequence",
             indicate_homologies = TRUE
             )












