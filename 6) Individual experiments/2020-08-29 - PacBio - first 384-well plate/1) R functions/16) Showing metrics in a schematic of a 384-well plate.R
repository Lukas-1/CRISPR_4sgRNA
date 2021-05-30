### 15th May 2021 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")




# Define general functions ------------------------------------------------

MakeEmptyPlot <- function(x_limits = c(0, 1), y_limits = c(0, 1)) {
  plot(1, type = "n", axes = FALSE, ann = FALSE,
       xlim = x_limits, ylim = y_limits, xaxs = "i", yaxs = "i"
       )
}


BeeBarSubPlot <- function(summary_df,
                          show_column,
                          groups_column,
                          color_scheme = "Blues"
                          ) {

  is_binary <- grepl("^[Bb]inary_", show_column)

  if (!(is_binary)) {
    set.seed(1)
  }

  numeric_vec <- summary_df[[show_column]]

  is_percentage <- grepl("^(Count|Num)_", show_column)
  if (show_column == "Longest_subsequence") {
    y_limits <- c(0, 20)
  } else if (is_percentage || is_binary) {
    y_limits <- c(0, 1)
    if (is_percentage) {
      numeric_vec <- numeric_vec / summary_df[["Count_total"]]
    }
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
  if (is_binary || is_percentage) {
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

  bar_heights <- vapply(split_list, mean, numeric(1), na.rm = TRUE)

  for (i in seq_len(num_groups)) {
    rect(xleft   = x_positions - 0.35,
         xright  = x_positions + 0.35,
         ybottom = 0,
         ytop    = bar_heights,
         col     = brewer.pal(9, color_scheme)[[5]],
         border  = NA
         )
  }

  if (!(is_binary)) {
    beeswarm_df <- beeswarm(split_list,
                            priority = "random",
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
  }

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
                              sg_df,
                              indicate_homologies   = FALSE,
                              number_wells          = FALSE,
                              gap_fraction          = 0.25,
                              show_low_read_numbers = FALSE,
                              outline_few_reads     = FALSE,
                              few_reads_cutoff      = 20L
                              ) {


  assign("delete_summary_df",  summary_df,  envir = globalenv())
  assign("delete_show_column", show_column, envir = globalenv())
  assign("delete_sg_df",       sg_df,       envir = globalenv())

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

  if (show_low_read_numbers && !(outline_few_reads)) {
    low_read_number_colors <- ifelse(summary_df[["Count_total"]] < few_reads_cutoff,
                                     "#FFDDCC",
                                     "#FFFFFF"
                                     )
  } else {
    low_read_number_colors <- "#FFFFFF"
  }
  assign("delete_summary_df", summary_df, envir = globalenv())

  empty_colors <- ifelse(sg_df[["Empty_well"]],
                         brewer.pal(9, "Pastel1")[[6]],
                         low_read_number_colors
                         )
  assign("delete_empty_colors", empty_colors, envir = globalenv())


  draw_outlines <- indicate_homologies || outline_few_reads

  if (draw_outlines) {
    outline_color <- brewer.pal(9, "Reds")[[7]]
    no_outline_color <- brewer.pal(9, "Greys")[[7]]
    if (indicate_homologies && ("Longest_subsequence" %in% names(sg_df))) {
      are_to_outline <- sg_df[["Longest_subsequence"]] >= 8
    } else if (outline_few_reads) {
      are_to_outline <- summary_df[["Count_total"]] < few_reads_cutoff
    } else {
      are_to_outline <- rep(FALSE, nrow(sg_df))
    }
    outline_colors_vec <- ifelse(are_to_outline, outline_color, no_outline_color)

    if (indicate_homologies) {

      x_coord <- 0.3
      y_coord <- -0.09
      text_y <- y_coord + (well_height * 0.25)

      rect(xleft   = x_coord,
           xright  = x_coord + (well_width / 2),
           ybottom = y_coord,
           ytop    = y_coord + (well_height / 2),
           border  = no_outline_color,
           lwd     = 1,
           col     = NA,
           xpd     = NA
           )

      rect(xleft   = x_coord + 0.1,
           xright  = x_coord + 0.1 + (well_width / 2),
           ybottom = y_coord,
           ytop    = y_coord + (well_height / 2),
           border  = outline_color,
           lwd     = 1,
           col     = NA,
           xpd     = NA
           )

      text(x      = x_coord - 0.22,
           y      = text_y,
           adj    = c(0, 0.5),
           labels = "Longest subsequence:",
           cex    = 0.9,
           xpd    = NA
           )

      text(x      = x_coord + 0.018,
           y      = text_y,
           adj    = c(0, 0.5),
           labels = expression("" <= "7 bp"),
           col    = no_outline_color,
           cex    = 0.9,
           xpd    = NA
           )

      text(x      = x_coord + 0.118,
           y      = text_y,
           adj    = c(0, 0.5),
           labels = expression("" > "8 bp"),
           col    = outline_color,
           cex    = 0.9,
           xpd    = NA
           )
    } else if (outline_few_reads) {

      x_coord <- 0.3
      y_coord <- -0.09
      text_y <- y_coord + (well_height * 0.25)

      rect(xleft   = x_coord,
           xright  = x_coord + (well_width / 2),
           ybottom = y_coord,
           ytop    = y_coord + (well_height / 2),
           border  = outline_color,
           lwd     = 1,
           col     = NA,
           xpd     = NA
           )

      rect(xleft   = x_coord + 0.14,
           xright  = x_coord + 0.14 + (well_width / 2),
           ybottom = y_coord,
           ytop    = y_coord + (well_height / 2),
           border  = no_outline_color,
           lwd     = 1,
           col     = NA,
           xpd     = NA
           )

      text(x      = x_coord - 0.22,
           y      = text_y,
           adj    = c(0, 0.5),
           labels = "Total number of reads:",
           cex    = 0.9,
           xpd    = NA
           )

      text(x      = x_coord + 0.018,
           y      = text_y,
           adj    = c(0, 0.5),
           labels = expression("" < "20 reads"),
           col    = outline_color,
           cex    = 0.9,
           xpd    = NA
           )

      text(x      = x_coord + 0.158,
           y      = text_y,
           adj    = c(0, 0.5),
           labels = expression("" >= "20 reads"),
           col    = no_outline_color,
           cex    = 0.9,
           xpd    = NA
           )

    }
  }



  bar_heights <- ifelse(is.na(numeric_vec), 0, numeric_vec) * well_height

  for (i in seq_len(nrow(summary_df))) {

    well_number <- summary_df[i, "Well_number"]

    are_this_well_mat <- well_number_mat == well_number
    x_coord <- x_starts[col_coords_mat[are_this_well_mat]]
    y_coord <- y_starts[row_coords_mat[are_this_well_mat]]

    border_fraction <- 0.1
    if (draw_outlines) {
      rect(xleft   = x_coord - (well_width * border_fraction),
           xright  = x_coord + well_width + (well_width * border_fraction),
           ybottom = y_coord - (well_height  * border_fraction),
           ytop    = y_coord + well_height + (well_height * border_fraction),
           border  = NA,
           col     = outline_colors_vec[[i]],
           )
    }

    bar_height <- bar_heights[[i]]
    if (bar_height > 0) { # This prevents zero-height rectangles being shown as a thin line on some zoom levels in the PDF viewer
      rect(xleft   = x_coord,
           xright  = x_coord + well_width,
           ybottom = y_coord,
           ytop    = y_coord + bar_height,
           border  = NA,
           col     = brewer.pal(9, "Blues")[[9]]
           )
    }

    if (bar_height < well_height) {
      rect(xleft   = x_coord,
           xright  = x_coord + well_width,
           ybottom = y_coord + bar_height,
           ytop    = y_coord + well_height,
           border  = NA,
           col     = empty_colors[[i]]
           )
    }

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
                         sg_df,
                         width_to_height_ratio = 2,
                         number_wells          = FALSE,
                         top_text              = "Original 384-well plate",
                         indicate_homologies   = FALSE,
                         show_low_read_numbers = FALSE,
                         outline_few_reads     = FALSE,
                         few_reads_cutoff      = 20L,
                         well_gap              = if (indicate_homologies || outline_few_reads) 0.4 else 0.25,
                         binary_cutoff         = 0.5
                         ) {

  stopifnot(all(c("titles_list", "barcode_combos_df") %in% ls(envir = globalenv())))
  stopifnot(identical(summary_df[, "Well_number"], sg_df[, "Well_number"]))

  if (show_column == "Longest_subsequence") {
    summary_df[["Longest_subsequence"]] <- sg_df[["Longest_subsequence"]]
  } else if (grepl("[Cc]ount_mean_sg1to4$", show_column)) {
    count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
    summary_df[["Count_mean_sg1to4"]] <- rowMeans(as.matrix(summary_df[, count_columns]))
  }

  if (grepl("^[Bb]inary_", show_column)) {
    if (show_column == "Binary_all_four_guides") {
      use_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
      fraction_correct_mat <- as.matrix(summary_df[, use_columns]) / summary_df[["Count_total"]]
      binary_vec <- rowSums(fraction_correct_mat >= binary_cutoff) == 4
    } else {
      original_column <- sub("^[Bb]inary_", "", show_column)
      original_column <- paste0(toupper(substr(original_column, 1, 1)),
                                substr(original_column, 2, nchar(original_column))
                                )
      stopifnot(grepl("^(Count|Num)_", original_column))
      binary_vec <- (summary_df[, original_column] / summary_df[, "Count_total"])>= binary_cutoff
    }
    summary_df[[show_column]] <- as.integer(binary_vec)
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
  use_titles_list <- titles_list
  names(use_titles_list) <- toupper(names(use_titles_list))
  use_title <- use_titles_list[[toupper(show_column)]]
  CenterText(as.expression(bquote(bold(.(use_title)))),
             y_position = 0.35,
             use_cex = 1.1
             )

  ## Draw the well layout
  for (i in 1:5) {
    MakeEmptyPlot()
  }

  text_y <- 0.4

  CenterText(top_text, y_position = text_y * (1 / space_1))
  WellLayoutBarplot(summary_df,
                    show_column,
                    sg_df,
                    gap_fraction          = well_gap,
                    indicate_homologies   = indicate_homologies,
                    number_wells          = number_wells,
                    show_low_read_numbers = show_low_read_numbers,
                    outline_few_reads     = outline_few_reads,
                    few_reads_cutoff      = few_reads_cutoff
                    )

  ## Draw the bar plots (with superimposed swarm plots)

  summary_df[["Row_groups"]]    <- factor(barcode_combos_df[["Row_number"]][summary_df[["Well_number"]]],
                                          levels = unique(barcode_combos_df[["Row_number"]])
                                          )
  summary_df[["Column_groups"]] <- factor(barcode_combos_df[["Column_number"]][summary_df[["Well_number"]]],
                                          levels = unique(barcode_combos_df[["Column_number"]])
                                          )

  summary_df <- summary_df[!(sg_df[["Empty_well"]]), ]

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




# Define functions used for the first two sequencing runs -----------------

DrawAllSchematicsForOnePlate <- function() {

  for (filter_category in c("original", "filtered reads", "filtered gRNAs")) {
    for (smrtlink_version in c(7, 9)) {
      for (highlight_homologies in c(FALSE, TRUE)) {
        for (show_well_numbers in c(FALSE, TRUE)) {

          if (show_well_numbers && !(highlight_homologies)) {
            next
          }

          file_prefix <- paste0("Schematics of a 384-well plate - SmrtLink ", smrtlink_version, " - ")

          if (filter_category == "original") {
            file_postfix <- " - original"
            df_name <- "original_summary_df"
          } else if (filter_category == "filtered reads") {
            file_postfix <- " - filtered"
            df_name <- "filtered_summary_df"
          } else if (filter_category == "filtered gRNAs") {
            file_postfix <- " - filtered gRNAs"
            df_name <- "filtered_gRNAs_df"
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
            ccs7_df_list <- sl7_ccs7_df_list
          } else if (smrtlink_version == 9) {
            ccs3_df_list <- sl9_ccs3_df_list
            ccs5_df_list <- sl9_ccs5_df_list
            ccs5_df_list <- sl9_ccs7_df_list
          }

          # Export plots ------------------------------------------------------------

          expansion_factor <- 5

          pdf(file = file.path(plots_output_directory,
                               paste0("SmrtLink ", smrtlink_version),
                               "Schematics of a 384-well plate",
                               paste0(file_prefix, "CCS3", file_postfix)
                               ),
              width  = 2 * expansion_factor,
              height = 1 * expansion_factor
              )
          for (i in seq_along(titles_list)) {
            BarPlotPanel(ccs3_df_list[[df_name]],
                         names(titles_list)[[i]],
                         sg_sequences_df,
                         indicate_homologies = highlight_homologies,
                         number_wells = show_well_numbers
                         )
          }
          dev.off()


          pdf(file = file.path(plots_output_directory,
                               paste0("SmrtLink ", smrtlink_version),
                               "Schematics of a 384-well plate",
                               paste0(file_prefix, "CCS5", file_postfix)
                               ),
              width  = 2 * expansion_factor,
              height = 1 * expansion_factor
              )
          for (i in seq_along(titles_list)) {
            BarPlotPanel(ccs5_df_list[[df_name]],
                         names(titles_list)[[i]],
                         sg_sequences_df,
                         indicate_homologies = highlight_homologies,
                         number_wells = show_well_numbers
                         )
          }
          dev.off()


          pdf(file = file.path(plots_output_directory,
                               paste0("SmrtLink ", smrtlink_version),
                               "Schematics of a 384-well plate",
                               paste0(file_prefix, "CCS7", file_postfix)
                               ),
              width  = 2 * expansion_factor,
              height = 1 * expansion_factor
              )
          for (i in seq_along(titles_list)) {
            BarPlotPanel(ccs5_df_list[[df_name]],
                         names(titles_list)[[i]],
                         sg_sequences_df,
                         indicate_homologies = highlight_homologies,
                         number_wells = show_well_numbers
                         )
          }
          dev.off()

        }
      }
    }
  }
}





