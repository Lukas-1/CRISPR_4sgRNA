### 15th May 2021 ###



# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")




# General functions -------------------------------------------------------

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

  if (!(is_binary) && (!(all(is.na(numeric_vec))))) {
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
                              indicate_homologies     = FALSE,
                              number_wells            = FALSE,
                              gap_fraction            = 0.25,
                              show_low_read_numbers   = FALSE,
                              outline_few_reads       = FALSE,
                              few_reads_cutoff        = 20L,
                              label_rows_with_letters = TRUE,
                              label_cex               = 0.5,
                              label_space_adjust      = 0,
                              label_font              = 2
                              ) {

  MakeEmptyPlot()

  numeric_vec <- summary_df[, show_column]
  if (show_column == "Longest_subsequence") {
    numeric_vec <- numeric_vec / 20
  } else if (grepl("^(Count|Num)_", show_column)) {
    numeric_vec <- numeric_vec / summary_df[["Count_total"]]
  } else if (!(all(is.na(numeric_vec)))) { # Avoid a warning message
    numeric_vec <- numeric_vec / max(numeric_vec, na.rm = TRUE)
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

  text(x      = -(x_gap + diff(grconvertX(c(0, label_space_adjust), from = "lines", to = "user"))),
       y      = y_starts[seq_len(16)] + (well_height / 2),
       labels = rev(if (label_rows_with_letters) LETTERS[1:16] else 1:16),
       cex    = label_cex,
       font   = label_font,
       col    = "gray30",
       adj    = c(1, 0.5),
       xpd    = NA

       )
  text(x      = x_starts[seq_len(24)] + (well_width / 2),
       y      = -(y_gap + diff(grconvertY(c(0, label_space_adjust), from = "lines", to = "user"))),
       labels = seq_len(24),
       cex    = label_cex,
       font   = label_font,
       col    = "gray30",
       adj    = c(0.5, 1),
       xpd    = NA
       )

  if (show_low_read_numbers && !(outline_few_reads)) {
    low_read_number_colors <- ifelse(summary_df[["Count_total"]] < few_reads_cutoff,
                                     "#FFDDCC",
                                     "#FFFFFF"
                                     )
  } else {
    low_read_number_colors <- "#FFFFFF"
  }

  empty_colors <- ifelse(sg_df[["Empty_well"]],
                         brewer.pal(9, "Pastel1")[[6]],
                         low_read_number_colors
                         )

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




CenterText <- function(show_text, y_position = 0.3, use_cex = 1, text_color = "black") {
  text(x      = 0.5,
       y      = y_position,
       labels = show_text,
       col    = text_color,
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
                         top_text_color        = "black",
                         indicate_homologies   = FALSE,
                         show_low_read_numbers = FALSE,
                         outline_few_reads     = FALSE,
                         few_reads_cutoff      = 20L,
                         well_gap              = if (indicate_homologies || outline_few_reads) 0.4 else 0.25,
                         binary_cutoff         = NULL
                         ) {

  stopifnot(all(c("titles_list", "barcode_combos_df") %in% ls(envir = globalenv())))
  stopifnot(identical(summary_df[, "Well_number"], sg_df[, "Well_number"]))

  use_titles_list <- titles_list
  names(use_titles_list) <- toupper(names(use_titles_list))
  use_title <- use_titles_list[[toupper(show_column)]]

  if (is.null(binary_cutoff)) {
    if (grepl("^[0-9]{1,3}%", use_title)) {
      percentage_used <- strsplit(use_title, "%", fixed = TRUE)[[1]][[1]]
      binary_cutoff <- as.numeric(percentage_used) / 100
    } else {
      binary_cutoff <- 0.5
    }
  }

  if (show_column == "Longest_subsequence") {
    summary_df[["Longest_subsequence"]] <- sg_df[["Longest_subsequence"]]
  } else if (grepl("[Cc]ount_mean_sg1to4$", show_column)) {
    count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
    summary_df[["Count_mean_sg1to4"]] <- rowMeans(as.matrix(summary_df[, count_columns]))
  } else if (grepl("[Cc]ount_mean_pr_sg1to4$", show_column)) {
    count_columns <- paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4)
    summary_df[["Count_mean_pr_sg1to4"]] <- rowMeans(as.matrix(summary_df[, count_columns]))
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
  CenterText(as.expression(bquote(bold(.(use_title)))),
             y_position = 0.35,
             use_cex = 1.1
             )

  ## Draw the well layout
  for (i in 1:5) {
    MakeEmptyPlot()
  }

  text_y <- 0.4
  CenterText(top_text,
             y_position = text_y * (1 / space_1),
             text_color = top_text_color
             )

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




# Functions used for the first two sequencing runs ------------------------

DrawAllSchematicsForOnePlate <- function() {

  show_metrics <- setdiff(names(titles_list), "Num_cross_plate_contaminated")

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
          for (show_metric in show_metrics) {
            BarPlotPanel(ccs3_df_list[[df_name]],
                         show_metric,
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
          for (show_metric in show_metrics) {
            BarPlotPanel(ccs5_df_list[[df_name]],
                         show_metric,
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
          for (show_metric in show_metrics) {
            BarPlotPanel(ccs5_df_list[[df_name]],
                         show_metric,
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




# Functions used for multi-plate experiments ------------------------------

DrawSchematicsForAllPlates <- function(export_PNGs = TRUE) {

  required_objects <- c("use_plate_numbers", "plates_df", "sg_sequences_df")
  stopifnot(all(required_objects %in% ls(envir = globalenv())))

  if (export_PNGs) {
    file_formats <- c("png", "pdf")
  } else {
    file_formats <- "pdf"
  }

  ccs_numbers <- c(3, 5, 7)
  accuracy_percentages <- c(99, 99.9, 99.99)
  df_list_names <- paste0("ccs", ccs_numbers, "_df_list")
  ccs_are_present <- df_list_names %in% ls(envir = globalenv())
  ccs_numbers <- ccs_numbers[ccs_are_present]
  accuracy_percentages <- accuracy_percentages[ccs_are_present]

  if ("Highlight_color" %in% names(plates_df)) {
    top_text_colors <- plates_df[, "Highlight_color"]
  } else {
    top_text_colors <- rep("black", nrow(plates_df))
  }

  count_metrics <- c(
    # "Count_sg1_cr1", "Count_sg2_cr2", "Count_sg3_cr3", "Count_sg4_cr4",
    "Count_mean_sg1to4",
    "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
    "Count_pr_at_least_1", "Count_pr_at_least_2", "Count_pr_at_least_3", "Count_pr_all_4",
    "Count_all_4_promoters", "Count_whole_plasmid"
  )

  percentages_metrics <- c(
    "Num_contaminated_reads", "Num_under_2kb",
    count_metrics,
    "Num_reads_with_deletions_exceeding_20bp",
    "Num_reads_with_sgRNA_deletion",
    "Num_reads_with_deletions_spanning_tracrRNAs",
    "Num_reads_with_deletions_spanning_promoters",
    "Num_cross_plate_contaminated",
    "Num_contaminated_reads_aligned"
  )

  binarized_metrics <- c(
    "Binary_all_four_guides",
    paste0("Binary_", tolower(substr(count_metrics, 1, 1)),
           substr(count_metrics, 2, nchar(count_metrics))
           ),
    "Binary_num_contaminated_reads",
    "Binary_num_reads_with_deletions_exceeding_20bp",
    "Binary_num_reads_with_deletions_spanning_tracrRNAs"
  )


  plate_labels <- paste0("Plate #", plates_df[["Plate_number"]], " \u2013 ", plates_df[["Plate_name"]])
  plate_number_width <- max(nchar(as.character(plates_df[["Plate_number"]])))

  for (file_format in file_formats) {
    for (i in seq_along(ccs_numbers)) {
      use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))

      filter_stages <- c("original_summary_df", "filtered_summary_df", "filtered_cross_plate_df")
      filter_labels <- c("i) unfiltered", "ii) filtered", "iii) cross-plate")
      df_are_present <- filter_stages %in% names(use_df_list)
      filter_stages <- filter_stages[df_are_present]
      filter_labels <- filter_labels[df_are_present]

      for (filter_stage in seq_along(filter_stages)) {
        df_name <- filter_stages[[filter_stage]] # "filtered_gRNAs_df"

        use_summary_df <- use_df_list[[df_name]]

        folder_name <- paste0("CCS", ccs_numbers[[i]],
                              " (", accuracy_percentages[[i]], ") - ",
                              filter_labels[[filter_stage]]
                              )
        if (file_format == "png") {
          folder_path <- file.path(PNGs_output_directory, folder_name)
        } else {
          folder_path <- file.path(plots_output_directory, folder_name)
        }

        dir.create(folder_path, showWarnings = FALSE)

        message(paste0("Exporting ", file_format, " images into the folder: ", folder_name, "..."))

        file_prefix <- folder_name

        for (number_wells in c(TRUE, FALSE)) {
          for (binarize in c(TRUE, FALSE)) {

            if (binarize) {
              current_metrics <- binarized_metrics
            } else {
              current_metrics <- percentages_metrics
            }

            for (j in seq_along(current_metrics)) {

              if (binarize) {
                sub_folder_path <- file.path(folder_path, "Binarized - ")
              } else {
                sub_folder_path <- file.path(folder_path, "Percentages - ")
              }
              if (number_wells) {
                sub_folder_path <- paste0(sub_folder_path, "numbered")
              } else {
                sub_folder_path <- paste0(sub_folder_path, "plain")
              }
              dir.create(sub_folder_path, showWarnings = FALSE)

              current_metric_name <- sub("^Num_reads_with_", "", current_metrics[[j]])
              current_metric_name <- sub("^(Num|Count)_", "", current_metric_name)
              current_metric_name <- sub("^Binary_num_reads_with_", "Binary_", current_metric_name)
              current_metric_name <- sub("^Binary_(num|count)_", "Binary_", current_metric_name)
              current_metric_name <- paste0(formatC(j, width = 2, flag = "0"),
                                            ") ",  current_metric_name
                                            )

              expansion_factor <- 5
              use_width <- 2 * expansion_factor
              use_height <- 1 * expansion_factor

              if (file_format == "pdf") {
                pdf_name <- paste0(file_prefix, " - ", current_metric_name, ".pdf")
                pdf(file   = file.path(sub_folder_path, pdf_name),
                    width  = use_width,
                    height = use_height
                    )
              } else if (file_format == "png") {
                metric_path <- file.path(sub_folder_path, current_metric_name)
                dir.create(metric_path, showWarnings = FALSE)
              }

              for (plate_number in use_plate_numbers) {
                if (file_format == "png") {
                  plate_name <- paste0(formatC(plate_number, width = plate_number_width, flag = "0"), ") - ",
                                       plates_df[["Plate_name"]][plates_df[["Plate_number"]] == plate_number]
                                       )
                  file_name <- paste0(file_prefix, " - ", plate_name, ".png")
                  png(filename = file.path(metric_path, file_name),
                      width    = use_width,
                      height   = use_height,
                      units    = "in",
                      res      = 600
                      )
                }
                summary_sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
                sg_sub_df <- sg_sequences_df[sg_sequences_df[["Plate_number"]] %in% plate_number, ]
                plate_index <- which(plates_df[["Plate_number"]] == plate_number)
                BarPlotPanel(summary_sub_df,
                             current_metrics[[j]],
                             sg_sub_df,
                             number_wells          = number_wells,
                             top_text              = plate_labels[[plate_index]],
                             top_text_color        = top_text_colors[[plate_index]],
                             show_low_read_numbers = TRUE,
                             outline_few_reads     = binarize
                             )
                if (file_format == "png") {
                  dev.off()
                }
              }
              if (file_format == "pdf") {
                dev.off()
              }
            }
          }
        }
      }
    }
  }
  return(invisible(NULL))
}





