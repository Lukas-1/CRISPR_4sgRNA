### 1st August 2021 ###




# Import packages and source code -----------------------------------------

library("vioplot")
library("png")




# Define constants --------------------------------------------------------

use_width <- 6.9
use_height <- 4.7

column_labels_list <- list(
  "Count_at_least_1"                            = expression("" >= "1 correct", "gRNA"),
  "Count_at_least_2"                            = expression("" >= "2 correct", "gRNAs"),
  "Count_at_least_3"                            = expression("" >= "3 correct", "gRNAs"),
  "Count_all_4"                                 = expression("All 4 gRNAs", "are correct"),

  "Count_pr_at_least_1"                         = expression("" >= "1 correct", "gRNA"),
  "Count_pr_at_least_2"                         = expression("" >= "2 correct", "gRNAs"),
  "Count_pr_at_least_3"                         = expression("" >= "3 correct", "gRNAs"),
  "Count_pr_all_4"                              = expression("All 4 gRNAs", "are correct"),

  "Count_no_contam_at_least_1"                  = expression("" >= "1 correct", "gRNA"),
  "Count_no_contam_at_least_2"                  = expression("" >= "2 correct", "gRNAs"),
  "Count_no_contam_at_least_3"                  = expression("" >= "3 correct", "gRNAs"),
  "Count_no_contam_all_4"                       = expression("All 4 gRNAs", "are correct"),

  "Correct_sg1_cr1"                             = expression("sg1 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),
  "Correct_sg2_cr2"                             = expression("sg2 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),
  "Correct_sg3_cr3"                             = expression("sg3 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),
  "Correct_sg4_cr4"                             = expression("sg4 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),

  "Correct_sg1"                                 = expression("sg1 is", "correct"),
  "Correct_sg2"                                 = expression("sg2 is", "correct"),
  "Correct_sg3"                                 = expression("sg3 is", "correct"),
  "Correct_sg4"                                 = expression("sg4 is", "correct"),

  "Count_sg1_cr1"                               = expression("sg1 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),
  "Count_sg2_cr2"                               = expression("sg2 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),
  "Count_sg3_cr3"                               = expression("sg3 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),
  "Count_sg4_cr4"                               = expression("sg4 (+" * scriptscriptstyle(" ") * "cr)", "is correct"),

  "Count_pr1_sg1_cr1"                           = expression("sg1 (+" * scriptscriptstyle(" ") * "pr/cr)", "is correct"),
  "Count_pr2_sg2_cr2"                           = expression("sg2 (+" * scriptscriptstyle(" ") * "pr/cr)", "is correct"),
  "Count_pr3_sg3_cr3"                           = expression("sg3 (+" * scriptscriptstyle(" ") * "pr/cr)", "is correct"),
  "Count_pr4_sg4_cr4"                           = expression("sg4 (+" * scriptscriptstyle(" ") * "pr/cr)", "is correct"),

  "Num_contaminated_reads"                      = expression("Contamination", "(full read)"),
  "Num_contaminated_reads_aligned"              = expression("Contamination", "(aligned read)"),
  "Num_cross_plate_contaminated"                = expression("Cross-plate", "contamination"),

  "Num_reads_with_sgRNA_deletion"               = expression("Dels of", "" >= "1 gRNA"),
  "Num_reads_with_deletions_exceeding_20bp"     = expression("Deletions", "" >= "20bp"),
  "Num_reads_with_deletions_spanning_tracrRNAs" = expression("... span", "tracRNAs"),
  "Num_reads_with_deletions_spanning_promoters" = expression("... span", "promoters"),
  "Num_reads_with_deletions_spanning_sg_cr"     = expression("... span", "sg" * scriptstyle(" ") * "+" * scriptstyle(" ") * "cr"),
  "Num_reads_with_deletions_spanning_sgRNAs"    = expression("... span", "gRNAs")
)


column_groups_list <- list(
  "At_least_num_guides"  = c(paste0("Count_at_least_", 1:3), "Count_all_4"),
  "No_contam_num_guides" = paste0("Count_no_contam_", c(paste0("at_least_", 1:3), "all_4")),
  "Promoter_num_guides"  = c(paste0("Count_pr_at_least_", 1:3), "Count_pr_all_4"),
  "Count_sg_cr"          = paste0("Count_sg", 1:4, "_cr", 1:4),
  "Count_pr_sg_cr"       = paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4),
  "Deletions"            = c("Num_reads_with_deletions_exceeding_20bp",
                             "Num_reads_with_sgRNA_deletion",
                             "Num_reads_with_deletions_spanning_tracrRNAs",
                             "Num_reads_with_deletions_spanning_promoters"
                             ),
  "Contaminations"       = c("Num_contaminated_reads", "Num_contaminated_reads_aligned", "Num_cross_plate_contaminated")
)



column_group_subtitles_list <- list(
  "At_least_num_guides"  = expression("The gRNA" * scriptscriptstyle(" ") *
                                      "+" * scriptscriptstyle(" ") *
                                      "tracrRNA sequences must be 100% correct."
                                      ),
  "No_contam_num_guides" = expression("The gRNA" * scriptscriptstyle(" ") *
                                      "+" * scriptscriptstyle(" ") *
                                      "tracrRNA must be 100% correct. Contaminations are excluded."
                                      ),
  "Promoter_num_guides"  = expression("The gRNA" * scriptscriptstyle(" ") *
                                      "+" * scriptscriptstyle(" ") *
                                      "tracrRNA must be 100% correct, and the promoter 95% correct."
                                      )
)


column_group_subtitles_list <- c(
  column_group_subtitles_list,
  list(
    "Count_sg_cr"        = column_group_subtitles_list[["At_least_num_guides"]],
    "Count_pr_sg_cr"     = column_group_subtitles_list[["Promoter_num_guides"]]
  )
)


count_columns <- c(
  paste0("Count_at_least_", 1:3), "Count_all_4",
  "Count_all_4_promoters", "Count_whole_plasmid",
  paste0("Count_sg", 1:4, "_cr", 1:4)
)

no_contam_count_columns <- sub("Count_", "Count_no_contam_", count_columns, fixed = TRUE)

count_columns <- c(
  count_columns,
  paste0("Count_pr_at_least_", 1:3), "Count_pr_all_4",
  paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4)
)


num_columns <- c(
  "Num_contaminated_reads",
  "Num_contaminated_reads_aligned",
  "Num_cross_plate_contaminated",
  "Num_reads_with_sgRNA_deletion",
  "Num_reads_with_deletions_exceeding_20bp",
  "Num_reads_with_deletions_spanning_tracrRNAs",
  "Num_reads_with_deletions_spanning_promoters",
  "Num_reads_with_deletions_spanning_sg_cr",
  "Num_reads_with_deletions_spanning_sgRNAs"
)

alignment_columns <- c(
  "Correct_TpR_DHFR", "Mutation_TpR_DHFR", "Deletion_TpR_DHFR",
  paste0("Contamination_sg", 1:4),
  paste0("Correct_sg",       1:4),
  paste0("Deletion_sg",      1:4),
  paste0("Mutation_sg",      1:4),
  paste0("Correct_sg",       1:4, "_cr", 1:4),
  paste0("Deletion_sg",      1:4, "_cr", 1:4),
  paste0("Mutation_sg",      1:4, "_cr", 1:4),
  paste0("Contamination_sg", 1:4, "_cr", 1:4)
)





# Helper functions for creating plots -------------------------------------

GetPlateSelection <- function(plate_names) {
  require_objects <- c("plate_selection_list", "plate_selection_titles_list")
  stopifnot(all(require_objects %in% ls(envir = globalenv())))
  if (is.null(plate_names)) {
    plate_names <- "All plates"
  }
  if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
    plate_names <- plate_selection_list[[plate_names]]
  } else if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_titles_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
  } else {
    main_title <- "PacBio sequencing"
  }
  results_list <- list("plate_names" = plate_names,
                       "title" = main_title
                       )
  return(results_list)
}



DrawGridlines <- function(y_limits, extra_grid_lines = TRUE) {
  y_range <- y_limits[[2]] - y_limits[[1]]
  divide_by <- 20
  grid_seq <- seq(from = y_limits[[1]], to = y_limits[[2]], by = y_range / divide_by)
  are_main <- rep(c(TRUE, FALSE), divide_by)[seq_along(grid_seq)]
  segments(x0   = par("usr")[[1]],
           x1   = par("usr")[[2]],
           y0   = grid_seq[are_main],
           col  = "gray88",
           lend = "butt",
           xpd  = NA
           )
  if (extra_grid_lines) {
    segments(x0   = par("usr")[[1]],
             x1   = par("usr")[[2]],
             y0   = grid_seq[!(are_main)],
             col  = "gray95",
             lend = "butt",
             xpd  = NA
             )
  }
  return(invisible(NULL))
}






# Functions for aggregating and tidying data ------------------------------

StandardizeCounts <- function(input_df) {
  stopifnot(all(c("count_columns", "alignment_columns") %in% ls(envir = globalenv())))

  assign("delete_input_df", input_df, envir = globalenv())
  count_mat <- as.matrix(input_df[, count_columns])
  count_mat <- sweep(count_mat, 1, input_df[["Count_total"]], "/")

  no_contam_count_mat <- as.matrix(input_df[, no_contam_count_columns])
  no_contam_count_mat <- sweep(no_contam_count_mat, 1, input_df[["Count_total_no_contam"]], "/")

  remaining_mat <- as.matrix(input_df[, c(num_columns, alignment_columns)])
  remaining_mat <- sweep(remaining_mat, 1, input_df[["Count_total"]], "/")

  results_mat <- cbind(count_mat, no_contam_count_mat, remaining_mat)
  return(results_mat)
}



PreparePlates <- function(input_df, plate_names) {
  stopifnot("plates_df" %in% ls(envir = globalenv()))
  plate_matches <- match(plate_names, plates_df[, "Plate_name"])
  plate_numbers <- plates_df[["Plate_number"]][plate_matches]
  input_df <- input_df[input_df[["Plate_number"]] %in% plate_numbers, ]
  row.names(input_df) <- NULL
  results_mat <- StandardizeCounts(input_df)
  return(results_mat)
}





# Functions for plotting summary statistics -------------------------------

SetUpPercentagePlot <- function(use_columns,
                                use_y_limits,
                                main_title,
                                point_cex        = 1.5,
                                legend_pch       = 21,
                                side_gap         = 0.5,
                                draw_grid        = TRUE,
                                extra_grid_lines = FALSE,
                                for_embedded_PNG = FALSE,
                                make_plot        = TRUE,
                                add_mean         = FALSE,
                                y_axis_label     = "Percentage of reads",
                                sub_title        = NULL,
                                draw_legend      = TRUE
                                ) {

  stopifnot("column_labels_list" %in% ls(envir = globalenv()))

  ## Determine group positions
  num_groups <- length(use_columns)
  group_positions <- seq_len(num_groups)
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  numeric_axis_pos <- pretty(use_y_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  numeric_axis_labels <- paste0(format(numeric_axis_pos * 100), "%")


  if (make_plot) {
    final_y_range <- numeric_limits[[2]] - numeric_limits[[1]]
    y_gap <- final_y_range * 0.0075
    final_y_limits <- c(numeric_limits[[1]] - y_gap, numeric_limits[[2]] + y_gap)
    ## Set up the plot canvas
    plot(1,
         xlim = group_limits,
         ylim = final_y_limits,
         xaxs = "i",
         yaxs = "i",
         type = "n",
         axes = FALSE,
         ann  = FALSE
         )
  }


  ## Draw the grid
  if (draw_grid) {
    DrawGridlines(numeric_limits,
                  extra_grid_lines = extra_grid_lines
                  )
  }

  if (for_embedded_PNG) {
    return(invisible(NULL))
  }

  ## Draw the axis and axis label
  axis(2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, 0.38, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1
       )

  are_no_contam <- grepl("_no_contam_", use_columns, fixed = TRUE)
  contams_excluded <- any(are_no_contam)
  if (contams_excluded) {
    stopifnot(all(are_no_contam))
    # y_axis_label <- paste0(y_axis_label, " (excluding contaminations)")
  }
  mtext(text = y_axis_label, side = 2, line = 3)


  ## Draw the title
  title(main_title, cex.main = 1.1, font.main = 1, line = 1.8)


  ## Draw the x axis labels
  large_gap_lines <- 2.5
  start_line <- 0.65

  if (!(is.null(sub_title))) {
    mtext(text = VerticalAdjust(sub_title),
          at   = mean(group_positions),
          line = start_line + large_gap_lines * 1.3,
          side = 1,
          cex  = 0.9,
          col  = "gray50"
          )
  } else {
    start_line <- start_line + 0.3
  }

  for (i in seq_along(group_positions)) {
    mtext(text = VerticalAdjust(column_labels_list[[use_columns[[i]]]][[1]]),
          at   = group_positions[[i]],
          line = start_line,
          side = 1
          )

    mtext(text = VerticalAdjust(column_labels_list[[use_columns[[i]]]][[2]]),
          at   = group_positions[[i]],
          line = start_line + (large_gap_lines / 2),
          side = 1
          )
  }



  ## Prepare for drawing the legend
  y_mid <- 0.5
  large_gap <- diff(grconvertY(c(0, large_gap_lines), from = "lines", to = "npc"))
  small_gap <- large_gap / 2

  labels_list <- list(
    c("Colony-",
      "picked",
      "controls"
    ),
    c("Polyclonal",
      "plasmid",
      "population"
    )
  )

  is_large_gap <- lapply(labels_list, function(x) {
    if (length(x) == 1) {
      FALSE
    } else {
      c(TRUE, rep(FALSE, length(x) - 1))
    }
  })

  gaps_vec <- ifelse(unlist(is_large_gap), large_gap, small_gap)
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")
  x_start <- 1 + diff(grconvertX(c(0, 1.25), from = "lines", to = "npc"))

  if (draw_legend) {
  ## Draw the legend
    text(x      = grconvertX(x = x_start, from = "npc", to = "user"),
         y      = y_pos,
         cex    = 1,
         labels = sapply(unlist(labels_list), VerticalAdjust),
         adj    = c(0, 0.5),
         xpd    = NA
         )

    colors_mat <- GetColorMat()
    points(x   = rep(grconvertX(x = x_start - 0.0, from = "npc", to = "user"), 2),
           y   = y_pos[c(1, 4)],
           cex = point_cex,
           pch = legend_pch,
           bg  = colors_mat[, "opaque_fill"],
           col = colors_mat[, "opaque_outline"],
           xpd = NA
           )
  }
  return(invisible(NULL))
}




GetColorMat <- function(hue_A = "Purples", hue_B = "Blues") {

  points_alpha <- 0.8
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

  fill_A    <- brewer.pal(9, hue_A)[[5]]
  fill_B    <- brewer.pal(9, hue_B)[[5]]
  outline_A <- brewer.pal(9, hue_A)[[9]]
  outline_B <- brewer.pal(9, hue_B)[[9]]

  original_fills    <- c(fill_A, fill_B)
  opaque_fills      <- vapply(original_fills, Palify, "", fraction_pale = 0.2, USE.NAMES = FALSE)
  transparent_fills <- paste0(original_fills, alpha_hex)

  original_outlines    <- c(outline_A, outline_B)
  opaque_outlines      <- vapply(original_outlines, Palify, "", fraction_pale = 0.2, USE.NAMES = FALSE)
  transparent_outlines <- paste0(original_outlines, alpha_hex)

  all_vars <- c("opaque_fills", "transparent_fills", "opaque_outlines", "transparent_outlines")
  results_mat <- do.call(cbind, sapply(all_vars, get, envir = sys.frame(sys.parent(0)), simplify = FALSE))
  colnames(results_mat) <- sub("s$", "", colnames(results_mat))

  return(results_mat)
}




LollipopPlot <- function(input_df,
                         plate_names       = NULL,
                         use_columns       = c("Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4"),
                         set_mar           = TRUE,
                         label_percentages = TRUE,
                         use_y_limits      = c(0, 1),
                         black_percentages = TRUE,
                         sub_title         = NULL
                         ) {

  plates_list <- GetPlateSelection(plate_names)
  plate_names <- plates_list[["plate_names"]]

  if (set_mar) {
    old_mar <- par("mar" = c(5.5, 5, 4, 8))
  }

  control_mat <- PreparePlates(input_df, plates_df[["Plate_name"]][plates_df[["Colony_picked"]]])
  selected_mat <- PreparePlates(input_df, plate_names)

  assign("delete_input_df", input_df, envir = globalenv())
  assign("delete_control_mat", control_mat, envir = globalenv())
  assign("delete_selected_mat", selected_mat, envir = globalenv())
  assign("delete_use_columns", use_columns, envir = globalenv())

  control_metrics  <- colMeans(control_mat[, use_columns,  drop = FALSE], na.rm = TRUE)
  selected_metrics <- colMeans(selected_mat[, use_columns, drop = FALSE], na.rm = TRUE)

  point_cex <- 1.5
  SetUpPercentagePlot(use_columns, use_y_limits, plates_list[["title"]], point_cex,
                      y_axis_label = "Mean percentage of reads",
                      sub_title = sub_title
                      )

  spaced_percent <- 2.5
  max_space <- (use_y_limits[[2]] - use_y_limits[[1]]) * (spaced_percent / 100)
  are_spaced <- (abs(control_metrics - selected_metrics) > max_space) %in% c(TRUE, NA)

  num_groups <- length(use_columns)
  group_positions <- seq_len(num_groups)
  if (any(are_spaced)) {
    segments(x0   = group_positions[are_spaced],
             x1   = group_positions[are_spaced],
             y0   = selected_metrics[are_spaced],
             y1   = control_metrics[are_spaced],
             col  = "gray75",
             lend = "butt",
             lwd  = 1.5,
             xpd  = NA
             )
  }

  colors_mat <- GetColorMat()

  for (i in 1:2) {
    points(x   = group_positions,
           y   = list(control_metrics, selected_metrics)[[i]],
           cex = point_cex,
           pch = 21,
           col = ifelse(are_spaced, colors_mat[i, "opaque_outline"], colors_mat[i, "transparent_outline"]),
           bg  = ifelse(are_spaced, colors_mat[i, "opaque_fill"], colors_mat[i, "transparent_fill"]),
           xpd = NA
           )
  }

  if (label_percentages) {
    if (black_percentages) {
      percentage_color <- "black"
    } else {
      percentage_color <- "gray70"
    }
    perc_x_gap <- diff(grconvertX(c(0, 1.4), from = "lines", to = "user"))
    text(x      = group_positions + perc_x_gap,
         y      = c(control_metrics, selected_metrics),
         labels = paste0(format(round(c(control_metrics, selected_metrics) * 100, digits = 1)), "%"),
         cex    = 0.7,
         pch    = 16,
         col    = percentage_color,
         xpd    = NA
         )
  }

  par(old_mar)
  return(invisible(NULL))
}



SummaryBoxPlot <- function(input_df,
                           plate_names       = NULL,
                           use_columns       = c("Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4"),
                           set_mar           = TRUE,
                           label_percentages = TRUE,
                           use_y_limits      = c(0, 1),
                           embed_PNG         = FALSE,
                           sub_title         = NULL
                           ) {

  set.seed(1) # For reproducible jitter

  plates_list <- GetPlateSelection(plate_names)
  plate_names <- plates_list[["plate_names"]]

  use_mar <- c(5.5, 5, 4, 8)
  if (set_mar) {
    old_mar <- par("mar" = use_mar)
  }

  control_mat <- PreparePlates(input_df, plates_df[["Plate_name"]][plates_df[["Colony_picked"]]])
  selected_mat <- PreparePlates(input_df, plate_names)

  control_list <- lapply(use_columns, function(x) control_mat[, x])
  selected_list <- lapply(use_columns, function(x) selected_mat[, x])

  control_unlisted <- unlist(control_list)
  selected_unlisted <- unlist(selected_list)

  point_cex <- 1.5

  ## Prepare the raster graphics device
  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(file_output_directory, "temp.png")
    cur_mai <- par("mar") * 0.2
    temp_width <- use_width - sum(cur_mai[c(2, 4)])
    temp_height <- use_height - sum(cur_mai[c(1, 3)])

    png(file   = temp_path,
        width  = temp_width,
        height = temp_height,
        units  = "in",
        res    = 900,
        bg     = "transparent"
        )

    par(mar = rep(0, 4))

    SetUpPercentagePlot(use_columns, use_y_limits, side_gap = 0.3,
                        for_embedded_PNG = TRUE, extra_grid_lines = TRUE
                        )
  } else {
    SetUpPercentagePlot(use_columns, use_y_limits, plates_list[["title"]], point_cex,
                        side_gap = 0.3, legend_pch = 22, extra_grid_lines = TRUE,
                        sub_title = sub_title
                        )
  }


  colors_mat <- GetColorMat()
  num_groups <- length(use_columns)
  group_positions <- seq_len(num_groups)

  use_wex <- 0.4

  use_gap <- 0.35
  control_pos  <- group_positions - (use_gap / 2)
  selected_pos <- group_positions + (use_gap / 2)

  RemoveAllZeros <- function(input_list) {
    lapply(input_list, function(x) {
      if (all(x == 0, na.rm = TRUE)) {
        0
      } else {
        x
      }
    })
  }

  vioplot(RemoveAllZeros(control_list),
          at       = control_pos,
          pchMed   = NA,
          drawRect = FALSE,
          col      = brewer.pal(9, "Purples")[[4]],
          border   = NA,
          wex      = use_wex,
          add      = TRUE,
          axes     = FALSE
          )

  vioplot(RemoveAllZeros(selected_list),
          at       = selected_pos,
          pchMed   = NA,
          drawRect = FALSE,
          col      = brewer.pal(9, "Blues")[[4]],
          border   = NA,
          wex      = use_wex,
          add      = TRUE,
          axes     = FALSE
          )

  control_jittered  <- control_pos[rep(seq_along(control_list), lengths(control_list))] +
                       rnorm(n = length(control_unlisted), mean = 0, sd = 0.03)

  selected_jittered <- selected_pos[rep(seq_along(selected_list), lengths(selected_list))] +
                       rnorm(n = length(selected_unlisted), mean = 0, sd = 0.03)

  points_alpha <- 0.2
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
  ## Draw the jittered points
  points(x   = control_jittered,
         y   = control_unlisted,
         cex = 0.4,
         col = paste0(brewer.pal(9, "Purples")[[8]], alpha_hex),
         pch = 16
         )
  points(x   = selected_jittered,
         y   = selected_unlisted,
         cex = 0.4,
         col = paste0(brewer.pal(9, "Blues")[[8]], alpha_hex),
         pch = 16
         )

  if (embed_PNG) {
    dev.off()
    raster_array <- readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)

    SetUpPercentagePlot(use_columns, use_y_limits, main_title = "",
                        side_gap = 0.3,
                        for_embedded_PNG = TRUE, draw_grid = FALSE
                        )

    rasterImage(raster_array,
                xleft = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                )

    SetUpPercentagePlot(use_columns, use_y_limits, plates_list[["title"]], point_cex,
                        legend_pch = 22, draw_grid = FALSE,
                        make_plot = FALSE, sub_title = sub_title
                        )
  }

  ## Draw the superimposed box plots
  boxplot(control_list,
          at         = control_pos,
          boxwex     = use_wex * 0.4,
          outline    = FALSE,
          names      = rep.int("", length(group_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * 3,
          col        = brewer.pal(9, "Purples")[[2]],
          border     = brewer.pal(9, "Purples")[[9]],
          add        = TRUE,
          axes       = FALSE,
          lwd        = 1
          )

  if (!(all(is.na(selected_mat)))) {
    boxplot(selected_list,
            at         = selected_pos,
            boxwex     = use_wex * 0.4,
            outline    = FALSE,
            names      = rep.int("", length(group_positions)),
            whisklty   = "blank",
            staplewex  = 0,
            whisklwd   = 0,
            staplelty  = 0,
            medlwd     = par("lwd") * 3,
            col        = brewer.pal(9, "Blues")[[2]],
            border     = brewer.pal(9, "Blues")[[9]],
            add        = TRUE,
            axes       = FALSE,
            lwd        = 1
            )
  }

  par(old_mar)
  return(invisible(NULL))
}





# Functions for plotting theoretical/expected proportions -----------------

TheoreticalAtLeastCounts <- function(black_percentages = TRUE,
                                     sub_title = column_group_subtitles_list[["At_least_num_guides"]],
                                     sg_cr_correct_rate = 0.965
                                     ) {

  use_mar <- c(5.5, 5, 4, 8)
  old_mar <- par("mar" = use_mar)

  SetUpPercentagePlot(c("Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4"),
                      c(0, 1), "Theoretical probabilities (gRNA error rate = 3.5%)",
                      draw_legend = FALSE, sub_title = sub_title
                      )

  error_rate <- (1 - sg_cr_correct_rate)
  four_probs <- rev(vapply(0:3, function(x) pbinom(x, size = 4, prob = error_rate), numeric(1)))
  points(x = 1:4, y = four_probs, pch = 16)

  colors_mat <- GetColorMat("Greens")

  group_positions <- 1:4
  point_cex <- 1.5

  points(x   = group_positions,
         y   = four_probs,
         cex = point_cex,
         pch = 21,
         col = colors_mat[1, "opaque_outline"],
         bg  = colors_mat[1, "opaque_fill"],
         xpd = NA
         )

  perc_x_gap <- diff(grconvertX(c(0, 1.8), from = "lines", to = "user"))
  text(x      = group_positions + perc_x_gap,
       y      = four_probs,
       labels = paste0(c(round(four_probs[[1]] * 100, digits = 4),
                         round(four_probs[[2]] * 100, digits = 2),
                         round(four_probs[[3]] * 100, digits = 1),
                         round(four_probs[[4]] * 100, digits = 1)
                         )
                       , "%"),
       cex    = 0.7,
       pch    = 16,
       col    = if (black_percentages) "black" else "gray70",
       xpd    = NA
       )
  par(old_mar)
  old_mar <- par("mar" = use_mar)
}






# Functions for writing all plots to disk ---------------------------------

DrawAllLollipopsAndViolins <- function(export_PNGs = TRUE) {


  message("Exporting plots for theoretical (expected) probabilities...")
  use_file_name <- "Theoretical_at_least_num_guides"
  pdf(file   = file.path(file_output_directory, "Theoretical", paste0(use_file_name, ".pdf")),
      width  = use_width,
      height = use_height
      )
  TheoreticalAtLeastCounts()
  dev.off()

  if (export_PNGs) {
    use_file_name <- "Theoretical_at_least_num_guides"
    png(file   = file.path(PNGs_output_directory, "Theoretical", paste0(use_file_name, ".png")),
        width  = use_width,
        height = use_height,
        units  = "in",
        res    = 600
        )
    TheoreticalAtLeastCounts()
    dev.off()
  }

  ccs_numbers <- c(3, 5, 7)
  accuracy_percentages <- c(99, 99.9, 99.99)

  if (export_PNGs) {
    file_formats <- c("png", "pdf")
  } else {
    file_formats <- "pdf"
  }

  for (file_format in file_formats) {

    for (plot_type in c("Lollipop", "Box")) {

      if (plot_type == "Box") {
        UseFunction <- function(...) SummaryBoxPlot(..., embed_PNG = TRUE)
      } else {
        UseFunction <- LollipopPlot
      }

      for (i in seq_along(ccs_numbers)) {

        use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))

        filter_stages <- "filtered_summary_df"
        filter_labels <- "ii) filtered"
        if ("filtered_cross_plate_df" %in% names(use_df_list)) {
          filter_stages <- c(filter_stages, "filtered_cross_plate_df")
          filter_labels <- c(filter_labels, "iii) filtered cross-plate")
        }

        for (filter_stage in seq_along(filter_stages)) {

          df_name <- filter_stages[[filter_stage]] # "filtered_gRNAs_df"

          use_df <- use_df_list[[df_name]]

          folder_name <- paste0("CCS", ccs_numbers[[i]],
                                " (", accuracy_percentages[[i]], ") - ",
                                filter_labels[[filter_stage]]
                                )

          message(paste0("Exporting ", file_format, " images of ",
                         tolower(plot_type), " plots for the following data: ",
                         folder_name, "..."
                         )
                  )

          if (file_format == "pdf") {
            folder_path <- file.path(file_output_directory, folder_name)
          } else if (file_format == "png") {
            folder_path <- file.path(PNGs_output_directory, folder_name)
          }
          dir.create(folder_path, showWarnings = FALSE)

          for (label_percentages in c(TRUE, FALSE)) {

            if (plot_type == "Box") {
              if (label_percentages) {
                next
              }
            } else {
              if (label_percentages) {
                message(paste0("  ... with labelled percentages..."))
              } else {
                message(paste0("  ... without labelled percentages..."))
              }
            }

            for (j in seq_along(column_groups_list)) {
              metric <- names(column_groups_list)[[j]]
              message(paste0("    ... for the metric: '", metric, "'"))
              file_name <- paste0(plot_type, " plot - ", j, ") ",
                                  sub("Count_", "", metric, fixed = TRUE),
                                  ".pdf"
                                  )
              sub_folder_name <- paste0(plot_type, " plots")
              if (plot_type == "Lollipop") {
                if (label_percentages) {
                  sub_folder_name <- paste0(sub_folder_name, " - with percentages")
                } else {
                  sub_folder_name <- paste0(sub_folder_name, " - plain")
                }
              }
              sub_folder_path <- file.path(folder_path, sub_folder_name)
              dir.create(sub_folder_path, showWarnings = FALSE)
              if (file_format == "pdf") {
                pdf(file = file.path(sub_folder_path, file_name),
                    width = use_width, height = use_height
                    )
              } else if (file_format == "png") {
                metric_path <- file.path(sub_folder_path, paste0(j, ") ", metric))
                dir.create(metric_path, showWarnings = FALSE)
              }

              for (k in seq_along(plate_selection_titles_list)) {
                plate_selection <- names(plate_selection_titles_list)[[k]]
                # message(paste0("          for the plates: ", plate_selection))
                if (file_format == "png") {
                  file_name <- paste0(plate_selection_prefixes[[k]],
                                      ") ", plate_selection, ".png"
                                      )
                  png(file   = file.path(metric_path, file_name),
                      width  = use_width,
                      height = use_height,
                      units  = "in",
                      res    = 600
                      )
                }
                UseFunction(use_df,
                            plate_names = plate_selection,
                            use_columns = column_groups_list[[metric]],
                            label_percentages = label_percentages,
                            sub_title = column_group_subtitles_list[[metric]]
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







