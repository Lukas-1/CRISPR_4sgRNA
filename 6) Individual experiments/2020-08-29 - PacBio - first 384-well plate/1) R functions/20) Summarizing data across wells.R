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
    "Count_sg_cr"    = column_group_subtitles_list[["At_least_num_guides"]],
    "Count_pr_sg_cr" = column_group_subtitles_list[["Promoter_num_guides"]]
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
    if (plate_names == "Colony-picked") {
      main_title <- "Colony-picked controls"
    } else {
      main_title <- plate_selection_titles_list[[plate_names]]
    }
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

GetMinQuality <- function(input_integer) {
  is_even <- (input_integer %% 2) == 0
  num_nines <- ceiling(input_integer / 2)
  digits_string <- paste0(rep("9", num_nines), collapse = "")
  if (is_even) {
    digits_string <- paste0(digits_string, "5")
  }
  final_number <- as.numeric(paste0("0.", digits_string))
  return(final_number)
}



SampleControlReads <- function(reads_df, plates_df, ccs_number = 7) {

  control_plates <- plates_df[plates_df[, "Colony_picked"], "Plate_number"]
  are_control_plates <- reads_df[, "Plate_number"] %in% control_plates
  exclude_columns <- "Pool"

  non_controls_df <- reads_df[!(are_control_plates), names(reads_df) != exclude_columns]
  median_read_count <- median(as.integer(table(non_controls_df[, "Combined_ID"])))

  controls_df <- reads_df[are_control_plates, names(reads_df) != exclude_columns]
  are_passing <- (controls_df[, "Passes_barcode_filters"] == 1) &
                 (controls_df[, "Read_quality"] >= GetMinQuality(ccs_number)) &
                 (controls_df[, "Num_full_passes"] >= ccs_number)
  if (!(all(are_passing))) {
    controls_df <- controls_df[are_passing, ]
  }

  controls_df_list <- split(controls_df, controls_df[, "Well_number"])
  set.seed(1)
  controls_df_list <- lapply(controls_df_list, function(x) {
    random_indices <- sample(seq_len(nrow(x)), median_read_count)
    x[random_indices, ]
  })
  sampled_controls_df <- do.call(rbind.data.frame,
                                 c(controls_df_list,
                                   list(stringsAsFactors = FALSE,
                                        make.row.names = FALSE
                                        )
                                   ))

  sampled_controls_df[, "Plate_number"] <- control_plates[[1]]
  sampled_controls_df[, "Combined_ID"] <- paste0(
    "Plate", formatC(control_plates[[1]], width = 3, flag = "0"),
    "_Well", formatC(sampled_controls_df[, "Well_number"], width = 3, flag = "0")
  )

  results_df <- rbind.data.frame(sampled_controls_df,
                                 non_controls_df,
                                 stringsAsFactors = FALSE,
                                 make.row.names = FALSE
                                 )
  return(results_df)
}




ConsolidateControls <- function(summary_df, plates_df) {

  ID_columns <- c("Combined_ID", "Plate_number", "Well_number")
  exclude_columns <- c(
    "Mean_read_quality", "Expected_from_close_wells",
    "Mean_distance", "Expected_distance", "Distance_p_value",
    grep("Perc", names(summary_df), fixed = TRUE, value = TRUE)
  )
  sum_columns <- setdiff(names(summary_df), c(ID_columns, exclude_columns))

  control_plates <- plates_df[plates_df[, "Colony_picked"], "Plate_number"]
  are_controls <- summary_df[, "Plate_number"] %in% control_plates
  controls_df <- summary_df[are_controls, ]
  not_controls_df <- summary_df[!(are_controls), ]

  sums_mat <- do.call(cbind, sapply(sum_columns, function(x) {
    tapply(controls_df[, x], controls_df[, "Well_number"], sum)
  }, simplify = FALSE))

  use_control_plate <- control_plates[[1]]
  sum_controls_df <- data.frame(
    controls_df[controls_df[, "Plate_number"] %in% use_control_plate, ID_columns],
    sums_mat,
    stringsAsFactors = FALSE
  )
  results_df <- rbind.data.frame(sum_controls_df,
                                 not_controls_df[, c(ID_columns, sum_columns)],
                                 stringsAsFactors = FALSE,
                                 make.row.names = FALSE
                                 )
  return(results_df)
}



StandardizeCounts <- function(input_df) {
  stopifnot(all(c("count_columns", "num_columns", "alignment_columns") %in% ls(envir = globalenv())))

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



HandleEmptyWells <- function(input_mat) {
  exclude_NA <- grepl("mut|del|contam", colnames(input_mat), ignore.case = TRUE)
  include_NA <- grepl("correct|count",  colnames(input_mat), ignore.case = TRUE)
  if (all(include_NA)) {
    are_NA <- is.na(input_mat)
    input_mat[are_NA] <- 0
  } else if (!(all(exclude_NA))) {
    stop("Unexpected column names!")
  }
  return(input_mat)
}





# Functions for plotting summary statistics -------------------------------

SetUpPercentagePlot <- function(use_columns,
                                use_y_limits,
                                main_title,
                                title_line         = 1.8,
                                title_cex          = 1.1,
                                point_cex          = 1.5,
                                legend_pch         = 21,
                                side_gap           = 0.5,
                                draw_grid          = TRUE,
                                extra_grid_lines   = FALSE,
                                for_embedded_PNG   = FALSE,
                                make_plot          = TRUE,
                                add_mean           = FALSE,
                                y_axis_label       = "Percentage of reads",
                                y_label_line       = 3,
                                sub_title          = NULL,
                                draw_legend        = TRUE,
                                bottom_labels_list = NULL,
                                bottom_start_y     = 0.65,
                                bottom_y_gap       = 2.5,
                                bottom_large_y     = 1,
                                side_y_gap         = bottom_y_gap,
                                side_labels_list   = list(c("Colony-", "picked", "controls"),
                                                          c("Polyclonal", "plasmid", "population")
                                                          ),
                                points_centered    = FALSE,
                                side_legend_x      = 1.25,
                                side_legend_x_gap  = 0
                                ) {

  stopifnot("column_labels_list" %in% ls(envir = globalenv()))

  ## Determine group positions
  assign("delete_group_positions", use_columns, envir = globalenv())
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
    y_gap <- final_y_range * 0.02
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
       las      = 1,
       lwd      = par("lwd")
       )

  are_no_contam <- grepl("_no_contam_", use_columns, fixed = TRUE)
  contams_excluded <- any(are_no_contam)
  if (contams_excluded) {
    stopifnot(all(are_no_contam))
    # y_axis_label <- paste0(y_axis_label, " (excluding contaminations)")
  }
  mtext(text = y_axis_label, side = 2, line = y_label_line, cex = par("cex"))


  ## Draw the title
  title(main_title, cex.main = title_cex, font.main = 1, line = title_line)


  ## Draw the x axis labels

  if (!(is.null(sub_title))) {
    mtext(text = VerticalAdjust(sub_title),
          at   = mean(group_positions),
          line = bottom_start_y + bottom_y_gap * 1.3,
          side = 1,
          cex  = 0.9 * par("cex"),
          col  = "gray50"
          )
  } else {
    bottom_start_y <- bottom_start_y + 0.3
  }

  if (is.null(bottom_labels_list)) {
    bottom_labels_list <- column_labels_list
  }

  for (i in seq_along(group_positions)) {
    mtext(text = VerticalAdjust(bottom_labels_list[[use_columns[[i]]]][[1]]),
          at   = group_positions[[i]],
          line = bottom_start_y,
          side = 1,
          cex  = par("cex")
          )

    mtext(text = VerticalAdjust(bottom_labels_list[[use_columns[[i]]]][[2]]),
          at   = group_positions[[i]],
          line = bottom_start_y + (bottom_y_gap / 2),
          side = 1,
          cex  = par("cex")
          )
  }


  ## Prepare for drawing the legend
  y_mid <- 0.5
  large_gap <- diff(grconvertY(c(0, side_y_gap), from = "lines", to = "npc"))
  small_gap <- large_gap / 2

  is_large_gap <- lapply(side_labels_list, function(x) {
    if (length(x) == 1) {
      FALSE
    } else {
      c(TRUE, rep(FALSE, length(x) - 1))
    }
  })

  gaps_vec <- ifelse(unlist(is_large_gap), large_gap * bottom_large_y, small_gap)
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")
  x_start <- 1 + diff(grconvertX(c(0, side_legend_x), from = "lines", to = "npc"))

  if (draw_legend) {
  ## Draw the legend
    text(x      = grconvertX(x = x_start, from = "npc", to = "user"),
         y      = y_pos,
         cex    = 1,
         labels = sapply(unlist(side_labels_list), VerticalAdjust),
         adj    = c(0, 0.5),
         xpd    = NA
         )

    colors_mat <- GetColorMat()
    side_legend_x_gap <- diff(grconvertX(c(0, side_legend_x_gap), from = "lines", to = "npc"))
    points(x   = rep(grconvertX(x = x_start - side_legend_x_gap, from = "npc", to = "user"), 2),
           y   = if (points_centered) y_pos[c(2, 5)] else y_pos[c(1, 4)],
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
                         custom_title      = NULL,
                         ...
                         ) {

  plates_list <- GetPlateSelection(plate_names)
  plate_names <- plates_list[["plate_names"]]

  input_df <- ConsolidateControls(input_df, plates_df)

  control_mat  <- PreparePlates(input_df, plates_df[plates_df[, "Colony_picked"], "Plate_name"])
  selected_mat <- PreparePlates(input_df, plate_names)[, use_columns, drop = FALSE]

  assign("delete_selected_mat_1", selected_mat, envir = globalenv())
  assign("delete_plate_names", plate_names, envir = globalenv())

  control_mat  <- HandleEmptyWells(control_mat[, use_columns, drop = FALSE])
  selected_mat <- HandleEmptyWells(selected_mat[, use_columns, drop = FALSE])

  assign("delete_input_df", input_df, envir = globalenv())
  assign("delete_control_mat", control_mat, envir = globalenv())
  assign("delete_selected_mat", selected_mat, envir = globalenv())
  assign("delete_use_columns", use_columns, envir = globalenv())

  control_metrics  <- colMeans(control_mat[, use_columns,  drop = FALSE], na.rm = TRUE)
  selected_metrics <- colMeans(selected_mat[, use_columns, drop = FALSE], na.rm = TRUE)

  point_cex <- 1.5
  if (is.null(custom_title)) {
    use_title <- plates_list[["title"]]
  } else {
    use_title <- custom_title
  }
  if (set_mar) {
    old_mar <- par("mar" = c(5.5, 5, 4, 8))
  }
  SetUpPercentagePlot(use_columns, use_y_limits, use_title, point_cex,
                      y_axis_label = "Mean percentage of reads", ...
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

  if (set_mar) {
    par(old_mar)
  }
  return(invisible(NULL))
}



CustomBoxPlot <- function(input_list,
                          at_positions,
                          use_brewer_pal,
                          use_wex       = 0.4,
                          draw_whiskers = TRUE,
                          median_lwd    = 3
                          ) {

  if (draw_whiskers) {
    quantile_mat <- t(sapply(input_list, quantile, probs = c(0.05, 0.95), na.rm = TRUE))
    whisker_color <- brewer.pal(9, use_brewer_pal)[[5]]
    segments(x0   = at_positions,
             y0   = quantile_mat[, 1],
             y1   = quantile_mat[, 2],
             col  = whisker_color,
             lwd  = par("lwd") * 1.5,
             xpd  = NA,
             lend = "butt"
             )
    # whisker_half_width <- 0.025
    # segments(x0   = at_positions - whisker_half_width,
    #          x1   = at_positions + whisker_half_width,
    #          y0   = quantile_mat[, 1],
    #          col  = whisker_color,
    #          lwd  = par("lwd"),
    #          xpd  = NA,
    #          lend = "butt"
    #          )
    # segments(x0   = at_positions - whisker_half_width,
    #          x1   = at_positions + whisker_half_width,
    #          y0   = quantile_mat[, 2],
    #          col  = whisker_color,
    #          lwd  = par("lwd"),
    #          xpd  = NA,
    #          lend = "butt"
    #          )
  }

  boxplot(input_list,
          at         = at_positions,
          boxwex     = use_wex * 0.4,
          outline    = FALSE,
          names      = rep.int("", length(at_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * median_lwd,
          col        = brewer.pal(9, use_brewer_pal)[[2]],
          border     = brewer.pal(9, use_brewer_pal)[[9]],
          add        = TRUE,
          axes       = FALSE,
          lwd        = par("lwd") * 1.5
          )

  return(invisible(NULL))

}


ViolinClippedBorders <- function(violin_list,
                                 at_positions  = seq_along(violin_list),
                                 border_colors = "black",
                                 fill          = FALSE,
                                 use_wex       = 0.4,
                                 use_lwd       = 0.75,
                                 clip_distance = 0.015,
                                 use_quantiles = c(0.005, 0.995)
                                 ) {

  assign("delete_violin_list", violin_list, envir = globalenv())
  assign("delete_at_positions", at_positions, envir = globalenv())

  if (length(border_colors) == 1) {
    border_colors <- rep(border_colors, length(violin_list))
  }
  violin_list <- RemoveIfOnlyZeros(violin_list)

  stopifnot(length(violin_list) == length(at_positions))
  for (i in seq_along(violin_list)) {
    use_borders <- vapply(use_quantiles, function(x) {
      use_border <- quantile(violin_list[[i]], probs = x, na.rm = TRUE, names = FALSE)
      return(grconvertY(use_border, from = "user", to = "npc"))
    }, numeric(1))
    bottom_border <- use_borders[[1]]
    top_border    <- use_borders[[2]]
    distance_to_top <- 1 - top_border
    distance_to_bottom <- bottom_border
    if (fill) {
      ## We only want to clip off the extreme ends of the violin,
      ## not the area close to the axis.

      if (distance_to_top > distance_to_bottom) {
        bottom_border <- 0
        if (top_border < 0.3) {
          top_border <- 0.3
        }
      } else {
        top_border <- 1
        if (bottom_border > 0.7) {
          bottom_border <- 0.7
        }
      }
    } else {
      if (distance_to_bottom > 0.5) {
        bottom_border <- bottom_border - 0.05
      } else if (distance_to_top > 0.5) {
        top_border <- top_border + 0.05
      }
      ## We do not want the border of the violin to extend along the bottom
      ## or top of the plot, hence the clipping at 0.02 and 99.8.
      bottom_border <- max(c(clip_distance, bottom_border))
      top_border    <- min(c(1 - clip_distance, top_border))
    }
    do.call(clip, as.list(c(par("usr")[c(1, 2)],
                            grconvertY(c(bottom_border, top_border), from = "npc", to = "user")
                            )))
    vioplot(violin_list[[i]],
            at       = at_positions[[i]],
            pchMed   = NA,
            drawRect = FALSE,
            col      = if (fill) border_colors[[i]] else "#00000000",
            border   = if (fill) border_colors[[i]] else adjustcolor(border_colors[[i]], alpha.f = 0.7),
            wex      = use_wex,
            lwd      = use_lwd,
            add      = TRUE,
            axes     = FALSE
            )
    do.call(clip, as.list(par("usr")))
  }
  return(invisible(NULL))
}



RemoveIfOnlyZeros <- function(input_list) {
  lapply(input_list, function(x) {
    if (all(x == 0, na.rm = TRUE)) {
      0
    } else {
      x
    }
  })
}


SummaryBoxPlot <- function(input_df,
                           plate_names       = NULL,
                           use_columns       = c("Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4"),
                           set_mar           = TRUE,
                           label_percentages = TRUE,
                           use_y_limits      = c(0, 1),
                           embed_PNG         = FALSE,
                           embed_PNG_res     = 900,
                           custom_title      = NULL,
                           draw_whiskers     = FALSE,
                           median_lwd        = 3,
                           box_wex           = 0.4,
                           violin_wex        = 0.4,
                           use_side_gap      = 0.3,
                           controls_x_gap    = 0.35,
                           ...
                           ) {

  set.seed(1) # For reproducible jitter

  plates_list <- GetPlateSelection(plate_names)
  plate_names <- plates_list[["plate_names"]]

  use_mar <- c(5.5, 5, 4, 8)
  if (set_mar) {
    old_mar <- par("mar" = use_mar)
  }

  input_df <- ConsolidateControls(input_df, plates_df)

  control_mat <- PreparePlates(input_df, plates_df[plates_df[, "Colony_picked"], "Plate_name"])
  selected_mat <- PreparePlates(input_df, plate_names)

  control_mat  <- HandleEmptyWells(control_mat[, use_columns, drop = FALSE])
  selected_mat <- HandleEmptyWells(selected_mat[, use_columns, drop = FALSE])

  control_list <- lapply(use_columns, function(x) control_mat[, x])
  selected_list <- lapply(use_columns, function(x) selected_mat[, x])

  control_unlisted <- unlist(control_list)
  selected_unlisted <- unlist(selected_list)

  point_cex <- 1.5

  if (is.null(custom_title)) {
    use_title <- plates_list[["title"]]
  } else {
    use_title <- custom_title
  }
  ## Prepare the raster graphics device
  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(file_output_directory, "temp.png")
    temp_width  <- par("pin")[[1]]
    temp_height <- par("pin")[[2]]
    current_par <- par(no.readonly = TRUE)

    png(filename = temp_path,
        width    = temp_width,
        height   = temp_height,
        units    = "in",
        res      = embed_PNG_res,
        bg       = "transparent"
        )

    par(lwd = current_par[["lwd"]])
    par(cex = current_par[["cex"]])
    par(mar = rep(0, 4))
    SetUpPercentagePlot(use_columns, use_y_limits, NULL, side_gap = use_side_gap,
                        for_embedded_PNG = TRUE, extra_grid_lines = TRUE
                        )
  } else {
    SetUpPercentagePlot(use_columns, use_y_limits, use_title, point_cex,
                        side_gap = use_side_gap, legend_pch = 22, extra_grid_lines = TRUE,
                        ...
                        )
  }


  colors_mat <- GetColorMat()
  num_groups <- length(use_columns)
  group_positions <- seq_len(num_groups)

  control_pos  <- group_positions - (controls_x_gap / 2)
  selected_pos <- group_positions + (controls_x_gap / 2)

  RemoveIfOnlyZeros <- function(input_list) {
    lapply(input_list, function(x) {
      if (all(x == 0, na.rm = TRUE)) {
        0
      } else {
        x
      }
    })
  }

  use_brewer_pals <- c("Purples", "Blues")

  if (draw_whiskers) {
    violin_colors <- vapply(use_brewer_pals, function(x) brewer.pal(9, x)[[3]], "")
  } else {
    violin_colors <- vapply(use_brewer_pals, function(x) brewer.pal(9, x)[[4]], "")
  }
  violin_border_width <- 1

  ViolinClippedBorders(control_list, control_pos, fill = TRUE,
                       border_colors = violin_colors[[1]],
                       use_wex = violin_wex, use_lwd = violin_border_width
                       )
  ViolinClippedBorders(selected_list, selected_pos, fill = TRUE,
                       border_colors = violin_colors[[2]],
                       use_wex = violin_wex, use_lwd = violin_border_width
                       )

  vioplot(RemoveIfOnlyZeros(selected_list),
          at       = selected_pos,
          pchMed   = NA,
          drawRect = FALSE,
          col      = violin_colors[[2]],
          border   = NA,
          wex      = violin_wex,
          add      = TRUE,
          axes     = FALSE
          )

  vioplot(RemoveIfOnlyZeros(control_list),
          at       = control_pos,
          pchMed   = NA,
          drawRect = FALSE,
          col      = violin_colors[[1]],
          border   = NA,
          wex      = violin_wex,
          add      = TRUE,
          axes     = FALSE
          )


  if (draw_whiskers) {
    jitter_sd <- 0.02
  } else {
    jitter_sd <- 0.03
  }

  control_jittered  <- control_pos[rep(seq_along(control_list), lengths(control_list))] +
                       rnorm(n = length(control_unlisted), mean = 0, sd = jitter_sd)

  selected_jittered <- selected_pos[rep(seq_along(selected_list), lengths(selected_list))] +
                       rnorm(n = length(selected_unlisted), mean = 0, sd = jitter_sd)

  if (draw_whiskers) {
    point_colors <- vapply(use_brewer_pals, function(x) {
      colorRampPalette(brewer.pal(9, x))(101)[[80]]
    }, "")
  } else {
    point_colors <- vapply(use_brewer_pals, function(x) brewer.pal(9, x)[[8]], "")
  }
  point_colors <- adjustcolor(point_colors, alpha.f = 0.2)

  ## Draw the jittered points
  points(x   = control_jittered,
         y   = control_unlisted,
         cex = 0.4,
         col = point_colors[[1]],
         pch = 16
         )
  points(x   = selected_jittered,
         y   = selected_unlisted,
         cex = 0.4,
         col = point_colors[[2]],
         pch = 16
         )

  if (embed_PNG) {
    dev.off()
    raster_array <- readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)

    SetUpPercentagePlot(use_columns, use_y_limits, main_title = "",
                        side_gap = use_side_gap,
                        for_embedded_PNG = TRUE, draw_grid = FALSE
                        )

    rasterImage(raster_array,
                xleft = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                )

    SetUpPercentagePlot(use_columns, use_y_limits, use_title, point_cex,
                        legend_pch = 22, draw_grid = FALSE,
                        make_plot = FALSE, ...
                        )
  }

  ## Draw the superimposed borders of the violin plots
  ViolinClippedBorders(control_list, control_pos,
                       border_colors = violin_colors[[1]],
                       use_wex = violin_wex, use_lwd = violin_border_width
                       )
  ViolinClippedBorders(selected_list, selected_pos,
                       border_colors = violin_colors[[2]],
                       use_wex = violin_wex, use_lwd = violin_border_width,
                       clip_distance = 0.15
                       )

  ## Draw the superimposed box plots
  assign("delete_control_list", control_list, envir = globalenv())
  CustomBoxPlot(control_list, control_pos, use_brewer_pals[[1]],
                use_wex = box_wex, draw_whiskers = draw_whiskers
                )

  if (!(all(is.na(selected_mat)))) {
    CustomBoxPlot(selected_list, selected_pos, use_brewer_pals[[2]],
                  use_wex = box_wex, draw_whiskers = draw_whiskers
                  )
  }

  if (set_mar) {
    par(old_mar)
  }
  return(invisible(NULL))
}



SummaryDfForPlates <- function(plate_names, use_summary_df) {
  plates_list <- GetPlateSelection(plate_names)
  plate_names <- plates_list[["plate_names"]]
  plate_matches <- match(plate_names, plates_df[, "Plate_name"])
  plate_numbers <- plates_df[["Plate_number"]][plate_matches]
  summary_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_numbers, ]
  results_list <- list(
    "summary_df" = summary_df,
    "title"      = plates_list[["title"]]
  )
  return(results_list)
}



SummaryStackedBars <- function(summary_df,
                               plate_names        = NULL,
                               consider_tracrRNAs = FALSE,
                               top_title          = NULL,
                               top_title_line     = 1.2,
                               top_title_cex      = 1.1,
                               top_title_font     = 2,
                               sub_title          = NULL,
                               set_mar            = TRUE,
                               y_label_line       = 3,
                               x_labels_line      = 0.5,
                               use_side_gap       = 0.5,
                               ...
                               ) {

  ## Filter the input data frame
  summary_df <- ConsolidateControls(summary_df, plates_df)
  assign("delete_consolidate_df", summary_df, envir = globalenv())
  if (!(is.null(plate_names))) {
    filtered_list <- SummaryDfForPlates(plate_names, summary_df)
    summary_df <- filtered_list[["summary_df"]]
    if (is.null(top_title)) {
      top_title <- filtered_list[["title"]]
    }
  }

  ## Prepare column names
  four_categories <- c("Correct", "Mutation", "Deletion", "Contamination")
  if (consider_tracrRNAs) {
    columns_list <- lapply(four_categories, function(x) {
      paste0(x, "_sg", 1:4, "_cr", 1:4)
    })
  } else {
    columns_list <- lapply(four_categories, function(x) {
      paste0(x, "_sg", 1:4)
    })
  }

  ## Extract data for columns
  counts_mat_list <- lapply(1:4, function(x) {
    use_columns <- sapply(columns_list, "[[", x)
    as.matrix(summary_df[, use_columns])
  })

  ## Check that all counts add up correctly
  counts_list <- c(list(summary_df[, "Count_total"]),
                   lapply(counts_mat_list, function(x) as.integer(rowSums(x)))
                   )
  assign("delete_counts_list", counts_list, envir = globalenv())
  stopifnot(length(unique(counts_list)) == 1)
  rm(counts_list)

  ## Convert counts to fractions. Ensure that NA values are counted.
  no_data_vec <- ifelse(summary_df[, "Count_total"] == 0L, 1L, 0L)
  pseudo_counts <- ifelse(summary_df[, "Count_total"] == 0L, 1L, summary_df[, "Count_total"])
  fractions_mat_list <- lapply(counts_mat_list, function(x) {
    x <- cbind(x, no_data_vec)
    apply(x, 2, function(x) x / pseudo_counts)
  })

  ## Summarize data across wells, then add the mean values
  fractions_mat <- do.call(cbind, lapply(fractions_mat_list, colMeans, na.rm = TRUE))
  dimnames(fractions_mat) <- list(c(four_categories, "No_data"), paste0("sg", 1:4))
  fractions_mat <- cbind(fractions_mat, "Mean" = rowMeans(fractions_mat))

  ## Determine group positions
  num_groups <- ncol(fractions_mat)
  group_positions <- seq_len(num_groups)
  side_gap <- use_side_gap
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  numeric_axis_pos <- pretty(c(0, 1))
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  numeric_axis_labels <- paste0(format(numeric_axis_pos * 100), "%")
  final_y_range <- numeric_limits[[2]] - numeric_limits[[1]]
  y_gap <- final_y_range * 0.02
  final_y_limits <- c(numeric_limits[[1]] - y_gap, numeric_limits[[2]] + y_gap)

  ## Set up the plot canvas
  if (set_mar) {
    old_mar <- par("mar" = c(5, 5, 4, 8))
  }
  plot(1,
       xlim = group_limits,
       ylim = final_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  DrawGridlines(numeric_limits, extra_grid_lines = TRUE)

  ## Draw the title and/or sub-title
  if (!(is.null(sub_title))) {
    mtext(sub_title, line = 1.8, cex = par("cex"), side = 1)
  }
  if (!(is.null(top_title) || isFALSE(top_title))) {
    mtext(top_title,
          line = top_title_line,
          cex  = par("cex") * top_title_cex,
          font = top_title_font
          )
  }

  ## Draw the axes and axis labels
  mtext(sapply(c(paste0("sg", 1:4), expression(italic("mean"))), VerticalAdjust),
        at   = seq_len(num_groups),
        side = 1,
        line = x_labels_line,
        cex  = par("cex")
        )
  axis(2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, 0.38, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1,
       lwd      = par("lwd")
       )
  mtext(text = "Percentage of reads", side = 2, line = y_label_line, cex = par("cex"))

  ## Prepare for drawing the bars
  four_colors <- c("#F9F4EC", "#DB678B", colorRampPalette(brewer.pal(9, "Blues")[c(4, 5)])(5)[[2]], "#601A4A", NA)

  num_categories <- nrow(fractions_mat)
  lower_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    if (x == 1) {
      rep(0, ncol(fractions_mat))
    } else {
      colSums(fractions_mat[seq_len(x - 1), , drop = FALSE])
    }
  })
  upper_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    colSums(fractions_mat[seq_len(x), , drop = FALSE])
  })

  ## Draw the bars
  bar_half_width <- 0.3
  for (i in seq_len(ncol(fractions_mat))) {
    for (j in 1:4) {
      rect(xleft   = group_positions[[i]] - bar_half_width,
           xright  = group_positions[[i]] + bar_half_width,
           ybottom = lower_borders_vec_list[[j]][[i]],
           ytop    = upper_borders_vec_list[[j]][[i]],
           col     = four_colors[[j]],
           border  = NA,
           xpd     = NA
           )
    }
  }

  ## Draw bar borders (for the "correct" category only)
  half_width <- bar_half_width - GetHalfLineWidth()
  segments(x0  = group_positions - half_width,
           y0  = 0,
           y1  = upper_borders_vec_list[[1]],
           col = Darken(four_colors[[1]], factor = 1.1),
           lend = "butt",
           xpd = NA
           )
  segments(x0  = group_positions + half_width,
           y0  = 0,
           y1  = upper_borders_vec_list[[1]],
           col = Darken(four_colors[[1]], factor = 1.1),
           lend = "butt",
           xpd = NA
           )
  ## Re-draw the x axis line
  segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]], y0 = 0, xpd = NA)

  ## Draw the side legend
  means_vec <- formatC(fractions_mat[, "Mean"] * 100, digits = 1, format = "f")
  labels_list <- lapply(1:4, function(x) {
    bquote_list <- list(bquote(.(four_categories[[x]])),
                        bquote(italic("(" * .(means_vec[[x]]) * "%)"))
                        )
    sapply(bquote_list, as.expression)
  })
  labels_order <- c(1, 4:2)
  DrawSideLegend(labels_list[labels_order],
                 four_colors[labels_order],
                 use_pch = 22, ...
                 )

  ## Final steps
  if (set_mar) {
    par(old_mar)
  }
  return(invisible(NULL))
}



ReadCountsBoxPlot <- function(summary_df,
                              plate_selections = list("CRISPRa", "CRISPRko"),
                              selection_labels = NULL,
                              x_labels_line    = 0.5,
                              y_label_line     = 3,
                              side_gap         = 0.5,
                              embed_PNG        = FALSE
                              ) {

  ## Prepare data on read counts
  selections_list <- lapply(plate_selections, function(x) {
    summary_results <- SummaryDfForPlates(x, summary_df)
    results_list <- c(list("count" = summary_results[["summary_df"]][, "Count_total"]),
                      summary_results["title"]
                      )
    return(results_list)
  })
  counts_list <- lapply(selections_list, function(x) x[["count"]])
  if (is.null(selection_labels)) {
    selection_labels <- vapply(selections_list, function(x) x[["title"]], "")
  }

  ## Determine group positions
  num_groups <- length(plate_selections)
  group_positions <- seq_len(num_groups)
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  y_span <- max(unlist(counts_list))
  y_space <- y_span * 0.02
  use_y_limits <- c(-y_space, y_span + y_space)

  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(file_output_directory, "temp.png")
    temp_width  <- par("pin")[[1]]
    temp_height <- par("pin")[[2]]
    current_par <- par(no.readonly = TRUE)

    png(filename = temp_path,
        width    = temp_width,
        height   = temp_height,
        units    = "in",
        res      = 900,
        bg       = "transparent"
        )

    par(lwd = current_par[["lwd"]])
    par(cex = current_par[["cex"]])
    par(mar = rep(0, 4))
  }

  ## Set up the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = use_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )

  ## Prepare the plot colors
  light_colors <- c("#FDE0EF", "#D1E5F0")
  dark_colors <- c("#E44499", "#2363B2")
  violin_colors <- vapply(seq_along(light_colors),
                          function(x) colorRampPalette(c(light_colors[[x]], dark_colors[[x]]))(5)[[2]],
                          ""
                          )

  ## Draw the violin plots (in the background)
  vioplot(counts_list,
          pchMed   = NA,
          drawRect = FALSE,
          col      = violin_colors,
          border   = violin_colors,
          wex      = 0.8,
          add      = TRUE,
          axes     = FALSE
          )

  ## Draw the x axis line
  segments(x0  = par("usr")[[1]],
           x1  = par("usr")[[2]] - ((par("usr")[[2]] - par("usr")[[1]]) * 0.005), # Prevent line end from being clipped (when embed_PNG is TRUE)
           y0  = 0,
           xpd = NA
           )

  ## Draw the jittered points
  set.seed(1)
  x_vec <- rep(group_positions, lengths(counts_list))
  x_vec <- x_vec + rnorm(n = length(x_vec), mean = 0, sd = 0.06)
  points(x   = x_vec,
         y   = unlist(counts_list),
         cex = 0.4,
         col = rep(adjustcolor(dark_colors, alpha.f = 0.2), lengths(counts_list)),
         pch = 16,
         xpd = NA
         )

  if (embed_PNG) {
    dev.off()
    raster_array <- readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)

    plot(1,
         xlim = group_limits,
         ylim = use_y_limits,
         xaxs = "i",
         yaxs = "i",
         type = "n",
         axes = FALSE,
         ann  = FALSE
         )

    rasterImage(raster_array,
                xleft   = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop   = par("usr")[[4]]
                )
  }

  ## Draw the violins' borders
  ViolinClippedBorders(counts_list,
                       group_positions,
                       border_colors = violin_colors,
                       use_wex       = 0.8,
                       use_quantiles = c(0, 1),
                       clip_distance = 0.03
                       )

  ## Draw the boxplot whiskers
  quantile_mat <- t(sapply(counts_list, quantile, probs = c(0.05, 0.95), na.rm = TRUE))
  segments(x0   = group_positions,
           y0   = quantile_mat[, 1],
           y1   = quantile_mat[, 2],
           col  = violin_colors,
           lwd  = par("lwd") * 1.5,
           lend = "butt",
           xpd  = NA
           )

  ## Draw the superimposed box plots
  boxplot(counts_list,
          at         = group_positions,
          boxwex     = 0.3,
          outline    = FALSE,
          names      = rep.int("", length(group_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * 2,
          col        = light_colors,
          border     = dark_colors,
          add        = TRUE,
          axes       = FALSE,
          lwd        = par("lwd") * 1.5
          )

  ## Draw the y axis and x and y axis labels
  axis(2,
       mgp      = c(3, 0.38, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1,
       lwd      = par("lwd")
       )
  mtext("HiFi reads per well", side = 2, line = y_label_line, cex = par("cex"))

  labels_splits <- strsplit(selection_labels, " ", fixed = TRUE)
  labels_top <- sapply(labels_splits, "[", 1)
  labels_bottom <- sapply(labels_splits, "[", 2)

  mtext(sapply(labels_top, VerticalAdjust),
        at   = seq_len(num_groups),
        side = 1,
        line = x_labels_line,
        cex  = par("cex")
        )

  mtext(ifelse(is.na(labels_bottom),
               NA,
               sapply(labels_bottom, VerticalAdjust)
               ),
        at   = seq_len(num_groups),
        side = 1,
        line = x_labels_line + 1,
        cex  = par("cex")
        )

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
    png(filename = file.path(PNGs_output_directory, "Theoretical", paste0(use_file_name, ".png")),
        width    = use_width,
        height   = use_height,
        units    = "in",
        res      = 600
        )
    TheoreticalAtLeastCounts()
    dev.off()
  }

  ccs_numbers <- c(3, 5, 7)
  accuracy_percentages <- c(99, 99.9, 99.99)

  df_list_names <- paste0("ccs", ccs_numbers, "_df_list")
  ccs_are_present <- df_list_names %in% ls(envir = globalenv())
  ccs_numbers <- ccs_numbers[ccs_are_present]
  accuracy_percentages <- accuracy_percentages[ccs_are_present]

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

        filter_stages <- c("original_summary_df", "filtered_summary_df", "filtered_cross_plate_df")
        filter_labels <- c("i) unfiltered", "ii) filtered", "iii) filtered cross-plate")
        df_are_present <- filter_stages %in% names(use_df_list)
        filter_stages <- filter_stages[df_are_present]
        filter_labels <- filter_labels[df_are_present]

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
                  png(filename = file.path(metric_path, file_name),
                      width    = use_width,
                      height   = use_height,
                      units    = "in",
                      res      = 600
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



DrawAllSummaryBarPlots <- function() {

  ccs_numbers <- c(3, 5, 7)
  accuracy_percentages <- c(99, 99.9, 99.99)

  df_list_names <- paste0("ccs", ccs_numbers, "_df_list")
  ccs_are_present <- df_list_names %in% ls(envir = globalenv())
  ccs_numbers <- ccs_numbers[ccs_are_present]
  accuracy_percentages <- accuracy_percentages[ccs_are_present]

  for (i in seq_along(ccs_numbers)) {

    use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))

    filter_stages <- c("original_summary_df", "filtered_summary_df", "filtered_cross_plate_df")
    filter_labels <- c("i) unfiltered", "ii) filtered", "iii) filtered cross-plate")
    df_are_present <- filter_stages %in% names(use_df_list)
    filter_stages <- filter_stages[df_are_present]
    filter_labels <- filter_labels[df_are_present]

    for (filter_stage in seq_along(filter_stages)) {

      df_name <- filter_stages[[filter_stage]] # "filtered_gRNAs_df"
      use_df <- use_df_list[[df_name]]
      folder_name <- paste0("CCS", ccs_numbers[[i]],
                            " (", accuracy_percentages[[i]], ") - ",
                            filter_labels[[filter_stage]]
                            )

      message(paste0("Exporting summary stacked bar plots for the following data: ",
                     folder_name, "..."
                     )
              )
      folder_path <- file.path(file_output_directory, folder_name)
      dir.create(folder_path, showWarnings = FALSE)

      sub_folder_path <- file.path(folder_path, "Summary bar plots")
      dir.create(sub_folder_path, showWarnings = FALSE)

      all_selections <- union(names(plate_selection_list),
                              names(plate_selection_titles_list)
                              )

      for (include_tracrRNAs in c(TRUE, FALSE)) {

        file_name <- paste0("Summary bar plots - ",
                            if (include_tracrRNAs) "including tracrRNA" else "only protospacer",
                            ".pdf"
                            )

        pdf(file = file.path(sub_folder_path, file_name),
            width = use_width, height = use_height
            )
        for (k in seq_along(all_selections)) {
          SummaryStackedBars(use_df, plate_names = all_selections[[k]],
                             consider_tracrRNAs = include_tracrRNAs
                             )
        }
        dev.off()
      }
    }
  }
  return(invisible(NULL))
}



