### 8th June 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("vioplot")

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions





# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output", "Figures", "Summaries across wells")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))
load(file.path(sql2_R_objects_directory, "17) Check for cross-plate contaminations.RData"))
load(file.path(sql2_R_objects_directory, "19) Identify and characterize large deletions.RData"))





# Define constants --------------------------------------------------------

column_labels_list <- list(
  "Count_at_least_1"                            = expression("" >= "1 correct", "gRNA"),
  "Count_at_least_2"                            = expression("" >= "2 correct", "gRNAs"),
  "Count_at_least_3"                            = expression("" >= "3 correct", "gRNAs"),
  "Count_all_4"                                 = expression("All 4 gRNAs", "are correct"),

  "Count_no_contam_at_least_1"                  = expression("" >= "1 correct", "gRNA"),
  "Count_no_contam_at_least_2"                  = expression("" >= "2 correct", "gRNAs"),
  "Count_no_contam_at_least_3"                  = expression("" >= "3 correct", "gRNAs"),
  "Count_no_contam_all_4"                       = expression("All 4 gRNAs", "are correct"),

  "Correct_sg1_cr1"                             = expression("sg1 (+" * scriptscriptstyle(" ") * "crRNA)", "is correct"),
  "Correct_sg2_cr2"                             = expression("sg2 (+" * scriptscriptstyle(" ") * "crRNA)", "is correct"),
  "Correct_sg3_cr3"                             = expression("sg3 (+" * scriptscriptstyle(" ") * "crRNA)", "is correct"),
  "Correct_sg4_cr4"                             = expression("sg4 (+" * scriptscriptstyle(" ") * "crRNA)", "is correct"),

  "Correct_sg1"                                 = expression("sg1 is", "correct"),
  "Correct_sg2"                                 = expression("sg2 is", "correct"),
  "Correct_sg3"                                 = expression("sg3 is", "correct"),
  "Correct_sg4"                                 = expression("sg4 is", "correct"),

  "Num_contaminated_reads"                      = expression("Contamination", "(full read)"),
  "Num_contaminated_reads_aligned"              = expression("Contamination", "(aligned read)"),
  "Num_cross_plate_contaminated"                = expression("Cross-plate", "contamination"),

  "Num_reads_with_deletions_exceeding_20bp"     = expression("Deletions", "" >= "20bp"),
  "Num_reads_with_deletions_spanning_tracrRNAs" = expression("Deletions that", "span tracRNAs"),
  "Num_reads_with_deletions_spanning_promoters" = expression("Deletions that", "span promoters"),
  "Num_reads_with_deletions_spanning_sg_cr"     = expression("Deletions that", "span sg" * scriptstyle(" ") * "+" * scriptstyle(" ") * "cr"),
  "Num_reads_with_deletions_spanning_sgRNAs"    = expression("Deletions that", "span sgRNAs")
)


column_groups_list <- list(
  "At_least_num_guides"  = c(paste0("Count_at_least_", 1:3), "Count_all_4"),
  "No_contam_num_guides" = paste0("Count_no_contam_", c(paste0("at_least_", 1:3), "all_4")),
  "Correct_sg_cr"        = paste0("Correct_sg", 1:4, "_cr", 1:4),
  "Deletions"            = c("Num_reads_with_deletions_exceeding_20bp", "Num_reads_with_deletions_spanning_tracrRNAs", "Num_reads_with_deletions_spanning_promoters"),
  "Contaminations"       = c("Num_contaminated_reads", "Num_contaminated_reads_aligned", "Num_cross_plate_contaminated")
)


plate_selection_titles_list <- list(
  "All plates"              = "All plates (purified using columns)",
  "Bead-purified"           = "Plates HA-11 and HO-1 (purified using beads)",
  "Matched column-purified" = "Plates HA-11 and HO-1 (purified using columns)",
  "Library plates"          = "4sg library plates (purified using columns)"
)


count_columns <- c(
  paste0("Count_at_least_", 1:3),
  "Count_all_4", "Count_all_4_promoters", "Count_whole_plasmid",
  paste0("Count_sg", 1:4, "_cr", 1:4)
)


no_contam_count_columns <- sub("Count_", "Count_no_contam_", count_columns, fixed = TRUE)

num_columns <- c(
  "Num_contaminated_reads",
  "Num_contaminated_reads_aligned",
  "Num_cross_plate_contaminated",
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






# Define functions for aggregating and tidying data -----------------------

StandardizeCounts <- function(input_df) {
  stopifnot(all(c("count_columns", "alignment_columns") %in% ls(envir = globalenv())))

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


AddContaminationSummary <- function(ccs_df_list, contamin_df) {

  pass_filters <- ccs_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
  use_zmws <- ccs_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]
  are_to_use <- contamin_df[["ZMW"]] %in% use_zmws

  cross_well_contam_df  <- contamin_df[are_to_use & !(contamin_df[["Are_cross_plate"]]), ]
  cross_plate_contam_df <- contamin_df[are_to_use & contamin_df[["Are_cross_plate"]], ]


  summary_df <- ccs_df_list[["filtered_summary_df"]]

  cross_well_contam_df[["Reference_ID"]]  <- factor(cross_well_contam_df[["Reference_ID"]],
                                                    levels = summary_df[["Combined_ID"]]
                                                    )
  cross_plate_contam_df[["Reference_ID"]] <- factor(cross_plate_contam_df[["Reference_ID"]],
                                                    levels = summary_df[["Combined_ID"]]
                                                    )

  well_contam_splits  <- split(cross_well_contam_df[["ZMW"]], cross_well_contam_df[["Reference_ID"]])
  plate_contam_splits <- split(cross_plate_contam_df[["ZMW"]], cross_plate_contam_df[["Reference_ID"]])

  well_contam_splits <- lapply(well_contam_splits, unique)
  plate_contam_splits <- lapply(plate_contam_splits, unique)

  summary_df[["Num_contaminated_reads_aligned"]] <- lengths(well_contam_splits)
  summary_df[["Num_cross_plate_contaminated"]] <- lengths(plate_contam_splits)


  for (column_name in c("Num_contaminated_reads_aligned", "Num_cross_plate_contaminated")) {
    summary_df[[column_name]] <- ifelse(summary_df[["Count_total"]] == 0,
                                        NA,
                                        summary_df[[column_name]]
                                        )
  }

  ccs_df_list[["filtered_summary_df"]] <- summary_df
  return(ccs_df_list)
}




AddDeletionSummary <- function(ccs_df_list, del_df) {

  pass_filters <- ccs_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
  use_zmws <- ccs_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]
  are_to_use <- del_df[["ZMW"]] %in% use_zmws

  del_df  <- del_df[are_to_use, ]

  summary_df <- ccs_df_list[["filtered_summary_df"]]

  del_df[["Combined_ID"]] <- factor(del_df[["Combined_ID"]],
                                    levels = summary_df[["Combined_ID"]]
                                    )

  del_splits <- split(del_df[["ZMW"]], del_df[["Combined_ID"]])
  del_splits <- lapply(del_splits, unique)
  summary_df[["Num_reads_with_deletions_exceeding_20bp"]] <- lengths(del_splits)

  span_columns <- c("Span_tracrRNAs", "Span_promoters", "Span_sg_cr", "Span_sgRNAs")
  span_list <- lapply(span_columns, function(x) {
    are_spanning <- del_df[[x]]
    span_splits <- split(del_df[["ZMW"]][are_spanning],
                         del_df[["Combined_ID"]][are_spanning]
                         )
    span_splits <- lapply(span_splits, unique)
    return(lengths(span_splits))
  })
  span_mat <- do.call(cbind, span_list)
  colnames(span_mat) <- sub("Span_", "Num_reads_with_deletions_spanning_", span_columns, fixed = TRUE)
  summary_df <- data.frame(summary_df, span_mat, stringsAsFactors = FALSE)
  ccs_df_list[["filtered_summary_df"]] <- summary_df
  return(ccs_df_list)
}



AddNoContamCounts <- function(ccs_df_list, contamin_df) {

  indiv_df <- ccs_df_list[["individual_reads_df"]]
  summary_df <- ccs_df_list[["filtered_summary_df"]]

  are_contaminated <- (indiv_df[["ZMW"]] %in% contamin_df[["ZMW"]]) |
                      ((indiv_df[["Num_contaminating_guides"]] >= 1) %in% TRUE)

  pass_filters <- indiv_df[["Passes_filters"]] == 1
  no_contam_df <- indiv_df[(!(are_contaminated)) & pass_filters, ]

  binary_columns <- c(paste0("sg", 1:4, "_cr", 1:4),
                      paste0("at_least_", 1:3),
                      "all_4", "all_4_promoters", "whole_plasmid"
                      )
  no_contam_mat <- as.matrix(no_contam_df[, binary_columns])

  wells_fac <- factor(no_contam_df[["Combined_ID"]],
                      levels = summary_df[["Combined_ID"]]
                      )

  summary_vec_list <- tapply(seq_len(nrow(no_contam_mat)),
                             wells_fac,
                             function(x) colSums(no_contam_mat[x, , drop = FALSE])
                             )
  are_empty <- vapply(summary_vec_list, is.null, logical(1))
  if (any(are_empty)) {
    summary_vec_list[are_empty] <- list(integer(ncol(no_contam_mat)))
  }

  summary_counts_mat <- do.call(rbind, summary_vec_list)
  mode(summary_counts_mat) <- "integer"

  colnames(summary_counts_mat) <- paste0("Count_no_contam_", colnames(summary_counts_mat))
  summary_df <- data.frame(summary_df,
                           "Count_total_no_contam" = tabulate(wells_fac),
                           summary_counts_mat,
                           stringsAsFactors = FALSE
                           )
  ccs_df_list[["filtered_summary_df"]] <- summary_df

  assign("delete_summary_df", summary_df, envir = globalenv())
  return(ccs_df_list)
}




# Define functions for plotting summary statistics ------------------------

SetUpPercentagePlot <- function(use_columns,
                                use_y_limits,
                                main_title,
                                point_cex        = 1.5,
                                legend_pch       = 21,
                                side_gap         = 0.5,
                                extra_grid_lines = FALSE
                                ) {

  stopifnot("column_labels_list" %in% ls(envir = globalenv()))

  are_no_contam <- grepl("_no_contam_", use_columns, fixed = TRUE)
  contams_excluded <- any(are_no_contam)

  if (contams_excluded) {
    stopifnot(all(are_no_contam))
    y_axis_label <- "Percentage of reads (excluding contaminations)"
  } else {
    y_axis_label <- "Percentage of reads"
  }

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


  ## Set up the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )

  ## Draw the grid and axis
  segments(x0   = par("usr")[[1]],
           x1   = par("usr")[[2]],
           y0   = seq(0, 1, by = 0.1),
           col  = "gray88",
           lend = "butt",
           xpd  = NA
           )

  if (extra_grid_lines) {
    segments(x0   = par("usr")[[1]],
             x1   = par("usr")[[2]],
             y0   = seq(0.05, 0.95, by = 0.1),
             col  = "gray95",
             lend = "butt",
             xpd  = NA
             )
  }

  axis(2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, 0.38, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1,
       lwd      = par("lwd")
       )

  mtext(text = y_axis_label, side = 2, line = 3)


  ## Draw the title
  title(main_title, cex.main = 1.1, font.main = 1, line = 2)


  ## Draw the x axis labels
  large_gap_lines <- 2.5
  start_line <- 1.25

  for (i in seq_along(group_positions)) {
    mtext(text = VerticalAdjust(column_labels_list[[use_columns[[i]]]][[1]]),
          at   = group_positions[[i]],
          line = 1.25,
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
  large_gap <- diff(grconvertY(c(0, large_gap_lines), from = "lines", to = "user"))
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
  x_start <- 1 + diff(grconvertX(c(0, 1.5), from = "lines", to = "npc"))

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

  return(invisible(NULL))
}




GetColorMat <- function() {

  points_alpha <- 0.8
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

  fill_A    <- brewer.pal(9, "Purples")[[5]]
  fill_B    <- brewer.pal(9, "Blues")[[5]]
  outline_A <- brewer.pal(9, "Purples")[[9]]
  outline_B <- brewer.pal(9, "Blues")[[9]]

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
                         set_par           = TRUE,
                         label_percentages = TRUE,
                         use_y_limits      = c(0, 1)
                         ) {

  require_objects <- c("plate_selection_list", "plate_selection_titles_list")

  stopifnot(all(require_objects %in% ls(envir = globalenv())))

  if (is.null(plate_names)) {
    plate_names <- "All plates"
  }
  if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
    plate_names <- plate_selection_list[[plate_names]]
  } else {
    main_title <- "PacBio sequencing"
  }

  if (set_par) {
    use_mar <- par("mar" = c(5, 5, 4.5, 8))
  }

  control_mat <- PreparePlates(input_df, "Intctrl")
  selected_mat <- PreparePlates(input_df, plate_names)

  control_metrics  <- colMeans(control_mat[, use_columns,  drop = FALSE], na.rm = TRUE)
  selected_metrics <- colMeans(selected_mat[, use_columns, drop = FALSE], na.rm = TRUE)

  point_cex <- 1.5
  SetUpPercentagePlot(use_columns, use_y_limits, main_title, point_cex)

  spaced_percent <- 2.5
  max_space <- (use_y_limits[[2]] - use_y_limits[[1]]) * (spaced_percent / 100)
  are_spaced <- abs(control_metrics - selected_metrics) > max_space

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

  points(x   = group_positions,
         y   = control_metrics,
         cex = point_cex,
         pch = 21,
         col = ifelse(are_spaced, colors_mat[1, "opaque_outline"], colors_mat[1, "transparent_outline"]),
         bg  = ifelse(are_spaced, colors_mat[1, "opaque_fill"], colors_mat[1, "transparent_fill"]),
         xpd = NA
         )

  points(x   = group_positions,
         y   = selected_metrics,
         cex = point_cex,
         pch = 21,
         col = ifelse(are_spaced, colors_mat[2, "opaque_outline"], colors_mat[2, "transparent_outline"]),
         bg  = ifelse(are_spaced, colors_mat[2, "opaque_fill"], colors_mat[2, "transparent_fill"]),
         xpd = NA
         )

  if (label_percentages) {
    perc_x_gap <- diff(grconvertX(c(0, 1.4), from = "lines", to = "user"))
    text(x      = group_positions + perc_x_gap,
         y      = control_metrics,
         labels = paste0(format(round(control_metrics * 100, digits = 1)), "%"),
         cex    = 0.7,
         pch    = 16,
         col    = "gray70",
         xpd    = NA
         )
    text(x      = group_positions + perc_x_gap,
         y      = selected_metrics,
         labels = paste0(format(round(selected_metrics * 100, digits = 1)), "%"),
         cex    = 0.7,
         pch    = 16,
         col    = "gray70",
         xpd    = NA
         )
  }

  par(use_mar)
  return(invisible(NULL))
}



SummaryBoxPlot <- function(input_df,
                           plate_names       = NULL,
                           use_columns       = c("Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4"),
                           set_par           = TRUE,
                           label_percentages = TRUE,
                           use_y_limits      = c(0, 1)
                           ) {

  set.seed(1) # For reproducible jitter

  require_objects <- c("plate_selection_list", "plate_selection_titles_list")
  stopifnot(all(require_objects %in% ls(envir = globalenv())))

  if (is.null(plate_names)) {
    plate_names <- "All plates"
  }
  if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
    plate_names <- plate_selection_list[[plate_names]]
  } else {
    main_title <- "PacBio sequencing"
  }

  if (set_par) {
    use_mar <- par("mar" = c(5, 5, 4.5, 8))
  }

  control_mat <- PreparePlates(input_df, "Intctrl")
  selected_mat <- PreparePlates(input_df, plate_names)

  control_list <- lapply(use_columns, function(x) control_mat[, x])
  selected_list <- lapply(use_columns, function(x) selected_mat[, x])

  control_unlisted <- unlist(control_list)
  selected_unlisted <- unlist(selected_list)

  point_cex <- 1.5
  SetUpPercentagePlot(use_columns, use_y_limits, main_title, point_cex,
                      side_gap = 0.3, legend_pch = 22, extra_grid_lines = TRUE
                      )

  colors_mat <- GetColorMat()
  num_groups <- length(use_columns)
  group_positions <- seq_len(num_groups)

  use_wex <- 0.4

  use_gap <- 0.35
  control_pos  <- group_positions - (use_gap / 2)
  selected_pos <- group_positions + (use_gap / 2)

  RemoveAllZeros <- function(input_list) {
    lapply(input_list, function(x) {
      if (all(x == 0)) {
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

  ## Draw the superimposed boxplots
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

  par(use_mar)
  return(invisible(NULL))
}





# Define plate selections -------------------------------------------------

are_beads <- grepl("-beads", plates_df[["Plate_name"]], fixed = TRUE)
are_controls <- plates_df[["Plate_name"]] == "Intctrl"

non_library_plates <- c("Vac-1", "PD_A", "PD_O")
are_columns <- !(are_beads | are_controls)

are_library <- are_columns & (!(plates_df[["Plate_name"]] %in% non_library_plates))

column_matched_plates <- sub("-beads", "", plates_df[["Plate_name"]][are_beads], fixed = TRUE)

plate_selection_list <- list(
  "All plates"              = plates_df[["Plate_name"]][are_columns],
  "Colony-picked"           = "Intctrl",
  "Bead-purified"           = plates_df[["Plate_name"]][are_beads],
  "Matched column-purified" = column_matched_plates,
  "Library plates"          = plates_df[["Plate_name"]][are_library]
)




# Explore stuff -----------------------------------------------------------

use_df <- ccs7_df_list[["filtered_summary_df"]]

are_bigger <- use_df[["Num_contaminated_reads_aligned"]] >
              use_df[["Num_contaminated_reads"]]


indiv_df <- ccs7_df_list[["individual_reads_df"]]

contam_matches <- match(contam_df[["ZMW"]], indiv_df[["ZMW"]])

contam_df[["Num_contaminating_guides"]] <- indiv_df[["Num_contaminating_guides"]][contam_matches]

are_bigger <- (contam_df[["Num_contaminating_guides"]] %in% 0) &
              (!(contam_df[["Are_cross_plate"]]))

head(contam_df[are_bigger, ])





# Export lollipop plots ---------------------------------------------------

ccs_numbers <- c(3, 5, 7)
accuracy_percentages <- c(99, 99.9, 99.99)


for (plot_type in c("Box", "Lollipop")) {

  if (plot_type == "Box") {
    UseFunction <- SummaryBoxPlot
  } else {
    UseFunction <- LollipopPlot
  }

  for (i in seq_along(ccs_numbers)) {

    use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))
    use_df_list <- AddContaminationSummary(use_df_list, contam_df)
    use_df_list <- AddDeletionSummary(use_df_list, deletions_df)
    use_df_list <- AddNoContamCounts(use_df_list, contam_df)

    for (filter_stage in 2) {

      df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"

      use_df <- use_df_list[[df_name]]

      folder_name <- paste0("CCS", ccs_numbers[[i]],
                            " (", accuracy_percentages[[i]], ") - ",
                            c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                            )
      folder_path <- file.path(file_output_directory, folder_name)
      dir.create(folder_path, showWarnings = FALSE)

      for (label_percentages in c(TRUE, FALSE)) {

        if (label_percentages && (plot_type == "Box")) {
          next
        }

        for (i in seq_along(column_groups_list)) {
          metric <- names(column_groups_list)[[i]]
          file_name <- paste0(plot_type, " plot - ", i, ") ", metric, ".pdf")
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
          pdf(file = file.path(sub_folder_path, file_name),
              width = 8, height = 6
              )
          for (plate_selection in names(plate_selection_titles_list)) {
            print(plate_selection)
            UseFunction(use_df,
                        plate_names = plate_selection,
                        use_columns = column_groups_list[[metric]],
                        label_percentages = label_percentages
                        )
          }
          dev.off()
        }
      }
    }
  }
}





# Try stuff ---------------------------------------------------------------

SummaryBoxPlot(use_df, "All plates")



LollipopPlot(use_df, "All plates")
LollipopPlot(use_df, "Bead-purified")




LollipopPlot(use_df,
             "All plates",
             paste0("Correct_sg", 1:4)
             )

LollipopPlot(use_df,
             "Bead-purified",
             paste0("Correct_sg", 1:4)
             )



LollipopPlot(use_df,
             "All plates",
             paste0("Correct_sg", 1:4, "_cr", 1:4)
             )

LollipopPlot(use_df,
             "Bead-purified",
             paste0("Correct_sg", 1:4, "_cr", 1:4)
             )





LollipopPlot(use_df,
             "Bead-purified",
             c("Num_contaminated_reads",
               "Num_contaminated_reads_aligned",
               "Num_cross_plate_contaminated"
               ),
             use_y_limits = c(0, 0.05)
             )






LollipopPlot(use_df,
             "Bead-purified",
             c("Num_contaminated_reads"),
             use_y_limits = c(0, 0.05)
             )



LollipopPlot(use_df,
             "All plates",
             c("Num_contaminated_reads",
               "Num_contaminated_reads_aligned",
               "Num_cross_plate_contaminated"
               )
             )


LollipopPlot(use_df,
             "All plates",
             c("Num_reads_with_deletions_exceeding_20bp",
               "Num_reads_with_deletions_spanning_tracrRNAs",
               "Num_reads_with_deletions_spanning_promoters"
               )
             )





