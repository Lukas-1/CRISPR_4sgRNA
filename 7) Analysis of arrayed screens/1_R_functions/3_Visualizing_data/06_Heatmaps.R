# 2021-12-27


# Load packages and source code -------------------------------------------

library("squash")
library("viridis")
library("RColorBrewer")



# Define functions --------------------------------------------------------

InvertTranspose <- function(input_mat) {
   t(apply(input_mat, 2, rev))
}

parula <- colorRampPalette(c("#352A87", "#0F5CDD", "#1481D6", "#06A4CA",
                             "#2EB7A4", "#87BF77", "#D1BB59", "#FEC832",
                             "#F9FB0E"
                             ))


MakeBreaks <- function(numeric_vec, num_breaks, trim = TRUE) {
  if (trim) {
    use_range <- quantile(numeric_vec, probs = c(0.02, 0.98), na.rm = TRUE)
  } else {
    use_range <- range(numeric_vec, na.rm = TRUE)
  }
  breaks_vec <- seq(from       = use_range[[1]],
                    to         = use_range[[2]],
                    length.out = num_breaks
                    )
  return(breaks_vec)
}


SigmoidalBreakPoints <- function(numeric_vec, use_for_log = 5, symm = TRUE) {

  numeric_vec <- as.numeric(numeric_vec)

  use_color_breaks <- exp(seq(log(use_for_log), log(100), length.out = 20)) - use_for_log
  use_color_breaks <- sort(unique(c(-use_color_breaks, use_color_breaks)))

  my_range <- range(numeric_vec, na.rm = TRUE)
  if (symm) {
    my_range <- c(-max(abs(my_range)), max(abs(my_range)))
  }
  my_range[1] <- my_range[1] - (diff(my_range) / 1000)
  my_range[2] <- my_range[2] + (diff(my_range) / 1000)

  use_color_breaks <- scales::rescale(use_color_breaks, to = my_range)
  use_color_breaks <- unique(round(use_color_breaks, digits = 3))
  return(use_color_breaks)
}



GetRepNumber <- function(column_name) {
  if (grepl("_rep", column_name)) {
    column_name_splits <- strsplit(column_name, "_", fixed = TRUE)[[1]]
    rep_string <- grep("^rep", column_name_splits, value = TRUE)
    rep_number <- as.integer(sub("rep", "", rep_string, fixed = TRUE))
    stopifnot(rep_number %in% 1:2)
  } else {
    rep_number <- NULL
  }
  return(rep_number)
}


BothRepColumns <- function(column_name) {
  rep_number <- GetRepNumber(column_name)
  if (rep_number == 1) {
    rep1_column <- column_name
    rep2_column <- sub("rep1", "rep2", column_name, fixed = TRUE)
  } else if (rep_number == 2) {
    rep1_column <- sub("rep2", "rep1", column_name, fixed = TRUE)
    rep2_column <- column_name
  }
  return(c(rep1_column, rep2_column))
}



BreaksForColumn <- function(input_df,
                            use_column,
                            num_empty_breaks          = 7,
                            num_other_breaks          = 50,
                            num_pos_breaks            = 10,
                            num_empty_to_other_breaks = 6,
                            num_other_to_pos_breaks   = 8,
                            use_custom_breaks         = NULL,
                            num_uniform_breaks        = NULL,
                            flatten_factor            = NULL,
                            weighting_for_controls    = NULL,
                            take_log2                 = FALSE,
                            use_vector                = NULL
                            ) {


  if (!(is.null(use_custom_breaks))) {
    return(list(breaks = use_custom_breaks, type = "custom"))
  }

  if (is.null(use_vector)) {
    mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
    are_empty <- input_df[, "Well_number_384"] %in% mat_384[, c(1, 24)]
    are_pos_ctrl <- input_df[, "Is_pos_ctrl"]

    rep_number <- GetRepNumber(use_column)

    if (is.null(rep_number)) {
      all_vec <- input_df[, use_column]
    } else {
      rep_columns <- BothRepColumns(use_column)
      all_vec <- c(input_df[, rep_columns[[1]]], input_df[, rep_columns[[2]]])
      are_empty <- rep(are_empty, 2)
      are_pos_ctrl <- rep(are_pos_ctrl, 2)
    }

  } else {
    all_vec <- use_vector
    weighting_for_controls <- FALSE
  }

  if (IsPValue(use_column)) {
    all_vec <- -log10(all_vec)
  } else if (take_log2) {
    all_vec <- log2(all_vec)
  }

  use_diverging <- any(all_vec < 0)
  if (use_diverging) {
    most_extreme_value <- abs(range(all_vec))
    limits_vec <- c(-most_extreme_value, most_extreme_value)
  } else {
    limits_vec <- range(all_vec)
  }

  if (is.null(num_uniform_breaks) && use_diverging) {
    if (is.null(flatten_factor)) {
      if (grepl("SSMD", use_column, fixed = TRUE)) {
        flatten_factor <- 0.5
      } else {
        flatten_factor <- 1
      }
    }
    sigmoidal_breaks <- SigmoidalBreakPoints(limits_vec, use_for_log = flatten_factor)
    return(list(breaks = sigmoidal_breaks, type = "sigmoidal"))
  }

  if (is.null(weighting_for_controls)) {
    weighting_for_controls <- NeedsWeighting(use_column)
  }
  if (!(weighting_for_controls)) {
    num_uniform_breaks <- 100
  }

  if (!(is.null(num_uniform_breaks))) {
    uniform_breaks <- MakeBreaks(limits_vec, trim = FALSE, num_breaks = num_uniform_breaks)
    return(list(breaks = uniform_breaks, type = "uniform"))
  }

  if (any(all_vec < 0)) {
    most_extreme_value <- abs(range(all_vec))
    limits_vec <- c(-most_extreme_value, most_extreme_value)
    uniform_breaks <- SigmoidalBreakPoints(1:1000, use_for_log = 0.5)
    return(uniform_breaks)
  }

  if (!(is.null(use_vector))) {
    stop(paste0("If 'use_vector' is specified, then 'weighting_for_controls' ",
                "must not be TRUE. Please check your input parameters."
                ))
  }

  empty_vec <- all_vec[are_empty]
  pos_vec   <- all_vec[are_pos_ctrl]
  other_vec <- all_vec[!(are_empty | are_pos_ctrl)]

  # stopifnot(max(empty_vec) < min(other_vec))
  # stopifnot(max(other_vec) < min(pos_vec))

  empty_breaks <- MakeBreaks(empty_vec, num_empty_breaks)
  empty_breaks <- c(min(all_vec), empty_breaks)

  other_breaks <- MakeBreaks(other_vec, num_other_breaks)

  pos_breaks <- MakeBreaks(pos_vec, num_pos_breaks)
  pos_breaks <- c(pos_breaks, max(all_vec))

  empty_to_other_breaks <- MakeBreaks(c(max(empty_breaks), min(other_breaks)),
                                      num_breaks = num_empty_to_other_breaks,
                                      trim = FALSE
                                      )

  other_to_pos_breaks   <- MakeBreaks(c(max(other_breaks), min(pos_breaks)),
                                      num_breaks = num_other_to_pos_breaks,
                                      trim = FALSE
                                      )

  all_breaks <- c(empty_breaks,
                  empty_to_other_breaks[2:(length(empty_to_other_breaks) - 1)],
                  other_breaks,
                  other_to_pos_breaks[2:(length(other_to_pos_breaks) - 1)],
                  pos_breaks
                  )
  all_breaks <- sort(unique(all_breaks))

  return(list(breaks = all_breaks, type = "weighted by controls"))
}


IsPValue <- function(column_name) {
  grepl("p_val", column_name, ignore.case = TRUE)
}


NeedsWeighting <- function(column_name) {
  grepl("_rep|Hit_strength", column_name)
}


PrettyRound <- function(numeric_vec) {
  old_scipen <- options("scipen" = 2)
  results_vec <- ifelse(abs(numeric_vec) > 100,
                        round(numeric_vec),
                        ifelse(abs(numeric_vec) > 10,
                               round(numeric_vec, digits = 1),
                               ifelse(abs(numeric_vec) > 1,
                                      round(numeric_vec, digits = 2),
                                      ifelse(abs(numeric_vec) > 0.01,
                                             signif(numeric_vec, digits = 2),
                                             signif(numeric_vec, digits = 1)
                                             )
                                      )
                               )
                        )
  results_vec <- as.character(results_vec)
  options(old_scipen)
  return(results_vec)
}


# lseq <- function(from, to, length.out) {
#   # logarithmic spaced sequence, from https://stackoverflow.com/a/29963530
#   exp(seq(log(from), log(to), length.out = length.out))
# }



PlateSchematic <- function(input_df,
                           plate_number,
                           main_title = NULL,
                           label_genes = FALSE,
                           color_by_source_plate = FALSE
                           ) {

  if (is.numeric(plate_number)) {
    plate_number <- as.character(as.roman(plate_number))
  }
  if (is.null(main_title)) {
    main_title <- paste0("Plate ", plate_number)
  }

  are_this_plate <- input_df[, "Plate_number_384"] == plate_number
  are_NT  <- input_df[, "Is_NT_ctrl"][are_this_plate]
  are_pos <- input_df[, "Is_pos_ctrl"][are_this_plate]
  are_gene <- !(is.na(input_df[, "Entrez_ID"][are_this_plate]))
  are_mCherry <- input_df[, "Target_flag"][are_this_plate] %in% "mCherry"
  if (label_genes) {
    gene_labels <- input_df[, "Gene_symbol"][are_this_plate]
    gene_splits <- strsplit(gene_labels[are_gene], "-", fixed = TRUE)
    was_split <- lengths(gene_splits) > 1
    gene_labels[are_gene][was_split] <- lapply(gene_splits[was_split],
                                               function(x) {
                                                 if (all(nchar(x) >= 4)) {
                                                   paste0(x, collapse = "-\n")
                                                 } else {
                                                   paste0(x, collapse = "-")
                                                 }
                                               })
  }

  type_colors <- c(
    "Empty"         = "gray95",
    "Untransfected" = "gray80",
    "Gene"          = brewer.pal(5, "Blues")[[3]],
    "NT"            = brewer.pal(5, "Blues")[[5]],
    "Pos"           = brewer.pal(5, "Reds")[[4]],
    "mCherry"       = brewer.pal(5, "Reds")[[2]]
  )
  labels_list <- list(
    "Empty"         = c("Empty", "well"),
    "Untransfected" = c("No", "virus"),
    "Gene"          = c("TF",
                        "gene"
                        ),
    "NT"            = c("NT",
                        "control"
                        ),
    "Pos"           = c("Pos.",
                        "control"
                        ),
    "mCherry"       = "mCherry"
  )

  color_vec <- rep(type_colors[["Untransfected"]], 384)
  color_vec[are_NT]   <- type_colors[["NT"]]
  color_vec[are_pos]  <- type_colors[["Pos"]]
  color_vec[are_gene] <- type_colors[["Gene"]]

  if (any(are_mCherry)) {
    color_vec[are_mCherry] <- type_colors[["mCherry"]]
  } else {
    type_colors <- type_colors[names(type_colors) != "mCherry"]
    labels_list <- labels_list[names(labels_list) != "mCherry"]
  }

  if (color_by_source_plate) {
    gene_colors <- brewer.pal(9, "Purples")[3:8]
    type_colors[["Gene"]] <- gene_colors[[4]]
    rounds_plates_384 <- split(1:12, rep(1:4, each = 3))
    all_plate_numbers_384 <- as.integer(as.roman(input_df[, "Plate_number_384"]))
    rounds_plates_96 <- lapply(rounds_plates_384, function(x) {
      results_vec <- unique(input_df[all_plate_numbers_384 %in% x, "Plate_number_96"])
      sort(results_vec[!(is.na(results_vec))])
    })
    colors_plates_96 <- unlist(lapply(rounds_plates_96, function(x) gene_colors[seq_along(x)]), use.names = FALSE)
    plate_numbers_96 <- input_df[, "Plate_number_96"][are_this_plate][are_gene]
    color_vec[are_gene] <- colors_plates_96[plate_numbers_96]
  }

  color_mat <- matrix(color_vec, nrow = 16, ncol = 24, byrow = TRUE)
  color_mat[, c(1, 24)] <- type_colors[["Empty"]]

  ## Plot heatmap
  old_mai <- par(mai = c(1.3, 0.7, 1.1, 0.7))
  cimage(zcol = InvertTranspose(color_mat),
         xlab    = "",
         ylab    = "",
         xlabels = rep("", ncol(color_mat)),
         ylabels = rep("", nrow(color_mat)),
         tcl     = 0,
         mgp     = c(3, 0.3, 0),
         bty     = "n",
         axes    = FALSE
         )
  horizontal_lines <- seq(from = 0.5, to = nrow(color_mat) + 0.5, by = 1)
  vertical_lines   <- seq(from = 0.5, to = ncol(color_mat) + 0.5, by = 1)
  grid_color <- "gray50"
  use_lwd <- 0.8
  segments(x0 = 0.5, x1 = ncol(color_mat) + 0.5, y0 = horizontal_lines,
           col = grid_color, lwd = use_lwd, xpd = NA
           )
  segments(x0 = vertical_lines, y0 = 0.5, y1 = nrow(color_mat) + 0.5,
           col = grid_color, lwd = use_lwd, xpd = NA
           )

  ## Label with gene names
  if (label_genes) {
    text(x      = rep(1:24, times = 16),
         y      = rev(rep(1:16, each = 24)),
         labels = gene_labels,
         cex    = 0.3,
         font   = 4,
         xpd    = NA,
         )
  }

  ## Label with row and column names
  text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 0.8), from = "lines", to = "user")),
       y      = seq_len(nrow(color_mat)),
       labels = rev(LETTERS[seq_len(nrow(color_mat))]),
       cex    = par("cex") * 0.8,
       col    = "black",
       xpd    = NA
       )
  text(x      = seq_len(ncol(color_mat)),
       y      = par("usr")[[4]] + diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
       labels = seq_len(ncol(color_mat)),
       cex    = par("cex") * 0.8,
       col    = "black",
       xpd    = NA
       )


  ## Draw the legend
  start_y <- par("usr")[[3]] - diff(grconvertY(c(0, 1.5), from = "lines", to = "user"))
  x_space <- diff(grconvertX(c(0, 2.8), from = "lines", to = "user"))
  x_positions <- seq(0, x_space * length(type_colors), length.out = length(type_colors))
  x_positions <- x_positions - (mean(x_positions) - grconvertX(0.5, from = "npc", to = "user"))
  rect(xleft   = x_positions - 0.5,
       xright  = x_positions + 0.5,
       ybottom = start_y - 1,
       ytop    = start_y,
       col     = type_colors,
       border  = grid_color,
       xpd     = NA
       )

  text(x      = x_positions,
       y      = start_y - 1 - diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
       labels = sapply(labels_list, function(x) VerticalAdjust(x[[1]])),
       cex    = 0.8,
       xpd    = NA
       )
  second_lines <- sapply(labels_list, function(x) {
    if (length(x) == 2) {
      VerticalAdjust(x[[2]])
    } else {
      NA
    }})
  text(x      = x_positions,
       y      = start_y - 1 - diff(grconvertY(c(0, 1.65), from = "lines", to = "user")),
       labels = second_lines,
       cex    = 0.8,
       xpd    = NA
       )
  title(main_title, line = 2.9)

  par(old_mai)
  return(invisible(NULL))
}



HeatMap384 <- function(numeric_vec,
                       use_breaks,
                       ColorFunction  = NULL,
                       main_title     = "",
                       use_subtext    = "",
                       label_values   = FALSE,
                       use_minuslog10 = FALSE,
                       take_log2      = FALSE,
                       uniform_legend = FALSE,
                       legend_cex     = 0.8,
                       well_label_cex = 0.4
                       ) {

  stopifnot(length(numeric_vec) == 384)

  ## Prepare data for plotting
  if (use_minuslog10) {
    final_vec <- -log10(numeric_vec)
  } else if (take_log2) {
    final_vec <- log2(numeric_vec)
  } else {
    final_vec <- numeric_vec
  }
  use_mat <- matrix(final_vec, nrow = 16, ncol = 24, byrow = TRUE)

  ## Create a matrix of colors for the heatmap
  if (is.null(ColorFunction)) {
    if (any(use_breaks < 0)) {
      ColorFunction <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
    } else {
      ColorFunction <- cividis
    }
  }
  use_mat <- matrix(final_vec, nrow = 16, ncol = 24, byrow = TRUE)
  use_cmap <- makecmap(use_mat, breaks = use_breaks, colFn = ColorFunction,
                       include.lowest = TRUE
                       )
  color_mat <- cmap(use_mat, map = use_cmap)

  ## Create a vector of colors for the legend
  if (uniform_legend) {
    num_steps <- 200
    legend_colors_seq <- seq(from       = min(use_breaks),
                             to         = max(use_breaks),
                             length.out = num_steps
                             )
    legend_colors_vec <- cmap(legend_colors_seq, map = use_cmap)
  } else {
    breaks_mat <- cbind(use_breaks[seq_len(length(use_breaks) - 1)],
                        use_breaks[seq_len(length(use_breaks) - 1) + 1L]
                        )
    legend_values_vec <- rowMeans(breaks_mat)
    legend_colors_vec <- use_cmap[["colors"]]
    num_steps <- length(legend_values_vec)
    stopifnot(length(legend_colors_vec) == num_steps)
  }

  ## Plot heatmap
  old_mai <- par(mai = c(1.3, 0.7, 1.1, 0.7))
  cimage(zcol = InvertTranspose(color_mat),
         xlab    = "",
         ylab    = "",
         xlabels = rep("", ncol(color_mat)),
         ylabels = rep("", nrow(color_mat)),
         tcl     = 0,
         mgp     = c(3, 0.3, 0),
         bty     = "n",
         axes    = FALSE
         )
  horizontal_lines <- seq(from = 0.5, to = nrow(use_mat) + 0.5, by = 1)
  vertical_lines   <- seq(from = 0.5, to = ncol(use_mat) + 0.5, by = 1)
  grid_color <- "white"
  use_lwd <- 0.2
  segments(x0 = 0.5, x1 = ncol(use_mat) + 0.5, y0 = horizontal_lines,
           col = grid_color, lwd = use_lwd, xpd = NA
           )
  segments(x0 = vertical_lines, y0 = 0.5, y1 = nrow(use_mat) + 0.5,
           col = grid_color, lwd = use_lwd, xpd = NA
           )

  ## Plot heatmap legend
  step_width <- (par("usr")[[2]] - par("usr")[[1]]) / num_steps
  start_y <- par("usr")[[3]] - diff(grconvertY(c(0, 2), from = "lines", to = "user"))
  legend_height <- diff(grconvertX(c(0, 0.8), from = "lines", to = "user"))
  step_left_vec <- seq(from = par("usr")[[1]],
                       to   = par("usr")[[2]] - step_width,
                       by   = step_width
                       )
  rect(xleft    = step_left_vec,
       xright   = step_left_vec + step_width,
       ybottom  = start_y,
       ytop     = start_y + legend_height,
       col      = legend_colors_vec,
       border   = NA,
       xpd      = NA
       )
  ## Prepare the axis labels for the heatmap legend
  if (uniform_legend) {
    legend_text_values <- pretty(legend_colors_seq, n = 10)
    use_tolerance <- diff(range(legend_colors_seq)) * 0.003
    are_within_limits <- (legend_text_values >= (min(legend_colors_seq) - use_tolerance)) &
                         (legend_text_values <= (max(legend_colors_seq) + use_tolerance))

    legend_text_values <- legend_text_values[are_within_limits]
    legend_text_vec <- format(legend_text_values, trim = TRUE)

    legend_text_range <- c(par("usr")[[1]] + (step_width / 2),
                           par("usr")[[2]] - (step_width / 2)
                           )
    legend_text_positions <- scales::rescale(legend_text_values,
                                             from = range(use_breaks),
                                             to   = legend_text_range
                                             )
  } else {
    use_indices <- seq(from = 1,
                       to = num_steps,
                       by = round((num_steps / 10) - 1)
                       )
    legend_text_values <- legend_values_vec[use_indices]
    legend_text_vec <- PrettyRound(legend_text_values)
    legend_text_positions <- (step_left_vec + (step_width / 2))[use_indices]
  }
  ## Draw the axis labels and tick marks for the heatmap legend
  text(x      = legend_text_positions,
       y      = start_y - diff(grconvertY(c(0, 1.0), from = "lines", to = "user")),
       labels = legend_text_vec,
       adj    = c(0.5, 0),
       cex    = legend_cex,
       xpd    = NA
       )
  segments(x0  = legend_text_positions,
           y0  = start_y,
           y1  = start_y - diff(grconvertY(c(0, 0.3), from = "lines", to = "user")),
           xpd = NA
           )
  ## Draw a box around the heatmap legend
  rect(xleft   = par("usr")[[1]],
       xright  = par("usr")[[2]],
       ybottom = start_y,
       ytop    = start_y + legend_height,
       xpd     = NA
       )

  ## Annotate the plot with row and column labels, title and subtext
  text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 0.8), from = "lines", to = "user")),
       y      = seq_len(nrow(use_mat)),
       labels = rev(LETTERS[seq_len(nrow(use_mat))]),
       cex    = par("cex") * 0.8,
       col    = "black",
       xpd    = NA
       )
  text(x      = seq_len(ncol(use_mat)),
       y      = par("usr")[[4]] + diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
       labels = seq_len(ncol(use_mat)),
       cex    = par("cex") * 0.8,
       col    = "black",
       xpd    = NA
       )
  title(main_title, line = 2.9)
  if (is.expression(use_subtext)) {
    use_subtext <- VerticalAdjust(use_subtext)
  } else {
    if (use_minuslog10) {
      use_subtext <- sub("p value", "-log10 p value", use_subtext, fixed = TRUE)
    }
    use_subtext <- FormatPlotMath(use_subtext)
  }
  mtext(use_subtext, side = 1, line = 4.2)

  ## Label wells with their numerical values
  if (label_values) {
    x_positions <- rep(1:24, times = 16)
    y_positions <- rep(16:1, each = 24)
    use_dark_text <- colMeans(col2rgb(as.character(t(color_mat))) / 255) > 0.5
    text(x      = x_positions,
         y      = y_positions,
         labels = PrettyRound(numeric_vec),
         col    = ifelse(use_dark_text, "gray25", "gray75"),
         cex    = well_label_cex,
         font   = 2,
         xpd    = NA
         )
  }
  par(old_mai)
  return(invisible(NULL))
}



AveragedHeatmap <- function(input_df,
                            use_column,
                            take_median            = FALSE,
                            both_replicates        = TRUE,
                            main_title             = "",
                            use_subtext            = "",
                            ColorFunction          = NULL,
                            label_values           = FALSE,
                            take_log2              = FALSE,
                            uniform_legend         = NULL,
                            use_one_scale          = TRUE,
                            weighting_for_controls = NULL,
                            ...
                            ) {

  ## Calculate mean values for display in the heatmap
  plates_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  has_replicates <- grepl("_rep", use_column, fixed = TRUE)
  if (has_replicates) {
    rep_columns <- BothRepColumns(use_column)
  }
  if (both_replicates && has_replicates) {
    vec_list <- c(split(input_df[, rep_columns[[1]]], plates_vec),
                  split(input_df[, rep_columns[[2]]], plates_vec)
                  )
  } else {
    vec_list <- split(input_df[, use_column], plates_vec)
  }
  all_mat <- do.call(cbind, vec_list)
  if (take_median) {
    average_vec <- apply(all_mat, 1, median)
  } else {
    average_vec <- rowMeans(all_mat)
  }


  ## Compute breakpoints for the color scale
  if (use_one_scale) {
    breaks_list <- BreaksForColumn(input_df, use_column, take_log2 = take_log2,
                                   weighting_for_controls = weighting_for_controls,
                                   ...
                                   )
  } else {
    mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
    are_empty <- input_df[, "Well_number_384"] %in% c(mat_384[, c(1, 24)])
    are_pos_ctrl <- input_df[, "Target_flag"] %in% "Pos. control"
    empty_layouts <- unique(split(are_empty, plates_vec))
    pos_layouts <- unique(split(are_pos_ctrl, plates_vec))
    if ((length(empty_layouts) == 1) && (length(pos_layouts) == 1)) {
      ## If positive controls and empty wells show exactly the same distribution
      ## on all plates, then we can respect the 'weighting_for_controls'
      ## parameter. Otherwise, it must be ignored.
      are_first_plate <- plates_vec == plates_vec[[1]]
      include_columns <- c("Well_number_384", "Is_pos_ctrl", use_column)
      if (has_replicates) {
        include_columns <- union(include_columns, rep_columns)
      }
      submit_df <- input_df[are_first_plate, include_columns]
      submit_df[, use_column] <- average_vec
      if (has_replicates) {
        for (column_name in rep_columns) {
          input_df[, column_name] <- input_df[, use_column]
        }
      }
      breaks_list <- BreaksForColumn(submit_df, use_column,
                                     take_log2 = take_log2,
                                     weighting_for_controls = weighting_for_controls,
                                     ...
                                     )
    } else {
      if (isTRUE(weighting_for_controls)) {
        message("If 'use_one_scale' is FALSE (i.e. values are to scaled to ",
                "the plate average only),", "\n",
                "then 'weighting_for_controls' is ",
                "ignored. The parameter is not applicable due to", "\n",
                "the different distribution of controls across plates."
                )
      }
      breaks_list <- BreaksForColumn(NULL, use_column,
                                     use_vector = average_vec,
                                     take_log2 = take_log2,
                                     weighting_for_controls = weighting_for_controls,
                                     ...
                                     )
    }
  }

  uniform_legend <- breaks_list[["type"]] == "uniform"
  HeatMap384(numeric_vec    = average_vec,
             use_breaks     = breaks_list[["breaks"]],
             ColorFunction  = ColorFunction,
             main_title     = main_title,
             use_subtext    = use_subtext,
             label_values   = label_values,
             use_minuslog10 = IsPValue(use_column),
             take_log2      = take_log2,
             uniform_legend = uniform_legend,
             legend_cex     = if (use_column == "CellTiterGlo_raw") 0.7 else 0.8,
             well_label_cex = if (use_column == "CellTiterGlo_raw") 0.3 else 0.4
             )

  return(invisible(NULL))
}



HeatmapForPlate <- function(input_df,
                            plate_number,
                            use_column,
                            main_title     = NULL,
                            use_subtext    = "",
                            ColorFunction  = NULL,
                            show_z_prime   = NULL,
                            label_values   = FALSE,
                            take_log2      = FALSE,
                            uniform_legend = NULL,
                            use_one_scale  = TRUE,
                            ...
                            ) {

  has_replicates <- grepl("_rep", use_column, fixed = TRUE)

  if (is.numeric(plate_number)) {
    plate_number <- as.character(as.roman(plate_number))
  }
  are_this_plate <- input_df[, "Plate_number_384"] == plate_number

  if (is.null(main_title)) {
    main_title <- paste0("Plate ", plate_number)
    if (has_replicates) {
      rep_number <- GetRepNumber(use_column)
      main_title <- paste0(main_title, ", replicate ", rep_number)
    }
  }

  if (is.null(show_z_prime)) {
    show_z_prime <- has_replicates
  }

  if (show_z_prime) {
    z_prime <- Calculate_Z_Prime(input_df[are_this_plate, ], use_column)
    z_prime_text <- paste0("Z' = ", formatC(z_prime, digits = 2, format = "f"))
    if (is.null(use_subtext) || (nchar(use_subtext) == 0)) {
      use_subtext <- z_prime_text
    } else {
      if (is.expression(use_subtext)) {
        z_prime_text <- as.expression(bquote(.(z_prime_text)))
        use_subtext <- ConcatenateExpressions(list(use_subtext, z_prime_text), my_sep = ", ")
      } else {
        use_subtext <- paste0(use_subtext, ", ", z_prime_text)
      }
    }
  }

  if (use_one_scale) {
    breaks_list <- BreaksForColumn(input_df, use_column, take_log2 = take_log2, ...)
  } else {
    include_columns <- c("Well_number_384", "Is_pos_ctrl", use_column)
    if (has_replicates) {
      rep_columns <- BothRepColumns(use_column)
      include_columns <- union(include_columns, rep_columns)
      for (column_name in rep_columns) {
        input_df[, column_name] <- input_df[, use_column]
      }
    }
    breaks_list <- BreaksForColumn(input_df[are_this_plate, include_columns],
                                   use_column, take_log2 = take_log2, ...
                                   )
  }
  if (is.null(uniform_legend)) {
    uniform_legend <- breaks_list[["type"]] == "uniform"
  }

  HeatMap384(numeric_vec    = input_df[are_this_plate, use_column],
             use_breaks     = breaks_list[["breaks"]],
             ColorFunction  = ColorFunction,
             main_title     = main_title,
             use_subtext    = use_subtext,
             label_values   = label_values,
             use_minuslog10 = IsPValue(use_column),
             take_log2      = take_log2,
             uniform_legend = uniform_legend,
             legend_cex     = if (use_column == "CellTiterGlo_raw") 0.7 else 0.8,
             well_label_cex = if (use_column == "CellTiterGlo_raw") 0.3 else 0.4
             )

  return(invisible(NULL))
}




ExportAllHeatmaps <- function(input_df, export_PNGs = TRUE, only_schematics = FALSE) {

  heatmap_width <- 8
  heatmap_height <- 6.5

  message("Exporting plate schematics...")
  for (label_genes in c(TRUE, FALSE)) {
    for (color_by_96wp in c(TRUE, FALSE)) {

      if (label_genes) {
        label_string <- "Heatmap schematic - genes labelled"
        PNG_folder_name <- "PNGs - heatmap schematic - genes labelled"
      } else {
        label_string <- "Heatmap schematic"
        PNG_folder_name <- "PNGs - heatmap schematic"
      }
      if (color_by_96wp) {
        label_string <- paste0(label_string, " - colored by source plate")
        PNG_folder_name <- paste0(PNG_folder_name, " - colored by source plate")
      }

      pdf(file = file.path(output_dir, "Figures", "Heatmap schematic", paste0(label_string, ".pdf")),
          width = heatmap_width, height = heatmap_height
          )
      for (plate_number in 1:12) {
        PlateSchematic(input_df, plate_number, label_genes = label_genes,
                       color_by_source_plate = color_by_96wp
                       )
      }
      dev.off()

      for (plate_number in 1:12) {
        PNG_folder_path <- file.path(output_dir, "Figures", "Heatmap schematic", PNG_folder_name)
        dir.create(PNG_folder_path, showWarnings = FALSE)
        file_name <- paste0(label_string, " - plate", plate_number, ".png")
        png(filename = file.path(output_dir, "Figures", "Heatmap schematic", PNG_folder_name, file_name),
            width = heatmap_width, height = heatmap_height, units = "in", res = 600
            )
        PlateSchematic(input_df, plate_number, label_genes = label_genes,
                       color_by_source_plate = color_by_96wp
                       )
        dev.off()
      }
    }
  }

  if (only_schematics) {
    return(invisible(NULL))
  }

  plate_average_text <- "Mean of all plates and replicates (well effect)"

  for (export_PDFs in c(TRUE, FALSE)) {

    if (export_PDFs) {
      message("Exporting PDF files...")
    } else {
      message("Exporting PNG files... ")
    }

    for (label_cells in c(FALSE, TRUE)) {

      if (label_cells) {
        heatmaps_folder <- "Heatmaps (labelled)"
        message("Exporting heatmaps with labels showing numerical values... ")
      } else {
        heatmaps_folder <- "Heatmaps"
        message("Exporting heatmaps without well labels... ")
      }

      for (scaled_per_plate in c(FALSE, TRUE)) {
        if (scaled_per_plate) {
          sub_folder_prefix <- "Scaled per plate"
          message("... Exporting heatmaps that are scaled per plate... ")
        } else {
          sub_folder_prefix <- "Scaled across plates"
          message("... Exporting heatmaps that share the same color scale across plates... ")
        }
        for (condense_controls in c(TRUE, FALSE)) {
          if (condense_controls) {
            sub_folder <- paste0(sub_folder_prefix, " - controls condensed")
            message("... Exporting heatmaps where the scale is 'squeezed' for empty wells and positive controls... ")
          } else {
            sub_folder <- paste0(sub_folder_prefix, " - uniform")
            message("... Exporting heatmaps using a uniform/linear color scale... ")
          }

          if (export_PDFs) {
            for (i in seq_along(column_file_names)) {

              current_column <- names(column_file_names)[[i]]
              message("...... for the metric:  ", current_column)
              has_replicates <- !(is.null(GetRepNumber(current_column)))

              file_name <- paste0("Heatmaps - ", i, ") ", column_file_names[[i]], ".pdf")
              column_subtext <- long_column_labels[[i]]

              pdf(file = file.path(output_dir, "Figures", heatmaps_folder, sub_folder, file_name),
                  width = heatmap_width, height = heatmap_height
                  )
              AveragedHeatmap(input_df, current_column,
                              main_title             = plate_average_text,
                              use_subtext            = column_subtext,
                              label_values           = label_cells,
                              use_one_scale          = !(scaled_per_plate),
                              weighting_for_controls = condense_controls
                              )
              AveragedHeatmap(input_df, current_column,
                              take_median            = TRUE,
                              main_title             = sub("Mean", "Median", plate_average_text, fixed = TRUE),
                              use_subtext            = column_subtext,
                              label_values           = label_cells,
                              use_one_scale          = !(scaled_per_plate),
                              weighting_for_controls = condense_controls
                              )
              for (plate_number in 1:12) {
                if (has_replicates) {
                  rep_columns <- BothRepColumns(current_column)
                  for (rep_column in rep_columns) {
                    HeatmapForPlate(input_df, plate_number, rep_column,
                                    use_subtext            = column_subtext,
                                    label_values           = label_cells,
                                    use_one_scale          = !(scaled_per_plate),
                                    weighting_for_controls = condense_controls
                                    )
                  }
                } else {
                  HeatmapForPlate(input_df, plate_number, current_column,
                                  use_subtext            = column_subtext,
                                  label_values           = label_cells,
                                  use_one_scale          = !(scaled_per_plate),
                                  weighting_for_controls = condense_controls
                                  )
                }

              }
              dev.off()
            }
          } else {
            # Export PNG files
            if (!(export_PNGs)) {
              next
            }
            for (i in seq_along(column_file_names)) {

              current_column <- names(column_file_names)[[i]]
              message("...... for the metric:  ", current_column)

              has_replicates <- !(is.null(GetRepNumber(current_column)))
              column_subtext <- long_column_labels[[i]]

              folder_name <- paste0(i, ") ", sub("_rep1", "", current_column, fixed = TRUE))
              folder_path <- file.path(output_dir, "Figures", heatmaps_folder, "PNGs", sub_folder, folder_name)
              dir.create(folder_path, showWarnings = FALSE)

              file_name <- paste0("Heatmap - ", column_file_names[[i]],
                                  " -- mean across plates.png"
                                  )
              png(filename = file.path(folder_path, file_name),
                  width = heatmap_width, height = heatmap_height,
                  units = "in", res = 600
                  )
              AveragedHeatmap(input_df, current_column,
                              main_title             = plate_average_text,
                              use_subtext            = column_subtext,
                              label_values           = label_cells,
                              use_one_scale          = !(scaled_per_plate),
                              weighting_for_controls = condense_controls
                              )
              dev.off()

              file_name <- paste0("Heatmap - ", column_file_names[[i]],
                                  " -- median across plates.png"
                                  )
              png(filename = file.path(folder_path, file_name),
                  width = heatmap_width, height = heatmap_height,
                  units = "in", res = 600
                  )
              AveragedHeatmap(input_df, current_column,
                              take_median            = TRUE,
                              main_title             = plate_average_text,
                              use_subtext            = column_subtext,
                              label_values           = label_cells,
                              use_one_scale          = !(scaled_per_plate),
                              weighting_for_controls = condense_controls
                              )
              dev.off()

              for (plate_number in 1:12) {
                if (has_replicates) {
                  rep_columns <- BothRepColumns(current_column)
                  for (j in 1:2) {
                    file_name <- paste0("Heatmap - ", column_file_names[[i]],
                                        " - plate ", plate_number, " rep ", j,
                                        ".png"
                                        )
                    png(filename = file.path(folder_path, file_name),
                        width = heatmap_width, height = heatmap_height,
                        units = "in", res = 600
                        )
                    HeatmapForPlate(input_df, plate_number, rep_columns[[j]],
                                    use_subtext            = column_subtext,
                                    label_values           = label_cells,
                                    use_one_scale          = !(scaled_per_plate),
                                    weighting_for_controls = condense_controls
                                    )
                    dev.off()
                  }
                } else {
                  file_name <- paste0("Heatmap - ", column_file_names[[i]],
                                      " - plate ", plate_number, ".png"
                                      )
                  png(filename = file.path(folder_path, file_name),
                      width = heatmap_width, height = heatmap_height,
                      units = "in", res = 600
                      )
                  HeatmapForPlate(input_df, plate_number, current_column,
                                  use_subtext            = column_subtext,
                                  label_values           = label_cells,
                                  use_one_scale          = !(scaled_per_plate),
                                  weighting_for_controls = condense_controls
                                  )
                  dev.off()
                }
              }
            }
          }
        }
      }
    }
  }
  return(invisible(NULL))
}


