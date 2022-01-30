### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("squash")
library("scales")





# Define plot dimensions --------------------------------------------------

use_height <- 7
use_width <- 6.5





# Choose default color function -------------------------------------------

DefaultColFun <- colorRampPalette(rev(viridis::magma(100)))
# DefaultColFun <- colorRampPalette(rev(sequential_hcl(n = 10, h = c(51, -341), c = c(60, NA, 0), l = c(12, 96), power = c(0.35, 1.3))))





# Define plot titles ------------------------------------------------------

ccs3_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "3 consensus reads "} *
                               "and " >= "99% accuracy)"
                               ))

ccs5_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "5 consensus reads "} *
                               "and " >= "99.9% accuracy)"
                               ))

ccs7_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "7 consensus reads "} *
                               "and " >= "99.99% accuracy)"
                               ))





# Helper functions --------------------------------------------------------

MakeEmptyPlot <- function(y_limits = c(0, 1), use_cex = NULL) {
  plot(1, type = "n", axes = FALSE, ann = FALSE,
       xlim = c(0, 1), ylim = y_limits, xaxs = "i", yaxs = "i"
       )
  if (!(is.null(use_cex))) {
    par("cex" = 1)
  }
}


DrawColorTrapezoid <- function(start_x     = 0,
                               start_y     = 0,
                               end_x       = 1,
                               end_y       = 1,
                               trapezoid_y = NULL,
                               use_colors  = colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
                               ) {

  my_x_positions <- seq(start_x, end_x, by = (end_x - start_x) / length(use_colors))
  num_positions <- length(my_x_positions)

  bottom_left_corners_x  <- my_x_positions[2:(num_positions - 1L)]
  top_left_corners_x     <- bottom_left_corners_x
  top_right_corners_x    <- my_x_positions[3:num_positions]
  bottom_right_corners_x <- top_right_corners_x

  draw_triangle <- !(!(is.null(trapezoid_y)) && (trapezoid_y == end_y))

  if (draw_triangle) {
    if (is.null(trapezoid_y)) {
      triangle_start_y <- start_y
    } else {
      triangle_start_y <- trapezoid_y
    }

    my_y_positions <- seq(triangle_start_y,
                          end_y,
                          by = (end_y - triangle_start_y) / length(use_colors)
                          )

    bottom_left_corners_y  <- rep_len(start_y, num_positions - 1L)
    top_left_corners_y     <- my_y_positions[2:(num_positions - 1L)]
    top_right_corners_y    <- my_y_positions[3:num_positions]
    bottom_right_corners_y <- rep_len(start_y, num_positions - 1L)
  }



  if (!(is.null(trapezoid_y))) {
    rect(xleft   = my_x_positions[[1]],
         xright  = my_x_positions[[2]],
         ybottom = start_y,
         ytop    = trapezoid_y,
         col     = use_colors[[1]],
         border  = NA,
         xpd     = NA
         )

    for (i in seq_len(length(use_colors) - 1L)) {
      my_color <- use_colors[[i + 1L]]
      rect(xleft   = bottom_left_corners_x[[i]],
           xright  = top_right_corners_x[[i]],
           ybottom = start_y,
           ytop    = trapezoid_y,
           col     = my_color,
           border  = NA,
           xpd     = NA
           )
    }


  }

  if (!(!(is.null(trapezoid_y)) && (trapezoid_y == end_y))) {
    polygon(x      = c(my_x_positions[[1]], my_x_positions[[2]], my_x_positions[[2]]),
            y      = c(my_y_positions[[1]], my_y_positions[[2]], my_y_positions[[1]]),
            col    = use_colors[[1]],
            border = NA,
            xpd    = NA
            )

    for (i in seq_len(length(use_colors) - 1L)) {
      my_color <- use_colors[[i + 1L]]
      polygon(x      = c(bottom_left_corners_x[[i]], top_left_corners_x[[i]],
                         top_right_corners_x[[i]], bottom_right_corners_x[[i]]
                         ),
              y      = c(bottom_left_corners_y[[i]], top_left_corners_y[[i]],
                         top_right_corners_y[[i]], bottom_right_corners_y[[i]]
                         ),
              col    = my_color,
              border = NA,
              xpd    = NA
              )
    }
  }
}





# Necessary for vertical bars and input
Transposify <- function(numeric_input) {
  t(apply(as.matrix(numeric_input), 2, rev))
}



Do_cimage <- function(color_input) {
# If a vector is provided, by default it will be turned into a HORIZONTAL matrix
  if (!(is.matrix(color_input))) {
    input_matrix <- as.matrix(color_input)
  } else if (nrow(color_input) == 1) {
    input_matrix <- t(color_input)
  } else {
    input_matrix <- Transposify(color_input)
  }
  n_col <- ncol(color_input)
  n_row <- nrow(color_input)
  cimage(x    = seq_len(n_col) / n_col - (1 / n_col * 0.5),
         y    = seq_len(n_row) / n_row - (1 / n_row * 0.5),
         zcol = input_matrix,
         axes = FALSE,
         xlab = "",
         ylab = "",
         add  = TRUE
         )
}



MakeInvisible <- function(expression_string,
                          keep_index,
                          keep_regex = paste0("color", keep_index, "(\""),
                          replace_regex = "color[0-9]+\\(\""
                          ) {
  results_string <- gsub(keep_regex, "textstyle(\"", expression_string, fixed = TRUE)
  results_string <- gsub(replace_regex, "phantom(\"", results_string)
  return(results_string)
}


ConcatenateExpressions <- function(expression_list, my_sep = "  \u2013  ") {
  literal_strings <- vapply(expression_list, StripExpression, "")
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}


VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("gh")), use_expression, expression(phantom("gh")))
  return(ConcatenateExpressions(my_list, my_sep = ""))
}


StripExpression <- function(my_expression) {
  if (is.character(my_expression)) {
    literal_string <- paste0("\"", capture.output(cat(my_expression)), "\"")
  } else {
    literal_string <- paste0(capture.output(my_expression), collapse = "")
    if (substr(literal_string, 1L, 11L) == "expression(") {
      literal_string <- substr(literal_string, 12L, nchar(literal_string) - 1L)
    }
  }
  return(literal_string)
}


Darken <- function(color, factor = 1.4) {
  # from https://gist.github.com/Jfortin1/72ef064469d1703c6b30
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}


palify_cache_101 <- list()
Palify <- function(myhex, fraction_pale = 0.5) {
  if (myhex %in% names(palify_cache_101)) {
    color_vec <- palify_cache_101[[myhex]]
  } else {
    color_vec <- colorRampPalette(c(myhex, "#FFFFFF"))(101)
    palify_cache_101[[myhex]] <- color_vec
    assign("palify_cache_101", palify_cache_101, envir = globalenv())
  }
  color_vec[[round(fraction_pale * 100) + 1]]
}



AddVerticalWhiteSpace <- function(color_matrix, groups_vec, whitespace_color = "#FFFFFF00", divider_width = 2L) {

  if (ncol(color_matrix) != length(groups_vec)) {
    stop("The number of columns of the color matrix must match the length of the group labels!")
  }
  num_groups <- length(unique(groups_vec))
  whitespace_bar <- as.matrix(rep(whitespace_color, nrow(color_matrix)))
  colnames(whitespace_bar) <- "divider"
  whitespace_bar <- whitespace_bar[, rep(1, divider_width)]
  indices_list <- lapply(unique(groups_vec), function(x) which(x == groups_vec))
  matrices_list <- do.call(c, lapply(seq_len(num_groups - 1),
                                     function(x) list(color_matrix[, indices_list[[x]], drop = FALSE], whitespace_bar)
                                     )
                           )
  matrices_list <- c(matrices_list, list(color_matrix[, indices_list[[num_groups]], drop = FALSE]))
  results_matrix <- do.call(cbind, matrices_list)
  return(results_matrix)
}



GappedPositionsVec <- function(groups_vec, gap_weight = 2L) {
  positions_vec <- rep(NA, length(groups_vec))
  current_index <- 0L
  current_block <- groups_vec[[1]]
  for (i in seq_along(groups_vec)) {
    if (current_block != groups_vec[[i]]) {
      current_block <- groups_vec[[i]]
      current_index <- current_index + gap_weight
    }
    current_index <- current_index + 1L
    positions_vec[[i]] <- current_index
  }
  return(positions_vec)
}



DrawSideLegend <- function(labels_list,
                           use_colors,
                           top_labels      = NULL,
                           use_pch         = NULL,
                           point_cex       = 1.5,
                           use_line_width  = 3,
                           lines_x_start   = 1,
                           lines_x_title   = lines_x_start - 1,
                           small_y_gap     = 1.25,
                           large_gap_ratio = 1.75,
                           segment_left    = lines_x_start - 0.3,
                           segment_right   = lines_x_start + 0.4,
                           point_x_start   = lines_x_start + 0.2
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing the legend
  y_mid <- 0.5
  small_gap <- diff(grconvertY(c(0, small_y_gap), from = "line", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_ratio
  have_multiline <- any(lengths(labels_list) > 1)

  if (have_multiline) {
    are_first <- unlist(lapply(labels_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  } else {
    gaps_vec <- rep(medium_gap, length(labels_list))
  }

  are_top_labels <- rep(FALSE, length(gaps_vec))
  if (!(is.null(top_labels))) {
    if (have_multiline) {
      are_first <- c(TRUE, rep(FALSE, length(top_labels) - 1))
      top_gaps_vec <- ifelse(are_first, large_gap, small_gap)
    } else {
      top_gaps_vec <- rep(medium_gap, length(top_labels))
    }
    are_top_labels <- c(rep(TRUE, length(top_labels)), are_top_labels)
    gaps_vec[[1]] <- (large_gap + 2 * medium_gap) / 3
    gaps_vec <- c(top_gaps_vec, gaps_vec)
    text_list <- c(as.list(top_labels), labels_list)
  } else {
    text_list <- labels_list
    are_top_labels <- rep(FALSE, length(gaps_vec))
  }
  gaps_vec[[1]] <- 0

  total_span <- sum(gaps_vec)

  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  LinesFromEdge <- function(num_lines) {
    par("usr")[[2]] + diff(grconvertX(c(0, num_lines), from = "lines", to = "user"))
  }

  ## Draw the legend
  text(x      = ifelse(are_top_labels,
                       LinesFromEdge(lines_x_title),
                       LinesFromEdge(lines_x_start)
                       ),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(text_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )

  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  assign("delete_groups_vec", groups_vec, envir = globalenv())
    assign("delete_y_pos", y_pos, envir = globalenv())
  assign("delete_are_top_labels", are_top_labels, envir = globalenv())

  groups_y_pos <- tapply(y_pos[!(are_top_labels)], groups_vec, mean)


  if (!(is.null(use_pch))) {
    points(x   = rep(LinesFromEdge(point_x_start), length(use_colors)),
           y   = groups_y_pos,
           col = "gray30",
           bg  = use_colors,
           pch = use_pch,
           cex = point_cex,
           xpd = NA
           )
  } else {
    segments(x0  = LinesFromEdge(segment_left),
             x1  = LinesFromEdge(segment_right),
             y0  = groups_y_pos,
             lwd = use_line_width * par("lwd"),
             col = use_colors,
             xpd = NA
             )
  }

  return(invisible(NULL))
}





# Plotting functions ------------------------------------------------------

PlotBarplotMat <- function(barplot_mat,
                           colors_vec,
                           positions_vec = seq_len(ncol(barplot_mat))
                           ) {

  num_categories <- nrow(barplot_mat)
  stopifnot(length(colors_vec) == num_categories)

  lower_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    if (x == 1) {
      rep(0, ncol(barplot_mat))
    } else {
      colSums(barplot_mat[seq_len(x - 1), , drop = FALSE])
    }
  })

  upper_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    colSums(barplot_mat[seq_len(x), , drop = FALSE])
  })

  for (i in seq_len(ncol(barplot_mat))) {
    for (j in seq_len(num_categories)) {
      rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
           xright  = (positions_vec[[i]] / max(positions_vec)),
           ybottom = lower_borders_vec_list[[j]][[i]],
           ytop    = upper_borders_vec_list[[j]][[i]],
           col     = colors_vec[[j]],
           border  = NA,
           xpd     = NA
           )
    }
  }
  return(invisible(NULL))
}



MakeColorBoxLegend <- function(labels_vec,
                               colors_vec,
                               y_pos,
                               x_pos,
                               aspect_ratio       = 1,
                               before_text        = NULL,
                               after_text         = NULL,
                               constant_multipler = 1,
                               use_constant_space = TRUE,
                               vertical_adjust    = 0.45,
                               x_space_adjust     = 3.3,
                               use_lwd            = 1
                               ){

  num_labels <- length(labels_vec)

  stopifnot(num_labels == length(colors_vec))

  base_width <- strwidth("a", units = "user")

  if (use_constant_space) {
    string_widths <- rep(base_width * constant_multipler, num_labels)
  } else {
    string_widths <- strwidth(labels_vec, units = "user")
  }

  rectangle_x_start  <- x_pos
  rectangle_x_spaces <- string_widths + (x_space_adjust * base_width)

  rectangle_x_mids <- c(rectangle_x_start,
                        rectangle_x_start +
                        cumsum(rectangle_x_spaces[seq_len(num_labels - 1)])
                        )
  rectangle_width  <- 1.3 * base_width
  rectangle_height <- rectangle_width * aspect_ratio
  rectangle_y_mid  <- y_pos

  rect(xleft   = rectangle_x_mids - (rectangle_width / 2),
       xright  = rectangle_x_mids + (rectangle_width / 2),
       ytop    = rectangle_y_mid + (rectangle_height / 2),
       ybottom = rectangle_y_mid - (rectangle_height / 2),
       xpd     = NA,
       col     = colors_vec,
       border  = "gray30",
       lwd     = use_lwd
       )

  text_gap <- rectangle_width * (-0.5)
  text_y_pos <- rectangle_y_mid - (vertical_adjust * base_width)

  text(x      = rectangle_x_mids + text_gap,
       y      = text_y_pos,
       labels = sapply(labels_vec, VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )
  side_space <- rectangle_width * 0.9

  text(x      = rectangle_x_start - (text_gap),
       y      = text_y_pos,
       labels = VerticalAdjust(before_text),
       xpd    = NA,
       adj    = c(1, 0.5)
       )
  text(x      = rectangle_x_mids[[length(rectangle_x_mids)]] +
                side_space,
       y      = text_y_pos,
       labels = VerticalAdjust(after_text),
       xpd    = NA,
       adj    = c(0, 0.5)
       )

}




DrawPercentCorrectBarplot <- function(summary_df,
                                      invert_barplot         = TRUE,
                                      add_gap                = FALSE,
                                      gap_weight             = 2L,
                                      show_correct_promoters = FALSE,
                                      ColorFunction          = DefaultColFun,
                                      show_zero_in_legend    = TRUE,
                                      prefix_text            = "% plasmids with ",
                                      postfix_text           = " correct gRNAs (including the tracrRNA)",
                                      text_at_bottom         = FALSE,
                                      narrow_lwd             = TRUE,
                                      all_text_bold          = FALSE,
                                      color_box_legend       = FALSE,
                                      aspect_ratio           = 2
                                      ) {

  if (add_gap) {
    positions_vec <- GappedPositionsVec(summary_df[, "Block"], gap_weight = gap_weight)
  } else {
    positions_vec <- seq_len(nrow(summary_df))
  }

  barplot_columns <- paste0("Perc_", c(paste0("at_least_", 1:3), "all_4"))
  if (show_correct_promoters) {
    barplot_columns <- c(barplot_columns, "Perc_all_4_promoters")
    barplot_colors <- ColorFunction(6)
    barplot_colors <- vapply(barplot_colors, Palify, fraction_pale = 0.3, "")
    BarPlotColorFun <- ColorFunction
  } else {
    BarPlotColorFun <- colorRampPalette(ColorFunction(100)[1:86])
    barplot_colors <- BarPlotColorFun(5)
    barplot_colors <- vapply(barplot_colors, Palify, fraction_pale = 0.2, "")
  }

  barplot_mat <- do.call(rbind, lapply(barplot_columns, function(x) summary_df[[x]]))
  barplot_mat <- barplot_mat / 100
  barplot_mat <- rbind(
    "none"      = 1 - barplot_mat[1, ],
    "exactly_1" = barplot_mat[1, ] - barplot_mat[2, ],
    "exactly_2" = barplot_mat[2, ] - barplot_mat[3, ],
    "exactly_3" = barplot_mat[3, ] - barplot_mat[4, ],
    "all_4"     = barplot_mat[4, ]
  )

  if (invert_barplot) {
    PlotBarplotMat(barplot_mat, barplot_colors, positions_vec)
  } else {
    PlotBarplotMat(barplot_mat[rev(seq_len(nrow(barplot_mat))), ],
                   rev(barplot_colors),
                   positions_vec
                   )
  }

  tick_locations <- axTicks(2)
  tick_labels <- paste0(tick_locations * 100, "%")
  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.3,
       lwd      = if (narrow_lwd) 0.75 else 1,
       cex.axis = 0.9,
       pos      = par("usr")[[1]] - GetHalfLineWidth()
       )
  if (!(add_gap)) {
    DrawOuterBox(use_lwd = if (narrow_lwd) 0.5 else 1)
  }

  if (color_box_legend) {

    MakeColorBoxLegend(labels_vec         = as.character(0:4),
                       colors_vec         = barplot_colors,
                       y_pos              = par("usr")[[1]] - diff(grconvertY(c(0, 1), from = "lines", to = "user")),
                       x_pos              = 0.092,
                       after_text         = "correct sgRNAs",
                       aspect_ratio       = aspect_ratio,
                       use_constant_space = TRUE
                       )


  } else {

    if (show_correct_promoters) {
      color_text_vec <- c(paste0(if (all_text_bold) 'bold' else 'plain',
                                 '(color1("', prefix_text, '") * color1("") *'
                                 ),
                          'scriptscriptstyle(" ") *',
                          'bold(color2("0,")) * ',
                          'color1(" ") * bold(color3("1,")) * ',
                          'color1(" ") * bold(color4("2,")) * ',
                          'color1(" ") * bold(color5("3,")) * ',
                          'color1(" ") * bold(color6("4")) * ',
                          paste0('color1(" correct gRNAs ") * '),
                          'color7("or an entirely correct construct"))'
                          )
    } else {
      color_text_vec <- c(paste0(if (all_text_bold) 'bold' else 'plain',
                                 '(color1("', prefix_text, '") * color1("") *'
                                 ),
                          'scriptscriptstyle(" ") *',
                          'bold(color2("0,")) * ',
                          'color1(" ") * bold(color3("1,")) * ',
                          'color1(" ") * bold(color4("2,")) * ',
                          'color1(" ") * bold(color5("3,")) * ',
                          'color1(" or ") * bold(color6("4")) * ',
                          'color1("', postfix_text, '"))'
                          )
    }

    text_colors <- c("#000000",
                     "#DDE100",
                     BarPlotColorFun(length(barplot_colors) + 1)[2:(length(barplot_colors))]
                     )

    color_indices <- seq_along(text_colors)
    text_indices <- seq_along(color_text_vec)
    if (!(show_zero_in_legend)) {
      are_zero <- grepl('"0,"', color_text_vec, fixed = TRUE)
      text_indices <- which(!(are_zero))
      color_indices <- color_indices[-2]
    }
    color_text <- paste0(color_text_vec[text_indices], collapse = "")

    for (j in color_indices) {
      text(x      = 0.5,
           y      = if (text_at_bottom) -0.18 else 1.05,
           labels = VerticalAdjust(parse(text = MakeInvisible(color_text, j))),
           adj    = c(0.5, 0),
           xpd    = NA,
           col    = text_colors[[j]]
           )
    }
  }
}



GetHalfLineWidth <- function(y_axis = FALSE) {
  if (y_axis) {
    diff(grconvertY(c(0, par("lwd") / 96), from = "inches", to = "user")) / 2
  } else {
    diff(grconvertX(c(0, par("lwd") / 96), from = "inches", to = "user")) / 2
  }
}



DrawOuterBox <- function(use_lwd = 1, fill = FALSE) {
  half_lwd_x <- GetHalfLineWidth() * use_lwd
  half_lwd_y <- GetHalfLineWidth(y_axis = TRUE) * use_lwd
  rect(xleft   = par("usr")[[1]] - half_lwd_x,
       xright  = par("usr")[[2]] + half_lwd_x,
       ybottom = par("usr")[[3]] - half_lwd_y,
       ytop    = par("usr")[[4]] + half_lwd_y,
       xpd     = NA,
       lwd     = use_lwd * par("lwd"),
       ljoin   = "mitre",
       lmitre  = 30,
       col     = if (fill) "black" else NA
       )
}


DrawHeatMap <- function(summary_df,
                        add_gap               = FALSE,
                        gap_weight            = 2L,
                        ColorFunction         = DefaultColFun,
                        trapezoid_on_new_plot = TRUE,
                        trapezoid_start_x     = 0.85,
                        trapezoid_end_x       = 1,
                        trapezoid_start_y     = 0.65,
                        trapezoid_end_y       = 0.9,
                        trapezoid_y           = 0.72,
                        add_percent           = TRUE,
                        accuracy_label        = expression(bold("Accuracy")),
                        accuracy_label_cex    = 0.9,
                        use_lwd               = 0.5,
                        label_y_factor        = 1,
                        bold_percentages      = TRUE
                        ) {

  ## Draw the heatmap

  percent_columns <- paste0("Perc_sg", 1:4, "_cr", 1:4)
  numeric_mat <- t(as.matrix(summary_df[, percent_columns])) / 100

  my_breaks <- seq(0, 1.01, by = 0.01)
  my_breaks[c(1, length(my_breaks))] <- c(0, ceiling(my_breaks[length(my_breaks)]))

  my_cmap <- makecmap(numeric_mat, colFn = ColorFunction, breaks = my_breaks)
  my_color_mat <- cmap(numeric_mat, my_cmap)

  if (add_gap) {
    my_color_mat <- AddVerticalWhiteSpace(my_color_mat,
                                          groups_vec = summary_df[, "Block"],
                                          divider_width = gap_weight
                                          )
  }

  MakeEmptyPlot()
  Do_cimage(my_color_mat)
  x_range <- par("usr")[[2]] - par("usr")[[1]]

  text(x      = par("usr")[[1]] - (x_range * 0.018),
       y      = seq(0.875, 0.125, by = -0.25),
       labels = paste0("sg", 1:4),
       xpd    = NA,
       adj    = c(1, 0.5)
       )


  if (trapezoid_on_new_plot) {
    MakeEmptyPlot()
  }

  ## Draw the color indicator (trapezoid)

  DrawColorTrapezoid(start_x     = trapezoid_start_x,
                     end_x       = trapezoid_end_x,
                     start_y     = trapezoid_start_y,
                     end_y       = trapezoid_end_y,
                     trapezoid_y = trapezoid_y,
                     use_colors  = ColorFunction(100)
                     )

  line_width <- GetHalfLineWidth() * use_lwd
  trapezoid_text_start_x <- trapezoid_start_x + line_width
  trapezoid_text_end_x <- trapezoid_end_x - line_width

  trapezoid_x_range <- trapezoid_text_end_x - trapezoid_text_start_x
  trapezoid_y_range <- trapezoid_end_y - trapezoid_start_y

  text(x      = trapezoid_text_start_x - (trapezoid_x_range * 0.07),
       y      = trapezoid_start_y + (trapezoid_y_range * 0.3),
       adj    = c(1, 0.5),
       labels = accuracy_label,
       xpd    = NA,
       cex    = accuracy_label_cex
       )

  trapezoid_seq <- seq(trapezoid_text_start_x,
                       trapezoid_text_end_x,
                       by = (trapezoid_x_range / 5)
                       )
  percent_labels <- seq(0, 100, by = 20)
  if (add_percent) {
    percent_labels <- paste0(percent_labels, "%")
  }

  text(x      = trapezoid_seq,
       y      = trapezoid_start_y - (0.6 * trapezoid_y_range * label_y_factor),
       labels = percent_labels,
       font   = if (bold_percentages) 2 else 1,
       cex    = 0.6,
       xpd    = NA
       )
  segments(x0   = trapezoid_seq,
           x1   = trapezoid_seq,
           y0   = trapezoid_start_y - (0.352 * trapezoid_y_range),
           y1   = trapezoid_start_y - (0 * trapezoid_y_range),
           xpd  = NA,
           lwd  = use_lwd,
           lend = "square"
           )
  segments(x0   = trapezoid_text_start_x,
           x1   = trapezoid_text_end_x,
           y0   = trapezoid_start_y - (0 * trapezoid_y_range),
           y1   = trapezoid_start_y - (0 * trapezoid_y_range),
           xpd  = NA,
           lwd  = use_lwd,
           lend = "square"
           )
  return(invisible(NULL))
}




PrepareSummaryDf <- function(summary_df, reorder_wells = FALSE) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  percent_columns <- paste0("Perc_sg", 1:4, "_cr", 1:4)
  if (reorder_wells) {
    mean_accuracies <- rowMeans(as.matrix(summary_df[, percent_columns]))
    new_order <- order(summary_df[["Perc_all_4"]],
                       summary_df[["Perc_at_least_3"]],
                       summary_df[["Perc_at_least_2"]],
                       summary_df[["Perc_at_least_1"]],
                       mean_accuracies,
                       na.last = FALSE
                       )
  } else {
    new_order <- seq_len(nrow(summary_df))
  }

  if ("Empty_well" %in% names(sg_sequences_df)) {
    are_to_include <- !(sg_sequences_df[["Empty_well"]])
  } else {
    are_to_include <- rep(TRUE, nrow(sg_sequences_df))
  }

  use_indices <- new_order[are_to_include[new_order]]
  summary_df <- summary_df[use_indices, ]

  if ("Block" %in% names(sg_sequences_df)) {
    summary_df[["Block"]] <- sg_sequences_df[["Block"]][use_indices]
  }

  summary_df[["Longest_subsequence"]] <- sg_sequences_df[["Longest_subsequence"]][use_indices]

  row.names(summary_df) <- NULL
  return(summary_df)
}



DrawAccuracyHeatmap <- function(summary_df,
                                main_title             = NULL,
                                title_color            = "black",
                                ColorFunction          = DefaultColFun, # colorRampPalette(brewer.pal(9, "YlGnBu")),
                                reorder_wells          = TRUE,
                                show_zero_in_legend    = TRUE,
                                invert_barplot         = TRUE,
                                gap_weight             = 2L,
                                show_correct_promoters = FALSE
                                ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  assign("delete_summary_df_3", summary_df, envir = globalenv())

  if (reorder_wells) {
    add_gap <- FALSE
  } else {
    add_gap <- "Block" %in% names(sg_sequences_df)
  }
  summary_df <- PrepareSummaryDf(summary_df, reorder_wells = reorder_wells)
  num_wells <- nrow(summary_df)


  ## Set up the plot layout

  space_height <- 1
  golden_ratio <- (1 + sqrt(5)) / 2
  layout(cbind(rep(1, 11),
               3:(11 + 3 - 1),
               rep(2, 11)
               ),
         widths  = c(0.1, 0.8, 0.1),
         heights = c(space_height * 2.7,
                     3.5,
                     space_height * 0.13,
                     space_height * 0.5,
                     space_height * 1.2,
                     3 * golden_ratio,
                     space_height * 2,
                     space_height * 2.2,
                     space_height * 0.5,
                     space_height * 0.5,
                     space_height * 1
                     )
         )

  old_par <- par(mar = rep(0, 4))

  for (i in 1:3) {
    MakeEmptyPlot()
  }
  if (!(is.null(main_title))) {
    text(x      = 0.5,
         y      = 0.65,
         labels = main_title,
         col    = title_color,
         cex    = 1.1,
         xpd    = NA
         )

  }
  MakeEmptyPlot()


  ## Draw a vertical barplot

  DrawPercentCorrectBarplot(summary_df,
                            ColorFunction          = ColorFunction,
                            invert_barplot         = invert_barplot,
                            add_gap                = add_gap,
                            gap_weight             = gap_weight,
                            show_correct_promoters = show_correct_promoters,
                            show_zero_in_legend    = show_zero_in_legend
                            )


  ## Draw the % accuracy horizontal barplot

  are_at_least_75 <- summary_df[["Perc_all_4"]] >= 75
  num_above_75 <- sum(are_at_least_75, na.rm = TRUE)
  over_75_fraction <- num_above_75 / num_wells

  if (add_gap) {
    positions_vec <- GappedPositionsVec(summary_df[, "Block"],
                                        gap_weight = gap_weight
                                        )
  } else {
    positions_vec <- seq_len(nrow(summary_df))
  }

  for (i in 1:2) {
    MakeEmptyPlot()
  }
  two_grey_colors <- c("gray77", "gray40")

  if (reorder_wells) {
    rect(xleft   = c(0, 1 - over_75_fraction),
         xright  = c(1 - over_75_fraction, 1),
         ybottom = 0,
         ytop    = 1,
         border  = NA,
         lwd     = 0.5,
         col     = two_grey_colors
         )
  } else {
    for (i in seq_len(num_wells)) {
      assign("delete_are_at_least_75", are_at_least_75, envir = globalenv())
      is_75 <- are_at_least_75[[i]]
      if (is.na(is_75)) {
        use_color <- "white"
      } else {
        use_color <- two_grey_colors[[as.integer(is_75) + 1]]
      }
      rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
           xright  = (positions_vec[[i]] / max(positions_vec)),
           ybottom = 0,
           ytop    = 1,
           col     = use_color,
           border  = NA
           )
    }
  }

  text(x      = 0.5,
       y      = -0.65,
       adj    = c(0.5, 0.5),
       labels = bquote(bold(.(as.character(num_above_75)) *  " / " *
                            .(as.character(num_wells)) * " genes " *
                              "had 4 correct gRNAs in " >= "75% of reads"
                            )),
       font   = 2,
       col    = "black",
       xpd    = NA
       )
  MakeEmptyPlot()


  ## Draw the heatmap

  DrawHeatMap(summary_df, ColorFunction = ColorFunction,
              add_gap = add_gap, gap_weight = gap_weight
              )

  MakeEmptyPlot()

  ## Draw a vertical barplot

  perc_truncated <- summary_df[["Num_under_2kb"]] / summary_df[["Count_total"]]

  two_colors <- c("gray88", "gray40")

  for (i in seq_len(num_wells)) {
    if (is.na(perc_truncated[[i]])) {
      use_color <- "white"
    } else {
      use_color <- two_colors[[1]]
    }
    rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
         xright  = (positions_vec[[i]] / max(positions_vec)),
         ybottom = 0,
         ytop    = 1,
         col     = use_color,
         border  = NA
         )
    rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
         xright  = (positions_vec[[i]] / max(positions_vec)),
         ybottom = 0,
         ytop    = perc_truncated[[i]],
         col     = two_colors[[2]],
         border  = NA
         )
  }

  text(x      = 0.5,
       y      = 1.18,
       adj    = c(0.5, 0.5),
       labels = expression(bold("% short plasmids ("  <  "2 kb)")),
       font   = 2,
       col    = "black",
       xpd    = NA
       )


  for (i in 1:2) {
    MakeEmptyPlot()
  }

  ## Draw the bottom strip (indicating homologies / wells with low read numbers)

  if ("Longest_subsequence" %in% names("Longest_subsequence")) {
    are_to_highlight <- summary_df[["Longest_subsequence"]] >= 8
    strip_label <- bquote(bold(.(as.character(sum(are_to_highlight))) * " / " *
                                 .(as.character(num_wells)) * " genes " *
                                 "had " >= "8bp homology"
                               ))
  } else {
    are_to_highlight <- summary_df[["Count_total"]] < 20
    strip_label <- bquote(bold(.(as.character(sum(are_to_highlight))) * " / " *
                                 .(as.character(num_wells)) * " genes " *
                                 "had " < "20 reads"
                               ))
  }

  for (i in seq_len(num_wells)) {
    rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
         xright  = (positions_vec[[i]] / max(positions_vec)),
         ybottom = 0,
         ytop    = 1,
         col     = two_grey_colors[[as.integer(are_to_highlight[[i]]) + 1]],
         border  = NA
         )
  }

  text(x      = 0.5,
       y      = -0.65,
       adj    = c(0.5, 0.5),
       labels = strip_label,
       font   = 2,
       col    = "black",
       xpd    = NA
       )

  MakeEmptyPlot()
  par(old_par)
  layout(1)

  return(invisible(NULL))
}




SparseAccuracyHeatmap <- function(summary_df,
                                  main_title   = NULL,
                                  use_cex      = 0.7,
                                  use_lwd      = 0.8,
                                  margin_ratio = 10,
                                  use_width    = 3.6,
                                  use_height   = 3.05,
                                  title_color  = "black"
                                  ) {

  layout_mat <- cbind(rep(1, 5), 3:7, rep(2, 5))
  old_par <- par(mar = rep(0, 4))
  layout(layout_mat,
         heights = c(2.98, margin_ratio, 0.01, margin_ratio, 0.01),
         widths = c(1, margin_ratio, 1)
         )
  MakeEmptyPlot()
  MakeEmptyPlot()
  MakeEmptyPlot()

  if (!(is.null(main_title))) {
    text(x      = 0.5,
         y      = 0.4,
         labels = main_title,
         col    = title_color,
         cex    = use_cex * 1.2,
         xpd    = NA
         )
  }

  par(mai = c(0.3, 0.4, 0.13, 0.13), lwd = use_lwd, cex = use_cex)
  MakeEmptyPlot()
  DrawPercentCorrectBarplot(summary_df,
                            prefix_text      = "",
                            postfix_text     = " correct gRNAs",
                            text_at_bottom   = TRUE,
                            narrow_lwd       = FALSE,
                            color_box_legend = TRUE,
                            aspect_ratio     = (use_width * (margin_ratio / (margin_ratio + 2)) - sum(par("mai")[c(2, 4)])) /
                                               (use_height * (margin_ratio / (margin_ratio * 2 + 3)) - sum(par("mai")[c(1, 3)]))
                            )
  par(mar = rep(0, 4))
  MakeEmptyPlot()

  par(mai = c(0.3, 0.4, 0.13, 0.13), lwd = use_lwd, cex = use_cex)
  DrawHeatMap(summary_df,
              trapezoid_on_new_plot = FALSE,
              trapezoid_start_y     = -0.16,
              trapezoid_y           = -0.125,
              trapezoid_end_y       = -0.05,
              trapezoid_start_x     = 0.75,
              trapezoid_end_x       = 1,
              add_percent           = FALSE,
              accuracy_label        = "Accuracy (%)",
              accuracy_label_cex    = 1,
              label_y_factor        = 1.35,
              use_lwd               = 0.75,
              bold_percentages      = FALSE
              )
  par(mar = rep(0, 4))
  MakeEmptyPlot()

  par(old_par)
  layout(1)

  return(invisible(NULL))
}



ExportFiguresForManuscript <- function(summary_df, use_prefix) {

  use_width <- 3.0
  use_height <- 1.5
  use_lwd <- 0.8
  use_cex <- 0.7

  for (use_PDF in TRUE) {

    file_name <- paste0(use_prefix,
                        " - stacked barplot - CCS7 (filtered) - SmrtLink 7",
                        ".pdf"
                        )

    if (use_PDF) {
      pdf(file = file.path(manuscript_directory, file_name),
          width = use_width, height = use_height
          )
    }
    par(mai = c(0.4, 0.4, 0.13, 0.13), lwd = use_lwd, cex = use_cex)
    MakeEmptyPlot()
    DrawPercentCorrectBarplot(summary_df,
                              prefix_text      = "",
                              postfix_text     = " correct gRNAs",
                              text_at_bottom   = TRUE,
                              narrow_lwd       = FALSE,
                              color_box_legend = TRUE,
                              aspect_ratio     = (use_width  - sum(par("mai")[c(2, 4)])) /
                                                 (use_height - sum(par("mai")[c(1, 3)]))
                              )

    if (use_PDF) {
      dev.off()
    }

    file_name <- paste0(use_prefix,
                        " - heatmap - CCS7 (filtered) - SmrtLink 7",
                        ".pdf"
                        )
    if (use_PDF) {
      pdf(file = file.path(manuscript_directory, file_name),
          width = use_width, height = use_height
          )
    }
    par(mai = c(0.4, 0.4, 0.13, 0.13), lwd = use_lwd, cex = use_cex)
    DrawHeatMap(summary_df,
                trapezoid_on_new_plot = FALSE,
                trapezoid_start_y     = -0.16,
                trapezoid_y           = -0.125,
                trapezoid_end_y       = -0.05,
                trapezoid_start_x     = 0.75,
                trapezoid_end_x       = 1,
                add_percent           = FALSE,
                accuracy_label        = "Accuracy (%)",
                accuracy_label_cex    = 1,
                label_y_factor        = 1.15,
                use_lwd               = 0.75,
                bold_percentages      = FALSE
                )
    if (use_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}



PlateFileName <- function(plate_number) {
  plate_number_width <- max(nchar(as.character(plates_df[["Plate_number"]])))
  paste0(formatC(plate_number, width = plate_number_width, flag = "0"), ") - ",
         plates_df[["Plate_name"]][plates_df[["Plate_number"]] == plate_number]
         )
}


ModifiedAlterationBarplot <- function(summary_df,
                                      reorder_wells = FALSE,
                                      main_title = NULL,
                                      title_color = "black"
                                      ) {

  DrawAlterationBarplot(summary_df,
                        reorder_wells        = reorder_wells,
                        main_title           = main_title,
                        title_color          = title_color,
                        show_color_legend    = TRUE,
                        show_color_text      = FALSE,
                        top_space            = 1.75,
                        bottom_space         = 1.25,
                        space_height         = 0.75,
                        sg_label_cex         = 1.1,
                        horizontal_y_lab_pos = -0.08,
                        color_legend_x_pos   = 0.15,
                        title_y_pos          = 0.47,
                        color_box_y_pos      = 0.4
                        )
}






# Functions used for the first two sequencing runs ------------------------

AllSmrtLinkVersionsForOnePlate <- function(PlotFunction,
                                           folder_prefix,
                                           file_prefix
                                           ) {

  for (smrtlink_version in c(7, 9)) {

    version_folder <- paste0("SmrtLink ", smrtlink_version)
    version_path <- file.path(plots_output_directory, version_folder)

    ccs_numbers <- c(3, 5, 7)
    use_df_list <- lapply(ccs_numbers, function(x) {
      df_list_name <- paste0("sl", smrtlink_version, "_ccs", x, "_df_list")
      return(get(df_list_name))
    })
    names(use_df_list) <- paste0("ccs", ccs_numbers)

    file_name_prefix <- paste0(file_prefix, " - SmrtLink ", smrtlink_version, " - ")

    DrawAllPlots(PlotFunction,
                 use_df_list,
                 version_path,
                 folder_prefix,
                 file_name_prefix
                 )

  }
}



DrawAllSubsampledPlotsForOnePlate <- function(PlotFunction,
                                              folder_prefix,
                                              PDF_prefix = NULL,
                                              use_num_reps = 3
                                              ) {

  for (smrtlink_version in c(7, 9)) {

    version_folder <- paste0("SmrtLink ", smrtlink_version, " - subsampled")
    version_path <- file.path(plots_output_directory, version_folder)

    if (smrtlink_version == 7) {
      subsampled_list <- sl7_subsampled_list
    } else {
      subsampled_list <- sl7_subsampled_list
    }

    for (i in seq_along(subsampled_list)) {
      for (rep in seq_len(min(length(subsampled_list[[i]]), use_num_reps))) {

        sampling_level <- names(subsampled_list)[[i]]
        is_all <- sampling_level == "100% sampled"
        if (is_all) {
          retained_string <- "All"
        } else {
          retained_string <- paste0("Sampled ",
                                    sub(" sampled", "", sampling_level, fixed = TRUE),
                                    " of"
                                    )
        }
        title_prefix <- paste0(retained_string, " reads",
                               if (is_all) "" else paste0(" \u2013 repetition #", rep)
                               )

        ccs3_title <- as.expression(bquote(plain({.(title_prefix) *
            " (" >= "3 consensus reads "} *
              "and " >= "99% accuracy)"
            )))

        ccs5_title <- as.expression(bquote(plain({.(title_prefix) *
            " (" >= "5 consensus reads "} *
              "and " >= "99.9% accuracy)"
            )))

        ccs7_title <- as.expression(bquote(plain({.(title_prefix) *
            " (" >= "7 consensus reads "} *
              "and " >= "99.99% accuracy)"
            )))



        file_name_postfix <- paste0(" - ", letters[[i]], ") ",
                                    sub("% sampled", " percent", names(subsampled_list)[[i]], fixed = TRUE),
                                    " of reads"
                                    )
        if (!(is_all)) {
          file_name_postfix <- paste0(file_name_postfix, " - ",
                                      names(subsampled_list[[i]])[[rep]]
                                      )
        }

        DrawAllPlots(PlotFunction,
                     subsampled_list[[i]][[rep]],
                     version_path,
                     folder_prefix,
                     file_prefix    = "",
                     file_postfix   = file_name_postfix,
                     draw_PDFs      = FALSE,
                     use_ccs3_title = ccs3_title,
                     use_ccs5_title = ccs5_title,
                     use_ccs7_title = ccs7_title
                     )
      }
    }

    ccs_numbers <- c(3, 5, 7)

    for (reorder_wells in c(TRUE, FALSE)) {

      order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
      order_folder <- paste0(folder_prefix, " - ", order_folder)
      order_path <- file.path(version_path, order_folder)

      for (i in seq_along(ccs_numbers)) {
        for (filter_stage in 1:3) {

          df_name <- c("original_summary_df", "filtered_summary_df", "filtered_gRNAs_df")[[filter_stage]]
          ccs_name <- paste0("ccs", ccs_numbers[[i]])

          file_name <- paste0(c("CCS3 (99)", "CCS5 (99.9)", "CCS7 (99.99)")[[i]],
                              " - ",
                              c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                              )
          if (!(is.null(PDF_prefix))) {
            file_name <- paste0(PDF_prefix, " - ", file_name)
          }
          pdf(file = file.path(order_path, paste0(file_name, ".pdf")),
              height = use_height,
              width  = use_width
              )

          for (sample_level_index in seq_along(subsampled_list)) {

            sampling_level <- names(subsampled_list)[[sample_level_index]]
            is_all <- sampling_level == "100% sampled"
            if (is_all) {
              retained_string <- "All"
            } else {
              retained_string <- paste0("Sampled ",
                                        sub(" sampled", "", sampling_level, fixed = TRUE),
                                        " of"
                                        )
            }

            for (rep in seq_len(min(length(subsampled_list[[sample_level_index]]), use_num_reps))) {

              use_title <- paste0(retained_string, " reads",
                                  if (is_all) "" else paste0(" (repetition #", rep, ")")
                                  )
              PlotFunction(subsampled_list[[sample_level_index]][[rep]][[ccs_name]][[df_name]],
                           main_title = use_title, reorder_wells = reorder_wells
                           )
            }
          }

          dev.off()
        }
      }
    }
  }
}



DrawAllPlotsForOnePlate <- function(PlotFunction,
                                    df_list_list,
                                    use_dir,
                                    folder_prefix,
                                    file_prefix,
                                    file_postfix = "",
                                    draw_PDFs = TRUE,
                                    use_ccs3_title = ccs3_title,
                                    use_ccs5_title = ccs5_title,
                                    use_ccs7_title = ccs7_title
                                    ) {

  for (reorder_wells in c(FALSE, TRUE)) {

    order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
    order_folder <- paste0(folder_prefix, " - ", order_folder)
    plots_dir <- file.path(use_dir, order_folder)

    ccs3_df_list <- df_list_list[["ccs3"]]
    ccs5_df_list <- df_list_list[["ccs5"]]
    ccs7_df_list <- df_list_list[["ccs7"]]


    ## Draw the accuracy plots in the console

    PlotFunction(ccs3_df_list[["original_summary_df"]],
                 main_title = use_ccs3_title, reorder_wells = reorder_wells
                 )
    PlotFunction(ccs3_df_list[["filtered_summary_df"]],
                 main_title = use_ccs3_title, reorder_wells = reorder_wells
                 )
    PlotFunction(ccs3_df_list[["filtered_gRNAs_df"]],
                 main_title = use_ccs5_title, reorder_wells = reorder_wells
                 )

    PlotFunction(ccs5_df_list[["original_summary_df"]],
                 main_title = use_ccs5_title, reorder_wells = reorder_wells
                 )
    PlotFunction(ccs5_df_list[["filtered_summary_df"]],
                 main_title = use_ccs5_title, reorder_wells = reorder_wells
                 )
    PlotFunction(ccs5_df_list[["filtered_gRNAs_df"]],
                 main_title = use_ccs5_title, reorder_wells = reorder_wells
                 )

    PlotFunction(ccs7_df_list[["original_summary_df"]],
                 main_title = use_ccs7_title, reorder_wells = reorder_wells
                 )
    PlotFunction(ccs7_df_list[["filtered_summary_df"]],
                 main_title = use_ccs7_title, reorder_wells = reorder_wells
                 )
    PlotFunction(ccs7_df_list[["filtered_gRNAs_df"]],
                 main_title = use_ccs7_title, reorder_wells = reorder_wells
                 )


    ## Produce the accuracy PNGs

    SavePNG <- function(summary_df, file_name, main_title) {
      full_path <- file.path(plots_dir, paste0(file_name, ".png"))
      assign("delete_full_path", full_path, envir = globalenv())
      png(filename = full_path,
          res      = 600,
          height   = use_height,
          width    = use_width,
          units    = "in"
          )
      PlotFunction(summary_df, main_title = main_title,
                   reorder_wells = reorder_wells
                   )
      dev.off()
    }


    SavePNG(ccs7_df_list[["original_summary_df"]],
            paste0(file_prefix, "CCS7 (99.99) - i) unfiltered", file_postfix),
            main_title = use_ccs7_title
            )
    SavePNG(ccs7_df_list[["filtered_summary_df"]],
            paste0(file_prefix, "CCS7 (99.99) - ii) filtered", file_postfix),
            main_title = use_ccs7_title
            )
    SavePNG(ccs7_df_list[["filtered_gRNAs_df"]],
            paste0(file_prefix, "CCS7 (99.99) - iii) filtered gRNAs", file_postfix),
            main_title = use_ccs7_title
            )

    SavePNG(ccs5_df_list[["original_summary_df"]],
            paste0(file_prefix, "CCS5 (99.9) - i) unfiltered", file_postfix),
            main_title = use_ccs5_title
            )
    SavePNG(ccs5_df_list[["filtered_summary_df"]],
            paste0(file_prefix, "CCS5 (99.9) - ii) filtered", file_postfix),
            main_title = use_ccs5_title
            )
    SavePNG(ccs5_df_list[["filtered_gRNAs_df"]],
            paste0(file_prefix, "CCS5 (99.9) - iii) filtered gRNAs", file_postfix),
            main_title = use_ccs5_title
            )

    SavePNG(ccs3_df_list[["original_summary_df"]],
            paste0(file_prefix, "CCS3 (99) - i) unfiltered", file_postfix),
            main_title = use_ccs3_title
            )
    SavePNG(ccs3_df_list[["filtered_summary_df"]],
            paste0(file_prefix, "CCS3 (99) - ii) filtered", file_postfix),
            main_title = use_ccs3_title
            )
    SavePNG(ccs3_df_list[["filtered_gRNAs_df"]],
            paste0(file_prefix, "CCS3 (99) - iii) filtered gRNAs", file_postfix),
            main_title = use_ccs3_title
            )


    ## Produce the accuracy PDF

    if (draw_PDFs) {
      pdf(file = file.path(plots_dir, paste0(file_prefix, "i) unfiltered", file_postfix, ".pdf")),
          height = use_height,
          width  = use_width
          )
      PlotFunction(ccs7_df_list[["original_summary_df"]],
                   main_title = use_ccs7_title, reorder_wells = reorder_wells
                   )
      PlotFunction(ccs5_df_list[["original_summary_df"]],
                   main_title = use_ccs5_title, reorder_wells = reorder_wells
                   )
      PlotFunction(ccs3_df_list[["original_summary_df"]],
                   main_title = use_ccs3_title, reorder_wells = reorder_wells
                   )
      dev.off()


      pdf(file = file.path(plots_dir, paste0(file_prefix, "ii) filtered", file_postfix, ".pdf")),
          height = use_height,
          width  = use_width
          )
      PlotFunction(ccs7_df_list[["filtered_summary_df"]],
                   main_title = use_ccs7_title, reorder_wells = reorder_wells
                   )
      PlotFunction(ccs5_df_list[["filtered_summary_df"]],
                   main_title = use_ccs5_title, reorder_wells = reorder_wells
                   )
      PlotFunction(ccs3_df_list[["filtered_summary_df"]],
                   main_title = use_ccs3_title, reorder_wells = reorder_wells
                   )
      dev.off()


      pdf(file = file.path(plots_dir, paste0(file_prefix, "iii) filtered gRNAs", file_postfix, ".pdf")),
          height = use_height,
          width  = use_width
          )
      PlotFunction(ccs7_df_list[["filtered_gRNAs_df"]],
                   main_title = use_ccs7_title, reorder_wells = reorder_wells
                   )
      PlotFunction(ccs5_df_list[["filtered_gRNAs_df"]],
                   main_title = use_ccs5_title, reorder_wells = reorder_wells
                   )
      PlotFunction(ccs3_df_list[["filtered_gRNAs_df"]],
                   main_title = use_ccs3_title, reorder_wells = reorder_wells
                   )
      dev.off()

    }
  }
  return(invisible(NULL))
}





# Functions used for multi-plate experiments ------------------------------


# Define constants for sparse graphics
margin_ratio <- 10
sparse_width <- 3.0 * (1 + (2 / margin_ratio))
sparse_height <- 1.5 * (2 + (3 / margin_ratio)) - (0.2 * 2)



DrawBarplotsAndHeatmapsForAllPlates <- function(export_PNGs = TRUE) {

  required_objects <- c("use_plate_numbers", "plates_df", "library_df")

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
    title_colors <- plates_df[, "Highlight_color"]
  } else {
    title_colors <- rep("black", nrow(plates_df))
  }

  for (file_format in file_formats) {
    for (reorder_wells in c(FALSE, TRUE)) {

      if (reorder_wells) {
        message("Showing re-ordered wells...")
      } else {
        message("Showing wells in their original order...")
      }

      order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
      heatmaps_folder   <- paste0("Heatmaps - ", order_folder)
      barplots_folder   <- paste0("Stacked barplots - ", order_folder)
      sand_charts_folder <- "Sand charts"
      sparse_folder <- "Heatmaps - original - sparse"

      for (i in seq_along(ccs_numbers)) {

        use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))

        filter_stages <- c("original_summary_df", "filtered_summary_df", "filtered_cross_plate_df")
        filter_labels <- c("i) unfiltered", "ii) filtered", "iii) filtered cross-plate")
        df_are_present <- filter_stages %in% names(use_df_list)
        filter_stages <- filter_stages[df_are_present]
        filter_labels <- filter_labels[df_are_present]

        titles_list <- lapply(plates_df[["Plate_number"]], function(x) {
          bquote({"Plate #" * .(as.character(x)) * " \u2013 " *
              bold(.(plates_df[plates_df[["Plate_number"]] == x, "Plate_name"])) *
                   " (" >= .(as.character(ccs_numbers[[i]])) * " consensus reads "} *
                   "and " >= .(as.character(accuracy_percentages[[i]])) * "% accuracy)"
                   )
        })

        for (filter_stage in seq_along(filter_stages)) {

          df_name <- filter_stages[[filter_stage]] # "filtered_gRNAs_df"
          use_summary_df <- use_df_list[[df_name]]

          sel_name <- paste0("CCS", ccs_numbers[[i]],
                             " (", accuracy_percentages[[i]], ") - ",
                             filter_labels[[filter_stage]]
                             )

          message(paste0("Exporting ", file_format, " images for the following data: ", sel_name, "..."))

          if (file_format == "pdf") {
            pdf(file = file.path(plots_output_directory, heatmaps_folder, paste0("Heatmaps - ", sel_name, ".pdf")),
                width = use_width, height = use_height
                )
          } else if (file_format == "png") {
            sub_folder_path <- file.path(PNGs_output_directory, heatmaps_folder, sel_name)
            dir.create(sub_folder_path, showWarnings = FALSE)
          }
          for (plate_number in use_plate_numbers) {
            sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
            sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
            assign("sg_sequences_df", sg_sequences_df, envir = globalenv())
            if (file_format == "png") {
              file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
              png(file   = file.path(sub_folder_path, file_name),
                  width  = use_width,
                  height = use_height,
                  units  = "in",
                  res    = 600
                  )
            }
            plate_index <- which(plates_df[["Plate_number"]] == plate_number)
            DrawAccuracyHeatmap(sub_df,
                                reorder_wells = reorder_wells,
                                main_title    = titles_list[[plate_index]],
                                title_color   = title_colors[[plate_index]]
                                )
            if (file_format == "png") {
              dev.off()
            }
          }
          if (file_format == "pdf") {
            dev.off()
          }

          if (reorder_wells) {
            if (file_format == "pdf") {
              pdf(file = file.path(plots_output_directory, sand_charts_folder, paste0("Sand chart - ", sel_name, ".pdf")),
                  width = 3.8, height = 6.7
                  )
            } else if (file_format == "png") {
              sub_folder_path <- file.path(PNGs_output_directory, sand_charts_folder, sel_name)
              dir.create(sub_folder_path, showWarnings = FALSE)
            }
            for (plate_number in use_plate_numbers) {
              sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
              sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
              assign("sg_sequences_df", sg_sequences_df, envir = globalenv())
              if (file_format == "png") {
                file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
                png(file   = file.path(sub_folder_path, file_name),
                    width  = 3.8,
                    height = 6.7,
                    units  = "in",
                    res    = 600
                    )
              }
              plate_index <- which(plates_df[["Plate_number"]] == plate_number)
              DrawReorderedSandPlots(sub_df,
                                     main_title  = titles_list[[plate_index]],
                                     title_color = title_colors[[plate_index]]
                                     )
              if (file_format == "png") {
                dev.off()
              }
            }
            if (file_format == "pdf") {
              dev.off()
            }
          } else {

            if (file_format == "pdf") {
              pdf(file = file.path(plots_output_directory, sparse_folder, paste0("Heatmap - ", sel_name, ".pdf")),
                  width = sparse_width, height = sparse_height
                  )
            } else if (file_format == "png") {
              sub_folder_path <- file.path(PNGs_output_directory, sparse_folder, sel_name)
              dir.create(sub_folder_path, showWarnings = FALSE)
            }
            for (plate_number in use_plate_numbers) {
              sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
              sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
              assign("sg_sequences_df", sg_sequences_df, envir = globalenv())
              if (file_format == "png") {
                file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
                png(file   = file.path(sub_folder_path, file_name),
                    width  = sparse_width,
                    height = sparse_height,
                    units  = "in",
                    res    = 600
                    )
              }
              plate_index <- which(plates_df[["Plate_number"]] == plate_number)
              SparseAccuracyHeatmap(sub_df,
                                    main_title  = titles_list[[plate_index]],
                                    title_color = title_colors[[plate_index]]
                                    )
              if (file_format == "png") {
                dev.off()
              }
            }
            if (file_format == "pdf") {
              dev.off()
            }
          }

          alterations_factor <- 0.7
          alterations_width <- use_width * alterations_factor
          alterations_height <- use_height * alterations_factor
          if (file_format == "pdf") {
            pdf(file = file.path(plots_output_directory, barplots_folder, paste0("Stacked barplots - ", sel_name, ".pdf")),
                height = alterations_height, width = alterations_width
                )
          } else if (file_format == "png") {
            sub_folder_path <- file.path(PNGs_output_directory, barplots_folder, sel_name)
            dir.create(sub_folder_path, showWarnings = FALSE)
          }
          for (plate_number in use_plate_numbers) {
            sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
            sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
            assign("sg_sequences_df", sg_sequences_df, envir = globalenv())
            if (file_format == "png") {
              file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
              png(file   = file.path(sub_folder_path, file_name),
                  width  = alterations_width,
                  height = alterations_height,
                  units  = "in",
                  res    = 600
                  )
            }
            plate_index <- which(plates_df[["Plate_number"]] == plate_number)
            ModifiedAlterationBarplot(sub_df,
                                      reorder_wells = reorder_wells,
                                      main_title    = titles_list[[plate_index]],
                                      title_color   = title_colors[[plate_index]]
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
      message("")
    }
  }
  return(invisible(NULL))
}












