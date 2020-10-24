### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("squash")
library("scales")






# Define plot dimensions --------------------------------------------------

use_height <- 7
use_width <- 6.5





# Define plot titles ------------------------------------------------------

ccs3_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "3 consensus reads "} *
                               "and " >= "99% accuracy)"
                               ))

ccs5_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "5 consensus reads "} *
                               "and " >= "99.9% accuracy)"
                               ))






# Define functions --------------------------------------------------------

MakeEmptyPlot <- function(y_limits = c(0, 1)) {
  plot(1, type = "n", axes = FALSE, ann = FALSE,
       xlim = c(0, 1), ylim = y_limits, xaxs = "i", yaxs = "i"
       )
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
      polygon(x      = c(bottom_left_corners_x[[i]], top_left_corners_x[[i]], top_right_corners_x[[i]], bottom_right_corners_x[[i]]),
              y      = c(bottom_left_corners_y[[i]], top_left_corners_y[[i]], top_right_corners_y[[i]], bottom_right_corners_y[[i]]),
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
  cimage(zcol = input_matrix, axes = FALSE, xlab = "", ylab = "")
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
  assign("delete_expression_list", expression_list, envir = globalenv())
  literal_strings <- vapply(expression_list, StripExpression, "")
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}

VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("g")), use_expression, expression(phantom("h")))
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
  assign("delete_color_matrix", color_matrix, envir = globalenv())

  if (ncol(color_matrix) != length(groups_vec)) {
    stop("The number of columns of the color matrix must match the length of the group labels!")
  }
  num_groups <- length(unique(groups_vec))
  whitespace_bar <- as.matrix(rep(whitespace_color, nrow(color_matrix)))
  colnames(whitespace_bar) <- "divider"
  whitespace_bar <- whitespace_bar[, rep(1, divider_width)]
  indices_list <- lapply(unique(groups_vec), function(x) which(x == groups_vec))
  assign("delete_groups_vecs", groups_vec, envir = globalenv())
  assign("delete_color_matrix", color_matrix, envir = globalenv())
  assign("delete_indices_list", indices_list, envir = globalenv())
  assign("delete_whitespace_bar", whitespace_bar, envir = globalenv())
  assign("delete_divider_width", divider_width, envir = globalenv())
  matrices_list <- do.call(c, lapply(seq_len(num_groups - 1), function(x) list(color_matrix[, indices_list[[x]], drop = FALSE], whitespace_bar)))
  matrices_list <- c(matrices_list, list(color_matrix[, indices_list[[num_groups]], drop = FALSE]))
  assign("delete_matrices_list", matrices_list, envir = globalenv())
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




DrawAccuracyHeatmap <- function(summary_df,
                                main_title             = NULL,
                                ColorFunction          = colorRampPalette(rev(viridis::magma(100))), # colorRampPalette(brewer.pal(9, "YlGnBu")),
                                reorder_wells          = TRUE,
                                show_zero_in_legend    = FALSE,
                                gap_weight             = 2L,
                                show_correct_promoters = FALSE
                                ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))


  percent_columns <- paste0("Perc_sg", 1:4, "_cr", 1:4)

  if (reorder_wells) {
    mean_accuracies <- rowMeans(as.matrix(summary_df[, percent_columns]))
    new_order <- order(summary_df[["Perc_all_4"]],
                       summary_df[["Perc_at_least_3"]],
                       summary_df[["Perc_at_least_2"]],
                       summary_df[["Perc_at_least_1"]],
                       mean_accuracies
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
  rownames(summary_df) <- NULL

  num_wells <- length(use_indices)

  ## Set up the layout

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

  par(mar = rep(0, 4))

  for (i in 1:3) {
    MakeEmptyPlot()
  }
  if (!(is.null(main_title))) {
    text(x = 0.5,
         y = 0.65,
         labels = main_title,
         cex = 1.1,
         xpd = NA
         )

  }
  MakeEmptyPlot()


  ## Draw a vertical barplot

  are_at_least_75 <- summary_df[["Perc_all_4"]] >= 75
  num_above_75 <- sum(are_at_least_75)
  over_75_fraction <- num_above_75 / num_wells


  # six_colors <- c(colorRampPalette(brewer.pal(11, "PuOr"))(21)[[11]],
  #                 colorRampPalette(brewer.pal(11, "PuOr"))(20)[c(8, 5, 3)],
  #                 colorRampPalette(brewer.pal(11, "PuOr"))(21)[[17]],
  #                 colorRampPalette(brewer.pal(11, "RdBu"))(21)[[20]]
  #                 )

  # lighter_color <- "#ECEBE9"
  # darker_color <- toupper("#15396d")
  # other_color <- "#6E1566"
  # other_color <- "#380030"
  # five_colors <- c(colorRampPalette(c(lighter_color, darker_color))(5)[2:5],
  #                  other_color #colorRampPalette(c(darker_color, other_color))(5)[[4]]
  #                  )


  add_gap <- (!(reorder_wells)) && ("Block" %in% names(sg_sequences_df))
  if (add_gap) {
    block_vec <- sg_sequences_df[["Block"]][are_to_include]
    positions_vec <- GappedPositionsVec(block_vec, gap_weight = gap_weight)
  } else {
    positions_vec <- seq_len(num_wells)
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


  for (i in seq_len(num_wells)) {
    rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
         xright  = (positions_vec[[i]] / max(positions_vec)),
         ybottom = 0,
         ytop    = 1,
         col     = barplot_colors[[1]],
         border  = NA
         )
    for (j in seq_along(barplot_columns)) {
      rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
           xright  = (positions_vec[[i]] / max(positions_vec)),
           ybottom = 0,
           ytop    = summary_df[[barplot_columns[[j]]]][[i]] / 100,
           col     = barplot_colors[2:length(barplot_colors)][[j]],
           border  = NA
           )
    }
  }


  tick_locations <- axTicks(2)
  tick_labels <- paste0(tick_locations * 100, "%")
  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.3,
       lwd      = 0.75,
       cex.axis = 0.9
       )
  box(lwd = 0.5, xpd = NA)


  if (show_correct_promoters) {
    color_text_vec <- c('bold(color1("% plasmids with ") * color1("" >= "") *',
                        'scriptscriptstyle(" ") *',
                        'color2("0,") * ',
                        'color1(" ") * color3("1,") * ',
                        'color1(" ") * color4("2,") * ',
                        'color1(" ") * color5("3,") * ',
                        'color1(" ") * color6("4") * ',
                        'color1(" correct gRNAs ") * ',
                        'color7("or an entirely correct construct"))'
                        )
  } else {
    color_text_vec <- c('bold(color1("% plasmids with ") * color1("" >= "") *',
                        'scriptscriptstyle(" ") *',
                        'color2("0,") * ',
                        'color1(" ") * color3("1,") * ',
                        'color1(" ") * color4("2,") * ',
                        'color1(" ") * color5("3,") * ',
                        'color1(" or ") * color6("4") * ',
                        'color1(" correct gRNAs (including the tracrRNA)"))'
                        )
  }


  text_colors <- c("#000000",
                   Darken(barplot_colors[[1]], factor = 1.3),
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
         y      = 1.05,
         labels = VerticalAdjust(parse(text = MakeInvisible(color_text, j))),
         adj    = c(0.5, 0),
         xpd    = NA,
         col    = text_colors[[j]]
         )
  }


  ## Draw the % accuracy horizontal barplot

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
      rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
           xright  = (positions_vec[[i]] / max(positions_vec)),
           ybottom = 0,
           ytop    = 1,
           col     = two_grey_colors[[as.integer(are_at_least_75[[i]]) + 1]],
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

  numeric_mat <- t(as.matrix(summary_df[, percent_columns])) / 100

  my_breaks <- seq(0, 1.01, by = 0.01)
  my_breaks[c(1, length(my_breaks))] <- c(0, ceiling(my_breaks[length(my_breaks)]))

  my_cmap <- makecmap(numeric_mat, colFn = ColorFunction, breaks = my_breaks)
  my_color_mat <- cmap(numeric_mat, my_cmap)


  assign("delete_my_color_mat", my_color_mat, envir = globalenv())

  if (add_gap) {
    my_color_mat <- AddVerticalWhiteSpace(my_color_mat,
                                          groups_vec = block_vec,
                                          divider_width = gap_weight
                                          )
  }

  Do_cimage(my_color_mat)
  x_range <- par("usr")[[2]] - par("usr")[[1]]

  text(x      = par("usr")[[1]] - (x_range * 0.018),
       y      = 4:1,
       labels = paste0("sg", 1:4),
       xpd    = NA,
       adj    = c(1, 0.5)
       )

  MakeEmptyPlot()

  ## Draw the color indicator (trapezoid)

  trapezoid_start_x <- 0.85
  trapezoid_start_y <- 0.65
  trapezoid_end_y   <- 0.9

  DrawColorTrapezoid(start_x     = trapezoid_start_x,
                     end_x       = 1,
                     start_y     = trapezoid_start_y,
                     end_y       = trapezoid_end_y,
                     trapezoid_y = 0.72,
                     use_colors  = ColorFunction(100)
                     )

  text(x      = trapezoid_start_x - 0.009,
       y      = trapezoid_start_y + ((trapezoid_end_y - trapezoid_start_y) * 0.3),
       adj    = c(1, 0.5),
       labels = "Accuracy",
       xpd    = NA,
       font   = 2,
       cex    = 0.9
       )

  trapezoid_seq <- seq(trapezoid_start_x, 1, by = ((1 - trapezoid_start_x) / 5))

  text(x      = trapezoid_seq,
       y      = trapezoid_start_y - 0.16,
       labels = paste0(seq(0, 100, by = 20), "%"),
       font   = 2,
       cex    = 0.6,
       xpd    = NA
       )
  segments(x0   = trapezoid_seq,
           x1   = trapezoid_seq,
           y0   = trapezoid_start_y - 0.088,
           y1   = trapezoid_start_y - 0.02,
           xpd  = NA,
           lwd  = 0.5,
           lend = "butt"
           )

  segments(x0   = trapezoid_start_x,
           x1   = 1,
           y0   = trapezoid_start_y - 0.02,
           y1   = trapezoid_start_y - 0.02,
           lwd  = 0.5,
           lend = "butt"
           )


  MakeEmptyPlot()

  ## Draw a vertical barplot

  perc_truncated <- summary_df[["Num_under_2kb"]] / summary_df[["Count_total"]]

  two_colors <- c("gray88", "gray40")

  for (i in seq_len(num_wells)) {
    rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
         xright  = (positions_vec[[i]] / max(positions_vec)),
         ybottom = 0,
         ytop    = 1,
         col     = two_colors[[1]],
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

  MakeEmptyPlot()




  MakeEmptyPlot()

  ## Draw the strip indicating genes with homologies >=8 bp

  longest_subsequence_vec <- sg_sequences_df[["Longest_subsequence"]][use_indices]
  have_homology <- longest_subsequence_vec >= 8

  homology_colors <- two_grey_colors # brewer.pal(8, "Paired")[c(3, 4)]

  for (i in seq_len(num_wells)) {
    rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
         xright  = (positions_vec[[i]] / max(positions_vec)),
         ybottom = 0,
         ytop    = 1,
         col     = homology_colors[[as.integer(have_homology[[i]]) + 1]],
         border  = NA
         )
  }

  text(x      = 0.5,
       y      = -0.65,
       adj    = c(0.5, 0.5),
       labels = bquote(bold(.(as.character(sum(have_homology))) * " / " *
                            .(as.character(num_wells)) * " genes " *
                              "had " >= "8bp homology"
                            )),
       font   = 2,
       col    = "black",
       xpd    = NA
       )

  MakeEmptyPlot()


  return(invisible(NULL))
}



