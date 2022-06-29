# 2022-01-03



# Helper functions for modifying plotmath expressions ---------------------

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


Embolden <- function(my_string) {
  parse(text = paste0("bold(", StripExpression(my_string), ")"))
}


FormatPlotMath <- function(char_vec) {

  if (any(grepl("\"", char_vec, fixed = TRUE))) {
    stop("The FormatPlotMath function does not support input strings with double quotation marks!")
  }

  replace_strings <- c(
    "[Ll]og2"  = "log\"[2] * \"",
    "[Ll]og10" = "log\"[10] * \"",
    "[Dd]elta" = "\" * Delta * \"",
    "PrPc"     = "PrP\"^C * \""
  )

  results_vec <- paste0("phantom(gh) * \"", char_vec, "\" * phantom(gh)")
  for (sub_string in names(replace_strings)) {
    results_vec <- gsub(sub_string, replace_strings[[sub_string]], results_vec)
  }
  results_vec <- sapply(results_vec, function(x) parse(text = x), USE.NAMES = FALSE)
  return(results_vec)
}




# Functions for modifying colors ------------------------------------------

Darken <- function(color, factor = 1.4) {
  # from https://gist.github.com/Jfortin1/72ef064469d1703c6b30
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}


Palify <- function(colors_vec, fraction_pale = 0.5) {
  adjustcolor(colors_vec,
              offset    = c(rep(fraction_pale, 3), 0),
              transform = diag(c(rep(1 - fraction_pale, 3), 1))
              )

}




# Functions for creating plot elements ------------------------------------

MakeEmptyPlot <- function() {
  plot(1, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       type = "n", axes = FALSE, ann = FALSE
       )
}


AbbreviateDataAxis <- function(side = 2, mgp = 0.38, tcl = -0.3) {
  tick_pos <- axTicks(side)
  if (all(tick_pos[-1] >= (5 * 10^5))) {
    axis_labels <- paste0(format(tick_pos / 10^6), "M")
  } else if (all(tick_pos[-1] >= (5 * 10^4))) {
    axis_labels <- paste0(format(tick_pos / 1000), "k")
  } else {
    axis_labels <- format(tick_pos)
  }
  axis(side, labels = axis_labels, at = tick_pos,
       mgp = c(3, mgp, 0), tcl = tcl, las = 1
       )
}


DrawSideLegend <- function(labels_list,
                           use_colors,
                           border_colors        = NULL,
                           use_pch              = 16,
                           use_point_size       = 1.2,
                           lines_x_start        = 0.75,
                           y_mid                = 0.5,
                           small_gap_size       = 1.25,
                           large_gap_multiplier = 1.75,
                           point_x_start        = 0.15
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing the legend
  small_gap <- diff(grconvertY(c(0, small_gap_size), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_multiplier

  if (all(lengths(labels_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(labels_list))
    are_first <- rep(TRUE, length(labels_list))
  } else {
    are_first <- unlist(lapply(labels_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  x_text  <- 1 + diff(grconvertX(c(0, lines_x_start), from = "lines", to = "npc"))
  x_point <- 1 + diff(grconvertX(c(0, lines_x_start + point_x_start), from = "lines", to = "npc"))

  ## Draw the legend
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(labels_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )

  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))

  points(x   = rep(grconvertX(x = x_point, from = "npc", to = "user"), length(labels_list)),
         y   = tapply(y_pos, groups_vec, mean),
         cex = use_point_size,
         pch = use_pch,
         col = if (!(is.null(border_colors))) border_colors else use_colors,
         bg  = use_colors,
         xpd = NA
         )

  return(invisible(NULL))
}


DataAxisLimits <- function(data_vec, space_fraction = 0.04) {

   # If the minimum and maximum are closer together than the distance
   # from the minimum to zero, then the y axis should not include zero
   far_from_zero <- diff(range(data_vec)) < min(data_vec)
   if (far_from_zero) {
      use_limits <- range(data_vec)
   } else {
      use_limits <- range(c(0, data_vec))
   }
   use_space <- diff(use_limits) * space_fraction
   use_limits <- use_limits + (use_space * c(-1, 1))

   # If the minimum value is not negative or close to zero, the y range should stop at zero
   if (!(far_from_zero) && (min(data_vec) >= use_space)) {
      use_limits[[1]] <- 0
   }
   return(use_limits)
}



