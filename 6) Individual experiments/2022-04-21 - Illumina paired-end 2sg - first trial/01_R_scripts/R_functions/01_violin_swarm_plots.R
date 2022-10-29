## 2022-06-04



# Load packages and source code -------------------------------------------

library("RColorBrewer")
library("beeswarm")
library("sm")



# General functions for plotting ------------------------------------------

SetUpBoxPlot <- function(num_groups,
                         data_range,
                         use_y_limits  = NULL,
                         draw_axis     = TRUE,
                         draw_box      = TRUE,
                         indicate_zero = FALSE
                         ) {

  ## Determine group positions
  group_positions <- seq_len(num_groups)
  side_gap <- 0.5
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  if (is.null(use_y_limits)) {
    y_space <- (data_range[[2]] - data_range[[1]]) * 0.02
    use_y_limits <- c(data_range[[1]] - y_space, data_range[[2]] + y_space)
    if (use_y_limits[[1]] > 0) {
      use_y_limits[[1]] <- 0
    }
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

  if (indicate_zero && (use_y_limits[[1]] < 0)) {
    abline(h = 0, lty = "dotted", col = "gray80")
  }

  if (draw_axis) {
    axis(2,
         mgp = c(3, 0.55, 0),
         tcl = -0.375,
         las = 1,
         lwd = par("lwd")
         )
  }

  if (draw_box) {
    box(bty = "l", lwd = par("lwd"))
  }

  return(invisible(NULL))
}



RepositionByGroups <- function(groups_vec, gap_ratio = 1.25) {
  stopifnot(gap_ratio > 1)
  group_lengths <- rle(groups_vec)[["lengths"]]
  rle_vec <- rep(seq_along(group_lengths), group_lengths)
  group_span <- length(groups_vec) - 1L
  groups_list <- split(seq_along(rle_vec), rle_vec)
  are_first <- unlist(lapply(groups_list, function(x) {
    c(TRUE, rep(FALSE, length(x) - 1))
  }), use.names = FALSE)
  are_large_gap <- are_first[-1]
  small_gap_size <- group_span / ((sum(are_large_gap) * gap_ratio) + sum(!(are_large_gap)))
  large_gap_size <- small_gap_size * gap_ratio
  gaps_vec <- ifelse(are_large_gap, large_gap_size, small_gap_size)
  positions_vec <- c(1, cumsum(gaps_vec) + 1)
  stopifnot(abs(positions_vec[[length(groups_vec)]] - length(groups_vec)) <
            (.Machine$double.eps * length(groups_vec))
            )
  return(positions_vec)
}


GetHalfLineWidth <- function(y_axis = FALSE) {
  if (y_axis) {
    diff(grconvertY(c(0, par("lwd") / 96), from = "inches", to = "user")) / 2
  } else {
    diff(grconvertX(c(0, par("lwd") / 96), from = "inches", to = "user")) / 2
  }
}





# Functions for drawing violin plots --------------------------------------

SingleViolin <- function(numeric_vec,
                         at,
                         wex            = 0.8,
                         show_quantiles = c(0.25, 0.5, 0.75),
                         quantiles_lty  = ifelse(show_quantiles == 0.5, "longdash", "dashed"),
                         trim           = TRUE,
                         show_quartiles = TRUE,
                         use_sm_density = FALSE,
                         adjust         = 1,
                         violin_color   = "#9ECAE1",
                         line_color     = "black",
                         border_color   = line_color,
                         draw_border    = FALSE
                         ) {

  if (length(quantiles_lty) == 1) {
    quantiles_lty <- rep(quantiles_lty, length(show_quantiles))
  }

  are_NA <- is.na(numeric_vec)
  numeric_vec <- numeric_vec[!(are_NA)]

  data_limits <- range(numeric_vec)
  draw_lines <- length(show_quantiles) != 0
  if (draw_lines) {
    lines_vec <- quantile(numeric_vec, probs = show_quantiles)
  }
  if (use_sm_density) {
    use_args <- list(x = numeric_vec, display = "none")
    if (trim) {
      use_args <- c(use_args, list(xlim = data_limits))
    }
    sm_output_polygon <- do.call(sm::sm.density, use_args)
    estimates_vec <- sm_output_polygon[["estimate"]]
    eval_points_vec <- sm_output_polygon[["eval.points"]]
    if (draw_lines) {
      sm_output_lines <- do.call(sm::sm.density, c(use_args, list(eval.points = lines_vec)))
      line_estimates_vec <- sm_output_lines[["estimate"]]
    }
  } else {
    density_output <- density(numeric_vec, adjust = adjust)
    estimates_vec <- density_output[["y"]]
    eval_points_vec <- density_output[["x"]]
    if (draw_lines) {
      dens_vec <- cumsum(estimates_vec) / sum(estimates_vec)
      ecdf_Function <- stats::approxfun(dens_vec, eval_points_vec, ties = "ordered")
      ys <- ecdf_Function(show_quantiles) # these are all the y-values for quantiles
      ys <- quantile(numeric_vec, show_quantiles)
      # Get the violin bounds for the requested quantiles.
      line_estimates_vec <- (stats::approxfun(eval_points_vec, estimates_vec))(ys)
    }
    if (trim) {
      are_within_bounds <- (eval_points_vec >= data_limits[[1]]) &
                           (eval_points_vec <= data_limits[[2]])
      eval_points_vec <- eval_points_vec[are_within_bounds]
      estimates_vec <- estimates_vec[are_within_bounds]
    }
  }

  h_scale <- 0.5 / max(estimates_vec) * wex
  heights_vec <- estimates_vec * h_scale
  x_vec <- c(at - heights_vec, rev(at + heights_vec))
  y_vec <- c(eval_points_vec, rev(eval_points_vec))
  polygon(x_vec, y_vec, col = violin_color,
          border = if (draw_border) border_color else NA, xpd = NA
          )

  if (draw_lines) {
    line_heights_vec <- line_estimates_vec * h_scale - GetHalfLineWidth()
    for (i in seq_along(show_quantiles)) {
      lines(x    = c(at - line_heights_vec[[i]], at + line_heights_vec[[i]]),
            y    = rep(lines_vec[[i]], 2),
            lty  = quantiles_lty[[i]],
            col  = line_color,
            lend = "butt",
            xpd  = NA
            )
    }
  }
  return(invisible(NULL))
}



BringWithinLimits <- function(numeric_vec, lower_bound = NULL, upper_bound = NULL) {

  stopifnot(!(anyNA(numeric_vec)))

  are_finite <- is.finite(numeric_vec)
  data_range <- range(numeric_vec[are_finite])
  if (any(!(are_finite))) {
    numeric_vec[!(are_finite)] <- ifelse(numeric_vec[!(are_finite)] < 0,
                                         if (is.null(lower_bound)) data_range[[1]] else lower_bound - 1,
                                         if (is.null(upper_bound)) data_range[[2]] else upper_bound + 1
                                         )
    message(paste0(sum(!(are_finite)), " values were not finite (i.e. Inf or ",
                   "-Inf) and were replaced by either the upper or lower ",
                   "bounds (if given) or the maximum or minimum values in the ",
                   "data, respectively."
                   ))
  }
  if ((!(is.null(upper_bound))) || (!(is.null(lower_bound)))) {
    if (is.null(upper_bound)) {
      upper_bound <- Inf
    }
    if (is.null(lower_bound)) {
      lower_bound <- -Inf
    }
    exceed_lower <- numeric_vec < lower_bound
    exceed_upper <- numeric_vec > upper_bound
    exceed_bounds <- exceed_lower | exceed_upper
    numeric_vec[exceed_bounds] <- ifelse(numeric_vec[exceed_bounds] < lower_bound,
                                         lower_bound, upper_bound
                                         )
    message(paste0(sum(exceed_bounds), " values lay outside of the given ",
                   "upper or lower bounds and were truncated."
                   ))
    lower_bound_enforced <- any(exceed_lower)
    upper_bound_enforced <- any(exceed_upper)
  } else {
    lower_bound_enforced <- FALSE
    upper_bound_enforced <- FALSE
  }

  results_list <- list(
    "curtailed_vec"        = numeric_vec,
    "lower_bound_enforced" = lower_bound_enforced,
    "upper_bound_enforced" = upper_bound_enforced
  )
  return(results_list)
}



CurtailedAxisLabels <- function(tick_positions, upper_bound, lower_bound,
                                lower_bound_enforced, upper_bound_enforced
                                ) {
  plain_axis_labels <- format(tick_positions, trim = TRUE)
  axis_labels <- lapply(plain_axis_labels, function(x) as.expression(bquote(""[.(x)])))
  axis_labels <- sapply(axis_labels, function(x) x)
  if (lower_bound_enforced && (abs(lower_bound - tick_positions[[1]]) < 1e-15)) {
    axis_labels[1] <- as.expression(bquote(""["" <= scriptscriptstyle(.(if (tick_positions[[1]] < 0) "" else " "))
                                              * .(plain_axis_labels[[1]])]
    ))
  }
  if (upper_bound_enforced) {
    num_ticks <- length(tick_positions)
    if (abs(upper_bound - tick_positions[[num_ticks]]) < 1e-15) {
      axis_labels[num_ticks] <- as.expression(bquote(""["" >= scriptscriptstyle(" ") *
                                                          .(plain_axis_labels[[num_ticks]])]
      ))
    }
  }
  return(axis_labels)
}



BeeViolinPlot <- function(input_list,
                          groups_vec      = NULL,
                          gap_ratio       = 1.25,
                          brewer_pals     = c("Blues", "Reds", "Purples",
                                              "Greens", "Oranges", "Greys"
                                              ),
                          point_colors    = NULL,
                          violin_colors   = NULL,
                          line_colors     = NULL,
                          border_colors   = NULL,
                          point_cex       = 0.4,
                          use_spacing     = 0.8,
                          y_limits        = NULL,
                          lower_bound     = NULL,
                          upper_bound     = NULL,
                          indicate_zero   = TRUE,
                          draw_groups_n   = TRUE,
                          draw_points     = TRUE,
                          use_swarm       = TRUE,
                          axis_cex_factor = 1 / 0.7,
                          cloud_alpha     = 0.2,
                          cloud_sd        = 0.04,
                          embed_PNG       = FALSE,
                          ...
                          ) {

  num_groups <- length(input_list)
  if (is.null(violin_colors)) {
    violin_colors <- vapply(brewer_pals, function(x) brewer.pal(9, x)[[3]], "")
  } else if (length(violin_colors) == 1) {
    violin_colors <- rep(violin_colors, num_groups)
  }
  if (is.null(point_colors)) {
    point_colors <- vapply(brewer_pals, function(x) brewer.pal(9, x)[[7]], "")
  } else if (length(point_colors) == 1) {
    point_colors <- rep(point_colors, num_groups)
  }
  if (is.null(line_colors)) {
    line_colors <- point_colors
  } else if (length(line_colors) == 1) {
    line_colors <- rep(line_colors, num_groups)
  }
  if (is.null(border_colors)) {
    border_colors <- line_colors
  } else if (length(border_colors) == 1) {
    border_colors <- rep(border_colors, num_groups)
  }

  group_indices <- rep(seq_along(input_list), lengths(input_list))

  numeric_vec <- unlist(input_list, use.names = FALSE)
  are_finite <- is.finite(numeric_vec)
  curtailed_list <- BringWithinLimits(numeric_vec,
                                      lower_bound = lower_bound,
                                      upper_bound = upper_bound
                                      )
  numeric_vec <- curtailed_list[["curtailed_vec"]]
  lower_bound_enforced <- curtailed_list[["lower_bound_enforced"]]
  upper_bound_enforced <- curtailed_list[["upper_bound_enforced"]]
  assign("delete_curtailed_list", curtailed_list, envir = globalenv())
  numeric_list <- split(numeric_vec, group_indices)
  are_finite_list <- split(are_finite, group_indices)

  if (is.null(groups_vec)) {
    x_positions <- seq_len(num_groups)
  } else {
    x_positions <- RepositionByGroups(groups_vec, gap_ratio = gap_ratio)
  }

  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(figures_dir, "temp.png")
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

  SetUpBoxPlot(num_groups,
               data_range    = range(numeric_vec),
               use_y_limits  = y_limits,
               draw_axis     = FALSE,
               draw_box      = !(embed_PNG),
               indicate_zero = indicate_zero
               )

  for (i in seq_along(numeric_list)) {
    SingleViolin(numeric_list[[i]][are_finite_list[[i]]],
                 at            = x_positions[[i]],
                 violin_color  = violin_colors[[i]],
                 line_color    = line_colors[[i]],
                 border_color  = border_colors[[i]],
                 quantiles_lty = c("21", "52", "21"),
                 ...
                 )
  }

  if (draw_points) {
    set.seed(1)
    if (use_swarm) {
      beeswarm_df <- beeswarm::beeswarm(numeric_list,
                                        at       = x_positions,
                                        priority = "random",
                                        spacing  = use_spacing,
                                        cex      = point_cex,
                                        col      = point_colors,
                                        do.plot  = FALSE
                                        )

      points(beeswarm_df[, "x"],
             beeswarm_df[, "y"],
             pch = 16,
             col = beeswarm_df[, "col"],
             cex = point_cex,
             xpd = NA
             )

    } else {
      x_vec <- rep(x_positions, lengths(numeric_list))
      x_vec <- x_vec + rnorm(n = length(x_vec), mean = 0, sd = cloud_sd)
      cloud_colors <- adjustcolor(point_colors, alpha.f = cloud_alpha)
      points(x = x_vec, y = unlist(numeric_list),
             col = cloud_colors[rep(seq_along(numeric_list), lengths(numeric_list))],
             cex = point_cex, pch = 16, xpd = NA
             )
    }
  }

  if (embed_PNG) {
    dev.off()
    raster_array <- png::readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)
    SetUpBoxPlot(num_groups,
                 data_range    = range(numeric_vec),
                 use_y_limits  = y_limits,
                 draw_axis     = FALSE,
                 indicate_zero = FALSE
                 )
    rasterImage(raster_array,
                xleft   = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop   = par("usr")[[4]]
                )
  }

  ## Annotate plot
  axis_ticks <- axTicks(2)
  axis_labels <- CurtailedAxisLabels(axis_ticks,
                                     lower_bound          = lower_bound,
                                     upper_bound          = upper_bound,
                                     lower_bound_enforced = lower_bound_enforced,
                                     upper_bound_enforced = upper_bound_enforced
                                     )
  axis(2,
       at       = axis_ticks,
       labels   = axis_labels,
       mgp      = c(3, 0.5, 0),
       tcl      = -0.35,
       las      = 1,
       lwd      = par("lwd"),
       cex.axis = par("cex") * axis_cex_factor
       )
  if (draw_groups_n) {
    mtext(lengths(numeric_list), at = x_positions, side = 1, line = -0.45,
          cex = par("cex") * 0.3, col = "gray70", padj = 0
          )
  }

  return(invisible(x_positions))
}



