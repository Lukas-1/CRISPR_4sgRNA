## 2022-06-04



# Load packages and source code -------------------------------------------

library("RColorBrewer")
library("beeswarm")
library("sm")



# General functions for plotting ------------------------------------------

Darken <- function(color, factor = 1.4) {
  # from https://gist.github.com/Jfortin1/72ef064469d1703c6b30
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}

SetUpBoxPlot <- function(num_groups,
                         data_range,
                         group_limits  = NULL,
                         use_y_limits  = NULL,
                         draw_axis     = TRUE,
                         draw_box      = TRUE,
                         draw_grid     = FALSE,
                         indicate_zero = FALSE,
                         zero_lty      = "dotted",
                         zero_lwd      = 1,
                         zero_color    = "gray80",
                         grid_color    = "gray91",
                         grid_lty      = "solid",
                         grid_lwd      = 1,
                         side_gap      = 0.5,
                         right_gap     = side_gap,
                         num_ticks     = 5,
                         GridFunction  = NULL
                         ) {

  ## Determine x axis limits
  if (is.null(group_limits)) {
    group_positions <- seq_len(num_groups)
    group_limits <- c((min(group_positions) - side_gap)  - (num_groups * 0.04),
                      max(group_positions) + right_gap  + (num_groups * 0.04)
                      )
  }

  ## Determine y axis limits
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

  if (!(is.null(GridFunction))) {
    GridFunction()
  }

  tick_locations <- pretty(axTicks(2), n = num_ticks)
  if (draw_grid) {
    grid_lines <- tick_locations
    if (indicate_zero) {
      segments(x0  = par("usr")[[1]],
               x1  = par("usr")[[2]],
               y0  = 0,
               col = zero_color,
               lwd = par("lwd") * zero_lwd,
               lty = zero_lty,
               xpd = NA
               )
      grid_lines <- grid_lines[grid_lines != 0]
    }
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[2]],
             y0  = grid_lines,
             col = grid_color,
             lwd = par("lwd") * grid_lwd,
             lty = grid_lty,
             xpd = NA
             )
  } else if (indicate_zero && (use_y_limits[[1]] < 0)) {
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[2]],
             y0  = 0,
             col = zero_color,
             lty = zero_lty,
             lwd = par("lwd") * zero_lwd,
             xpd = NA
             )
  }

  if (draw_axis) {
    axis(2,
         at     = tick_locations,
         labels = format(tick_locations, trim = TRUE),
         mgp    = c(3, 0.55, 0),
         tcl    = -0.375,
         las    = 1,
         lwd    = par("lwd")
         )
  }

  if (draw_box) {
    box(bty = "l")
  }
  return(invisible(NULL))
}



RepositionByGroups <- function(groups_vec, gap_ratio = 1.25) {
  stopifnot(gap_ratio >= 1)
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



StartEmbedPNG <- function(temp_dir,
                          png_res           = 900,
                          use_cairo         = FALSE,
                          transparent_bg    = FALSE,
                          add_padding       = FALSE,
                          padding_in_inches = 0.01
                          ) {
  previous_device <- dev.cur()
  current_par <- par(no.readonly = TRUE)
  png_width <- par("pin")[[1]]
  png_height <- par("pin")[[2]]
  if (add_padding) {
    png_width <- png_width + (padding_in_inches * 2)
    png_height <- png_height + (padding_in_inches * 2)
  }
  png(filename = file.path(temp_dir, "temp.png"),
      width    = png_width,
      height   = png_height,
      units    = "in",
      res      = png_res,
      bg       = if (transparent_bg) "transparent" else "white",
      type     = if (use_cairo) "cairo-png" else "windows"
      )
  par(lwd = current_par[["lwd"]])
  par(cex = current_par[["cex"]])
  if (add_padding) {
    par(mai = rep(padding_in_inches, 4))
  } else {
    par(mai = rep(0, 4))
  }
  return(previous_device)
}



StopEmbedPNG <- function(current_device,
                         temp_dir,
                         draw_image        = TRUE,
                         add_padding       = FALSE,
                         padding_in_inches = 0.01,
                         make_empty_plot   = TRUE
                         ) {
  previous_usr <- par("usr")
  dev.off()
  temp_path <- file.path(temp_dir, "temp.png")
  raster_array <- png::readPNG(temp_path)
  file.remove(temp_path)
  dev.set(current_device)
  if (draw_image) {
    if (make_empty_plot) {
      plot(NA, xlim = previous_usr[c(1, 2)], ylim = previous_usr[c(3, 4)],
           xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE
           )
    }
    if (add_padding) {
      x_padding <- diff(grconvertX(c(0, padding_in_inches), from = "inches", to = "user"))
      y_padding <- diff(grconvertY(c(0, padding_in_inches), from = "inches", to = "user"))
      rasterImage(raster_array,
                  xleft   = par("usr")[[1]] - x_padding,
                  xright  = par("usr")[[2]] + x_padding,
                  ybottom = par("usr")[[3]] - y_padding,
                  ytop    = par("usr")[[4]] + y_padding,
                  xpd     = NA
                  )
    } else {
      rasterImage(raster_array, xleft = par("usr")[[1]], xright = par("usr")[[2]],
                  ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                  )
    }
    return(invisible(NULL))
  } else {
    return(invisible(raster_array))
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
                         draw_border    = FALSE,
                         use_lwd        = 1,
                         use_lend       = "butt"
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
      line_estimates_vec <- stats::approxfun(eval_points_vec, estimates_vec)(lines_vec) # get the violin bounds for the requested quantiles
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
          border = if (draw_border) border_color else NA,
          lwd = use_lwd * par("lwd"),
          xpd = NA
          )

  if (draw_lines) {
    line_heights_vec <- line_estimates_vec * h_scale - GetHalfLineWidth()
    for (i in seq_along(show_quantiles)) {
      lines(x    = c(at - line_heights_vec[[i]], at + line_heights_vec[[i]]),
            y    = rep(lines_vec[[i]], 2),
            lty  = quantiles_lty[[i]],
            col  = line_color,
            lend = use_lend,
            lwd  = use_lwd * par("lwd"),
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
                                lower_bound_enforced, upper_bound_enforced,
                                show_axis_truncation = TRUE
                                ) {
  plain_axis_labels <- format(tick_positions, trim = TRUE)
  plain_no_minus <- sub("-", "", plain_axis_labels, fixed = TRUE)
  axis_labels <- lapply(seq_along(tick_positions), function(x) {
    if (tick_positions[[x]] < 0) {
      as.expression(bquote(""["" - .(plain_no_minus[[x]])]))
    } else {
      as.expression(bquote(""[.(plain_axis_labels[[x]])]))
    }
  })
  axis_labels <- sapply(axis_labels, function(x) x)
  if (show_axis_truncation && lower_bound_enforced && (abs(lower_bound - tick_positions[[1]]) < 1e-15)) {
    if (tick_positions[[1]] < 0) {
      axis_labels[1] <- as.expression(bquote(""["" <= "" - "" * .(plain_no_minus[[1]])]))
    } else {
      axis_labels[1] <- as.expression(bquote(""["" <= scriptscriptstyle(" ") * .(plain_axis_labels[[1]])]))
    }
  }
  if (show_axis_truncation && upper_bound_enforced) {
    num_ticks <- length(tick_positions)
    if (abs(upper_bound - tick_positions[[num_ticks]]) < 1e-15) {
      if (tick_positions[[num_ticks]] < 0) {
        axis_labels[num_ticks] <- as.expression(bquote(""["" >= "" - "" * .(plain_no_minus[[num_ticks]])]))
      } else {
        axis_labels[num_ticks] <- as.expression(bquote(""["" >= scriptscriptstyle(" ") * .(plain_axis_labels[[num_ticks]])]))
      }
    }
  }
  return(axis_labels)
}




BeeViolinPlot <- function(input_list,
                          groups_vec      = NULL,
                          gap_ratio       = 1.25,
                          side_gap        = 0.5,
                          right_gap       = side_gap,
                          brewer_pals     = rep_len(c("Blues", "Reds", "Purples",
                                                      "Greens", "Oranges", "Greys"
                                                      ), length.out = length(input_list)
                                                    ),
                          point_colors    = NULL,
                          outline_colors  = NULL,
                          violin_colors   = NULL,
                          line_colors     = NULL,
                          border_colors   = NULL,
                          quantiles_lty   = c("21", "52", "21"),
                          point_cex       = 0.4,
                          use_spacing     = 0.8,
                          custom_spacing  = NULL,
                          swarm_method    = "compactswarm",
                          outline_points  = FALSE,
                          sina_plot       = FALSE,
                          sina_wex_factor = 1,
                          sina_jitter     = 0.01,
                          draw_violins    = (!(sina_plot)) || (length(input_list[[1]]) < 100),
                          sina_pt_colors  = draw_violins,
                          mini_box        = FALSE,
                          x_limits        = NULL,
                          y_limits        = NULL,
                          lower_bound     = NULL,
                          upper_bound     = NULL,
                          show_truncation = TRUE,
                          draw_grid       = FALSE,
                          indicate_zero   = TRUE,
                          zero_lty        = "dotted",
                          zero_lwd        = 1,
                          zero_color      = "gray80",
                          grid_color      = "gray91",
                          grid_lty        = "solid",
                          grid_lwd        = 1,
                          draw_groups_n   = TRUE,
                          draw_points     = TRUE,
                          use_swarm       = TRUE,
                          points_alpha    = if (sina_plot) {if (draw_violins) 0.8 else 0.4} else if (use_swarm) 1 else 0.2,
                          outlines_alpha  = points_alpha,
                          cloud_sd        = 0.04,
                          embed_PNG       = FALSE,
                          png_res         = 900,
                          show_x_axis     = TRUE,
                          show_y_axis     = TRUE,
                          text_vec        = NULL,
                          wex             = 0.8,
                          x_positions     = NULL,
                          num_ticks       = 5,
                          GridFunction    = NULL,
                          ...
                          ) {

  num_groups <- length(input_list)
  if (length(brewer_pals) == 1) {
    brewer_pals <- rep(brewer_pals, num_groups)
  }
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
  if (is.null(outline_colors)) {
    outline_colors <- violin_colors
  } else if (length(outline_colors) == 1) {
    outline_colors <- rep(outline_colors, num_groups)
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
  numeric_list <- split(numeric_vec, group_indices)
  are_finite_list <- split(are_finite, group_indices)

  if (is.null(x_positions)) {
    if (is.null(groups_vec)) {
      x_positions <- seq_len(num_groups)
    } else {
      x_positions <- RepositionByGroups(groups_vec, gap_ratio = gap_ratio)
    }
  }

  if (embed_PNG) {
    current_device <- StartEmbedPNG(figures_dir)
  }

  SetUpBoxPlot(num_groups,
               data_range    = range(numeric_vec),
               group_limits  = x_limits,
               use_y_limits  = y_limits,
               side_gap      = side_gap,
               right_gap     = right_gap,
               draw_axis     = FALSE,
               draw_box      = (!(embed_PNG)) && show_x_axis,
               draw_grid     = draw_grid,
               indicate_zero = indicate_zero,
               zero_lty      = zero_lty,
               zero_lwd      = zero_lwd,
               zero_color    = zero_color,
               grid_color    = grid_color,
               grid_lty      = grid_lty,
               grid_lwd      = grid_lwd,
               num_ticks     = num_ticks,
               GridFunction  = GridFunction
               )

  if (draw_violins) {
    for (i in seq_along(numeric_list)) {
      SingleViolin(numeric_list[[i]][are_finite_list[[i]]],
                   at            = x_positions[[i]],
                   violin_color  = violin_colors[[i]],
                   line_color    = line_colors[[i]],
                   border_color  = border_colors[[i]],
                   quantiles_lty = quantiles_lty,
                   wex           = wex,
                   ...
                   )
    }
  }

  if (draw_points) {
    set.seed(1)
    if (sina_plot) {
      sina_df <- sinaplot::sinaplot(numeric_list,
                                    col       = adjustcolor(if (sina_pt_colors) point_colors else violin_colors, alpha.f = points_alpha),
                                    plot      = FALSE,
                                    scale     = FALSE,
                                    maxwidth  = wex * sina_wex_factor,
                                    bin_limit = 1
                                    )
      x_vec <- rep(x_positions, lengths(numeric_list))
      if (sina_jitter > 0) {
        x_vec <- x_vec + rnorm(n = length(x_vec), mean = 0, sd = sina_jitter)
      }

      if (outline_points) {
        point_col_vec <- rep(adjustcolor(outline_colors, alpha.f = outlines_alpha), lengths(numeric_list))
      } else {
        point_col_vec <- sina_df[, "col"]
      }
      points(x_vec + (sina_df[, "scaled"] - sina_df[, "x"]),
             sina_df[, "y"],
             pch = if (outline_points) 21 else 16,
             col = point_col_vec,
             bg  = if (outline_points) sina_df[, "col"] else NA,
             cex = point_cex,
             lwd = par("lwd") * point_cex,
             xpd = NA
             )
    } else if (use_swarm) {
      for (i in seq_along(numeric_list)) {
        if (names(input_list)[[i]] %in% names(custom_spacing)) {
          this_spacing <- custom_spacing[[names(input_list)[[i]]]]
        } else {
          this_spacing <- use_spacing
        }
        beeswarm_df <- beeswarm::beeswarm(numeric_list[[i]],
                                          at       = x_positions[[i]],
                                          method   = swarm_method,
                                          priority = "random",
                                          spacing  = this_spacing,
                                          cex      = point_cex,
                                          col      = adjustcolor(point_colors[[i]], alpha.f = points_alpha),
                                          do.plot  = FALSE
                                          )

        if (outline_points) {
          point_col_vec <- rep(adjustcolor(outline_colors[[i]], alpha.f = outlines_alpha), length(numeric_list[[i]]))
        } else {
          point_col_vec <- beeswarm_df[, "col"]
        }

        points(beeswarm_df[, "x"],
               beeswarm_df[, "y"],
               pch = if (outline_points) 21 else 16,
               col = point_col_vec,
               bg  = if (outline_points) beeswarm_df[, "col"] else NA,
               cex = point_cex,
               lwd = par("lwd") * point_cex,
               xpd = NA
               )
        if (!(is.null(text_vec))) {
          this_text_vec <- text_vec[rep(seq_along(numeric_list), lengths(numeric_list)) == i]
          stopifnot(length(this_text_vec) == nrow(beeswarm_df))
          old_lheight <- par("lheight" = 1.3)
          text(x      = beeswarm_df[, "x"],
               y      = beeswarm_df[, "y"],
               labels = this_text_vec,
               col    = Darken(beeswarm_df[, "col"], factor = 1.2),
               cex    = 0.05,
               font   = 3,
               xpd    = NA
               )
          par(old_lheight)
        }
      }

    } else {
      x_vec <- rep(x_positions, lengths(numeric_list))
      x_vec <- x_vec + rnorm(n = length(x_vec), mean = 0, sd = cloud_sd)

      assign("delete_numeric_list", numeric_list, envir = globalenv())
      assign("delete_foo", adjustcolor(point_colors, alpha.f = points_alpha), envir = globalenv())

      fills_vec <- rep(adjustcolor(point_colors, alpha.f = points_alpha), lengths(numeric_list))
      if (outline_points) {
        point_col_vec <- rep(adjustcolor(outline_colors, alpha.f = outlines_alpha), lengths(numeric_list))
      } else {
        point_col_vec <- fills_vec
      }
      points(x   = x_vec,
             y   = unlist(numeric_list),
             col = point_col_vec,
             bg  = if (outline_points) fills_vec else NA,
             cex = point_cex,
             pch = if (outline_points) 21 else 16,
             lwd = par("lwd") * point_cex,
             xpd = NA
             )
    }
  }


  if (embed_PNG) {
    StopEmbedPNG(current_device, figures_dir)
    if (show_x_axis) {
      box(bty = "l")
    } else {
      segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]],
               y0 = par("usr")[[3]], xpd = NA
               )
    }
  }

  if (sina_plot) {
    if (mini_box) {
      quantile_mat <- t(sapply(input_list, quantile, probs = c(0.05, 0.95), na.rm = TRUE))
      segments(x0   = x_positions,
               y0   = quantile_mat[, 1],
               y1   = quantile_mat[, 2],
               col  = line_colors,
               lend = "butt",
               lwd  = par("lwd"),
               xpd  = NA
               )
      point_width <- (par("cin")[[2]] / par("pin")[[1]]) * (par("usr")[[2]] - par("usr")[[1]]) * 0.5 * par("cex") * 0.375
      boxplot(input_list,
              at         = x_positions,
              boxwex     = point_width * 2.5,
              outline    = FALSE,
              names      = rep.int("", length(x_positions)),
              whisklty   = "blank",
              staplewex  = 0,
              whisklwd   = 0,
              staplelty  = 0,
              medlwd     = NA,
              col        = line_colors,
              border     = NA,
              add        = TRUE,
              axes       = FALSE
              )
      points(x    = x_positions,
             y    = vapply(input_list, median, numeric(1)),
             col  = line_colors,
             bg   = if (identical(point_colors, line_colors)) violin_colors else point_colors,
             pch  = 23,
             cex  = 0.5,
             lwd  = par("lwd") * 0.5,
             xpd  = NA
             )
    } else if (!(draw_violins)) {
      for (i in seq_along(numeric_list)) {
        SingleViolin(numeric_list[[i]][are_finite_list[[i]]],
                     at            = x_positions[[i]],
                     violin_color  = NA,
                     line_color    = line_colors[[i]],
                     border_color  = NA,
                     quantiles_lty = quantiles_lty,
                     wex           = wex,
                     ...
                     )
      }
    }
  }

  ## Annotate plot
  axis_ticks <- pretty(range(axTicks(2)), n = num_ticks)
  axis_labels <- CurtailedAxisLabels(axis_ticks,
                                     lower_bound          = lower_bound,
                                     upper_bound          = upper_bound,
                                     lower_bound_enforced = lower_bound_enforced,
                                     upper_bound_enforced = upper_bound_enforced,
                                     show_axis_truncation = show_truncation
                                     )
  if (show_y_axis) {
    axis(2,
         at       = axis_ticks,
         labels   = axis_labels,
         mgp      = c(3, 0.5, 0),
         tcl      = -0.35,
         las      = 1,
         lwd      = par("lwd"),
         cex.axis = 1 / 0.7
         )
  } else if (show_x_axis) {
    segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]],
             y0 = par("usr")[[3]], xpd = NA
             )
  }
  if (draw_groups_n) {
    mtext(lengths(numeric_list), at = x_positions, side = 1, line = -0.45,
          cex = par("cex") * 0.3, col = "gray70", padj = 0
          )
  }

  return(invisible(x_positions))
}



