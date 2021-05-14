### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")





# Define constants --------------------------------------------------------

pdf_height <- 5.4
pdf_width <- 7.2
use_mar <- c(4.5, 4, 5.7, 2.5)

exclude_vars <- "Longest_subsequence"
exclude_vars <- c(
  exclude_vars,
  "Num_contam_wells",
  "Num_low_barcode_scores",
  "Num_low_quality_scores",
  "Num_low_bc_or_qual",
  "Mean_read_quality"
)

use_vars <-  setdiff(names(titles_list), exclude_vars)




# Define functions --------------------------------------------------------

MakeCorrTitle <- function(numeric_vec_1, numeric_vec_2, use_line = 0.75, bold_corr = TRUE, text_cex = 0.9) {
  corr_results <- cor.test(numeric_vec_1, numeric_vec_2)
  corr_value <- as.character(signif(corr_results[["estimate"]], digits = 2))
  p_value <- format(signif(corr_results[["p.value"]], digits = 1), scientific = 8)
  if (bold_corr) {
    title_expression <- as.expression(bquote(
      bold("Pearson's " * bolditalic("r") * " = " * .(corr_value)) *
           " ("  * italic("p") * " = " * .(p_value) * ")"
      ))
  } else {
    title_expression <- as.expression(bquote(
      "Pearson's " * italic("r") * " = " * .(corr_value) *
           " ("  * italic("p") * " = " * .(p_value) * ")"
      ))
  }

  text(x      = par("usr")[[1]] + (par("usr")[[2]] - par("usr")[[1]]) / 2,
       y      = par("usr")[[4]] + diff(grconvertY(c(0, use_line), from = "lines", to = "user")),
       labels = title_expression,
       cex    = text_cex,
       xpd    = NA
       )
  return(invisible(corr_results))
}




PlotBySharedSubsequence <- function(summary_df,
                                    show_column,
                                    use_spacing     = 0.6,
                                    use_boxwex      = 0.7,
                                    use_title       = NULL,
                                    y_axis_label    = "",
                                    text_cex        = 0.9,
                                    grid_light_gray = "gray98",
                                    grid_dark_gray  = "gray94",
                                    grid_lwd        = 0.75,
                                    corr_line       = 0.75,
                                    bold_corr       = TRUE,
                                    x_axis_label    = "Length of shared subsequence"
                                    ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  if (show_column == "Count_mean_sg1to4") {
    count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
    summary_df[["Count_mean_sg1to4"]] <- rowMeans(as.matrix(summary_df[, count_columns]))
  }

  if ("Empty_well" %in% names(sg_sequences_df)) {
    are_to_include <- !(sg_sequences_df[["Empty_well"]])
  } else {
    are_to_include <- rep(TRUE, nrow(sg_sequences_df))
  }

  if (is.null(use_title)) {
    use_title <- titles_list[[show_column]]
  }

  shared_bp_vec <- sg_sequences_df[["Longest_subsequence"]][are_to_include]
  numeric_vec <- summary_df[[show_column]][are_to_include]

  light_color <- brewer.pal(9, "Blues")[[2]]
  dark_color <- brewer.pal(9, "Blues")[[7]]


  is_percentage <- grepl("^(Count|Num)_", show_column)

  if (is_percentage) {
    numeric_vec <- numeric_vec / summary_df[["Count_total"]][are_to_include]
    y_limits <- c(0, 1)
  } else {
    numeric_vec <- numeric_vec
    y_limits <- c(0, max(numeric_vec))
  }

  plot(1,
       xlim = c(3.3, 19.7),
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       ann  = FALSE,
       axes = FALSE
       )

  x_mid <- par("usr")[[1]] + (par("usr")[[2]] - par("usr")[[1]]) / 2

  text(x       = x_mid,
       y       = par("usr")[[4]] + diff(grconvertY(c(0, 3.15), from = "lines", to = "user")),
       labels  = use_title,
       xpd     = NA,
       cex     = text_cex,
       font    = 1
       )

  tick_locations <- axTicks(2)
  if (is_percentage) {
    tick_labels <- paste0(tick_locations * 100, "%")
  } else {
    tick_labels <- TRUE
  }

  grid_locations <- seq(0, 1, 0.1)
  grid_colors <- ifelse(grid_locations %in% tick_locations, grid_dark_gray, grid_light_gray)
  segments(x0   = par("usr")[[1]],
           x1   = par("usr")[[2]],
           y0   = grid_locations,
           col  = grid_colors,
           lend = "butt",
           xpd  = NA
           )

  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.38, 0),
       tcl      = -0.3,
       lwd      = par("lwd"),
       cex.axis = text_cex
       )

  corr_results <- MakeCorrTitle(numeric_vec, shared_bp_vec,
                                use_line = corr_line, bold_corr = bold_corr,
                                text_cex = text_cex
                                )

  box(bty = "l")

  shared_bp_fac <- factor(shared_bp_vec, levels = 1:20)

  boxplot_fac <- shared_bp_fac
  boxplot_fac[table(boxplot_fac)[boxplot_fac] == 1] <- NA # I don't want the mean to be shown if there is only one point

  boxplot(numeric_vec ~ boxplot_fac,
          ylim      = c(0.5, 1),
          boxwex    = use_boxwex,
          outline   = FALSE,
          names     = rep("", 20),
          whisklty  = "blank",
          staplewex = 0,
          axes      = FALSE,
          whisklwd  = 0,
          staplelty = 0,
          col       = light_color,
          boxlwd    = 0.75,
          medlwd    = par("lwd") * 2,
          add       = TRUE
          )

  set.seed(1)
  point_cex <- 0.4
  beeswarm_df <- beeswarm(numeric_vec ~ shared_bp_fac,
                          priority = "random",
                          spacing  = use_spacing,
                          cex      = point_cex,
                          do.plot  = FALSE
                          )

  points(beeswarm_df[["x"]],
         beeswarm_df[["y"]],
         pch = 16,
         cex = point_cex,
         col = dark_color,
         xpd = NA
         )

  present_lengths <- seq(from = min(shared_bp_vec), to = max(shared_bp_vec))

  text(x      = present_lengths,
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
       labels = present_lengths,
       font   = 1,
       cex    = text_cex,
       xpd    = NA
       )
  text(x      = x_mid,
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 2.3), from = "lines", to = "user")),
       labels = x_axis_label,
       font   = 1,
       cex    = text_cex,
       xpd    = NA
       )
  text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 2.53), from = "lines", to = "user")),
       y      = grconvertY(0.5, from = "npc", to = "user"),
       labels = VerticalAdjust(y_axis_label),
       srt    = 90,
       xpd    = NA,
       adj    = c(0.5, 0)
       )
  return(invisible(corr_results))
}



