### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")






# Define functions --------------------------------------------------------

MakeCorrTitle <- function(numeric_vec_1, numeric_vec_2) {
  corr_results <- cor.test(numeric_vec_1, numeric_vec_2)
  title_expression <- as.expression(bquote(
    bold("Pearson's " * bolditalic("r") * " = " *
           .(as.character(signif(corr_results[["estimate"]], digits = 2)))) *
      " ("  * italic("p") * " = " *
      .(format(signif(corr_results[["p.value"]], digits = 1), scientific = 8)) * ")"
  ))
  text(x      = par("usr")[[1]] + (par("usr")[[2]] - par("usr")[[1]]) / 2,
       y      = par("usr")[[4]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.042),
       labels = title_expression,
       cex    = 0.8,
       xpd    = NA
       )
  return(invisible(corr_results))
}




PlotBySharedSubsequence <- function(summary_df, show_column) {

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
       xlim = c(3, 20),
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       ann  = FALSE,
       type = "n",
       axes = FALSE,
       mgp  = c(2.7, 1, 0)
       )

  x_mid <- par("usr")[[1]] + (par("usr")[[2]] - par("usr")[[1]]) / 2


  text(x       = x_mid,
       y       = par("usr")[[4]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.187),
       labels  = titles_list[[show_column]],
       xpd     = NA,
       cex     = 0.9,
       font    = 1
       )


  tick_locations <- axTicks(2)
  if (is_percentage) {
    tick_labels <- paste0(tick_locations * 100, "%")
  } else {
    tick_labels <- TRUE
  }

  if (is_percentage) {
    abline(h = seq(0.1, 0.9, 0.2), col = "gray98", lwd = 0.75)
  }
  abline(h = tick_locations, col = "gray94", lwd = 0.75)

  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.35,
       lwd      = 0.75,
       cex.axis = 0.9
       )

  corr_results <- MakeCorrTitle(numeric_vec, shared_bp_vec)

  box(bty = "l", lwd = 0.75)

  shared_bp_fac <- factor(shared_bp_vec, levels = 1:20)

  boxplot_fac <- shared_bp_fac
  boxplot_fac[table(boxplot_fac)[boxplot_fac] == 1] <- NA # I don't want the mean to be shown if there is only one point

  boxplot(numeric_vec ~ boxplot_fac,
          ylim      = c(0.5, 1),
          boxwex    = 0.7,
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
                          spacing  = 0.6,
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
       y      = par("usr")[[3]] - ((par("usr")[[4]] - par("usr")[[3]]) * 0.04),
       labels = present_lengths,
       font   = 1,
       cex    = 0.8,
       xpd    = NA
       )

  text(x      = x_mid,
       y      = par("usr")[[3]] - ((par("usr")[[4]] - par("usr")[[3]]) * 0.12),
       labels = "Length of shared subsequence",
       font   = 1,
       cex    = 0.9,
       xpd    = NA
       )

  return(invisible(corr_results))
}


