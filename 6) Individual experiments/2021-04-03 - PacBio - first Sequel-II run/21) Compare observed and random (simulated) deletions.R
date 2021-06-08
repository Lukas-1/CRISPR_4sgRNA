### 7th June 2021 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output", "Figures", "Compare simulated and observed deletions")





# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "20) Simulate random deletions (to compare with observed deletions).RData"))





# Define functions --------------------------------------------------------

MakeWhiskers <- function(numeric_vec, two_SE_whiskers = TRUE, number_of_SEs = if (two_SE_whiskers) 2 else 1) {
  my_mean <- mean(numeric_vec, na.rm = TRUE)
  my_sd <- sd(numeric_vec, na.rm = TRUE)
  group_n <- sum(!(is.na(numeric_vec)))
  group_se_mean <- my_sd / sqrt(group_n)
  lower_bound <- my_mean - (number_of_SEs * group_se_mean)
  upper_bound <- my_mean + (number_of_SEs * group_se_mean)
  results_vec <- c("Mean" = my_mean, "Lower_bound" = lower_bound, "Upper_bound" = upper_bound)
  return(results_vec)
}



SimVsObsDeletions <- function(use_simul_mat, use_del_df) {

  stopifnot("labels_list" %in% ls(envir = globalenv()))

  set.seed(1) # For reproducible jitter

  ## Prepare the relevant metrics
  new_column_order <- order(match(colnames(use_simul_mat), names(labels_list)))
  use_simul_mat <- use_simul_mat[, new_column_order]
  simul_means <- colMeans(use_simul_mat)
  real_means <- colSums(use_del_df[, colnames(use_simul_mat)]) / nrow(use_del_df)


  ## Determine group positions
  num_groups <- ncol(use_simul_mat)
  group_positions <- seq_len(num_groups)
  group_limits <- c((min(group_positions) - 0.5) - (num_groups * 0.04),
                    max(group_positions) + 0.5 + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  use_numeric_limits <- c(0, 1)
  numeric_axis_pos <- pretty(use_numeric_limits)
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

  segments(x0   = par("usr")[[1]],
           x1   = par("usr")[[2]],
           y0   = seq(0.05, 0.95, by = 0.1),
           col  = "gray95",
           lend = "butt",
           xpd  = NA
           )

  segments(x0   = par("usr")[[1]],
           x1   = par("usr")[[2]],
           y0   = seq(0, 1, by = 0.1),
           col  = "gray85",
           lend = "butt",
           xpd  = NA
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

  mtext(text = "Reads containing deletions of 20 bp or more",
        side = 2,
        line = 3.2
        )


  ## Plot the bars
  bar_width <- 0.3
  final_bar_width <- bar_width * ((max(group_positions) - min(group_positions)) / (num_groups - 1))
  x_gap <- 0.36

  simulated_color <- brewer.pal(9, "Blues")[[4]]
  real_color <- brewer.pal(9, "Purples")[[7]]

  rect(xleft   = group_positions - (x_gap / 2) - (final_bar_width / 2),
       xright  = group_positions - (x_gap / 2) + (final_bar_width / 2),
       ybottom = 0,
       ytop    = simul_means,
       col     = simulated_color,
       border  = NA,
       xpd     = NA
       )

  rect(xleft   = group_positions + (x_gap / 2) - (final_bar_width / 2),
       xright  = group_positions + (x_gap / 2) + (final_bar_width / 2),
       ybottom = 0,
       ytop    = real_means,
       col     = real_color,
       border  = NA,
       xpd     = NA
       )



  ## Plot the points (for the simulations)
  points_alpha <- 0.2
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
  for (i in seq_len(num_groups)) {
    jitter_vec <- rnorm(n = nrow(use_simul_mat), mean = 0, sd = 0.025)
    points(x   = (group_positions[[i]] - (x_gap / 2)) + jitter_vec,
           y   = use_simul_mat[, i],
           cex = 0.5,
           pch = 16,
           col = paste0(brewer.pal(9, "Blues")[[8]], alpha_hex),
           xpd = NA
           )
  }

  box(bty = "l")


  ## Draw the group labels
  sim_obs_y_pos <- par("usr")[[3]] - diff(grconvertY(c(0, 1), from = "lines", to = "user"))
  sim_obs_cex <- 0.9

  text(x       = group_positions - (x_gap / 2),
       y       = sim_obs_y_pos,
       labels  = "Sim",
       cex     = sim_obs_cex,
       xpd     = NA
       )
  text(x       = group_positions + (x_gap / 2),
       y       = sim_obs_y_pos,
       labels  = "Obs",
       cex     = sim_obs_cex,
       xpd     = NA
       )

  text(x       = group_positions,
       y       = par("usr")[[3]] - diff(grconvertY(c(0, 2.75), from = "lines", to = "user")),
       labels  = unlist(labels_list[colnames(use_simul_mat)]),
       xpd     = NA
       )

  title("Simulated (random) vs. observed deletions", cex.main = 1.1)

  return(invisible(NULL))
}




# Define labels -----------------------------------------------------------

labels_list <- list(
  "Span_sgRNAs"    = "Span sgRNAs",
  "Span_tracrRNAs" = "Span tracrRNAs",
  "Span_sg_cr"     = "Span sg+cr",
  "Span_promoters" = "Span promoters"
)





# Explore the metrics -----------------------------------------------------

colMeans(one_read_simul_mat) * 100
colSums(one_read_del_df[, colnames(one_read_simul_mat)])  / nrow(one_read_del_df) * 100

colMeans(all_reads_simul_mat) * 100
colSums(all_reads_del_df[, colnames(one_read_simul_mat)]) / nrow(all_reads_del_df) * 100





# Draw plots --------------------------------------------------------------

SimVsObsDeletions(one_read_simul_mat, one_read_del_df)
SimVsObsDeletions(all_reads_simul_mat, all_reads_del_df)


SimVsObsDeletions(one_read_simul_mat[, c(1, 2, 4)], one_read_del_df)




pdf(file = file.path(file_output_directory, "Compare simulated and observed deletions - 4 categories.pdf"),
    width = 7.25, height = 5.5
    )

par(mar = c(5, 5, 4, 3))

SimVsObsDeletions(one_read_simul_mat, one_read_del_df)

dev.off()





pdf(file = file.path(file_output_directory, "Compare simulated and observed deletions - 3 categories.pdf"),
    width = 6.5, height = 5.5
    )

par(mar = c(5, 5, 4, 3))

SimVsObsDeletions(one_read_simul_mat[, c(1, 2, 4)], one_read_del_df)

dev.off()






