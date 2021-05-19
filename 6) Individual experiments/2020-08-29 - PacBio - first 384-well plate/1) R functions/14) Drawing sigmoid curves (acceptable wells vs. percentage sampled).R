### 9th February 2021 ###




# Import packages and source code -----------------------------------------

library("beeswarm")
library("RColorBrewer")






# Define constants --------------------------------------------------------

ccs3_title <- expression(plain({"" >= "3 consensus reads "} *
                                 "and " >= "99% accuracy"
                               ))

ccs5_title <- expression(plain({"" >= "5 consensus reads "} *
                                 "and " >= "99.9% accuracy"
                               ))

ccs7_title <- expression(plain({"" >= "7 consensus reads "} *
                                 "and " >= "99.99% accuracy"
                               ))



# Define functions --------------------------------------------------------

DrawAllSigmoidCurves <- function() {

  summary_df_names <- c("original_summary_df", "filtered_summary_df", "filtered_gRNAs_df")
  use_colors <- c("Blues", "Reds", "Purples", "Greys")
  num_reads <- c(10, 20, 30, 50)
  legend_labels <- lapply(as.character(num_reads),
                          function(x) bquote(plain("" >= .(x) ~ "reads"))
                          )
  legend_labels <- sapply(legend_labels, as.expression)

  use_height <- 7
  use_width <- 7

  for (smrtlink_version in c(7, 9)) {

    version_folder <- paste0("SmrtLink ", smrtlink_version)
    version_path <- file.path(plots_output_directory, version_folder)

    use_sample_list <- get(paste0("sl", smrtlink_version, "_subsampled_list"))
    for (df_i in seq_along(summary_df_names)) {
      file_name <- paste0("Sigmoid curves - ",
                          "SmrtLink ", smrtlink_version, " - ",
                          c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[df_i]]
                          )

      pdf(file = file.path(version_path, "Sigmoid curves", paste0(file_name, ".pdf")),
          height = use_height,
          width  = use_width
          )
      old_mar <- par(mar = c(5, 5.5, 5, 2.5))

      for (ccs_number in c(3, 5, 7)) {
        for (i in seq_along(use_colors)) {
          DrawSigmoidCurve(use_sample_list,
                           summary_df_name  = summary_df_names[[df_i]],
                           use_ccs          = ccs_number,
                           add_curve        = if (i == 1) FALSE else TRUE,
                           threshold_number = num_reads[[i]],
                           use_color        = use_colors[[i]]
                           )
          legend(x        = 65,
                 y        = 20,
                 lwd      = 2,
                 legend   = legend_labels,
                 col      = vapply(use_colors, function(x) brewer.pal(9, x)[[5]], ""),
                 text.col = vapply(use_colors, function(x) brewer.pal(9, x)[[7]], ""),
                 bty      = "n",
                 seg.len  = 1,
                 adj      = c(0.12, 0.5),
                 inset    = 0.02,
                 xpd      = NA
                 )
        }
      }
      par(old_mar)
      dev.off()
    }
  }
  return(invisible(NULL))
}






DrawGridLines <- function(line_positions,
                          use_color,
                          horizontal  = TRUE,
                          grid_lwd    = 0.5,
                          lower_limit = NULL,
                          upper_limit = 100
                          ) {
  if (horizontal) {
    segments(x0   = if (is.null(lower_limit)) par("usr")[[1]] else lower_limit,
             x1   = upper_limit,
             y0   = line_positions,
             y1   = line_positions,
             col  = use_color,
             lwd  = grid_lwd,
             lend = "butt",
             xpd  = NA
             )
  } else {
    segments(y0   = if (is.null(lower_limit)) par("usr")[[3]] else lower_limit,
             y1   = upper_limit,
             x0   = line_positions,
             x1   = line_positions,
             col  = use_color,
             lwd  = grid_lwd,
             xpd  = NA,
             lend = "butt"
             )
  }
  return(invisible(NULL))
}




DrawGrid <- function(line_positions,
                     main_color    = "gray86",
                     between_color = "gray94",
                     horizontal    = TRUE,
                     lower_limit   = NULL,
                     upper_limit   = 100
                     ) {
  are_main <- (line_positions %% 20) == 0
  DrawGridLines(line_positions[are_main],
                main_color,
                horizontal = horizontal,
                lower_limit = lower_limit,
                upper_limit = upper_limit
                )
  DrawGridLines(line_positions[!(are_main)],
                between_color,
                horizontal = horizontal,
                lower_limit = lower_limit,
                upper_limit = upper_limit
                )
  return(invisible(NULL))
}




SetUpPlot <- function(main_title = "") {

  plot(1,
       type = "n",
       xlim = c(-5, 105),
       ylim = c(-2, 102),
       xaxs = "i",
       yaxs = "i",
       tcl  = -0.35,
       las  = 1,
       mgp  = c(2.7, 0.6, 0),
       ylab = "Wells with sufficient numbers of reads",
       xlab = "Reads sampled",
       main = main_title,
       axes = FALSE
       )
  grid_lwd <- 0.5
  main_color <- "gray86"
  between_color <- "gray94"

  DrawGrid(c(0, seq(20, 100, by = 10)),
           horizontal = TRUE
           )
  DrawGrid(seq(0, 20, by = 10),
           horizontal = TRUE,
           upper_limit = 60
           )
  DrawGrid(c(seq(0, 60, by = 10), 100),
           horizontal = FALSE
           )
  DrawGrid(seq(60, 100, by = 10),
           horizontal = FALSE,
           lower_limit = 20
           )

  y_axis_pos <- axTicks(2)
  x_axis_pos <- axTicks(1)
  axis(2,
       at     = y_axis_pos,
       labels = paste0(y_axis_pos, "%"),
       tcl    = -0.35,
       las    = 1,
       mgp    = c(2.5, 0.45, 0)
       )
  axis(1,
       at     = x_axis_pos,
       labels = paste0(x_axis_pos, "%"),
       tcl    = -0.35,
       las    = 1,
       mgp    = c(2.5, 0.45, 0)
       )
  segments(x0  = par("usr")[[1]],
           x1  = 100,
           y0  = par("usr")[[3]],
           y1  = par("usr")[[3]],
           xpd = NA
           )
  segments(x0  = par("usr")[[1]],
           x1  = par("usr")[[1]],
           y0  = par("usr")[[3]],
           y1  = 100,
           xpd = NA
           )
  return(invisible(NULL))
}




MakeTransparent <- function(color_string, fraction = 0.5) {
  rgb_val <- col2rgb(color_string)
  new_col <- rgb(rgb_val[1], rgb_val[2], rgb_val[3],
                 max = 255,
                 alpha = (1 - fraction) * 255
                 )
  return(new_col)
}



DrawSigmoidCurve <- function(sampled_list,
                             use_ccs = 7L,
                             threshold_number = 20L,
                             add_curve = FALSE,
                             use_color = "Blues",
                             summary_df_name = "filtered_summary_df"
                             ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  if ("Empty_well" %in% names(sg_sequences_df)) {
    are_empty <- sg_sequences_df[["Empty_well"]]
  } else {
    are_empty <- rep(FALSE, nrow(sg_sequences_df))
  }

  ccs_name <- paste0("ccs", use_ccs)

  use_title <- get(paste0(ccs_name, "_title"))

  num_meet_target <- lapply(seq_along(sampled_list), function(x) {
    lapply(seq_along(sampled_list[[x]]), function(y) {
      sum(sampled_list[[x]][[y]][[ccs_name]][[summary_df_name]][["Count_total"]][!(are_empty)] > threshold_number)
    })
  })
  num_meet_target <- lapply(sampled_list, function(x) {
    lapply(x, function(y) {
      sum(y[[ccs_name]][[summary_df_name]][["Count_total"]][!(are_empty)] >= threshold_number)
    })
  })

  rep_mat_list <- lapply(num_meet_target, function(x) {
    cbind("Rep_number"  = as.integer(sub("rep", "", names(x), fixed = TRUE)),
          "Num_passing" = unlist(x, use.names = FALSE)
          )
  })
  sample_mat_list <- lapply(seq_along(rep_mat_list), function(x) {
    cbind("Fraction_sampled" = as.numeric(sub("% sampled", "", names(rep_mat_list)[[x]])) / 100,
          rep_mat_list[[x]]
          )
   })

  sample_mat <- do.call(rbind, sample_mat_list)
  sample_mat <- cbind(sample_mat,
                      "Fraction_passing" = sample_mat[, "Num_passing"] / sum(!(are_empty))
                      )

  if (!(add_curve)) {
    SetUpPlot(main_title = use_title)
  }

  sampled_list <- split(sample_mat[, "Fraction_passing"] * 100,
                        sample_mat[, "Fraction_sampled"]
                        )

  x_vec <- sample_mat[, "Fraction_sampled"] * 100
  y_vec <- sample_mat[, "Fraction_passing"] * 100

  fit <- try(nls(y ~ SSlogis(x, Asym, xmid, scal),
                 data = data.frame(x = x_vec, y = y_vec)
                 ))
  if (class(fit) == "try-error") {
    message("A sigmoid curve could not be fitted successfully. Some jitter was added.")
    set.seed(1)
    fit <- try(nls(y ~ SSlogis(x, Asym, xmid, scal),
                   data = data.frame(x = jitter(x_vec), y = jitter(y_vec))
                   ))
  }

  use_seq <- seq(0, 100, length.out = 100)
  lines(use_seq,
        predict(fit, newdata = data.frame(x = use_seq)),
        lwd = 3,
        col = brewer.pal(9, use_color)[[5]]
        )

  x_positions <- sort(unique(sample_mat[, "Fraction_sampled"]) * 100)
  group_means <- vapply(sampled_list, mean, numeric(1))
  mean_length <- 3
  segments(x0   = x_positions - (mean_length / 2),
           x1   = x_positions + (mean_length / 2),
           y0   = group_means,
           y1   = group_means,
           xpd  = NA,
           lwd  = 2,
           col  = MakeTransparent(brewer.pal(9, use_color)[[8]], 0.3),
           lend = "butt"
           )
  set.seed(1)
  beeswarm(sampled_list,
           at       = x_positions,
           add      = TRUE,
           pch      = 16,
           col      = MakeTransparent(brewer.pal(9, use_color)[[9]]),
           cex      = 0.4,
           spacing  = 0.8,
           priority = "random"
           )
  return(invisible(NULL))
}



