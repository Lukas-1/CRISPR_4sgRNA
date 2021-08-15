### 2nd August 2021 ###




# Import packages and source code -----------------------------------------

library("beeswarm")
library("RColorBrewer")





# Define functions --------------------------------------------------------

ComparePlates <- function(summary_df,
                          show_column,
                          use_cex          = 0.175,
                          beeswarm_spacing = 0.7,
                          beeswarm_corral  = "none",
                          side_space       = 0.2,
                          order_by_rank    = TRUE
                          ) {

  stopifnot("plates_df" %in% ls(envir = globalenv()))

  if (show_column == "Count_mean_sg1to4") {
    count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
    numeric_vec <- rowMeans(as.matrix(summary_df[, count_columns]))
  } else if (show_column == "Count_mean_pr_sg1to4") {
    count_columns <- paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4)
    numeric_vec <- rowMeans(as.matrix(summary_df[, count_columns]))
  } else {
    numeric_vec <- summary_df[[show_column]]
  }
  is_percentage <- grepl("^(Count|Num)_", show_column) &&
                   (!(show_column == "Count_total"))
  if (is_percentage) {
    numeric_vec <- numeric_vec / summary_df[["Count_total"]]
    y_limits <- c(0, 1)
  } else {
    y_max <-  max(numeric_vec, na.rm = TRUE)
    if (show_column == "Mean_read_quality") {
      y_min <- min(numeric_vec, na.rm = TRUE)
      y_span <- y_max - y_min
      y_range <- c(y_min - (y_span * 0.02), y_max + (y_span * 0.02))
      y_limits <- range(pretty(y_range))
    } else {
      y_limits <- c(0, y_max * 1.02)
    }
  }

  if (order_by_rank && ("Plate_rank" %in% names(plates_df))) {
    plates_order <- order(plates_df[, "Plate_rank"])
  } else {
    plates_order <- seq_len(nrow(plates_df))
  }

  groups_fac <- factor(summary_df[["Plate_number"]],
                       levels = plates_df[["Plate_number"]][plates_order]
                       )
  num_groups <- nlevels(groups_fac)
  stopifnot(num_groups == nrow(plates_df))

  group_limits <- c((1 - side_space) - (num_groups * 0.04),
                    (num_groups + side_space) + (num_groups * 0.04)
                    )
  plot(1,
       xlim = group_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       ann  = FALSE,
       axes = FALSE
       )
  title(titles_list[[show_column]], cex.main = 1)

  tick_locations <- axTicks(2)
  if (is_percentage) {
    tick_labels <- paste0(tick_locations * 100, "%")
  } else {
    tick_labels <- TRUE
  }
  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       mgp      = c(3, 0.45, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1,
       lwd      = par("lwd")
       )
  box(bty = "l")

  suppressWarnings(boxplot(numeric_vec ~ groups_fac, # Suppress the warning: no non-missing arguments to min / max
                           boxwex    = 0.7,
                           outline   = FALSE,
                           names     = rep("", nlevels(groups_fac)),
                           whisklty  = "blank",
                           staplewex = 0,
                           axes      = FALSE,
                           whisklwd  = 0,
                           staplelty = 0,
                           col       = brewer.pal(9, "Blues")[[2]],
                           boxlwd    = 0.75,
                           medlwd    = par("lwd") * 2,
                           add       = TRUE,
                           xpd       = NA
                           ))

  set.seed(1)
  beeswarm_df <- beeswarm(numeric_vec ~ groups_fac,
                          spacing  = beeswarm_spacing,
                          priority = "random",
                          corral   = beeswarm_corral,
                          cex      = use_cex,
                          do.plot  = FALSE
                          )
  points(beeswarm_df[["x"]],
         beeswarm_df[["y"]],
         pch = 16,
         cex = use_cex,
         col = brewer.pal(9, "Blues")[[8]],
         xpd = NA
         )

  plate_labels <- paste0("#", plates_df[["Plate_number"]], " ", plates_df[["Plate_name"]])

  text(x      = seq_len(num_groups),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.4), from = "lines", to = "user")),
       labels = plate_labels[plates_order],
       srt    = 45,
       adj    = c(1, 0.5),
       cex    = 0.6,
       xpd    = NA
       )
  return(invisible(NULL))
}





DrawAllPlateComparisons <- function(export_PNGs      = TRUE,
                                    use_cex          = 0.175,
                                    side_space       = 0.2,
                                    use_width        = 9.75,
                                    use_height       = 6.5,
                                    beeswarm_spacing = 0.7,
                                    beeswarm_corral  = "none",
                                    order_by_rank    = TRUE
                                    ) {

  stopifnot("titles_list" %in% ls(envir = globalenv()))

  if (export_PNGs) {
    file_formats <- c("png", "pdf")
  } else {
    file_formats <- "pdf"
  }

  ccs_numbers <- c(3, 5, 7)
  accuracy_percentages <- c(99, 99.9, 99.99)

  use_metrics <- setdiff(names(titles_list), "Longest_subsequence")
  use_metrics <- grep("^Binary_", use_metrics, value = TRUE, invert = TRUE)


  for (file_format in file_formats) {
    for (i in seq_along(ccs_numbers)) {
      use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))
      filter_stages <- c("original_summary_df", "filtered_summary_df")
      filter_labels <- c("i) unfiltered", "ii) filtered")
      if ("filtered_cross_plate_df" %in% names(use_df_list)) {
        filter_stages <- c(filter_stages, "filtered_cross_plate_df")
        filter_labels <- c(filter_labels, "iii) filtered cross-plate")
      }
      for (filter_stage in seq_along(filter_stages)) {
        df_name <- filter_stages[[filter_stage]] # "filtered_gRNAs_df"

        use_summary_df <- use_df_list[[df_name]]

        sel_name <- paste0("CCS", ccs_numbers[[i]],
                           " (", accuracy_percentages[[i]], ") - ",
                           c("i) unfiltered", "ii) filtered", "iii) filtered cross-plate")[[filter_stage]]
                           )
        message(paste0("Exporting ", file_format, " images into the folder: ", sel_name, "..."))
        if (file_format == "pdf") {
          pdf(file   = file.path(plots_output_directory, paste0(sel_name, ".pdf")),
              width  = use_width,
              height = use_height
              )
          for (use_metric in use_metrics) {
            message(paste0("   metric: '", use_metric, "'"))
            ComparePlates(use_summary_df,
                          use_metric,
                          use_cex          = use_cex,
                          beeswarm_spacing = beeswarm_spacing,
                          beeswarm_corral  = beeswarm_corral,
                          side_space       = side_space,
                          order_by_rank    = order_by_rank
                          )
          }
          dev.off()
        } else if (file_format == "png") {
          sub_folder_path <- file.path(PNGs_output_directory, sel_name)
          dir.create(sub_folder_path, showWarnings = FALSE)
          for (j in seq_along(use_metrics)) {
            message(paste0("   metric: '", use_metrics[[j]], "'"))
            file_name <- paste0(formatC(j, width = 2, flag = "0"), ") ",
                                use_metrics[[j]], ".png"
                                )
            png(file   = file.path(sub_folder_path, file_name),
                width  = use_width,
                height = use_height,
                units  = "in",
                res    = 600
                )
            ComparePlates(use_summary_df,
                          use_metrics[[j]],
                          use_cex          = use_cex,
                          beeswarm_spacing = beeswarm_spacing,
                          beeswarm_corral  = beeswarm_corral,
                          side_space       = side_space,
                          order_by_rank    = order_by_rank
                          )
            dev.off()
          }
        }
        message("")
      }
    }
  }
  return(invisible(NULL))
}












