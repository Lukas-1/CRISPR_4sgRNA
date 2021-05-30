### 28 May 2020 ###



# Import packages and source code -----------------------------------------

library("beeswarm")
library("RColorBrewer")

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Compare plates")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))




# Prepare titles ----------------------------------------------------------

titles_list <- c(list("Count_total" = "Number of reads per well"),
                 titles_list
                 )




# Define functions --------------------------------------------------------

ComparePlates <- function(summary_df, show_column) {

  stopifnot("plates_df" %in% ls(envir = globalenv()))

  if (show_column == "Count_mean_sg1to4") {
    count_columns <- paste0("Count_sg", 1:4, "_cr", 1:4)
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
    numeric_vec <- numeric_vec
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

  groups_fac <- factor(summary_df[["Plate_number"]],
                       levels = plates_df[["Plate_number"]][order(plates_df[["Plate_rank"]])]
                       )
  num_groups <- nlevels(groups_fac)
  stopifnot(num_groups == nrow(plates_df))

  use_space <- 0.3
  group_limits <- c((1 - use_space) - (num_groups * 0.04),
                    (num_groups + use_space) + (num_groups * 0.04)
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

  boxplot(numeric_vec ~ groups_fac,
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
          )

  set.seed(1)
  use_cex <- 0.15
  beeswarm_df <- beeswarm(numeric_vec ~ groups_fac,
                          spacing  = 0.8,
                          priority = "random",
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

  assign("delete_beeswarm_df", beeswarm_df, envir = globalenv())

  plate_labels <- paste0("#", plates_df[["Plate_number"]], " ", plates_df[["Plate_name"]])

  text(x      = seq_len(num_groups),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.4), from = "lines", to = "user")),
       labels = plate_labels[order(plates_df[["Plate_rank"]])],
       srt    = 45,
       adj    = c(1, 0.5),
       cex    = 0.6,
       xpd    = NA
       )
  return(invisible(NULL))
}




# Export graphics ---------------------------------------------------------

ccs_numbers <- c(3, 5, 7)
accuracy_percentages <- c(99, 99.9, 99.99)

use_metrics <- setdiff(names(titles_list), "Longest_subsequence")

for (i in seq_along(ccs_numbers)) {
  use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))
  for (filter_stage in 1:2) {
    df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"

    use_summary_df <- use_df_list[[df_name]]

    file_name <- paste0("CCS", ccs_numbers[[i]],
                        " (", accuracy_percentages[[i]], ") - ",
                        c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                        )

    pdf(file   = file.path(plots_output_directory, paste0("Compare plates - ", file_name, ".pdf")),
        width  = 3.9 * 2.5,
        height = 2.6 * 2.5
        )
    for (use_metric in use_metrics) {
      ComparePlates(use_summary_df, use_metric)
    }
    dev.off()
  }
}






