### 5th August 2021 ###



# Define maps -------------------------------------------------------------

use_width <- 8
use_height <- 6

count_titles_list <- list(
  "Count_sg1_cr1"         = "Reads where sg1 (+ tracrRNA) is correct",
  "Count_sg2_cr2"         = "Reads where sg2 (+ tracrRNA) is correct",
  "Count_sg3_cr3"         = "Reads where sg3 (+ tracrRNA) is correct",
  "Count_sg4_cr4"         = "Reads where sg3 (+ tracrRNA) is correct",

  "Count_at_least_1"      = "Reads with at least 1 correct gRNA",
  "Count_at_least_2"      = "Reads with at least 2 correct gRNAs",
  "Count_at_least_3"      = "Reads with at least 3 correct gRNAs",
  "Count_all_4"           = "Reads where all 4 gRNAs are correct",

  "Count_all_4_promoters" = "Reads with 4 correct gRNAs (including promoters)",
  "Count_whole_plasmid"   = "Reads where the entire plasmid sequence is correct"
)

count_no_contam_titles_list <- count_titles_list
count_no_contam_titles_list <- paste0(count_no_contam_titles_list,
                                      " (excluding contaminations)"
                                      )
names(count_no_contam_titles_list) <- sub("Count_",
                                          "Count_no_contam_",
                                          names(count_titles_list),
                                          fixed = TRUE
                                          )

titles_list <- c(
  count_titles_list,
  count_no_contam_titles_list,
  "Num_contaminated_reads"                      = "Reads that match gRNAs from other wells",
  "Num_contaminated_reads_aligned"              = "Aligned reads that match gRNAs from other wells",
  "Num_reads_with_deletions_exceeding_20bp"     = "Reads with deletions (of 20 base pairs or more)",
  "Num_reads_with_deletions_spanning_tracrRNAs" = "Reads with deletions that span tracrRNAs",
  "Num_reads_with_deletions_spanning_promoters" = "Reads with deletions that span promoters",
  "Num_reads_with_deletions_spanning_sg_cr"     = "Reads with deletions that span sg+cr sequences",
  "Num_reads_with_deletions_spanning_sgRNAs"    = "Reads with deletions that span gRNAs"
)


export_metrics <- c("Num_contaminated_reads",
                    "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
                    "Count_all_4_promoters", "Count_whole_plasmid",
                    "Num_reads_with_deletions_exceeding_20bp",
                    "Num_reads_with_deletions_spanning_tracrRNAs",
                    "Count_no_contam_all_4",
                    "Count_sg1_cr1", "Count_sg2_cr2", "Count_sg3_cr3", "Count_sg4_cr4"
                    )




# Define functions --------------------------------------------------------

GetMinQuality <- function(input_integer) {
  is_even <- (input_integer %% 2) == 0
  num_nines <- ceiling(input_integer / 2)
  digits_string <- paste0(rep("9", num_nines), collapse = "")
  if (is_even) {
    digits_string <- paste0(digits_string, "5")
  }
  final_number <- as.numeric(paste0("0.", digits_string))
  return(final_number)
}



CustomSummarizeWells <- function(analysis_list, input_ccs_number) {
  required_objects <- c("sg_sequences_df", "contam_df", "deletions_df")
  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))
  min_quality <- GetMinQuality(input_ccs_number)
  reads_df <- analysis_list[["individual_reads_df"]]
  are_passing <- (reads_df[, "Passes_barcode_filters"] == 1) &
                 (reads_df[, "Read_quality"] >= min_quality) &
                 (reads_df[, "Num_full_passes"] >= input_ccs_number)
  selected_zmws <- reads_df[["ZMW"]][are_passing]
  summary_df_list <- SummarizeWells(plates_analysis_list,
                                    use_zmws          = selected_zmws,
                                    ID_column         = "Combined_ID",
                                    unique_IDs        = sg_sequences_df[["Combined_ID"]],
                                    deletions_df      = deletions_df,
                                    aligned_contam_df = contam_df
                                    )
  return(summary_df_list)
}



SummarizeWellsForAllCutoffs <- function(analysis_list) {
  num_passes <- seq_len(11)
  min_qualities <- vapply(num_passes, GetMinQuality, numeric(1))
  full_names <- paste0("CCS", num_passes, " / ", min_qualities * 100, "%")
  results_list_list <- lapply(seq_along(num_passes), function(x) {
    message(paste0("Calculating per-well summaries for ", full_names[[x]], "..."))
    return(CustomSummarizeWells(analysis_list, num_passes[[x]]))
  })
  names(results_list_list) <- full_names
  return(results_list_list)
}



ExtractMetricForAllCutoffs <- function(all_ccs_list,
                                       use_column,
                                       plate_names,
                                       filter_mean_quality = FALSE
                                       ) {

  stopifnot("plates_df" %in% ls(envir = globalenv()))

  plate_matches <- match(plate_names, plates_df[, "Plate_name"])
  plate_numbers <- plates_df[["Plate_number"]][plate_matches]

  if (filter_mean_quality) {
    use_df_name <- "filtered_summary_df"
  } else {
    use_df_name <- "original_summary_df"
  }
  if (grepl("Count_no_contam_", use_column, fixed = TRUE)) {
    total_column <- "Count_total_no_contam"
  } else {
    total_column <- "Count_total"
  }

  results_vec_list <- lapply(all_ccs_list, function(x) {
    are_these_plates <- x[[use_df_name]][, "Plate_number"] %in% plate_numbers
    stopifnot(any(are_these_plates))
    metric_vec <- x[[use_df_name]][are_these_plates, use_column]
    total_vec <- x[[use_df_name]][are_these_plates, total_column]
    stopifnot(identical("integer", class(metric_vec)))
    stopifnot(all(total_vec >= metric_vec, na.rm = TRUE))
    results_vec <- metric_vec / total_vec
    return(results_vec)
  })
  return(results_vec_list)
}




ViolinBoxAllCutoffs <- function(all_ccs_list,
                                use_column,
                                plate_names,
                                filter_mean_quality = FALSE,
                                use_y_limits        = c(0, 1),
                                set_mar             = TRUE,
                                embed_PNG           = FALSE
                                ) {

  set.seed(1) # For reproducible jitter

  plates_list <- GetPlateSelection(plate_names)
  plate_names <- plates_list[["plate_names"]]

  metric_list <- ExtractMetricForAllCutoffs(all_ccs_list,
                                            use_column,
                                            plate_names,
                                            filter_mean_quality = filter_mean_quality
                                            )
  if (filter_mean_quality) {
    df_name <- "filtered_summary_df"
  } else {
    df_name <- "original_summary_df"
  }
  if (grepl("Count_no_contam_", use_column, fixed = TRUE)) {
    total_column <- "Count_total_no_contam"
  } else {
    total_column <- "Count_total"
  }
  group_n <- vapply(all_ccs_list, function(x) sum(x[[df_name]][, total_column]), integer(1))

  name_splits <- strsplit(names(all_ccs_list), " / ", fixed = TRUE)
  ccs_labels <- sapply(name_splits, "[[", 1)
  ccs_labels <- sub("CCS", "", ccs_labels, fixed = TRUE)

  ccs_labels <- sapply(ccs_labels,
                       function(x) as.expression(bquote("CCS" * scriptscriptstyle(" ") * bold(.(x)))),
                       USE.NAMES = FALSE
                       )
  quality_labels <- sapply(name_splits, "[[", 2)

  if (set_mar) {
    old_mar <- par("mar" = c(5, 5, 4.5, 2))
  }

  ## Determine group positions
  num_groups <- length(metric_list)
  group_positions <- seq_len(num_groups)
  group_limits <- c((min(group_positions) - 0.3) - (num_groups * 0.04),
                     max(group_positions) + 0.3  + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  numeric_axis_pos <- pretty(use_y_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  numeric_axis_labels <- paste0(format(numeric_axis_pos * 100), "%")

  final_y_range <- numeric_limits[[2]] - numeric_limits[[1]]
  y_gap <- final_y_range * 0.0075
  final_y_limits <- c(numeric_limits[[1]] - y_gap, numeric_limits[[2]] + y_gap)


  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(file_output_directory, "temp.png")
    cur_mai <- par("mar") * 0.2
    temp_width <- use_width - sum(cur_mai[c(2, 4)])
    temp_height <- use_height - sum(cur_mai[c(1, 3)])
    png(file   = temp_path,
        width  = temp_width,
        height = temp_height,
        units  = "in",
        res    = 900,
        bg     = "transparent"
        )
    par(mar = rep(0, 4))
  }


  ## Set up the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = final_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  DrawGridlines(use_y_limits)

  ## Draw the violin
  use_wex <- 0.8
  vioplot(metric_list,
          at       = group_positions,
          pchMed   = NA,
          drawRect = FALSE,
          col      = brewer.pal(9, "Blues")[[4]],
          border   = NA,
          wex      = use_wex,
          add      = TRUE,
          axes     = FALSE
          )


  ## Draw the jittered points
  metric_unlisted <- unlist(metric_list, use.names = FALSE)
  jittered_vec  <- group_positions[rep(seq_along(metric_list), lengths(metric_list))] +
                   rnorm(n = length(metric_unlisted), mean = 0, sd = 0.05)
  points_alpha <- 0.2
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
  points(x   = jittered_vec,
         y   = metric_unlisted,
         cex = 0.4,
         col = paste0(brewer.pal(9, "Blues")[[8]], alpha_hex),
         pch = 16
         )


  if (embed_PNG) {
    dev.off()
    raster_array <- readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)
    plot(1,
       xlim = group_limits,
       ylim = final_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
    rasterImage(raster_array,
                xleft = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                )
  }


  ## Draw the superimposed boxplots
  boxplot(metric_list,
          at         = group_positions,
          boxwex     = use_wex * 0.4,
          outline    = FALSE,
          names      = rep.int("", length(group_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * 3,
          col        = brewer.pal(9, "Blues")[[2]],
          border     = brewer.pal(9, "Blues")[[9]],
          add        = TRUE,
          axes       = FALSE,
          lwd        = 1
          )


  ## Draw the axis and axis label
  axis(2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, 0.38, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1
       )

  y_axis_label <- titles_list[[use_column]]
  mtext(text = y_axis_label, side = 2, line = 3)


  ## Draw the title
  title(plates_list[["title"]], cex.main = 1.1, font.main = 1, line = 2)


  ## Draw the x axis labels
  large_gap_lines <- 2
  start_line <- 0.7

  group_numbers <- formatC(group_n, big.mark = "'")

  for (i in seq_along(group_positions)) {
    x_label_cex <- 0.7
    mtext(text = VerticalAdjust(ccs_labels[[i]]),
          at   = group_positions[[i]],
          line = start_line,
          side = 1,
          cex  = x_label_cex
          )

    mtext(text = VerticalAdjust(quality_labels[[i]]),
          at   = group_positions[[i]],
          line = start_line + (large_gap_lines / 2),
          side = 1,
          cex  = x_label_cex
          )

    mtext(text = VerticalAdjust(group_numbers[[i]]),
          at   = group_positions[[i]],
          line = start_line + (large_gap_lines),
          side = 1,
          cex  = 0.5,
          col  = "gray60"
          )
  }

  par(old_mar)
  return(invisible(NULL))
}




DrawAllCutoffComparisons <- function(export_PNGs = TRUE) {

  stopifnot("export_metrics" %in% ls(envir = globalenv()))


  metrics_seq <- seq_along(export_metrics)
  file_numbers <- formatC(seq_along(export_metrics),
                          flag = "0",
                          width = max(nchar(as.character(metrics_seq)))
                          )


  if (export_PNGs) {
    file_formats <- c("png", "pdf")
  } else {
    file_formats <- "pdf"
  }

  for (file_format in file_formats) {
    for (filter_stage in 1:2) {
      sub_folder_name <- c("i) unfiltered", "ii) filtered")[[filter_stage]]

      message(paste0("Exporting ", file_format, " images into the following folder: ", sub_folder_name, "..."))

      if (file_format == "pdf") {
        folder_path <- file.path(file_output_directory, sub_folder_name)
      } else if (file_format == "png") {
        folder_path <- file.path(PNGs_output_directory, sub_folder_name)
      }
      for (i in metrics_seq) {
        metric <- export_metrics[[i]]

        message(paste0(" ... for the metric: '", metric, "'"))

        metric_name <- sub("^Num_reads_with_", "", metric)
        metric_name <- sub("^(Num|Count)_", "", metric_name)
        metric_name <- paste0(toupper(substr(metric_name, 1, 1)),
                              substr(metric_name, 2, nchar(metric_name))
                              )
        metric_name <- paste0(file_numbers[[i]], ") ", metric_name)
        if (file_format == "pdf") {
          pdf(file = file.path(folder_path, paste0("Quality cut-offs - ", metric_name, ".pdf")),
              width = use_width, height = use_height
              )
        } else if (file_format == "png") {
          sub_folder_path <- file.path(folder_path, metric_name)
          dir.create(sub_folder_path, showWarnings = FALSE)
        }
        for (j in seq_along(plate_selection_titles_list)) {
          plate_selection <- names(plate_selection_titles_list)[[j]]

          message(paste0("    ... and for the plate selection: '", plate_selection, "'"))

          if ((metric == "Num_contaminated_reads") && (plate_selection == "Colony-picked")) {
            custom_y_limits <- c(0, 0.1)
          } else {
            custom_y_limits <- c(0, 1)
          }
          if (file_format == "png") {
            file_name <- paste0(plate_selection_prefixes[[j]],
                                ") ", plate_selection, ".png"
                                )
            png(file   = file.path(sub_folder_path, file_name),
                width  = use_width,
                height = use_height,
                units  = "in",
                res    = 600
                )
          }
          ViolinBoxAllCutoffs(summary_list_list,
                              use_column          = metric,
                              plate_names         = plate_selection,
                              filter_mean_quality = filter_stage == 2,
                              use_y_limits        = custom_y_limits,
                              embed_PNG           = TRUE
                              )
          if (file_format == "png") {
            dev.off()
          }
        }
        if (file_format == "pdf") {
          dev.off()
        }
      }
    }
  }
  return(invisible(NULL))
}






