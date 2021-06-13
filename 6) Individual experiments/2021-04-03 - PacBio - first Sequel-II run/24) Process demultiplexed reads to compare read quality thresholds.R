### 12th June 2021 ###




# Import packages and source code -----------------------------------------

library("vioplot")
library("png")

CRISPR_root_directory      <- "~/CRISPR"
plate1_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
sql2_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
R_functions_directory      <- file.path(plate1_directory, "1) R functions")
sql2_R_functions_directory <- file.path(sql2_directory, "1) R functions")

source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))
source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions
source(file.path(sql2_R_functions_directory, "03) Summarizing data across wells.R"))




# Define folder paths -----------------------------------------------------

sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output", "Figures", "Explore read quality cut-offs")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads - plates_analysis_list.RData"))
load(file.path(sql2_R_objects_directory, "22) Summarize data across wells - plate selections.RData"))





# Define maps -------------------------------------------------------------

use_width <- 8
use_height <- 6

titles_list <- list(
  "Num_contaminated_reads" = "Reads that match gRNAs from other wells",

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


export_metrics <- c("Num_contaminated_reads",
                    "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
                    "Count_all_4_promoters", "Count_whole_plasmid",
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
  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))
  min_quality <- GetMinQuality(input_ccs_number)
  reads_df <- analysis_list[["individual_reads_df"]]
  are_passing <- (reads_df[, "Passes_barcode_filters"] == 1) &
                 (reads_df[, "Read_quality"] >= min_quality) &
                 (reads_df[, "Num_full_passes"] >= input_ccs_number)
  selected_zmws <- reads_df[["ZMW"]][are_passing]
  summary_df_list <- SummarizeWells(plates_analysis_list,
                                    use_zmws = selected_zmws,
                                    ID_column = "Combined_ID",
                                    unique_IDs = sg_sequences_df[["Combined_ID"]]
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

  results_vec_list <- lapply(all_ccs_list, function(x) {
    are_these_plates <- x[[use_df_name]][, "Plate_number"] %in% plate_numbers
    stopifnot(any(are_these_plates))
    metric_vec <- x[[use_df_name]][are_these_plates, use_column]
    total_vec <- x[[use_df_name]][are_these_plates, "Count_total"]
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

  require_objects <- c("plate_selection_list", "plate_selection_titles_list")
  stopifnot(all(require_objects %in% ls(envir = globalenv())))
  if (is.null(plate_names)) {
    plate_names <- "All plates"
  }
  if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
    plate_names <- plate_selection_list[[plate_names]]
  } else {
    main_title <- "PacBio sequencing"
  }

  metric_list <- ExtractMetricForAllCutoffs(all_ccs_list,
                                            use_column,
                                            plate_names,
                                            filter_mean_quality = filter_mean_quality
                                            )


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
  title(main_title, cex.main = 1.1, font.main = 1, line = 2)


  ## Draw the x axis labels
  large_gap_lines <- 2
  start_line <- 0.75

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
  }

  par(old_mar)

  return(invisible(NULL))
}




# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)





# Summarize wells for different read quality cut-offs ---------------------

summary_list_list <- SummarizeWellsForAllCutoffs(plates_analysis_list)





# Modify plate selections -------------------------------------------------

plate_selection_titles_list <- c(
  list("Colony-picked" = "Single-colony picked controls"),
  plate_selection_titles_list
)





# Produce example plots ---------------------------------------------------

ViolinBoxAllCutoffs(summary_list_list, "Count_all_4", "All plates")

ViolinBoxAllCutoffs(summary_list_list, "Count_all_4", "All plates",
                    filter_mean_quality = TRUE
                    )

ViolinBoxAllCutoffs(summary_list_list, "Num_contaminated_reads", "All plates")


ViolinBoxAllCutoffs(summary_list_list, "Num_contaminated_reads",
                    "Colony-picked", use_y_limits = c(0, 0.1)
                    )





# Export violin/box plots -------------------------------------------------

for (filter_stage in 1:2) {
  sub_folder_name <- c("i) unfiltered", "ii) filtered")[[filter_stage]]
  folder_path <- file.path(file_output_directory, sub_folder_name)

  for (i in seq_along(export_metrics)) {
    metric <- export_metrics[[i]]
    file_name <- paste0("Quality cut-offs - ", i, ") ", metric, ".pdf")
    pdf(file = file.path(folder_path, file_name),
        width = use_width, height = use_height
        )
    for (plate_selection in names(plate_selection_titles_list)) {
      if ((metric == "Num_contaminated_reads") && (plate_selection == "Colony-picked")) {
        custom_y_limits <- c(0, 0.1)
      } else {
        custom_y_limits <- c(0, 1)
      }
      print(plate_selection)
      ViolinBoxAllCutoffs(summary_list_list,
                          use_column          = metric,
                          plate_names         = plate_selection,
                          filter_mean_quality = filter_stage == 2,
                          use_y_limits        = custom_y_limits,
                          embed_PNG           = TRUE
                          )
    }
    dev.off()
  }
}






