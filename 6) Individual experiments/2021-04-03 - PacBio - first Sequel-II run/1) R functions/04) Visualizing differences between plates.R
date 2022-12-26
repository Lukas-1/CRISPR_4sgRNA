### 2nd August 2021 ###




# Import packages and source code -----------------------------------------

library("beeswarm")
library("RColorBrewer")





# Functions for drawing swarm plots ---------------------------------------

ComparePlates <- function(summary_df,
                          show_column,
                          use_cex          = 0.175,
                          beeswarm_spacing = 0.7,
                          beeswarm_corral  = "none",
                          side_space       = 0.2,
                          order_by_rank    = TRUE,
                          exclude_controls = FALSE
                          ) {

  stopifnot("plates_df" %in% ls(envir = globalenv()))

  if (exclude_controls) {
    control_plates <- plates_df[plates_df[, "Colony_picked"], "Plate_number"]
    summary_df <- summary_df[!(summary_df[, "Plate_number"]) %in% control_plates, ]
    plates_df <- plates_df[!(plates_df[, "Plate_number"]) %in% control_plates, ]
  }

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
         col = brewer.pal(9, "Blues")[[7]],
         xpd = NA
         )

  plate_labels <- paste0("#", plates_df[["Plate_number"]], " ", plates_df[["Plate_name"]])

  if ("Highlight_color" %in% names(plates_df)) {
    label_colors <- plates_df[, "Highlight_color"]
  } else {
    label_colors <- "black"
  }

  text(x      = seq_len(num_groups),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.4), from = "lines", to = "user")),
       labels = plate_labels[plates_order],
       srt    = 45,
       adj    = c(1, 0.5),
       cex    = 0.6,
       xpd    = NA,
       col    = label_colors
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
                                    order_by_rank    = TRUE,
                                    exclude_CCS3     = FALSE,
                                    exclude_controls = FALSE
                                    ) {

  stopifnot("titles_list" %in% ls(envir = globalenv()))

  if (export_PNGs) {
    file_formats <- c("png", "pdf")
  } else {
    file_formats <- "pdf"
  }

  ccs_numbers <- c(3, 5, 7)
  accuracy_percentages <- c(99, 99.9, 99.99)
  df_list_names <- paste0("ccs", ccs_numbers, "_df_list")
  ccs_are_present <- df_list_names %in% ls(envir = globalenv())
  ccs_numbers <- ccs_numbers[ccs_are_present]
  accuracy_percentages <- accuracy_percentages[ccs_are_present]

  if (exclude_CCS3) {
    ccs_numbers <- ccs_numbers[-1]
    accuracy_percentages <- accuracy_percentages[-1]
  }

  use_metrics <- setdiff(names(titles_list), "Longest_subsequence")
  use_metrics <- grep("^Binary_", use_metrics, value = TRUE, invert = TRUE)

  for (file_format in file_formats) {
    for (i in seq_along(ccs_numbers)) {
      use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))

      filter_stages <- c("original_summary_df", "filtered_summary_df", "filtered_cross_plate_df")
      filter_labels <- c("i) unfiltered", "ii) filtered", "iii) filtered cross-plate")
      df_are_present <- filter_stages %in% names(use_df_list)
      filter_stages <- filter_stages[df_are_present]
      filter_labels <- filter_labels[df_are_present]

      for (filter_stage in seq_along(filter_stages)) {
        df_name <- filter_stages[[filter_stage]] # "filtered_gRNAs_df"

        use_summary_df <- use_df_list[[df_name]]

        sel_name <- paste0("CCS", ccs_numbers[[i]],
                           " (", accuracy_percentages[[i]], ") - ",
                           filter_labels[[filter_stage]]
                           )
        message(paste0("Exporting ", file_format, " images for the following data: ", sel_name, "..."))
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
                          order_by_rank    = order_by_rank,
                          exclude_controls = exclude_controls
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
                          order_by_rank    = order_by_rank,
                          exclude_controls = exclude_controls
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




# Functions for drawing ridgeline plots -----------------------------------

MakeDensityPolygons <- function(x_vec, y_vec) {
  are_empty <- y_vec == 0
  group_lengths <- rle(are_empty)[["lengths"]]
  rle_vec <- rep(seq_along(group_lengths), group_lengths)
  polygon_df <- data.frame(
    "x"        = x_vec,
    "y"        = y_vec,
    "Is_empty" = are_empty,
    "Group"    = rle_vec,
    stringsAsFactors = FALSE
  )
  polygon_df <- polygon_df[!(are_empty), ]
  row.names(polygon_df) <- NULL
  polygon_df[, "Group"] <- match(polygon_df[, "Group"], unique(polygon_df[, "Group"]))
  split_df_list <- split(polygon_df, polygon_df[, "Group"])
  mat_list <- lapply(split_df_list, function(x) {
    cbind(x = c(x[, "x"][[1]], x[, "x"], x[, "x"][[nrow(x)]]),
          y = c(0, x[, "y"], 0)
          )
  })
  return(mat_list)
}



PlotDensityList <- function(density_list,
                            use_colors,
                            y_positions   = rev(seq_along(density_list)),
                            min_y         = 0,
                            ridge_height  = 3,
                            border_colors = "white",
                            use_alpha     = 0.8
                            ) {
  fill_colors <- adjustcolor(use_colors, alpha.f = use_alpha)
  x_list <- lapply(density_list, "[[", "x")
  y_list <- lapply(density_list, "[[", "y")
  y_list <- lapply(y_list, function(x) (x - min(x)) / max(x - min(x)))
  original_y_list <- y_list
  if (min_y > 0) {
    y_list <- lapply(y_list, function(x) ifelse(x >= min_y, x, 0))
  }
  y_list <- lapply(y_list, "*", ridge_height)
  polygon_mat_list_list <- mapply(MakeDensityPolygons, x_list, y_list, SIMPLIFY = FALSE)
  y_list <- Map("+", y_list, y_positions)
  polygon_mat_list_list <- lapply(seq_along(density_list), function(x) {
    mat_list <- polygon_mat_list_list[[x]]
    mat_list <- lapply(mat_list, function(y) {
      y[, "y"] <- y[, "y"] + y_positions[[x]]
      return(y)
    })
    return(mat_list)
  })
  fills_vec <- rep(fill_colors, length.out = length(density_list))
  borders_vec <- rep(border_colors, length.out = length(density_list))
  lapply(seq_along(density_list), function(x) {
    for (polygon_mat in polygon_mat_list_list[[x]]) {
      polygon(x      = polygon_mat[, "x"],
              y      = polygon_mat[, "y"],
              col    = fills_vec[[x]],
              border = NA,
              xpd    = NA
              )
    }
    lines(x      = x_list[[x]],
          y      = ifelse(original_y_list[[x]] > 0.05, y_list[[x]], NA),
          col    = borders_vec[[x]],
          xpd    = NA
          )

  })
  return(invisible(NULL))
}



PlateRidgelines <- function(summary_df,
                            num_slots       = NULL,
                            fill_colors     = c("#452e84", "#baace2"),
                            border_colors   = c("#ddd5f0", "#795dc6"),
                            show_title      = NULL,
                            embed_PNG       = FALSE,
                            PNG_lwd         = 1,
                            ridge_height    = 3,
                            zebra_pattern   = FALSE,
                            zebra_palify    = 0.4,
                            plate_prefix    = FALSE,
                            label_on_right  = FALSE,
                            plate_label_cex = 0.4,
                            grid_lwd        = 1,
                            reorder_plates  = FALSE,
                            for_powerpoint  = FALSE,
                            ...
                            ) {

  stopifnot("plates_df" %in% ls(envir = globalenv()))

  ## Prepare plate labels and re-order plates
  control_plates <- plates_df[plates_df[, "Colony_picked"], "Plate_number"]
  summary_df <- summary_df[!(summary_df[, "Plate_number"]) %in% control_plates, ]

  internal_plate_numbers <- unique(summary_df[, "Plate_number"])
  matches_vec <- match(internal_plate_numbers, plates_df[, "Plate_number"])
  plate_names <- plates_df[, "Plate_name"][matches_vec]
  plate_names <- sub("-beads", "", plate_names, fixed = TRUE)
  plate_names <- sub("plus", "+", plate_names, fixed = TRUE)
  plate_names <- sub("and", "&", plate_names, fixed = TRUE)
  first_plates <- sapply(strsplit(plate_names, " ", fixed = TRUE), "[[", 1)
  plate_name_numbers <- as.integer(gsub("[^0-9]", "", first_plates))
  new_order <- order(plate_name_numbers)
  plate_names <- plate_names[new_order]
  internal_plate_numbers <- internal_plate_numbers[new_order]
  if (!(plate_prefix)) {
    plate_names <- paste0("Plate ", gsub("H[AO]_", "", plate_names))
    plate_names[plate_names == "Plate 5 & 53"] <- "Plates 5 & 53"
  }
  if (for_powerpoint) {
    plate_names[plate_names == "Plate 5+"] <- "Plate 05"
    plate_names <- gsub(" ", "_", plate_names, fixed = TRUE)
    plate_list <- strsplit(plate_names, "_", fixed = TRUE)
    plate_names <- vapply(plate_list, function(x) paste0(rev(x), collapse = "_"), "")
  }

  if (reorder_plates) {
    plates_fac <- factor(summary_df[, "Plate_number"], levels = unique(summary_df[, "Plate_number"]))
    mean_percentages <- tapply(summary_df[, "Count_all_4"] / summary_df[, "Count_total"],
                               plates_fac,
                               mean,
                               na.rm = TRUE
                               )
    mean_percentages <- sort(mean_percentages, decreasing = TRUE)
    new_order <- order(match(summary_df[, "Plate_number"], names(mean_percentages)))
    plate_names <- plate_names[order(match(internal_plate_numbers, names(mean_percentages)))]
  } else {
    new_order <- order(match(summary_df[, "Plate_number"], internal_plate_numbers))
  }

  summary_df <- summary_df[new_order, ]
  plates_fac <- factor(summary_df[, "Plate_number"], levels = unique(summary_df[, "Plate_number"]))

  ## Compute densities
  DensitiesForColumn <- function(use_column, omit_NA = TRUE) {
    integer_vec <- summary_df[, use_column]
    if (omit_NA) {
      fraction_vec <- (integer_vec / summary_df[, "Count_total"])
      are_NA <- is.na(fraction_vec)
      fraction_vec <- fraction_vec[!(are_NA)]
      use_fac <- plates_fac[!(are_NA)]
    } else {
      use_fac <- plates_fac
      are_NA <- is.na(integer_vec)
      integer_vec[are_NA] <- 0
      fraction_vec <- integer_vec / summary_df[, "Count_total"]
      fraction_vec[is.nan(fraction_vec)] <- 0
    }
    tapply(fraction_vec, use_fac, density, from = 0, to = 1)
  }

  # summary_df[, "Count_fewer_than_two"] <- summary_df[, "Count_total"] - summary_df[, "Count_at_least_2"]
  # less_two_density_list <- DensitiesForColumn("Count_fewer_than_two", omit_NA = TRUE)
  # deletion_density_list <- DensitiesForColumn("Num_reads_with_sgRNA_deletion", omit_NA = TRUE)
  # atleast3_density_list <- DensitiesForColumn("Count_at_least_3", omit_NA = FALSE)

  all4_density_list     <- DensitiesForColumn("Count_all_4", omit_NA = FALSE)
  atleast1_density_list <- DensitiesForColumn("Count_at_least_1", omit_NA = FALSE)

  ## Prepare x and y axes
  if (is.null(num_slots)) {
    num_slots <- nlevels(plates_fac)
  }
  y_positions <- rev(seq_len(nlevels(plates_fac))) + (num_slots - nlevels(plates_fac))
  bottom_gap <- 1
  y_top <- num_slots + 2
  y_limits <- c(1 - bottom_gap, num_slots + ridge_height + 0.5)

  ## Set up the plot region
  if (embed_PNG) {
    current_device <- StartEmbedPNG(plots_output_directory, add_padding = TRUE,
                                    use_cairo = TRUE, png_res = 600
                                    )
  }

  plot.new()
  plot.window(xlim = c(0, 1), ylim = y_limits, xaxs = "i", yaxs = "i")

  segments(x0  = seq(0, 1, by = 0.1),
           y0  = par("usr")[[3]],
           y1  = y_top,
           col = "gray80",
           lwd = par("lwd") * grid_lwd,
           xpd = NA
           )
  segments(x0  = seq(0.05, 0.95, by = 0.1),
           y0  = par("usr")[[3]],
           y1  = y_top,
           col = "gray90",,
           lwd = par("lwd") * grid_lwd,
           xpd = NA
           )

  ## Draw densities
  PlotDensityList(all4_density_list,
                  if (zebra_pattern) c(Palify(fill_colors[[1]], 0.2), fill_colors[[1]]) else fill_colors[[1]],
                  border_colors[[1]],
                  min_y        = 0.02,
                  y_positions  = y_positions,
                  ridge_height = ridge_height,
                  ...
                  )
  PlotDensityList(atleast1_density_list,
                  if (zebra_pattern) c(Palify(fill_colors[[2]], zebra_palify), fill_colors[[2]]) else fill_colors[[2]],
                  border_colors[[2]],
                  min_y        = 0.04,
                  y_positions  = y_positions,
                  ridge_height = ridge_height,
                  ...
                  )


  if (embed_PNG) {
    current_device <- StopEmbedPNG(current_device, plots_output_directory, add_padding = TRUE)
  }

  ## Annotate plot
  mtext(plate_names,
        side = if (label_on_right) 4 else 2,
        line = if (label_on_right) 0.5 else 0.3,
        at   = y_positions,
        cex  = par("cex") * plate_label_cex,
        col  = if (zebra_pattern) c("gray67", "gray50") else "gray50",
        las  = 2,
        padj = 0
        )

  x_ticks <- seq(0, 1, by = 0.1)
  axis(1, at = x_ticks, labels = x_ticks * 100, tcl = -0.325, mgp = c(3, 0.35, 0),
       gap.axis = 0.25, lwd = par("lwd")
       )
  mtext("Percentage of reads", side = 1, line = 1.5, cex = par("cex"))


  ## Draw legend
  original_text_vec <- c(expression(scriptstyle(" ") * "4 correct gRNAs"), expression("" >= "1 correct gRNA"))
  text_vec <- sapply(original_text_vec, VerticalAdjust)
  pch1_pos <- 0
  y_pos <- par("usr")[[3]] - diff(grconvertY(c(0, 3.7), from = "lines", to = "user"))
  point_y_pos <- y_pos + ((par("cxy")[2] / pi) * par("cex") * 0.2)
  point_radius <- (par("cxy")[1] / pi) * par("cex")
  pch_gap <- point_radius * 1.7
  text1_pos <- pch1_pos + pch_gap
  pch2_pos <- text1_pos + strwidth(original_text_vec[[1]]) + diff(grconvertX(c(0, 1.25), from = "char", to = "user"))
  text2_pos <- pch2_pos + pch_gap
  x_span <- text2_pos + strwidth(original_text_vec[[2]])
  min_x <- grconvertX(0.5, from = "npc", to = "user") - (x_span / 2)
  max_x <- grconvertX(0.5, from = "npc", to = "user") + (x_span / 2)
  points(x   = rescale(pch1_pos, from = c(0, x_span), to = c(min_x, max_x)),
         y   = point_y_pos,
         pch = 22,
         col = if (sum(col2rgb(border_colors[[1]])) < sum(col2rgb(fill_colors[[1]]))) border_colors[[1]] else "black",
         bg  = fill_colors[[1]],
         cex = 1.4,
         xpd = NA
         )
  text(x      = rescale(text1_pos - strwidth("gh"), from = c(0, x_span), to = c(min_x, max_x)),
       y      = y_pos,
       labels = text_vec[[1]],
       adj    = c(0, 0.5),
       xpd    = NA
       )
  points(x   = rescale(pch2_pos, from = c(0, x_span), to = c(min_x, max_x)),
         y   = point_y_pos,
         pch = 22,
         col = if (sum(col2rgb(border_colors[[2]])) < sum(col2rgb(fill_colors[[2]]))) border_colors[[2]] else "black",
         bg  = fill_colors[[2]],
         cex = 1.4,
         xpd = NA
         )
  text(x      = rescale(text2_pos - strwidth("gh"), from = c(0, x_span), to = c(min_x, max_x)),
       y      = y_pos,
       labels = text_vec[[2]],
       adj    = c(0, 0.5),
       xpd    = NA
       )

  if (!(is.null(show_title))) {
    text(x      = grconvertX(0.5, from = "npc", to = "user"),
         y      = y_top + diff(grconvertY(c(0, 0.7), from = "lines", to = "user")),
         labels = show_title,
         adj    = c(0.5, 0),
         xpd    = NA
         )
  }

  return(invisible(NULL))
}

