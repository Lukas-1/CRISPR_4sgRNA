# 2022-01-17


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define labels -----------------------------------------------------------

controls_labels <- list(
  "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
  "NT"   = c("Non-targeting", "controls"),
  "Gene" = c("Genes in ", "CRISPRa", "library")
)



# Define functions --------------------------------------------------------

AggregateWells <- function(input_df, use_column, by_row = FALSE) {

  are_NT    <- input_df[, "Is_NT_ctrl"]
  are_pos   <- input_df[, "Is_pos_ctrl"]
  are_gene  <- !(is.na(input_df[, "Entrez_ID"]))
  are_valid <- are_NT | are_pos | are_gene

  category_vec <- rep(NA, nrow(input_df))
  category_vec[are_gene] <- "gene"
  category_vec[are_NT]   <- "non-targeting control"
  category_vec[are_pos]  <- "positive control"

  if (!("Plate_ID" %in% names(input_df))) {
    input_df[, "Plate_ID"] <- input_df[, "Plate_number_384"]
  }

  use_columns <- c("Plate_ID", "Well_number_384", "x_position",
                   "Is_NT_ctrl", "Is_pos_ctrl", use_column
                   )
  use_df <- input_df[are_valid, use_columns]
  use_df[, "Category"] <- category_vec[are_valid]

  if (by_row) {
    rows_vec <- rep(1:16, each = 24)
    use_df[, "Row_number"] <- rows_vec[use_df[, "Well_number_384"]]
    division_column <- "Row_number"
  } else {
    columns_vec <- rep(1:24, times = 16)
    use_df[, "Column_number"] <- columns_vec[use_df[, "Well_number_384"]]
    division_column <- "Column_number"
  }

  groups_fac <- interaction(factor(use_df[, "Plate_ID"],
                                   levels = unique(use_df[, "Plate_ID"])
                                   ),
                            use_df[, division_column],
                            use_df[, "Category"],
                            sep = "_",
                            lex.order = TRUE,
                            drop = TRUE
                            )
  assign("delete_groups_fac", groups_fac, envir = globalenv())

  split_df_list <- split(use_df, groups_fac)

  mean_columns <- c("x_position", use_column)
  split_df_list <- lapply(split_df_list, function(x) {
    means_vec <- vapply(mean_columns, function(y) mean(x[, y]), numeric(1))
    x <- x[1, ]
    for (i in seq_along(mean_columns)) {
      x[, mean_columns[[i]]] <- means_vec[[i]]
    }
    return(x)
  })
  results_df <- do.call(rbind.data.frame, c(split_df_list, stringsAsFactors = FALSE, make.row.names = FALSE))
  results_df <- results_df[, names(results_df) != "Well_number_384"]
  results_df[, "Entrez_ID"] <- ifelse(results_df[, "Category"] == "gene", TRUE, NA)
  return(results_df)
}



StackReplicatePlates <- function(input_df) {

  ## Lengthen ("stack") the data frame so that replicate plates follow each other

  rep1_columns <- grep("_rep1", names(input_df), value = TRUE, fixed = TRUE)
  rep2_columns <- grep("_rep2", names(input_df), value = TRUE, fixed = TRUE)
  stripped_columns <- sub("_rep1", "", rep1_columns)
  stopifnot(identical(stripped_columns, sub("_rep2", "", rep2_columns)))

  rep1_df <- input_df[, rep1_columns, drop = FALSE]
  rep2_df <- input_df[, rep2_columns, drop = FALSE]
  names(rep1_df) <- stripped_columns
  names(rep2_df) <- stripped_columns
  reps_df <- rbind.data.frame(rep1_df, rep2_df, make.row.names = FALSE,
                              stringsAsFactors = FALSE
                              )

  are_rep_columns <- grepl("_rep[12]", names(input_df))
  num_rows <- nrow(input_df)
  non_rep_df <- input_df[rep(seq_len(num_rows), 2), !(are_rep_columns)]
  non_rep_df[, "Replicate_number"] <- c(rep(1L, num_rows), rep(2L, num_rows))
  results_df <- data.frame(non_rep_df, reps_df, stringsAsFactors = FALSE)
  new_order <- order(match(results_df[, "Plate_number_384"],
                           results_df[, "Plate_number_384"]
                           ))
  results_df <- results_df[new_order, ]
  row.names(results_df) <- NULL
  return(results_df)
}



PlateWellPlot <- function(input_df,
                          use_column      = "Raw_rep1",
                          show_title      = NULL,
                          y_axis_label    = NULL,
                          point_size      = 0.6,
                          order_by_column = FALSE,
                          aggregate_wells = FALSE,
                          emphasize_NT    = FALSE
                          ) {

  ## Prepare axis and plot labels
  if (is.null(y_axis_label)) {
    if (("short_column_labels" %in% ls(envir = globalenv())) &&
        (use_column %in% names(short_column_labels))
    ) {
      y_axis_label <- FormatPlotMath(short_column_labels[[use_column]])
    } else {
      y_axis_label <- use_column
    }
  }
  if (is.null(show_title)) {
    show_title <- "Plate-well series"
    if (aggregate_wells) {
      if (order_by_column) {
        show_title <- paste0(show_title, ", column means")
      } else {
        show_title <- paste0(show_title, ", row means")
      }
    } else {
      if (order_by_column) {
        show_title <- paste0(show_title, ", ordered by column")
      } else {
        show_title <- paste0(show_title, ", ordered by row")
      }
    }
  }

  ## Re-order wells
  if (order_by_column) {
    columns_vec <- rep(1:24, times = 16)
    long_columns_vec <- columns_vec[input_df[, "Well_number_384"]]
    new_order <- order(match(input_df[, "Plate_number_384"], input_df[, "Plate_number_384"]),
                       long_columns_vec
                       )
    input_df <- input_df[new_order, ]
    row.names(input_df) <- NULL
  }

  ## Lengthen ("stack") the data frame so that replicate plates follow each other
  has_replicates <- grepl("_rep", use_column, fixed = TRUE)
  if (has_replicates) {
    if (!(grepl("_rep1", use_column, fixed = TRUE))) {
      stop("Please supply the column name for replicate 1! Both replicates will be displayed.")
    }
    rep2_column <- sub("_rep1", "_rep2", use_column, fixed = TRUE)
    stack_columns <- c("Plate_number_384", "Well_number_384",
                       "Is_NT_ctrl", "Is_pos_ctrl", "Entrez_ID",
                       use_column, rep2_column
                       )
    input_df <- StackReplicatePlates(input_df[, stack_columns])
  }

  ## Prepare data for plotting
  plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
  if (has_replicates) {
    replicate_numbers_vec <- input_df[, "Replicate_number"]
    plate_reps_vec <- paste0(plate_numbers_vec, "-rep", replicate_numbers_vec)
    input_df[, "Plate_ID"] <- plate_reps_vec
    use_column <- sub("_rep1", "", use_column, fixed = TRUE)
  }
  x_position_vec <- seq_len(nrow(input_df)) # Note: The well number includes empty wells, etc.
  input_df[, "x_position"] <- x_position_vec
  x_limits <- range(x_position_vec)

  if (aggregate_wells) {
    input_df <- AggregateWells(input_df, use_column, by_row = !(order_by_column))
  }

  are_NT      <- input_df[, "Is_NT_ctrl"]
  are_posctrl <- input_df[, "Is_pos_ctrl"]
  are_gene    <- !(is.na(input_df[, "Entrez_ID"]))
  are_valid   <- are_NT | are_posctrl | are_gene

  numeric_vec <- input_df[, use_column]
  valid_vec <- numeric_vec[are_valid]

  ## Prepare axis limits
  x_limits <- x_limits + (diff(x_limits) * 0.04 * c(-1, 1))
  y_limits <- DataAxisLimits(valid_vec)

  ## Prepare graphical parameters
  old_mar <- par(mar = c(7, 4.6, 3.8, 7.5))
  use_mgp <- c(3, 0.65, 0)
  use_tcl <- -0.45

  ## Set up the plot region
  plot(1,
       xlim = x_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       type = "n",
       ann  = FALSE
       )

  ## Annotate the plot
  title(show_title, cex.main = 1.1)
  AbbreviateDataAxis(2, mgp = use_mgp[[2]], tcl = use_tcl)
  axis(1, mgp = use_mgp, cex.axis = par("cex"), tcl = use_tcl)
  mtext(y_axis_label, side = 2, line = 3)
  label_x_pos <- par("usr")[[1]] + diff(grconvertX(c(0, 0.2), from = "lines", to = "user"))
  mtext("Well #:",
        side = 1,
        line = use_mgp[[2]],
        at   = label_x_pos,
        adj  = 1
        )

  ## Draw indicator lines in the plot region
  if (y_limits[[1]] != 0) {
    abline(h = 0, col = "gray75", lty = "dotted")
  }
  if (grepl("PercActivation|[Ff]oldNT|p_value", use_column)) {
    abline(h = 1, col = "gray75", lty = "dotted")
  }

  ## Prepare for drawing the plate number indicator bar
  plate_x_starts <- tapply(x_position_vec, plate_numbers_vec, function(x) x[[1]]) - 0.5
  plate_x_ends   <- tapply(x_position_vec, plate_numbers_vec, function(x) x[[length(x)]]) + 0.5

  rectangle_height <- diff(grconvertY(c(0, 0.9), from = "lines", to = "user"))
  start_y <- par("usr")[[3]] - diff(grconvertY(c(0, 3.2), from = "lines", to = "user"))
  text_colors <- c("gray25", "white")

  ## Draw the plate number indicator bar
  text(x      = label_x_pos,
       y      = start_y - (rectangle_height / 2),
       labels = "Plate #:", adj = c(1, 0.5),
       xpd    = NA
       )
  rect(xleft   = plate_x_starts,
       xright  = plate_x_ends,
       ybottom = start_y - rectangle_height,
       ytop    = start_y,
       col     = brewer.pal(9, "Blues")[c(3, 7)],
       xpd     = NA
       )
  text(x      = rowMeans(cbind(plate_x_starts, plate_x_ends)),
       y      = start_y - (rectangle_height / 2),
       labels = as.character(as.roman(unique(plate_numbers_vec))),
       cex    = par("cex") * 0.8,
       col    = text_colors,
       xpd    = NA
       )
  abline(v = c(plate_x_starts, plate_x_ends[[length(plate_x_ends)]]),
         col = "gray90"
         )

  if (has_replicates) {

    rep_x_starts <- tapply(x_position_vec, plate_reps_vec, function(x) x[[1]]) - 0.5
    rep_x_ends   <- tapply(x_position_vec, plate_reps_vec, function(x) x[[length(x)]]) + 0.5

    ## Draw the replicate number indicator bar
    rep_y_gap <- diff(grconvertY(c(0, 1.8), from = "lines", to = "user"))
    text(x      = label_x_pos,
         y      = start_y - rep_y_gap - (rectangle_height / 2),
         labels = "Replicate #:", adj = c(1, 0.5),
         xpd    = NA
         )
    rect(xleft   = rep_x_starts,
         xright  = rep_x_ends,
         ybottom = start_y - rep_y_gap - rectangle_height,
         ytop    = start_y - rep_y_gap,
         col     = brewer.pal(9, "Purples")[c(3, 7)],
         xpd     = NA
         )
    text(x      = rowMeans(cbind(rep_x_starts, rep_x_ends)),
         y      = start_y - rep_y_gap - (rectangle_height / 2),
         labels = rle(replicate_numbers_vec)[["values"]],
         cex    = par("cex") * 0.8,
         col    = text_colors,
         xpd    = NA
         )

    abline(v = c(rep_x_starts, rep_x_ends[[length(rep_x_ends)]]),
           col = "gray80", lty = "dotted"
           )
  }
  box()

  ## Draw the points of the plate-well series

  pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
  NT_ctrl_color <- colorRampPalette(brewer.pal(5, "Blues")[c(4, 5)])(3)[[2]]
  if (emphasize_NT) {
    gene_point_color <- adjustcolor("gray40", alpha.f = 0.1)
  } else {
    gene_point_color <- adjustcolor("black", alpha.f = 0.3)
  }

  points(input_df[are_gene, "x_position"],
         numeric_vec[are_gene],
         pch = 16,
         col = gene_point_color,
         cex = point_size
         )


  points(input_df[are_posctrl, "x_position"],
         numeric_vec[are_posctrl],
         pch = 16,
         col = adjustcolor(pos_ctrl_color, alpha.f = 0.5),
         cex = point_size
         )
  points(input_df[are_NT, "x_position"],
         numeric_vec[are_NT],
         pch = 16,
         col = adjustcolor(NT_ctrl_color, alpha.f = 0.5),
         cex = point_size
         )

  ## Draw a legend for the points
  DrawSideLegend(labels_list = controls_labels,
                 use_colors = c(pos_ctrl_color, NT_ctrl_color, "black")
                 )

  par(old_mar)

  return(invisible(NULL))
}



ExportAllPlateSeriesPlots <- function(input_df,
                                      top_folder,
                                      plot_width = 8,
                                      plot_height = 5
                                      ) {

  for (export_PDF in c(TRUE, FALSE)) {
    if (export_PDF) {
      message("Exporting PDF files...")
    } else {
      message("Exporting PNG files...")
    }
    for (emphasize_NT in c(FALSE, TRUE)) {
      if (emphasize_NT) {
        plots_path <- file.path(top_folder, "b) NT controls emphasized")
        message("Exporting the version where non-targeting controls are emphasized...")
      } else {
        plots_path <- file.path(top_folder, "a) Standard version")
        message("Exporting the standard version (where NT controls may be obscured by genes)...")
      }
      for (by_row in c(TRUE, FALSE)) {
        for (take_means in c(FALSE, TRUE)) {
          file_name <- "Plate well series plots"
          if (by_row) {
            if (take_means) {
              file_name <- paste0(file_name, ", row means")
              sub_folder_name <- "Row means"
              message("... row means are shown...")
            } else {
              file_name <- paste0(file_name, ", by row")
              sub_folder_name <- "Wells ordered by row"
              message("... wells are ordered by row...")
            }
          } else {
            if (take_means) {
              file_name <- paste0(file_name, ", column means")
              sub_folder_name <- "Column means"
              message("... column means are shown...")
            } else {
              file_name <- paste0(file_name, ", by column")
              sub_folder_name <- "Wells ordered by column"
              message("... wells are ordered by column...")
            }
          }
          if (export_PDF) {
            pdf(file = file.path(plots_path, paste0(file_name, ".pdf")),
                width = plot_width, height = plot_height
                )
          }
          for (i in seq_along(column_file_names)) {
            use_column <- names(column_file_names)[[i]]
            if (!(export_PDF)) {
              PNG_name <- paste0(file_name, " - ", i, ") ",
                                 sub("_rep1", "", use_column, fixed = TRUE),
                                 ".png"
                                 )
              png(filename = file.path(plots_path, "PNGs", sub_folder_name, PNG_name),
                  width = plot_width, height = plot_height, res = 600, units = "in"
                  )
            }
            PlateWellPlot(input_df,
                          use_column,
                          order_by_column = !(by_row),
                          aggregate_wells = take_means,
                          emphasize_NT    = emphasize_NT
                          )
            if (!(export_PDF)) {
              dev.off()
            }
          }
          if (export_PDF) {
            dev.off()
          }
        }
      }
    }
  }
  return(invisible(NULL))
}



