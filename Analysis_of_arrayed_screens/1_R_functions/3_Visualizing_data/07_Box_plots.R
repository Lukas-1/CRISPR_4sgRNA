# 2022-02-01


# Load packages and source code -------------------------------------------

library("RColorBrewer")
library("beeswarm")
library("vioplot")



# Define functions --------------------------------------------------------

BeeBox <- function(numeric_vec,
                   groups_fac,
                   group_colors,
                   use_swarm           = TRUE,
                   use_y_limits        = NULL,
                   y_axis_label        = "",
                   group_positions     = NULL,
                   group_labels        = NULL,
                   group_labels_cex    = 0.9,
                   group_labels_line   = 0.5,
                   group_label_lheight = 0.9,
                   draw_group_labels   = TRUE,
                   use_spacing         = 0.5,
                   point_cex           = 0.8,
                   horiz_lines         = NULL,
                   indicate_zero       = TRUE,
                   indicate_n          = TRUE,
                   points_alpha        = 0.4,
                   opaque_box          = NULL,
                   violin_wex          = 0.9
                   ) {

  assign("delete_numeric_vec", numeric_vec, envir = globalenv())
  assign("delete_groups_fac", groups_fac, envir = globalenv())
  assign("delete_group_colors", group_colors, envir = globalenv())

  stopifnot(length(numeric_vec) == length(groups_fac))

  if (is.null(group_labels)) {
    group_labels <- levels(groups_fac)
  }
  num_groups <- nlevels(groups_fac)
  if (is.null(group_positions)) {
    group_positions <- seq_len(num_groups)
  }

  numeric_list <- split(numeric_vec, groups_fac)

  ## Determine group positions
  side_gap <- 0.5
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  if (is.null(use_y_limits)) {
    use_y_limits <- range(numeric_vec)
  }
  y_space <- (use_y_limits[[2]] - use_y_limits[[1]]) * 0.02
  use_y_limits <- c(use_y_limits[[1]] - y_space, use_y_limits[[2]] + y_space)
  if (use_y_limits[[1]] > 0) {
    use_y_limits[[1]] <- 0
  }

  ## Set up the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = use_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  if (indicate_zero && (use_y_limits[[1]] < 0)) {
    abline(h = 0, col = "gray75", lty = "dotted")
  }
  if (!(is.null(horiz_lines))) {
    for (horiz_line in horiz_lines) {
      abline(h = horiz_line, col = "gray75", lty = "dotted")
    }
  }

  ## Prepare the plot colors
  light_colors  <- vapply(group_colors, function(x) brewer.pal(9, x)[[2]], "")
  dark_colors   <- vapply(group_colors, function(x) brewer.pal(9, x)[[9]], "")
  violin_colors <- vapply(group_colors, function(x) colorRampPalette(brewer.pal(9, x)[3:4])(3)[[2]], "")
  point_colors  <- vapply(group_colors, function(x) brewer.pal(9, x)[[7]], "")

  ## Draw the violin plots (in the background)
  vioplot(split(numeric_vec, groups_fac),
          at       = group_positions,
          pchMed   = NA,
          drawRect = FALSE,
          col      = violin_colors,
          border   = NA,
          wex      = violin_wex,
          add      = TRUE,
          axes     = FALSE
          )

  ## Draw the x axis line
  box(bty = "l")

  ## Draw the jittered points
  set.seed(1)
  if (use_swarm) {
    beeswarm_df <- beeswarm(numeric_list, at = group_positions, do.plot = FALSE,
                            cex = point_cex, spacing = use_spacing,
                            priority = "random"
                            )
    x_vec <- beeswarm_df[, "x"]
    x_pos <- rep(group_positions, lengths(numeric_list))
    x_deviations_list <- split(x_pos - x_vec,
                               match(beeswarm_df[, "x.orig"], levels(groups_fac))
                               )
    if (is.null(opaque_box)) {
      are_opaque <- vapply(x_deviations_list, function(x) any(abs(x) > 0.11), logical(1))
    } else {
      are_opaque <- opaque_box
    }

    light_colors <- ifelse(are_opaque,
                           adjustcolor(light_colors, alpha.f = 0.7),
                           adjustcolor(light_colors, alpha.f = 0.4)
                           )
  } else {
    x_vec <- rep(group_positions, lengths(numeric_list))
    x_vec <- x_vec + rnorm(n = length(x_vec), mean = 0, sd = 0.06)
  }
  points(x   = x_vec,
         y   = unlist(numeric_list),
         cex = point_cex,
         pch = 16,
         col = rep(adjustcolor(point_colors, alpha.f = points_alpha), lengths(numeric_list)),
         xpd = NA
         )

  ## Draw the superimposed box plots
  boxplot(numeric_list,
          at        = group_positions,
          boxwex    = 0.3,
          outline   = FALSE,
          names     = rep.int("", length(group_positions)),
          whisklty  = "blank",
          staplewex = 0,
          whisklwd  = 0,
          staplelty = 0,
          medlwd    = par("lwd") * 2,
          col       = light_colors,
          border    = dark_colors,
          add       = TRUE,
          axes      = FALSE,
          lwd       = par("lwd") * 1.5
          )

  ## Draw the y axis and x and y axis labels
  AbbreviateDataAxis(2)
  mtext(y_axis_label, side = 2, line = 2.8, cex = par("cex"))

  labels_splits <- strsplit(group_labels, " ", fixed = TRUE)
  labels_top <- sapply(labels_splits, "[", 1)
  labels_bottom <- sapply(labels_splits, "[", 2)

  if (indicate_n) {
    group_labels_line <- group_labels_line + 0.6
    mtext(lengths(numeric_list),
          at   = group_positions,
          side = 1,
          line = -0.3,
          cex  = 0.5,
          font = 2,
          col  = "gray80"
          )
  }

  if (draw_group_labels) {
    mtext(sapply(labels_top, VerticalAdjust),
          at   = group_positions,
          side = 1,
          line = group_labels_line,
          cex  = par("cex") * group_labels_cex
          )

    mtext(ifelse(is.na(labels_bottom),
                 NA,
                 sapply(labels_bottom, VerticalAdjust)
                 ),
          at   = group_positions,
          side = 1,
          line = group_labels_line + group_label_lheight,
          cex  = par("cex") * group_labels_cex
          )
  }
  return(invisible(NULL))
}



PositionsForSubgroups <- function(subgroup_levels_list, space_ratio = 0.7) {
  num_subgroups <- length(unlist(subgroup_levels_list, use.names = FALSE))
  total_space <- num_subgroups - 1L
  large_space_after <- unlist(lapply(subgroup_levels_list, function(x) {
    if (length(x) == 1) {
      TRUE
    } else {
      c(rep(FALSE, length(x) - 1), TRUE)
    }
  }))
  are_large_space <- large_space_after[-length(large_space_after)]
  large_space_size <- total_space / (sum(are_large_space) + (sum(!(are_large_space)) * space_ratio))
  positions_vec <- rep(NA, num_subgroups)
  positions_vec[c(1, length(positions_vec))] <- c(1, length(positions_vec))
  for (i in seq_len(length(positions_vec) - 2) + 1) {
    new_pos <- positions_vec[[i - 1]]
    if (large_space_after[[i - 1]]) {
      new_pos <- new_pos + large_space_size
    } else {
      new_pos <- new_pos + large_space_size * space_ratio
    }
    positions_vec[[i]] <- new_pos
  }
  return(positions_vec)
}



BeeBoxPlates <- function(input_df,
                         use_column,
                         show_subgroups = FALSE,
                         split_NT       = FALSE,
                         compare_group  = NULL,
                         plate_number   = NULL,
                         show_96wp      = FALSE,
                         common_scale   = TRUE,
                         point_cex      = 0.7,
                         indicate_n     = TRUE
                         ) {

  if (show_subgroups && split_NT) {
    stop("Please choose either 'show_subgroups' or 'split_NT' to be TRUE!")
  }

  if (is.null(plate_number)) {
    if (show_96wp) {
      use_title <- "All 96-well source plates"
    } else {
      use_title <- "All plates"
    }
  } else {
    if (show_96wp) {
      use_title <- paste0("96-well plate #", plate_number)
    } else {
      use_title <- paste0("Plate ", plate_number)
    }
  }

  ## Summarize replicates, if appropriate
  has_replicates <- grepl("_rep", use_column, fixed = TRUE)
  if (has_replicates) {
    if (!(grepl("_rep1", use_column, fixed = TRUE))) {
      stop("Please supply the column name for replicate 1! Both replicates will be displayed.")
    }
    rep2_column <- sub("_rep1", "_rep2", use_column, fixed = TRUE)
    numeric_vec <- rowMeans(input_df[, c(use_column, rep2_column)])
  } else {
    numeric_vec <- input_df[, use_column]
  }

  ## Determine the y axis label
  y_axis_label <- FormatPlotMath(short_column_labels[[use_column]])

  ## Assign each well to a group
  mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
  groups_vec <- ifelse(input_df[, "Is_NT_ctrl"],
                       "NT control",
                       ifelse(input_df[, "Is_pos_ctrl"],
                              "Pos. control",
                              ifelse(input_df[, "Well_number_384"] %in% mat_384[, c(1, 24)],
                                     "Empty well",
                                     ifelse(!(is.na(input_df[, "Entrez_ID"])),
                                            "Gene",
                                            ifelse(input_df[, "Target_flag"] %in% "mCherry",
                                                   "mCherry",
                                                   "No virus"
                                                   )
                                            )
                                     )
                              )
                       )

  ## Exclude wells that shouldn't be used for calculating the y axis limits
  are_to_exclude <- (groups_vec == "Empty well")
  if ("Is_problematic" %in% names(input_df)) {
    are_to_exclude <- are_to_exclude | input_df[, "Is_problematic"]
  }
  if (max(input_df[, "CellTiterGlo_raw"]) > 10^7) {
    glo_cutoff <- 10^6
  } else {
    glo_cutoff <- 20000
  }
  were_killed <- (as.integer(as.roman(input_df[, "Plate_number_384"])) >= 10) &
                 (input_df[, "Well_number_384"] %in% mat_384[c(1, 16), ]) &
                 (input_df[, "CellTiterGlo_raw"] < glo_cutoff)
  are_to_exclude <- are_to_exclude | were_killed
  use_y_limits <- range(numeric_vec[!(are_to_exclude)])


  ## Exclude additional wells
  are_to_exclude <- are_to_exclude |
                    ((groups_vec == "No virus") & (!(is.na(input_df[, "Target_ID"]))))
  if (!(is.null(plate_number))) {
    if ((!(show_96wp)) && (!(is.numeric(plate_number)))) {
      plate_number <- as.integer(as.roman(plate_number))
    }
    if (show_96wp) {
      plate_numbers_vec <- input_df[, "Plate_number_96"]
    } else {
      if (!(is.numeric(plate_number))) {
        plate_number <- as.integer(as.roman(plate_number))
      }
      plate_numbers_vec <- as.integer(as.roman(input_df[, "Plate_number_384"]))
    }
    are_to_exclude <- are_to_exclude | (!(plate_numbers_vec %in% plate_number))
  } else if (show_96wp) {
    are_to_exclude <- are_to_exclude | is.na(input_df[, "Plate_number_96"])
  }

  if (show_96wp) {
    well_numbers_vec <- input_df[, "Well_number_96"]
  } else {
    well_numbers_vec <- input_df[, "Well_number_384"]
  }
  if (split_NT) {
    are_to_exclude <- are_to_exclude | (groups_vec == "mCherry")
  }

  numeric_vec <- numeric_vec[!(are_to_exclude)]
  groups_vec <- groups_vec[!(are_to_exclude)]
  well_numbers_vec <- well_numbers_vec[!(are_to_exclude)]
  if (length(numeric_vec) > 384) {
    use_spacing <- 0.09
  } else {
    use_spacing <- 0.15
  }

  ## Create a factor that groups wells by the target type

  if (split_NT) {
    group_colors <- c(
      "Gene"         = "Purples",
      "Tubingen NT"  = "Blues",
      "Own NT"       = "Blues",
      "Pos. control" = "Reds",
      "No virus"     = "Greys"
    )
    are_NT <- groups_vec == "NT control"
    wells_384_vec <- input_df[, "Well_number_384"][!(are_to_exclude)][are_NT]
    are_own_NT <- wells_384_vec %in% mat_384[, c(2, 22)]
    groups_vec[are_NT] <- ifelse(are_own_NT, "Own NT", "Tubingen NT")
  } else {
    group_colors <- c(
      "Gene"         = "Purples",
      "NT control"   = "Blues",
      "Pos. control" = "Reds",
      "No virus"     = "Greys",
      "mCherry"      = "Oranges"
    )

  }
  groups_fac <- droplevels(factor(groups_vec, levels = names(group_colors)))

  ## Prepare plot parameters
  if (grepl("PercActivation|[Ff]oldNT|p_value", use_column)) {
    horiz_lines <- 1
  } else {
    horiz_lines <- NULL
  }
  old_par <- par(mar = c(3.8, 4.1, 3.4, 2))


  if (!(is.null(compare_group))) {
    wells_384_vec <- input_df[, "Well_number_384"][!(are_to_exclude)]
    are_own_NT <- wells_384_vec %in% mat_384[, c(2, 22)]
    if (compare_group == "Genes and NT") {
      if (show_96wp) {
        stop("When the value of the 'compare_group' argument is 'Genes and NT'",
             ", the 'show_96wp' argument is not supported!"
             )
      }
      are_this_group <- (groups_vec == "Gene") |
                        ((groups_vec == "NT control") & !(are_own_NT))
      plate_numbers_fac <- factor(input_df[, "Plate_number_384"][!(are_to_exclude)][are_this_group],
                                  levels = as.character(as.roman(1:12))
                                  )
      these_groups_fac <- factor(groups_vec[are_this_group])
      plate_groups_fac <- interaction(plate_numbers_fac, these_groups_fac,
                                      sep = " / ", lex.order = TRUE, drop = TRUE
                                      )
      subgroup_levels_list <- lapply(split(as.character(these_groups_fac), plate_numbers_fac), unique)
      subgroup_labels <- unlist(subgroup_levels_list, use.names = FALSE)
      use_level_colors <- group_colors[levels(droplevels(groups_fac[are_this_group]))][subgroup_labels]
      subgroup_labels <- c("Gene"       = "G",
                           "NT control" = "NT"
                           )[subgroup_labels]
      positions_vec <- PositionsForSubgroups(subgroup_levels_list, space_ratio = 0.55)
      group_labels_cex <- 0.7
      positions_list <- split(positions_vec,
                              rep(seq_along(subgroup_levels_list), lengths(subgroup_levels_list))
                              )
      group_labels_line <- 0.2
      violin_wex <- 0.7
      assign("delete_positions_list", positions_list, envir = globalenv())
    } else {
      if (compare_group == "Own NT control") {
        compare_group <- "NT control"
        are_this_group <- (groups_vec == compare_group) & are_own_NT
      } else {
        are_this_group <- groups_vec == compare_group
      }
      if (show_96wp) {
        group_labels_cex <- 0.75
        group_labels_line <- 0.4
        plate_groups_fac <- factor(input_df[, "Plate_number_96"][!(are_to_exclude)][are_this_group])
      } else {
        group_labels_cex <- 0.9
        group_labels_line <- 0.5
        plate_groups_fac <- factor(input_df[, "Plate_number_384"][!(are_to_exclude)][are_this_group],
                                   levels = as.character(as.roman(1:12))
                                   )
      }
      use_level_colors <- rep(group_colors[[compare_group]], nlevels(plate_groups_fac))
      subgroup_labels <- levels(plate_groups_fac)
      positions_vec <- NULL
      violin_wex <- 0.9
    }
    numeric_vec <- numeric_vec[are_this_group]

    titles_vec <- c(
      "Own NT control" = "In-house non-targeting controls",
      "Genes and NT"   = "Genes and interspersed NT controls",
      "Gene"           = "Genes",
      "NT control"     = "Non-targeting controls",
      "Pos. control"   = "Positive control"
    )
    if (compare_group %in% names(titles_vec)) {
      use_title <- titles_vec[[compare_group]]
    } else {
      use_title <- compare_group
    }
    if (show_96wp) {
      use_title <- paste(use_title, "across 96-well plates")
    } else {
      use_title <- paste(use_title, "across plates")
    }
    assign("delete_are_this_group", are_this_group, envir = globalenv())
    assign("delete_numeric_vec", numeric_vec, envir = globalenv())
    assign("delete_plate_numbers_fac", plate_groups_fac, envir = globalenv())

    BeeBox(numeric_vec, plate_groups_fac, use_level_colors,
           group_labels = subgroup_labels,
           y_axis_label = y_axis_label,
           use_spacing = use_spacing, use_y_limits = use_y_limits,
           point_cex = point_cex, horiz_lines = horiz_lines,
           indicate_n = indicate_n,
           group_labels_cex = group_labels_cex,
           group_positions = positions_vec, opaque_box = TRUE,
           group_labels_line = group_labels_line,
           violin_wex = violin_wex
           )
    if (compare_group == "Genes and NT") {
      mtext(levels(plate_numbers_fac),
            side = 1,
            at   = vapply(positions_list, mean, numeric(1)),
            line = 1.2 + if (indicate_n) 0.6 else 0,
            cex  = par("cex") * 0.8
            )
    } else {
      if (show_96wp) {
        x_axis_label <- "96-well plate number"
      } else {
        x_axis_label <- "Plate number"
      }
      mtext(x_axis_label,
            side = 1,
            line = (if (show_96wp) 1.55 else 1.65) + (if (indicate_n) 0.6 else 0),
            cex  = par("cex")
            )
    }
  } else if (show_subgroups) { # If show_subgroups is TRUE, wells are sub-divided by their physical location

    ## Split the control wells into subgroups
    subgroups_fac <- rep(NA, length(groups_fac))
    are_NT <- groups_fac == "NT control"
    subgroups_fac[are_NT] <- ifelse(well_numbers_vec[are_NT] %in% mat_384[, 2],
                                    "L",
                                    ifelse(well_numbers_vec[are_NT] %in% mat_384[, 22],
                                           "R", if (show_96wp) " " else "Mid"
                                           )
                                    )

    are_pos <- groups_fac == "Pos. control"
    subgroups_fac[are_pos] <- ifelse(well_numbers_vec[are_pos] %in% mat_384[, 4], "L", "R")

    ## Split the gene wells into subgroups
    if (show_96wp) {
      are_to_divide <- groups_fac == "Gene"
      mat_96 <- matrix(seq_len(96), nrow = 8, ncol = 12, byrow = TRUE)
      subgroups_fac[are_to_divide] <- ifelse(well_numbers_vec[are_to_divide] %in% mat_96[1:4, 1:6],
                                             "TL",
                                             ifelse(well_numbers_vec[are_to_divide] %in% mat_96[5:8, 1:6],
                                                    "BL",
                                                    ifelse(well_numbers_vec[are_to_divide] %in% mat_96[1:4, 7:12],
                                                           "TR",
                                                           "BR"
                                                           )
                                                    )
                                             )
    } else {
      are_to_divide <- !(are_pos | are_NT | groups_fac == "mCherry")
      subgroups_fac[are_to_divide] <- ifelse(well_numbers_vec[are_to_divide] %in% mat_384[1:8, 1:12],
                                             "TL",
                                             ifelse(well_numbers_vec[are_to_divide] %in% mat_384[9:16, 1:12],
                                                    "BL",
                                                    ifelse(well_numbers_vec[are_to_divide] %in% mat_384[1:8, 13:24],
                                                           "TR",
                                                           "BR"
                                                           )
                                                    )
                                             )
    }

    ## Create the final subgroups factor
    subgroups_fac[is.na(subgroups_fac)] <- " "
    subgroups_fac <- droplevels(factor(subgroups_fac, levels = c("TL", "TR", "BL", "BR", "L", "Mid", "R", " ")))
    stopifnot(!(anyNA(subgroups_fac)))
    interaction_fac <- interaction(groups_fac, subgroups_fac, sep = " / ",
                                   lex.order = TRUE, drop = TRUE
                                   )

    ## Assign labels and colors to the subgroups
    subgroup_levels_list <- lapply(split(as.character(subgroups_fac), groups_fac), unique)
    subgroup_levels_list <- lapply(subgroup_levels_list, function(x) x[order(match(x, levels(subgroups_fac)))])
    subgroup_labels <- unlist(subgroup_levels_list, use.names = FALSE)
    subgroup_colors <- rep(group_colors[levels(groups_fac)], lengths(subgroup_levels_list))

    ## Determine the x axis locations for the subgroups
    assign("delete_interaction_fac", interaction_fac, envir = globalenv())
    assign("delete_subgroup_levels_list", subgroup_levels_list, envir = globalenv())
    positions_vec <- PositionsForSubgroups(subgroup_levels_list)

    ## Create the bee/violin/box plot for the subgroups
    if (show_96wp) {
      group_labels_cex <- 0.9
    } else {
      group_labels_cex <- 0.8
    }
    BeeBox(numeric_vec, interaction_fac, subgroup_colors,
           y_axis_label = y_axis_label, group_labels = subgroup_labels,
           group_labels_cex = group_labels_cex, group_positions = positions_vec,
           group_label_lheight = 0.7, group_labels_line = 0.4,
           use_spacing = use_spacing, use_y_limits = use_y_limits,
           indicate_n = indicate_n, point_cex = point_cex,
           horiz_lines = horiz_lines
           )

    ## Create x axis labels for the super-groups
    positions_list <- split(positions_vec,
                            rep(seq_len(nlevels(groups_fac)), lengths(subgroup_levels_list))
                            )
    for (i in seq_along(positions_list)) {
      if (length(positions_list[[i]]) > 1) {
        segments(x0  = positions_list[[i]][[1]],
                 x1  = positions_list[[i]][[length(positions_list[[i]])]],
                 y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.45 + if (indicate_n) 0.6 else 0),
                                                         from = "lines", to = "user"
                                                         )),
                 xpd = NA,
                 col = "gray65"
                 )
      }
    }
    group_labels <- sub(" control", " ctrl", levels(groups_fac), fixed = TRUE)
    group_labels <- sub("Pos. ", "Pos ", group_labels, fixed = TRUE)
    mtext(group_labels,
          side = 1,
          at   = vapply(positions_list, mean, numeric(1)),
          line = 1.5 + if (indicate_n) 0.6 else 0,
          cex  = par("cex") * group_labels_cex
          )

  } else {
    ## Create the bee/violin/box plot for wells combined by target type
    if (split_NT && is.null(horiz_lines)) {
      horiz_lines <- median(numeric_vec[groups_vec == "Gene"])
    }
    BeeBox(numeric_vec, groups_fac,
           group_colors[levels(groups_fac)], y_axis_label = y_axis_label,
           use_spacing = use_spacing, use_y_limits = use_y_limits,
           point_cex = point_cex, horiz_lines = horiz_lines,
           indicate_n = indicate_n
           )
  }

  ## Draw the plot title
  title(use_title, cex.main = 1, line = 1.35)
  par(old_par)
  return(invisible(NULL))
}




ExportAllBoxPlots <- function(input_df, top_folder) {

  use_height <- 4
  narrow_width <- 5
  wide_width <- 6.2

  compare_groups <- c(
    "Gene", "Own NT control", "Pos. control", "Genes and NT"
  )

  for (export_PDF in c(TRUE, FALSE)) {

    if (export_PDF) {
      message("Exporting PDF files...")
    } else {
      message("Exporting PNG files...")
    }

    if (export_PDF) {
      figures_folder <- file.path(top_folder, "PDFs")
      dir.create(figures_folder, showWarnings = FALSE)
    } else {
      figures_folder <- file.path(top_folder, "PNGs")
      dir.create(figures_folder, showWarnings = FALSE)
    }

    for (show_subgroups in c(FALSE, TRUE)) {

      if (show_subgroups) {
        message("... wells are divided into subgroups by their physical location...")
      } else {
        message("... wells are combined according to their target type...")
      }

      if (show_subgroups) {
        use_width <- wide_width
      } else {
        use_width <- narrow_width
      }

      if (show_subgroups) {
        sub_folder <- file.path(figures_folder, "Subgroups")
      } else {
        sub_folder <- file.path(figures_folder, "Combined")
      }
      dir.create(sub_folder, showWarnings = FALSE)

      for (i in seq_along(column_file_names)) {
        use_column <- names(column_file_names)[[i]]
        folder_name <- paste0("Box plots - ", i, ") ", column_file_names[[i]])
        if (export_PDF) {
          pdf(file = file.path(sub_folder, paste0(folder_name, ".pdf")),
              width = use_width, height = use_height
              )
          BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups)
          if (show_subgroups) {
            for (compare_group in compare_groups) {
              BeeBoxPlates(input_df, use_column, compare_group = compare_group)
            }
            for (compare_group in c("Gene", "NT control")) {
              BeeBoxPlates(input_df, use_column, compare_group = compare_group,
                           show_96wp = TRUE
                           )
            }
          } else {
            BeeBoxPlates(input_df, use_column, split_NT = TRUE)
          }
          for (plate_number in as.character(as.roman(1:12))) {
            BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups,
                         plate_number = plate_number
                         )
          }
          BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups,
                       show_96wp = TRUE
                       )

          for (i in 1:22) {
            BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups,
                         plate_number = i, show_96wp = TRUE
                         )
          }
          dev.off()
        } else {
          sub_sub_folder <- file.path(sub_folder, folder_name)
          dir.create(sub_sub_folder, showWarnings = FALSE)
          png(filename = file.path(sub_sub_folder, "01) Box plot - all plates.png"),
              width = use_width, height = use_height, units = "in", res = 600
              )
          BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups)
          dev.off()

          if (show_subgroups) {
            for (i in seq_along(compare_groups)) {
              compare_group <- compare_groups[[i]]
              png(filename = file.path(sub_sub_folder,
                                       paste0(formatC(i + 1, flag = "0", width = 2),
                                              ") Box plot - ", tolower(compare_group),
                                              " across plates.png"
                                              )
                                       ),
                  width = use_width, height = use_height, units = "in", res = 600
                  )
              BeeBoxPlates(input_df, use_column, compare_group = compare_group)
              dev.off()
            }
            for (i in 1:2) {
              compare_group <- c("Gene", "NT control")[[i]]
              png(filename = file.path(sub_sub_folder,
                                       paste0(formatC(i + 1 + length(compare_groups), flag = "0", width = 2),
                                              ") Box plot - ", tolower(compare_group),
                                              " across 96wp.png"
                                              )
                                       ),
                  width = use_width, height = use_height, units = "in", res = 600
                  )
              BeeBoxPlates(input_df, use_column, compare_group = compare_group,
                           show_96wp = TRUE
                           )
              dev.off()
            }
            current_i <- 4 + length(compare_groups)
          } else {
            png(filename = file.path(sub_sub_folder,
                                     "02) Box plots - NT controls split.png"
                                     ),
                width = use_width, height = use_height, units = "in", res = 600
                )
            BeeBoxPlates(input_df, use_column, compare_group = compare_group,
                         show_96wp = TRUE
                         )
            dev.off()
            current_i <- 3
          }

          for (i in current_i:(current_i + 12 - 1)) {
            roman_number <- as.character(as.roman(i - current_i + 1))
            png(filename = file.path(sub_sub_folder,
                                     paste0(formatC(i, flag = "0", width = 2),
                                            ") Box plot - plate ", roman_number,
                                            ".png"
                                            )
                                     ),
                width = use_width, height = use_height, units = "in", res = 600
                )
            BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups,
                         plate_number = roman_number
                         )
            dev.off()
          }
          png(filename = file.path(sub_sub_folder,
                                   paste0(current_i + 12, ") Box plot - all 96-well plates.png")
                                   ),
              width = use_width, height = use_height, units = "in", res = 600
              )
          BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups,
                       show_96wp = TRUE
                       )
          dev.off()
          for (i in (current_i + 13):(current_i + 12 + 22)) {
            plate_number <- i - (current_i + 12)
            png(filename = file.path(sub_sub_folder,
                                     paste0(formatC(i, flag = "0", width = 2),
                                            ") Box plot - plate ", plate_number,
                                            ".png"
                                            )
                                     ),
                width = use_width, height = use_height, units = "in", res = 600
                )
            BeeBoxPlates(input_df, use_column, show_subgroups = show_subgroups,
                         plate_number = plate_number, show_96wp = TRUE
                         )
            dev.off()
          }
        }
      }
    }
  }
  return(invisible(NULL))
}





