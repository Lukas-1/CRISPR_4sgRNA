### 12th March 2022 ###



# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-03-11 - HNRNPK pooled screen")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Read in data, and rank hit genes.RData"))
load(file.path(R_objects_directory, "03) Find non-essential hit genes.RData"))
load(file.path(R_objects_directory, "06) Rank hit genes.RData"))





# Define global variables -------------------------------------------------

use_margin <-  c(4.25, 4, 3.5, 8.25)





# Define functions --------------------------------------------------------

TruncateValues <- function(numeric_vec, max_value) {
  ifelse(numeric_vec < 0,
         ifelse(numeric_vec < -(max_value), -(max_value), numeric_vec),
         ifelse(numeric_vec > max_value, max_value, numeric_vec)
         )
}



PrettyAxis <- function(use_side, axis_limits) {
  tick_positions <- pretty(axis_limits, n = 5)
  axis(use_side,
       mgp    = c(3, 0.55, 0),
       tcl    = -0.35,
       las    = 1,
       at     = tick_positions,
       labels = ifelse(tick_positions %in% pretty(axis_limits),
                       tick_positions,
                       ""
                       ),
       gap.axis = 0.5
       )
  return(invisible(NULL))
}



HitStrengthScatter <- function(input_df_1,
                               input_df_2,
                               use_max_value = 55,
                               selection = "Non-targeting",
                               x_label = "Hit strength (dataset A)",
                               y_label = "Hit strength (dataset B)"
                               ) {

  stopifnot("use_margin" %in% ls(envir = globalenv()))

  ## Perform checks

  required_columns <- c("Hit_strength", "gene_id", "Gene_symbol")
  stopifnot(all(required_columns %in% names(input_df_1)))
  stopifnot(all(required_columns %in% names(input_df_2)))


  ## Prepare the numeric data

  matches_vec <- match(input_df_1[, "gene_id"], input_df_2[, "gene_id"])
  stopifnot(!(anyNA(matches_vec)))
  x_vec <- TruncateValues(input_df_1[, "Hit_strength"], max_value = use_max_value)
  y_vec <- TruncateValues(input_df_2[matches_vec, "Hit_strength"], max_value = use_max_value)


  ## Set up the plot

  old_par <- par(mar = use_margin)

  axis_limits <- c(-(use_max_value - 0.5), (use_max_value + 0.5))

  plot(1,
       type = "n",
       xlim = axis_limits,
       ylim = axis_limits,
       xlab = x_label,
       ylab = y_label,
       axes = FALSE
       )

  PrettyAxis(1, axis_limits)
  PrettyAxis(2, axis_limits)

  abline(h = 0, v = 0, col = "gray70")


  ## Add points

  AddPoints(x_vec, y_vec, input_df_1[, "Gene_symbol"], selection = selection)

  box()

  par(old_par)

  return(invisible(NULL))
}





VolcanoPlot <- function(input_df, main_title, selection = "Non-targeting") {

  stopifnot("use_margin" %in% ls(envir = globalenv()))

  ## Perform checks

  required_columns <- c("log2 Ratio", "pValue", "Gene_symbol")
  stopifnot(all(required_columns %in% names(input_df)))


  ## Prepare the numeric data

  x_vec <- input_df[, "log2 Ratio"]
  y_vec <- -log10(input_df[, "pValue"])


  ## Set up the plot

  old_par <- par(mar = use_margin)

  plot(1,
       type = "n",
       xlim = c(-5.9, 5.9),
       ylim = c(0, 100),
       xlab = expression("log" * scriptscriptstyle(" ") * "2" *
                         scriptstyle(" ") * "(" * scriptscriptstyle(" ") *
                         "fold change" * scriptscriptstyle(" ") * ")"
                         ),
       ylab = expression("" - "log" * scriptscriptstyle(" ") * "10" *
                         scriptstyle(" ") * scriptscriptstyle(" ") * "(" *
                         italic("p") * " value" * scriptscriptstyle(" ") * ")"
                         ),
       mgp  = c(2.6, 0.55, 0),
       tcl  = -0.35,
       las  = 1
       )

  title(main_title, line = 1.1, cex.main = 1.1 * par("cex"))
  abline(h = 0, v = 0, col = "gray70")


  ## Add points

  AddPoints(x_vec, TruncateValues(y_vec, max_value = 100),
            input_df[, "Gene_symbol"], selection = selection
            )

  box()

  par(old_par)

  return(invisible(NULL))
}




AddPoints <- function(x_vec, y_vec, symbols_vec, selection, random_order = TRUE) {

  stopifnot(identical(length(x_vec), length(y_vec)))
  stopifnot(identical(length(x_vec), length(symbols_vec)))


  ## Set up the graphical parameters

  points_alpha <- 0.4
  bg_alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

  points_alpha <- 0.7
  fg_alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

  use_random_order <- FALSE
  legend_y <- 0.5

  bg_point_size <- 0.8
  fg_point_size <- 0.8

  border_colors_map <- NULL

  if (selection == "Essential") {

    legend_y <- 0.558

    points_alpha <- 0.4
    bg_alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

    points_alpha <- 0.4
    fg_alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

    gene_sets <- list(
      "Essential" = essential_df[essential_df[, "Category"] %in% "Essential", "Gene_symbol"]
    )
    fill_colors_map <- c(
      "Essential" = brewer.pal(9, "Reds")[[6]]
    )
    labels_list <- list(
      "Essential" = c("Essential genes", "in humans", "(n = 1833)")
    )

    use_random_order <- TRUE


  } else if (selection == "Non-targeting") {

    legend_y <- 0.52

    gene_sets <- list(
      "Non-targeting" = grep("non-targeting", symbols_vec, fixed = TRUE, value = TRUE)
    )
    fill_colors_map <- c(
      "Non-targeting" = brewer.pal(9, "Blues")[[5]]
    )
    labels_list <- list(
      "Non-targeting" = c("Non-targeting", "sgRNAs")
    )

  } else if (selection == "Top selected") {

    fg_point_size <- 1.1

    gene_sets <- list(
      "Top hits"          = Stefano_top_genes,
      "Top non-essential" = Stefano_non_essential_genes
    )

    fill_colors_map <- c(
      "Top hits"          = brewer.pal(9, "Reds")[[6]],
      "Top non-essential" = brewer.pal(9, "Blues")[[6]]
    )

    border_colors_map <- c(
      "Top hits"          = brewer.pal(9, "Reds")[[8]],
      "Top non-essential" = brewer.pal(9, "Blues")[[8]]
    )

    labels_list <- list(
      "Top hits"          = c("Top genes", "(overall, by ", "hit strength)"),
      "Top non-essential" = c("Top genes", "(filtered by", "essentiality)")
    )

  } else if (selection == "Top 2") {

    fg_point_size <- 1.2

    gene_sets <- list(
      "GTF3C6"  = "GTF3C6",
      "PYROXD1" = "PYROXD1"
    )

    fill_colors_map <- c(
      "GTF3C6" = brewer.pal(9, "Oranges")[[6]],
      "PYROXD1" = brewer.pal(9, "Blues")[[5]]
    )

    border_colors_map <- c(
      "GTF3C6" = brewer.pal(9, "Oranges")[[8]],
      "PYROXD1" = brewer.pal(9, "Blues")[[8]]
    )

    labels_list <- list(
      "GTF3C6" = expression(bolditalic("GTF3C6")),
      "PYROXD1" = expression(bolditalic("PYROXD1"))
    )

  }

  if (is.null(border_colors_map)) {
    border_colors_map <- fill_colors_map
  }

  if (use_random_order) {
    fill_colors_vec <- rep(brewer.pal(9, "Greys")[[7]], length(x_vec))
    border_colors_vec <- fill_colors_vec
    for (gene_set in names(gene_sets)) {
      are_this_set <- symbols_vec %in% gene_sets[[gene_set]]
      fill_colors_vec[are_this_set] <- fill_colors_map[[gene_set]]
      border_colors_vec[are_this_set] <- border_colors_map[[gene_set]]
    }
    set.seed(1)
    use_order <- sample(seq_along(x_vec))
    x_vec <- x_vec[use_order]
    y_vec <- y_vec[use_order]
    fill_colors_vec <- fill_colors_vec[use_order]
    border_colors_vec <- border_colors_vec[use_order]
    points(x_vec,
           y_vec,
           pch  = 16,
           cex  = bg_point_size,
           col  = paste0(fill_colors_vec, bg_alpha_hex)
           )

  } else {
    are_not_highlighted <- !(symbols_vec %in% unlist(gene_sets))
    points(x_vec[are_not_highlighted],
           y_vec[are_not_highlighted],
           pch  = 16,
           cex  = bg_point_size,
           col  = paste0(brewer.pal(9, "Greys")[[7]], bg_alpha_hex),
           )
    for (gene_set in names(gene_sets)) {
      are_this_set <- symbols_vec %in% gene_sets[[gene_set]]
      points(x_vec[are_this_set],
             y_vec[are_this_set],
             pch  = 21,
             cex  = fg_point_size,
             col  = paste0(border_colors_map[[gene_set]], fg_alpha_hex),
             bg   = paste0(fill_colors_map[[gene_set]], fg_alpha_hex),
             )
    }
  }

  DrawSideLegend(labels_list    = labels_list,
                 fill_colors    = fill_colors_map,
                 border_colors  = border_colors_map,
                 use_point_size = max(1.1, fg_point_size)
                 )

  return(invisible(NULL))
}



ConcatenateExpressions <- function(expression_list, my_sep = "  \u2013  ") {
  literal_strings <- vapply(expression_list, StripExpression, "")
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}


VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("gh")), use_expression, expression(phantom("gh")))
  return(ConcatenateExpressions(my_list, my_sep = ""))
}


StripExpression <- function(my_expression) {
  if (is.character(my_expression)) {
    literal_string <- paste0("\"", capture.output(cat(my_expression)), "\"")
  } else {
    literal_string <- paste0(capture.output(my_expression), collapse = "")
    if (substr(literal_string, 1L, 11L) == "expression(") {
      literal_string <- substr(literal_string, 12L, nchar(literal_string) - 1L)
    }
  }
  return(literal_string)
}


Embolden <- function(my_string) {
  my_string <- StripExpression(my_string)
  my_string <- sub("italic(", "bolditalic(", my_string, fixed = TRUE)
  parse(text = paste0("bold(", my_string, ")"))
}




DrawSideLegend <- function(labels_list,
                           fill_colors,
                           border_colors,
                           use_pch = 21,
                           use_point_size = 1
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(fill_colors)))
  stopifnot(identical(length(labels_list), length(border_colors)))

  ## Prepare for drawing the legend
  y_mid <- 0.5
  small_gap <- diff(grconvertY(c(0, 1.25), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * 1.75

  if (all(lengths(labels_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(labels_list))
    are_first <- rep(TRUE, length(labels_list))
  } else {
    are_first <- unlist(lapply(labels_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  x_text  <- 1 + diff(grconvertX(c(0, 0.75), from = "lines", to = "npc"))
  x_point <- 1 + diff(grconvertX(c(0, 0.9), from = "lines", to = "npc"))

  ## Draw the legend
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(labels_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )

  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))

  points(x   = rep(grconvertX(x = x_point, from = "npc", to = "user"), length(labels_list)),
         y   = tapply(y_pos, groups_vec, mean),
         cex = use_point_size,
         pch = use_pch,
         bg  = fill_colors,
         col = border_colors,
         xpd = NA
         )

  return(invisible(NULL))
}



AllHitStrengthScatters <- function(input_df_1,
                                   input_df_2,
                                   axis_label_1,
                                   axis_label_2,
                                   folder_name
                                   ) {

  folder_path <- file.path(file_output_directory, "Plots", folder_name)
  dir.create(folder_path, showWarnings = FALSE)

  use_width  <- 3 + (sum(use_margin[c(2, 4)]) / 5)
  use_height <- 3 + (sum(use_margin[c(1, 3)]) / 5)

  selection_names <- c("Non-targeting", "Essential", "Top selected", "Top 2")

  for (use_device in c("none", "png", "pdf")) {

    if (use_device == "pdf") {
      pdf(file.path(folder_path, paste0(folder_name, ".pdf")),
          width = use_width, height = use_height
          )
    }

    for (i in seq_along(selection_names)) {
      if (use_device == "png") {
        file_name <- paste0(folder_name, " - ", letters[[i]], ") ", selection_names[[i]])
        png(filename = file.path(folder_path, paste0(file_name, ".png")),
            width    = use_width,
            height   = use_height,
            units    = "in",
            res      = 600
            )
      }
      HitStrengthScatter(input_df_1, input_df_2,
                         x_label   = axis_label_1,
                         y_label   = axis_label_2,
                         selection = selection_names[[i]]
                         )
      if (use_device == "png") {
        dev.off()
      }
    }

    if (use_device == "pdf") {
      dev.off()
    }
  }

  return(invisible(NULL))
}




AllVolcanoPlots <- function(input_df, folder_name, main_title) {

  folder_path <- file.path(file_output_directory, "Plots", folder_name)
  dir.create(folder_path, showWarnings = FALSE)

  use_width  <- 3 + (sum(use_margin[c(2, 4)]) / 5)
  use_height <- 3 + (sum(use_margin[c(1, 3)]) / 5)

  selection_names <- c("Non-targeting", "Essential", "Top selected", "Top 2")

  for (use_device in c("none", "pdf", "png")) {

    if (use_device == "pdf") {
      pdf(file.path(folder_path, paste0(folder_name, ".pdf")),
          width = use_width, height = use_height
          )
    }

    for (i in seq_along(selection_names)) {
      if (use_device == "png") {
        file_name <- paste0(folder_name, " - ", letters[[i]], ") ", selection_names[[i]])
        png(filename = file.path(folder_path, paste0(file_name, ".png")),
            width    = use_width,
            height   = use_height,
            units    = "in",
            res      = 600
            )
      }
      VolcanoPlot(input_df, main_title, selection = selection_names[[i]])
      if (use_device == "png") {
        dev.off()
      }
    }

    if (use_device == "pdf") {
      dev.off()
    }
  }

  return(invisible(NULL))
}





HitStrengthBeeBox <- function(data_df,
                              main_title = "",
                              symbols_list = NULL,
                              super_groups_brewers = NULL,
                              small_label_symbols = NULL,
                              log2fc = FALSE
                              ) {

  stopifnot(!(anyNA(data_df[, "Gene_symbol"])))

  if (is.null(symbols_list)) {
    symbols_list <- list("Top 5 overall"       = top_10_HS[1:5],
                         "Top 5 non-essential" = top_10_NE[1:5]
                         )
  }

  if (is.null(super_groups_brewers)) {
    super_groups_brewers = c("Greens", "Reds", "Blues", "Purples")
  }

  if (is.null(small_label_symbols)) {
    small_label_symbols <- c()
  }

  ## Prepare the data
  NT_string <- "non-targeting"
  are_NT <- grepl(NT_string, data_df[, "Gene_symbol"], fixed = TRUE)


  groups_vec <- ifelse(are_NT, NT_string, data_df[, "Gene_symbol"])

  all_symbols <- c(NT_string, unlist(symbols_list))

  are_included <- groups_vec %in% all_symbols

  use_indices <- which(are_included)
  new_order <- order(match(groups_vec[use_indices], all_symbols))
  use_indices <- use_indices[new_order]
  groups_vec <- groups_vec[use_indices]
  groups_fac <- factor(groups_vec, levels = all_symbols)
  if (log2fc) {
    use_column <- "log2 Ratio"
  } else {
    use_column <- "Hit_strength"
  }
  numeric_vec <- data_df[use_indices, use_column]

  super_groups <- c(1L, rep(seq_along(symbols_list) + 1L, lengths(symbols_list)))
  groups_list <- split(levels(groups_fac), super_groups)
  super_groups_brewers <- super_groups_brewers[seq_along(groups_list)]


  group_labels <- c("NT ctrls", sapply(unlist(symbols_list), function(x) as.expression(bquote(italic(.(x))))))
  many_genes <- length(all_symbols) >= 20
  if (!(many_genes)) {
    group_labels <- sapply(group_labels, function(x) VerticalAdjust(x))
  }

  ## Determine group positions
  small_gap <- 0.75
  medium_gap <- 1
  large_gap <- 1
  if (all(lengths(groups_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(groups_list))
    are_first <- rep(TRUE, length(groups_list))
  } else {
    are_first <- unlist(lapply(groups_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }), use.names = FALSE)
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_position <- 1
  group_positions <- start_position + cumsum(gaps_vec)

  width <- 0.6 #2 / 3
  final_width <- width * ((max(group_positions) - min(group_positions)) / (nlevels(groups_fac) - 1))
  side_gap <- final_width
  group_limits <- c(group_positions[[1]] - side_gap,
                    group_positions[[length(group_positions)]] + side_gap
                    )


  ## Prepare the data axis
  use_numeric_limits <- c(min(numeric_vec) * 1.02, max(numeric_vec) * 1.02)
  numeric_axis_pos <- pretty(use_numeric_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  numeric_axis_labels <- format(numeric_axis_pos)


  ## Draw the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       xlab = "",
       ylab = if (log2fc) expression("Log"["2"] * " fold change") else "Hit strength",
       mgp = c(2.6, 1, 0)
       )
  abline(h = 0, col = "gray50")

  title(main_title, line = 1.1, cex.main = 1.1 * par("cex"))

  if (many_genes) {
    text(x      = group_positions,
         y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.3), from = "lines", to = "user")),
         labels = group_labels,
         srt    = 45,
         xpd    = NA,
         cex    = 0.8,
         adj    = c(1, 1)
         )
    group_label_line <- 3.7
  } else {
    mtext(group_labels, at = group_positions, side = 1, line = 0.6, padj = 0.5,
          cex = ifelse(levels(groups_fac) %in% small_label_symbols,
                       par("cex") * 0.7,
                       par("cex") * 0.9
                       )
          )
    group_label_line <- 2.2
  }


  for (i in seq_along(symbols_list)) {
    are_this_group <- levels(groups_fac) %in% symbols_list[[i]]
    if (sum(are_this_group) >= 3) {
      mtext(VerticalAdjust(names(symbols_list)[[i]]),
            at   = mean(group_positions[are_this_group]),
            side = 1,
            line = group_label_line,
            cex  = par("cex") * 0.9
            )
      segments(x0  = min(group_positions[are_this_group]),
               x1  = max(group_positions[are_this_group]),
               y0  = par("usr")[[3]] - diff(grconvertY(c(0, group_label_line - (if (many_genes) 0.3 else 0.4)), from = "lines", to = "user")),
               col = "gray70",
               lwd = par("lwd"),
               xpd = NA
               )
    }
  }


  ## Draw the boxplots
  box_colors <- sapply(super_groups_brewers, function(x) brewer.pal(8, x)[[2]])
  box_margins <- sapply(super_groups_brewers, function(x) brewer.pal(8, x)[[8]])
  boxplot(numeric_vec ~ groups_fac,
          at        = group_positions,
          boxwex    = c(0.7, rep(0.5, nlevels(groups_fac) - 1)),
          outline   = FALSE,
          names     = rep("", nlevels(groups_fac)),
          whisklty  = "blank",
          staplewex = 0,
          axes      = FALSE,
          whisklwd  = 0,
          staplelty = 0,
          col       = rep(box_colors, lengths(groups_list)),
          border    = rep(box_margins, lengths(groups_list)),
          boxlwd    = c(par("lwd") * 2, rep(par("lwd"), nlevels(groups_fac) - 1)),
          medlwd    = par("lwd") * 2,
          add       = TRUE,
          lend      = "butt"
          )

  ## Draw the points
  are_NT <- groups_vec == "non-targeting"
  set.seed(1)
  jittered_vec  <- group_positions[[1]] + rnorm(n = sum(are_NT), mean = 0, sd = 0.055)

  points_alpha <- 0.3
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
  if (many_genes) {
    point_cex <- 0.7
  } else {
    point_cex <- 1
  }

  point_colors <- sapply(super_groups_brewers, function(x) colorRampPalette(brewer.pal(9, x)[c(7, 8)])(5)[[2]])
  points(jittered_vec,
         numeric_vec[are_NT],
         pch = 16,
         col = paste0(point_colors[[1]], alpha_hex),
         cex = point_cex,
         xpd = NA
         )

  beeswarm_df <- beeswarm(numeric_vec[!(are_NT)] ~ droplevels(groups_fac[!(are_NT)]),
                          at       = group_positions[-1],
                          col      = rep(point_colors[-1], lengths(groups_list)[-1]),
                          priority = "density",
                          do.plot  = FALSE,
                          cex      = point_cex
                          )

  points(beeswarm_df[, "x"],
         beeswarm_df[, "y"],
         pch = 16,
         col = beeswarm_df[, "col"],
         cex = point_cex,
         xpd = NA
         )

  axis(2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, 0.54, 0),
       gap.axis = 0,
       tcl      = -0.4,
       las      = 1,
       lwd      = par("lwd"),
       cex      = par("cex"),
       xpd      = NA
       )

  box(bty = "o", lwd = par("lwd"))
  return(invisible(NULL))
}





All3BeeBoxes <- function(add_to_file_name, symbols_list = NULL, brewers_vec = NULL, small_label_symbols = NULL) {

  half_prefix <- paste0(" box plot - ", add_to_file_name, " - ")

  for (log2fc in c(FALSE, TRUE)) {

    if (log2fc) {
      sub_path <- file.path(beebox_path, "Log2FC")
      prefix <- paste0("Log2FC", half_prefix)
    } else {
      sub_path <- file.path(beebox_path, "Hit strength")
      prefix <- paste0("Log2FC", half_prefix)
    }

    png(file.path(sub_path, paste0(prefix, "1) NT vs. HNRNPK", ".png")),
        width = beebox_width, height = beebox_height, units = "in", res = 600
        )
    HitStrengthBeeBox(NT_v_HNRNPK_df, "NT vs. HNRNPK",
                      symbols_list = symbols_list,
                      super_groups_brewers = brewers_vec,
                      small_label_symbols = small_label_symbols,
                      log2fc = log2fc
                      )
    dev.off()


    png(file.path(sub_path, paste0(prefix, "2) Baseline vs. NT", ".png")),
        width = beebox_width, height = beebox_height, units = "in", res = 600
        )
    HitStrengthBeeBox(base_v_NT_df, "Baseline vs. NT",
                      symbols_list = symbols_list,
                      super_groups_brewers = brewers_vec,
                      small_label_symbols = small_label_symbols,
                      log2fc = log2fc
                      )
    dev.off()


    png(file.path(sub_path, paste0(prefix, "3) Baseline vs. HNRNPK", ".png")),
        width = beebox_width, height = beebox_height, units = "in", res = 600
        )
    HitStrengthBeeBox(base_v_HNRNPK_df, "Baseline vs. HNRNPK",
                      symbols_list = symbols_list,
                      super_groups_brewers = brewers_vec,
                      small_label_symbols = small_label_symbols,
                      log2fc = log2fc
                      )
    dev.off()

  }

  return(invisible(NULL))
}




# Check that the gene names are consistent --------------------------------

matches_vec <- match(base_v_NT_df[, "gene_id"], base_v_HNRNPK_df[, "gene_id"])
stopifnot(!(anyNA(matches_vec)))
stopifnot(identical(base_v_NT_df[, "Gene_symbol"], base_v_HNRNPK_df[, "Gene_symbol"][matches_vec]))




# Re-order by sgRNA sequence ----------------------------------------------

ReorderBySequence <- function(input_df) {
  input_df <- input_df[order(input_df[, "gene_id"]), ]
  row.names(input_df) <- NULL
  return(input_df)
}

base_v_HNRNPK_df <- ReorderBySequence(base_v_HNRNPK_df)
NT_v_HNRNPK_df <- ReorderBySequence(NT_v_HNRNPK_df)
base_v_NT_df  <- ReorderBySequence(base_v_NT_df)

for (df_name in c("base_v_HNRNPK_df", "NT_v_HNRNPK_df", "base_v_NT_df")) {
  use_df <- ReorderBySequence(get(df_name))
  are_NT <- grepl("Control-", use_df[, "gene_id"], fixed = TRUE)
  use_df[, "Gene_symbol"][are_NT] <- "non-targeting"
  assign(df_name, use_df, envir = globalenv())
}






# Calculate hit strengths -------------------------------------------------

base_v_HNRNPK_df[["Hit_strength"]] <- base_v_HNRNPK_df[, "log2 Ratio"] * (-log10(base_v_HNRNPK_df[, "pValue"]))
NT_v_HNRNPK_df[["Hit_strength"]] <- NT_v_HNRNPK_df[, "log2 Ratio"] * (-log10(NT_v_HNRNPK_df[, "pValue"]))
base_v_NT_df[["Hit_strength"]]  <- base_v_NT_df[, "log2 Ratio"] * (-log10(base_v_NT_df[, "pValue"]))




# Count the number of essential genes -------------------------------------

are_nontargeting <- grepl("non-targeting", base_v_NT_df[, "Gene_symbol"], fixed = TRUE)
unique_symbols <- unique(base_v_NT_df[!(are_nontargeting), "Gene_symbol"])
matches_vec <- match(unique_symbols, essential_df[, "Gene_symbol"])

categories_vec <- as.character(essential_df[matches_vec, "Category"])
categories_vec[is.na(categories_vec)] <- "No data"
use_levels <- c(levels(essential_df[, "Category"]), "No data")
categories_fac <- factor(categories_vec, levels = use_levels)







# Generate all plots ------------------------------------------------------

AllVolcanoPlots(NT_v_HNRNPK_df,
                folder_name = "1) Volcano plots - NT control vs. HNRNPK",
                main_title = "NT vs. HNRNPK"
                )

AllVolcanoPlots(base_v_NT_df,
                folder_name = "2) Volcano plots - Baseline vs. NT",
                main_title = "Baseline vs. NT"
                )

AllVolcanoPlots(base_v_HNRNPK_df,
                folder_name = "3) Volcano plots - Baseline vs. HNRNPK",
                main_title = "Baseline vs. HNRNPK"
                )


AllHitStrengthScatters(base_v_NT_df, NT_v_HNRNPK_df,
                       axis_label_1 = "Hit strength \u2013 Baseline vs. NT",
                       axis_label_2 = "Hit strength \u2013 NT vs. HNRNPK",
                       "4) Hit strength - NT-HNRNPK against Baseline-NT"
                       )

AllHitStrengthScatters(base_v_HNRNPK_df, NT_v_HNRNPK_df,
                       axis_label_1 = "Hit strength \u2013 Baseline vs. HNRNPK",
                       axis_label_2 = "Hit strength \u2013 NT vs. HNRNPK",
                       "5) Hit strength - NT-HNRNPK against Baseline-HNRNPK"
                       )

AllHitStrengthScatters(base_v_HNRNPK_df, base_v_NT_df,
                       axis_label_1 = "Hit strength \u2013 Baseline vs. HNRNPK",
                       axis_label_2 = "Hit strength \u2013 Baseline vs. NT",
                       "6) Hit strength - Baseline-NT against Baseline-HNRNPK"
                       )





# Show individual example plots -------------------------------------------

HitStrengthScatter(base_v_NT_df, NT_v_HNRNPK_df,
                   x_label   = "Hit strength \u2013 Baseline vs. NT",
                   y_label   = "Hit strength \u2013 NT vs. HNRNPK",
                   selection = "Non-targeting"
                   )

VolcanoPlot(base_v_NT_df, "Baseline vs. NT", selection = "Non-targeting")
VolcanoPlot(base_v_NT_df, "Baseline vs. NT", selection = "Essential")
VolcanoPlot(base_v_NT_df, "Baseline vs. NT", selection = "Top selected")
VolcanoPlot(base_v_NT_df, "Baseline vs. NT", selection = "Top 2")





# Show hit strength box plots ---------------------------------------------

use_small_labels <- c("GTF3C6", "PYROXD1")

HitStrengthBeeBox(NT_v_HNRNPK_df, "NT vs. HNRNPK", small_label_symbols = use_small_labels)
HitStrengthBeeBox(base_v_NT_df,  "Baseline vs. NT", small_label_symbols = use_small_labels)
HitStrengthBeeBox(base_v_HNRNPK_df, "Baseline vs. HNRNPK", small_label_symbols = use_small_labels)



beebox_width <- 7.5
beebox_height <- 4.2
beebox_path <- file.path(file_output_directory, "Plots", "9) sgRNA box plots")


All3BeeBoxes("c) top genes", small_label_symbols = use_small_labels)


All3BeeBoxes("a) Stefano's top hits", small_label_symbols = use_small_labels,
             symbols_list = list("Top hits (including essential genes)" = Stefano_top_genes)
             )

All3BeeBoxes("b) Stefano's non-essential hit genes", small_label_symbols = use_small_labels,
             symbols_list = list("Non-essential hit genes" = Stefano_non_essential_genes)
             )




# Save data ---------------------------------------------------------------

save(list = c("base_v_HNRNPK_df", "NT_v_HNRNPK_df", "base_v_NT_df"),
     file = file.path(R_objects_directory, "07) Show differential enrichment or depletion.RData")
     )







