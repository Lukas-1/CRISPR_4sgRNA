### 2022-10-31


# Load packages and source code -------------------------------------------

library("devEMF")

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))



# Define paths ------------------------------------------------------------

first_rdata_dir   <- file.path(first_illumina_trial_dir, "03_R_objects")
off_2sg_rdata_dir <- file.path(experiments_directory,
                               "2022-06-21 - Illumina paired-end 2sg - correct reference",
                               "03_R_objects"
                               )
off_4sg_rdata_dir <- file.path(experiments_directory,
                               "2022-09-02 - Illumina 4sg sequencing",
                               "03_R_objects"
                               )
project_dir       <- file.path(experiments_directory, "2022-10-31 - compare CRISPRoff screens")
output_dir        <- file.path(project_dir, "02_output_data")
figures_dir       <- output_dir
manuscript_dir    <- file.path(figures_dir, "Manuscript")
thesis_dir        <- file.path(figures_dir, "Thesis")



# Load data ---------------------------------------------------------------

load(file.path(off_2sg_rdata_dir, "09_create_figures_from_count_data.RData"))
load(file.path(off_2sg_rdata_dir, "10_recreate_figures_of_Nunez_et_al.RData"))
load(file.path(off_4sg_rdata_dir, "09_create_figures_from_count_data.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))

load(file.path(CRISPR_root_directory, "3) RData files", "1) General",
               "24) Enumerate pairs of genes at bidirectional promoters.RData"
               ))



# Define functions --------------------------------------------------------

ChooseOnePerEntrez <- function(input_df) {
  are_to_include <- (!(duplicated(input_df[, "Entrez_ID"]))) &
                    (!(input_df[, "Is_NT"]))
  input_df <- input_df[are_to_include, ]
  row.names(input_df) <- NULL
  return(input_df)
}


ScatterInputDf <- function(logfc_1_df, logfc_2_df, choose_rep = NULL) {
  if (is.null(choose_rep)) {
    logfc_column <- "Mean_log2FC"
  } else {
    logfc_column <- paste0("Log2FC_rep", choose_rep)
  }
  common_entrezs <- intersect(
    logfc_1_df[, "Entrez_ID"][(!(is.nan(logfc_1_df[, logfc_column])))],
    logfc_2_df[, "Entrez_ID"][(!(is.nan(logfc_2_df[, logfc_column])))]
  )
  common_entrezs <- sort(common_entrezs)
  logfc_1_df <- logfc_1_df[match(common_entrezs, logfc_1_df[, "Entrez_ID"]), ]
  logfc_2_df <- logfc_2_df[match(common_entrezs, logfc_2_df[, "Entrez_ID"]), ]
  results_df <- data.frame(
    "Entrez_ID" = common_entrezs,
    "Is_NT"     = FALSE,
    "Rep1_data" = logfc_1_df[, logfc_column],
    "Rep2_data" = logfc_2_df[, logfc_column]
  )
  return(results_df)
}



ThreeScatterPlots <- function(logfc_1_df,
                              logfc_2_df,
                              library_1_label,
                              library_2_label,
                              add_rep_to_label     = FALSE,
                              choose_rep           = NULL,
                              show_phenotype_score = TRUE,
                              num_cell_divisions   = 10L,
                              lower_bound          = -6 * (if (show_phenotype_score) 0.1 else 1),
                              upper_bound          = 6  * (if (show_phenotype_score) 0.1 else 1),
                              embed_PNG            = FALSE,
                              show_remaining_genes = FALSE,
                              layout_widths        = c(0.11, 0.26, 0.03, 0.26, 0.03, 0.26, 0.02),
                              layout_heights       = c(0.14, 0.66, 0.2),
                              x_axis_label_line    = 1.8,
                              y_axis_label_line    = 2.2,
                              label_gene_sets      = TRUE,
                              show_empty_ticks     = TRUE
                              ) {

  required_objects <- c("essentials_2020Q2_df", "non_essentials_2020Q2_df")
  stopifnot(all(required_objects %in% ls(envir = globalenv())))

  ## Prepare data
  scatter_df <- ScatterInputDf(logfc_1_df, logfc_2_df, choose_rep = choose_rep)
  xy_mat <- as.matrix(scatter_df[, c("Rep1_data", "Rep2_data")])
  if (show_phenotype_score) {
    xy_mat <- xy_mat / num_cell_divisions
  }

  ## Categorize genes
  are_essential <- scatter_df[, "Entrez_ID"] %in% essentials_2020Q2_df[, "Entrez_ID"]
  are_non_essential <- scatter_df[, "Entrez_ID"] %in% non_essentials_2020Q2_df[, "Entrez_ID"]

  ## Adjust axis limits
  xy_list <- lapply(1:2, function(x) {
    BringWithinLimits(xy_mat[, x], lower_bound = lower_bound, upper_bound = upper_bound)
  })
  xy_mat <- cbind(xy_list[[1]][["curtailed_vec"]], xy_list[[2]][["curtailed_vec"]])
  xy_range <- c(lower_bound, upper_bound)
  xy_space <- (xy_range[[2]] - xy_range[[1]]) * 0.025
  xy_lim <- c(xy_range[[1]] - xy_space, xy_range[[2]] + xy_space)

  ## Prepare axis tick labels
  tick_locations <- pretty(xy_range, n = 6)
  tick_labels_list <- lapply(1:2, function(i) {
    CurtailedAxisLabels(tick_locations,
                        lower_bound          = lower_bound,
                        upper_bound          = upper_bound,
                        lower_bound_enforced = xy_list[[i]][["lower_bound_enforced"]],
                        upper_bound_enforced = xy_list[[i]][["lower_bound_enforced"]]
                        )
  })

  ## Prepare axis labels
  library_1_label <- paste0(library_1_label, " (")
  library_2_label <- paste0(library_2_label, " (")
  if ((!(is.null(choose_rep))) && add_rep_to_label) {
    library_1_label <- paste0(library_1_label, "R", choose_rep, " ")
    library_2_label <- paste0(library_2_label, "R", choose_rep, " ")
  }
  if (show_phenotype_score) {
    axis_labels_list <- list(
      ConcatenateExpressions(list(library_1_label, expression("" * gamma * ")")), my_sep = ""),
      ConcatenateExpressions(list(library_2_label, expression("" * gamma * ")")), my_sep = "")
    )
  } else {
    axis_labels_list <- list(
      ConcatenateExpressions(list(library_1_label, expression("log"[2] * "FC")), my_sep = ""),
      ConcatenateExpressions(list(library_2_label, expression("log"[2] * "FC")), my_sep = ""),
    )
  }

  ## Prepare the 3 scatter plots
  are_selected_mat <- cbind(
    "Is_gene"          = rep(TRUE, nrow(xy_mat)),
    "Is_essential"     = are_essential,
    "Is_non_essential" = are_non_essential
  )
  point_colors <- c("black", "#6244ab", "#32853c") # "#745dac", "#45874d"
  point_colors[[1]] <- adjustcolor(point_colors[[1]], alpha.f = 0.3)
  point_colors[-1] <- adjustcolor(point_colors[-1], alpha.f = 0.4)
  top_labels <- c("All genes", "Essential genes", "Non-essential genes")
  if (show_remaining_genes) {
    are_selected_mat <- cbind(
      are_selected_mat[, 2:3],
      "Is_remaining" = !(are_essential | are_non_essential)
    )
    point_colors <- point_colors[c(2, 3, 1)]
    top_labels <- c(top_labels[c(2, 3)], "Remaining genes")
  }

  ## Set up the multi-plot layout
  layout_mat <- rbind(c(1, 2, 2, 2, 2, 2, 3),
                      5:11,
                      rep(4, 7)
                      )
  original_cex <- par("cex")
  layout(layout_mat,
         widths = layout_widths,
         heights = layout_heights
         )
  old_par <- par(mar = rep(0, 4), cex = original_cex)

  for (i in 1:5) {
    MakeEmptyPlot()
  }

  ## Draw the 3 scatter plots

  for (i in 1:3) {

    MakeEmptyPlot(xy_lim, xy_lim)

    if (embed_PNG) {
      PDF_mar <- par("mar")
      PDF_device <- dev.cur()
      temp_path <- file.path(output_dir, "temp.png")
      temp_width  <- par("pin")[[1]]
      temp_height <- par("pin")[[2]]
      current_par <- par(no.readonly = TRUE)
      png(filename = temp_path,
          width    = temp_width,
          height   = temp_height,
          units    = "in",
          res      = 900,
          bg       = "white"
          )
      par(lwd = current_par[["lwd"]])
      par(cex = current_par[["cex"]])
      par(mar = rep(0, 4))
      MakeEmptyPlot(xy_lim, xy_lim)
    }

    ## Set up plot canvas
    abline(v = 0, h = 0,  col = "gray85", lend = "butt")
    abline(a = 0, b = 1,  col = "gray85", lend = "butt")
    abline(a = 0, b = -1, col = "gray85", lend = "butt")

    ## Draw points
    points(xy_mat[are_selected_mat[, i], ],
           pch = 16,
           cex = 0.5,
           col = point_colors[[i]],
           xpd = NA
           )

    if (embed_PNG) {
      dev.off()
      raster_array <- png::readPNG(temp_path)
      file.remove(temp_path)
      dev.set(PDF_device)
      par(PDF_mar)
      rasterImage(raster_array,
                  xleft   = par("usr")[[1]], xright = par("usr")[[2]],
                  ybottom = par("usr")[[3]], ytop   = par("usr")[[4]]
                  )
    }

    ## Annotate plot
    if (label_gene_sets) {
      mtext(VerticalAdjust(top_labels[[i]]), line = 0.2, cex = par("cex"))
    }
    mtext(VerticalAdjust(axis_labels_list[[1]]), side = 1,
          line = x_axis_label_line, cex = par("cex")
          )
    if (i == 1) {
      mtext(VerticalAdjust(axis_labels_list[[2]]), side = 2,
            line = y_axis_label_line, cex = par("cex")
            )
    }
    ## Draw axes
    for (j in 1:2) {
      if (j == 1) {
        if ((which(tick_locations == 0) %% 2) == 0) {
          rep_vec <- c(FALSE, TRUE)
        } else {
          rep_vec <- c(TRUE, FALSE)
        }
        are_to_show <- rep(rep_vec, length.out = length(tick_locations))
        use_tick_labels <- ifelse(are_to_show, tick_labels_list[[j]], "")
      } else if ((i != 1) && (j == 2)) {
        use_tick_labels <- rep("", length(tick_locations))
      } else {
        use_tick_labels <- tick_labels_list[[j]]
      }
      if (show_empty_ticks || (j == 1) || (i == 1)) {
        axis(j,
             at       = tick_locations,
             labels   = use_tick_labels,
             mgp      = c(3, if (j == 1) 0.35 else 0.5, 0),
             tcl      = -0.3,
             las      = 1,
             lwd      = par("lwd"),
             cex.axis = 1 / 0.7
             )
      }
    }
    box()
    MakeEmptyPlot()
  }

  par(old_par)
  layout(1)
  return(invisible(NULL))
}



CompareScreenBars <- function(bars_vec,
                              short_labels      = FALSE,
                              use_tcl           = 0.375,
                              y_axis_label_line = 2.25,
                              y_axis_mgp        = 0.55,
                              gini_index        = TRUE,
                              y_axis_label      = if (gini_index) "Gini index" else NULL,
                              use_title         = if (gini_index) "Plasmid count heterogeneity at baseline" else NULL,
                              y_upper_limit     = if (gini_index) 0.5 else NULL,
                              lollipop          = FALSE,
                              lollipop_mode     = "standard",
                              lollipop_pch      = 16,
                              stem_color        = brewer.pal(9, "Blues")[[2]],
                              bar_color         = brewer.pal(9, "Blues")[[8]],
                              interrupt_x_axis  = FALSE,
                              draw_grid         = TRUE,
                              grid_lwd          = 0.75,
                              y_axis_n          = 5,
                              bar_width         = 0.45,
                              gap_ratio         = 1.6,
                              side_gap          = 0.5,
                              point_size_factor = 1,
                              label_libraries   = TRUE,
                              brackets_y_line   = 1.375,
                              group_labels_y    = 1.7
                              ) {

  stopifnot(length(bars_vec) %in% c(4, 6))
  stopifnot(lollipop_mode %in% c("standard", "bullseye", "bisected"))

  ## Determine bar positions
  groups_vec <- rep(seq_len(length(bars_vec) / 2), each = 2)
  bar_positions <- RepositionByGroups(groups_vec, gap_ratio = gap_ratio)
  num_bars <- length(bar_positions)
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - side_gap) - (num_bars * 0.04),
                    (max(bar_positions) + side_gap) + (num_bars * 0.04)
                     )

  ## Prepare the data axis
  use_numeric_limits <- c(0, max(bars_vec))
  if (is.null(y_upper_limit)) {
    y_upper_limit <- max(pretty(c(0, max(bars_vec))))
  }
  numeric_limits <- c(0, y_upper_limit)

  ## Draw the bars
  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)
  if (draw_grid) {
    if (all(numeric_limits[-1] >= 10)) {
      grid_pos <- pretty(numeric_limits, n = 15)
      are_major <- (grid_pos %% 10) < (10^-12)
    } else if (all(numeric_limits <= 0.5)) {
      grid_pos <- pretty(numeric_limits, n = 30)
      are_major <- ((grid_pos * 100) %% 10) < (10^-12)
    } else {
      grid_pos <- pretty(numeric_limits, n = 15)
      are_major <- ((grid_pos * 10) %% 5) < (10^-12)
    }
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[2]],
             y0  = grid_pos,
             col = ifelse(are_major, "gray88", "gray95"),
             lwd = par("lwd") * grid_lwd,
             xpd = NA
             )
  }
  if (lollipop) {
    segments(x0   = bar_positions,
             y0   = 0,
             y1   = par("usr")[[4]],
             col  = stem_color,
             lwd  = par("lwd") * 2,
             lend = "butt",
             xpd  = NA
             )
    if (lollipop_mode == "standard") {
      points(x   = bar_positions,
             y   = bars_vec,
             col = bar_color,
             cex = par("cex") * 1.5 * point_size_factor,
             pch = lollipop_pch,
             xpd = NA
             )
    } else if (lollipop_mode == "bisected") {
      points(x   = bar_positions,
             y   = bars_vec,
             col = bar_color,
             bg  = stem_color,
             cex = par("cex") * 1.7 * point_size_factor,
             pch = 21,
             xpd = NA
             )
      point_radius <- (par("cxy")[2] / pi) * par("cex")
      segments(x0  = bar_positions - point_radius,
               x1  = bar_positions + point_radius,
               y0  = bars_vec,
               col = bar_color,
               xpd = NA
               )
    } else if (lollipop_mode == "bullseye") {
      points(x   = bar_positions,
             y   = bars_vec,
             col = bar_color,
             cex = par("cex") * 1.7 * point_size_factor,
             pch = 16,
             xpd = NA
             )
      points(x   = bar_positions,
             y   = bars_vec,
             col = stem_color,
             cex = par("cex") * 0.25 * point_size_factor,
             pch = 16,
             xpd = NA
             )
    }
  } else {
    PlotBarplotMat(t(bars_vec),
                   colors_vec    = bar_color,
                   positions_vec = bar_positions,
                   bar_width     = bar_width
                   )
  }

  ## Draw the y axis
  tick_locations <- pretty(numeric_limits, n = y_axis_n)
  axis(2,
       at     = tick_locations,
       labels = format(tick_locations),
       las    = 2,
       mgp    = c(3, y_axis_mgp, 0),
       tcl    = -(use_tcl),
       lwd    = par("lwd")
       )
  if (!(is.null(y_axis_label))) {
    mtext(VerticalAdjust(y_axis_label),
          side = 2,
          line = y_axis_label_line,
          cex  = par("cex")
          )
  }

  ## Annotate bars
  mtext(text = rep(paste0("R", 1:2), 2),
        at = bar_positions, side = 1, line = 0.25, cex = par("cex")
        )
  if (label_libraries) {
    segments(x0  = tapply(bar_positions, groups_vec, min) - 0.2,
             x1  = tapply(bar_positions, groups_vec, max) + 0.2,
             y0  = par("usr")[[3]] - diff(grconvertY(c(0, brackets_y_line), from = "lines", to = "user")),
             col = "gray50",
             xpd = NA
             )
    group_labels <- c("CRISPRoff", "T.gonfio")
    if (!(short_labels)) {
      group_labels <- paste0(group_labels, " library")
    }
    if (length(bars_vec) == 6) {
      group_labels <- c(if (short_labels) "Re-analysis" else "Data re-analysis", group_labels)
    }
    mtext(text = group_labels,
          at   = tapply(bar_positions, groups_vec, mean),
          side = 1,
          line = group_labels_y,
          cex  = par("cex")
          )
  }
  if (!(is.null(use_title))) {
    mtext(use_title, side = 3, line = 1, cex = par("cex"))
  }
  box(bty = "l")
  if (interrupt_x_axis && (length(bars_vec) == 4)) {
    x_gap <- 0.15
    x_mid <- mean(bar_positions[2:3])
    rect(xleft   = x_mid - (x_gap / 2),
         xright  = x_mid + (x_gap / 2),
         ybottom = par("usr")[[3]] - diff(grconvertY(c(0, 0.3), from = "lines", to = "user")),
         ytop    = par("usr")[[4]] + diff(grconvertY(c(0, 0.3), from = "lines", to = "user")),
         col     = "white",
         border  = NA,
         xpd     = NA
         )
  }
  return(invisible(bar_positions))
}



DrawBottomLegend <- function(labels_list,
                             use_colors,
                             border_colors        = NULL,
                             use_pch              = 16,
                             use_point_size       = 1.2,
                             lines_x_start        = 0.7,
                             lines_y_start        = 1.3,
                             y_mid                = 0.5,
                             small_gap_size       = 1.25,
                             large_gap_multiplier = 1.4,
                             use_lwd              = 2,
                             line_x_distance      = -0.2
                             ) {

  lines_x_start <- lines_x_start - diff(grconvertX(c(0, strwidth(expression(""^"gh"))), from = "user", to = "lines"))

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing legend
  small_gap <- diff(grconvertY(c(0, small_gap_size), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_multiplier

  text_list <- labels_list
  if (all(lengths(text_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(text_list))
    are_first <- rep(TRUE, length(text_list))
  } else {
    are_first <- unlist(lapply(text_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- total_span + diff(grconvertY(c(0, lines_y_start), from = "lines", to = "npc"))
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  all_expressions <- sapply(unlist(labels_list), VerticalAdjust)

  x_text  <- par("usr")[[2]] -
             diff(grconvertX(c(0, lines_x_start), from = "lines", to = "user")) -
             max(strwidth(all_expressions))

  ## Draw legend
  x_user <- grconvertX(x = x_text, from = "npc", to = "user")
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       labels = all_expressions,
       adj    = c(0, 0.5),
       xpd    = NA
       )
  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  length_in_lines <- 0.5
  line_x_start <- x_user + diff(grconvertX(c(0, line_x_distance), from = "lines", to = "user"))
  segments(x0  = line_x_start,
           x1  = line_x_start + diff(grconvertX(c(0, length_in_lines), from = "lines", to = "user")),
           y0  = tapply(y_pos, groups_vec, mean),
           col = use_colors,
           lwd = par("lwd") * use_lwd,
           xpd = NA
           )

  return(invisible(NULL))
}



ThreeLinesROC <- function(ROC_df_list,
                          flip          = TRUE,
                          embed_PNG     = FALSE,
                          transparency  = TRUE,
                          use_lwd       = 2,
                          legend_lwd    = use_lwd * 1.25,
                          legend_order  = seq_along(ROC_df_list),
                          axis_limits   = c(0, 1),
                          legend_inside = TRUE,
                          middle_line   = !(legend_inside),
                          ...
                          ) {

  use_colors <- c(brewer.pal(9, "Greys")[[7]],
                  brewer.pal(9, "Blues")[[6]],
                  "#B8363F"
                  )
  if (transparency) {
    legend_colors <- c(Palify("black", fraction_pale = 1 - 0.55),
                       Palify(use_colors[-1], fraction_pale = 1 - 0.8)
                       )
    line_colors <- c(adjustcolor("black", alpha.f = 0.55),
                     adjustcolor(use_colors[-1], alpha.f = 0.8)
                     )
  } else {
    legend_colors <- use_colors
    line_colors <- use_colors
  }

  ROC_mat_list <- lapply(ROC_df_list, GetROCMat)
  AUC_vec <- vapply(ROC_mat_list, GetAUC, numeric(1))

  if (flip) {
    for (i in seq_along(ROC_mat_list)) {
      sens_vec <- ROC_mat_list[[i]][, "Sensitivity"]
      spec_vec <- ROC_mat_list[[i]][, "Specificity"]
      ROC_mat_list[[i]][, "Sensitivity"] <- spec_vec
      ROC_mat_list[[i]][, "Specificity"] <- sens_vec
    }
  }

  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(figures_dir, "temp.png")
    temp_width  <- par("pin")[[1]]
    temp_height <- par("pin")[[2]]
    current_par <- par(no.readonly = TRUE)
    png(filename = temp_path,
        width    = temp_width,
        height   = temp_height,
        units    = "in",
        res      = 900,
        bg       = "white",
        type     = "cairo-png"
        )
    par(lwd = current_par[["lwd"]])
    par(cex = current_par[["cex"]])
    par(mar = rep(0, 4))
  }

  MakeEmptyPlot(axis_limits, axis_limits)
  if (middle_line) {
    abline(a = 0, b = 1, col = "gray78", lty = "dashed")
  }
  if (!(embed_PNG)) {
    use_tcl <- -0.36
    axis(1, mgp = c(3, 0.4, 0), tcl = use_tcl, lwd = par("lwd"))
    mtext("False positive rate", side = 1, line = 1.75, cex = par("cex"))
    axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1, lwd = par("lwd"))
    mtext("True positive rate", side = 2, line = 2.3, cex = par("cex"))
    box()
  }
  for (i in seq_along(ROC_mat_list)) {
    lines(x   = 1 - ROC_mat_list[[i]][, "Specificity"],
          y   = ROC_mat_list[[i]][, "Sensitivity"],
          lwd = use_lwd,
          col = line_colors[[i]],
          xpd = NA
          )
  }

  if (embed_PNG) {
    dev.off()
    raster_array <- png::readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)
    MakeEmptyPlot(axis_limits, axis_limits)
    rasterImage(raster_array,
                xleft   = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop   = par("usr")[[4]]
                )
  }

  if (embed_PNG) {
    use_tcl <- -0.36
    axis(1, mgp = c(3, 0.4, 0), tcl = use_tcl, lwd = par("lwd"))
    mtext("False positive rate", side = 1, line = 2, cex = par("cex"))
    axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1, lwd = par("lwd"))
    mtext("True positive rate", side = 2, line = 2.3, cex = par("cex"))
    box()
  }

  ## Draw legend
  legend_vec <- c("Re-analysis", "CRISPRoff", "T.gonfio")
  AUC_legend_vec <- format(round(AUC_vec, digits = 2), nsmall = 2)
  AUC_legend_vec <- sapply(AUC_legend_vec, function(x) {
    as.expression(bquote("(AUC" * scriptscriptstyle(" ") * "=" * scriptscriptstyle(" ") * .(x) * ")"))
  })
  labels_list <- lapply(seq_along(legend_vec), function(x) c(as.expression(legend_vec[[x]]), AUC_legend_vec[[x]]))

  if (legend_inside) {
    DrawBottomLegend(labels_list = labels_list[legend_order],
                     use_colors = legend_colors[legend_order],
                     use_lwd = legend_lwd,
                     ...
                     )
  } else {
    DrawSideLegend(labels_list = labels_list[legend_order],
                   use_colors = legend_colors[legend_order],
                   use_lwd = legend_lwd,
                   draw_lines = TRUE,
                   ...
                   )
  }

  return(invisible(NULL))
}




# Compare Gini indices at baseline ----------------------------------------

CompareScreenBars(c(gini_indices_CRISPRoff[c(1, 2)], gini_indices_4sg[c(1, 2)]),
                  lollipop = FALSE
                  )

CompareScreenBars(c(gini_indices_CRISPRoff[c(1, 2)], gini_indices_4sg[c(1, 2)]),
                  lollipop = TRUE
                  )

pdf(file.path(manuscript_dir, "Gini index comparison.pdf"),
    width = 2, height = 1.75
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
CompareScreenBars(c(gini_indices_CRISPRoff[c(1, 2)], gini_indices_4sg[c(1, 2)]),
                  short_labels = TRUE, use_title = "",
                  use_tcl = 0.35, y_axis_label_line = 2.2, y_axis_mgp = 0.525
                  )
title("Count heterogeneity", cex.main = 1, font.main = 1)
dev.off()



devEMF::emf(file.path(thesis_dir, "2B) Gini index comparison.emf"),
            width = 1.55, height = 1.75, emfPlus = FALSE, coordDPI = 1500
            )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.7)
CompareScreenBars(c(gini_indices_CRISPRoff[c(1, 2)], gini_indices_4sg[c(1, 2)]),
                  short_labels = TRUE, use_title = "",
                  use_tcl = 0.35, y_axis_label_line = 1.8, y_axis_mgp = 0.525,
                  lollipop = TRUE, lollipop_pch = 18,
                  gap_ratio = 1.7, side_gap = 0.6, point_size_factor = 1.25,
                  grid_lwd = 0.7
                  )
title("Count heterogeneity", cex.main = 1, font.main = 1)
dev.off()




# Plot the numbers of missing plasmids at baseline ------------------------

devEMF::emf(file.path(thesis_dir, "2C) Missing plasmids at baseline.emf"),
            width = 1.55, height = 1.75, emfPlus = FALSE, coordDPI = 1500
            )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.7)
CompareScreenBars(c(num_missing_CRISPRoff[1:2], num_missing_4sg[1:2]),
                  short_labels = TRUE,
                  use_tcl = 0.35, y_axis_label_line = 1.8, y_axis_mgp = 0.525,
                  y_axis_label = "Number of missed plasmids",
                  use_title = "", gini_index = FALSE, bar_width = 0.4,
                  y_upper_limit = 100, bar_color = brewer.pal(9, "Blues")[[7]],
                  gap_ratio = 1.7, side_gap = 0.6, grid_lwd = 0.7
                  )
title("Absent at baseline", cex.main = 1, font.main = 1)
dev.off()




# Plot the separation between essential and non-essential genes -----------

CompareScreenBars(abs(c(separation_CRISPRoff_mat["SSMD", ], separation_4sg_mat["SSMD", ])),
                  gini_index = FALSE, y_axis_label = "SSMD"
                  )
CompareScreenBars(abs(c(separation_CRISPRoff_mat["Robust SSMD", ], separation_4sg_mat["Robust SSMD", ])),
                  gini_index = FALSE, y_axis_label = "SSMD*"
                  )
CompareScreenBars(abs(c(separation_CRISPRoff_mat["SSMD finite", ], separation_4sg_mat["SSMD finite", ])),
                  gini_index = FALSE, y_axis_label = "SSMD"
                  )

SSMD_methods <- c(
  "SSMD"                = expression("SSMD (\u2012infinity" %->% "\u20120.6)"),
  "SSMD finite"         = expression("SSMD (\u2012infinity excluded)"),
  "Robust SSMD"         = expression("SSMD* (\u2012infinity" %->% "\u20120.6)"),
  "Robust SSMD finite"  = expression("SSMD* (\u2012infinity excluded)")
)


pdf(file.path(output_dir, "SSMD - bar chart.pdf"),
    width = 2, height = 1.75
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
for (use_method in names(SSMD_methods)) {
  this_bars_vec <- abs(c(separation_CRISPRoff_mat[use_method, ],
                         separation_4sg_mat[use_method, ]
                         ))
  custom_y <- all(this_bars_vec < 1.5)
  CompareScreenBars(this_bars_vec,
                    short_labels = TRUE,
                    use_tcl = 0.35, y_axis_label_line = 2.2, y_axis_mgp = 0.525,
                    gini_index = FALSE, y_axis_label = SSMD_methods[[use_method]],
                    y_upper_limit = if (custom_y) 1.5 else NULL,
                    y_axis_n = if (custom_y) 3 else 5
                    )
}
dev.off()



for (lollipop_style in c("standard", "bullseye", "bisected")) {
  pdf(file.path(output_dir, paste0("SSMD - dot chart - ", lollipop_style, ".pdf")),
      width = 2, height = 1.75
      )
  old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
  for (use_method in names(SSMD_methods)) {
    this_bars_vec <- abs(c(separation_CRISPRoff_mat[use_method, ],
                           separation_4sg_mat[use_method, ]
                           ))
    custom_y <- all(this_bars_vec < 1.5)
    CompareScreenBars(this_bars_vec,
                      lollipop = TRUE, lollipop_mode = lollipop_style,
                      short_labels = TRUE,
                      use_tcl = 0.35, y_axis_label_line = 2.2, y_axis_mgp = 0.525,
                      gini_index = FALSE, y_axis_label = SSMD_methods[[use_method]],
                      y_upper_limit = if (custom_y) 1.5 else NULL,
                      y_axis_n = if (custom_y) 3 else 5
                      )
  }
  dev.off()
}



pdf(file.path(manuscript_dir, "Figure 7H - SSMD - bar chart.pdf"),
    width = 2, height = 2
    )
old_par <- par(mar = c(3, 4, 4, 1), cex = 0.6, lwd = 0.8)
this_bars_vec <- abs(c(separation_CRISPRoff_mat["Robust SSMD", ],
                       separation_4sg_mat["Robust SSMD", ]
                       ))
CompareScreenBars(this_bars_vec,
                  short_labels = TRUE,
                  use_tcl = 0.35, y_axis_label_line = 2.1, y_axis_mgp = 0.525,
                  gini_index = FALSE, y_axis_label = "SSMD*",
                  y_upper_limit = 1.5, draw_grid = FALSE,
                  y_axis_n = 3, bar_color = brewer.pal(9, "Blues")[[7]], #6385aa",
                  gap_ratio = 1.7, side_gap = 0.6
                  )
dev.off()



devEMF::emf(file.path(thesis_dir, "3C) SSMD - bar chart.emf"),
            width = 2.3, height = 2, emfPlus = FALSE, coordDPI = 3000
            )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.7)
this_bars_vec <- abs(c(separation_original_mat["Robust SSMD", ],
                       separation_CRISPRoff_mat["Robust SSMD", ],
                       separation_4sg_mat["Robust SSMD", ]
                       ))
bar_positions <- CompareScreenBars(
  this_bars_vec,
  short_labels = TRUE,
  use_tcl = 0.35, y_axis_label_line = 2, y_axis_mgp = 0.525,
  gini_index = FALSE, y_axis_label = "SSMD*",
  y_upper_limit = 1.5, draw_grid = TRUE,
  y_axis_n = 3, bar_color = "#7b98b7", stem_color = "#e8edf2",
  gap_ratio = 1.7, side_gap = 0.6,
  bar_width = 0.4, lollipop = TRUE, point_size_factor = 1.25,
  group_labels_y = 1.75
)
dev.off()



# Prioritize plasmids included in the new CRISPRoff 2sg library -----------

intersect_plasmids <- intersect(logfc_CRISPRoff_df[, "sgID"], logfc_original_df[, "sgID"])
matches_vec <- match(logfc_original_df[, "sgID"], logfc_CRISPRoff_df[, "sgID"])
logfc_original_df <- logfc_original_df[order(matches_vec), ]
row.names(logfc_original_df) <- NULL



# Choose one plasmid for each Entrez gene ID ------------------------------

logfc_4sg_df       <- ChooseOnePerEntrez(logfc_4sg_df)
logfc_CRISPRoff_df <- ChooseOnePerEntrez(logfc_CRISPRoff_df)
logfc_original_df  <- ChooseOnePerEntrez(logfc_original_df)

common_genes <- intersect(logfc_CRISPRoff_df[, "Entrez_ID"], logfc_4sg_df[, "Entrez_ID"])
common_genes <- intersect(logfc_original_df[, "Entrez_ID"], common_genes)

non_NA_common_genes <- intersect(common_genes,
                                 logfc_CRISPRoff_df[, "Entrez_ID"][!(is.na(logfc_CRISPRoff_df[, "Mean_log2FC"]))]
                                 )
non_NA_common_genes <- intersect(non_NA_common_genes,
                                 logfc_4sg_df[, "Entrez_ID"][!(is.na(logfc_4sg_df[, "Mean_log2FC"]))]
                                 )



# Create scatter plots ----------------------------------------------------

ReplicateScatterPlot(ScatterInputDf(logfc_4sg_df, logfc_CRISPRoff_df),
                     highlight_NT = FALSE,
                     axis_labels_list = list(expression("T.gonfio" ~ "(" * gamma * ")"),
                                             expression("CRISPRoff" ~ "(" * gamma * ")")
                                             )
                     )

ReplicateScatterPlot(ScatterInputDf(logfc_4sg_df, logfc_original_df),
                     highlight_NT = FALSE,
                     axis_labels_list = list(expression("T.gonfio" ~ "(" * gamma * ")"),
                                             expression("Nu\u00f1ez et al." ~ "(" * gamma * ")")
                                             )
                     )

ReplicateScatterPlot(ScatterInputDf(logfc_CRISPRoff_df, logfc_original_df),
                     highlight_NT = FALSE,
                     axis_labels_list = list(expression("CRISPRoff" ~ "(" * gamma * ")"),
                                             expression("Nu\u00f1ez et al." ~ "(" * gamma * ")")
                                             )
                     )

ThreeScatterPlots(logfc_4sg_df,
                  logfc_CRISPRoff_df,
                  "T.gonfio",
                  "CRISPRoff"
                  )
ThreeScatterPlots(logfc_4sg_df,
                  logfc_original_df,
                  "T.gonfio",
                  "Nu\u00f1ez et al."
                  )
ThreeScatterPlots(logfc_CRISPRoff_df,
                  logfc_original_df,
                  "CRISPRoff",
                  "Nu\u00f1ez et al."
                  )



# Export scatter plots ----------------------------------------------------

PDF_height <- 3
use_heights <- c(0.14, 0.66, 0.2)
use_widths <- c(0.11, 0.26, 0.03, 0.26, 0.03, 0.26, 0.02)
plot_height <- PDF_height * (use_heights[[2]] / sum(use_heights))
PDF_width <- (plot_height * 3) + ((sum(use_widths[c(1, 3, 5, 7)]) / use_widths[[2]]) * plot_height)

for (make_PDF in c(FALSE, TRUE)) {
  for (choose_rep in list(NULL, 1, 2)) {
    for (show_other_genes in c(FALSE, TRUE)) {

      file_name <- "Three plots - "
      if (is.null(choose_rep)) {
        file_name <- paste0(file_name, "both replicates - ")
      } else {
        file_name <- paste0(file_name, "replicate ", choose_rep, " - ")
      }
      if (show_other_genes) {
        file_name <- paste0(file_name, "remaining genes")
      } else {
        file_name <- paste0(file_name, "all genes")
      }

      if (make_PDF) {
        pdf(file = file.path(output_dir, paste0(file_name, ".pdf")),
            width = PDF_width, height = PDF_height
            )
        ThreeScatterPlots(logfc_4sg_df,
                          logfc_CRISPRoff_df,
                          "T.gonfio",
                          "CRISPRoff",
                          embed_PNG = TRUE,
                          show_remaining_genes = show_other_genes,
                          choose_rep = choose_rep
                          )
        ThreeScatterPlots(logfc_4sg_df,
                          logfc_original_df,
                          "T.gonfio",
                          "Nu\u00f1ez et al.",
                          embed_PNG = TRUE,
                          show_remaining_genes = show_other_genes,
                          choose_rep = choose_rep
                          )
        ThreeScatterPlots(logfc_CRISPRoff_df,
                          logfc_original_df,
                          "CRISPRoff",
                          "Nu\u00f1ez et al.",
                          embed_PNG = TRUE,
                          show_remaining_genes = show_other_genes,
                          choose_rep = choose_rep
                          )
      }

      use_mai <- c(0.75, 0.75, 0.4, 1.50)
      plot_height <- 2.5
      single_width <- plot_height + sum(use_mai[c(2, 4)])
      single_height <- plot_height + sum(use_mai[c(1, 3)])

      common_args <- list(
        highlight_NT         = FALSE,
        lower_bound          = -6,
        upper_bound          = 6,
        show_phenotype_score = TRUE,
        use_mar              = use_mai * 5,
        embed_PNG            = TRUE
      )

      if (make_PDF) {
        dev.off()
        file_name <- "Single plots - "
        if (is.null(choose_rep)) {
          file_name <- paste0(file_name, "both replicates")
        } else {
          file_name <- paste0(file_name, "replicate ", choose_rep)
        }
        pdf(file = file.path(output_dir, paste0(file_name, ".pdf")),
            width = single_width, height = single_height
            )
      }

      args_list <- c(list(ScatterInputDf(logfc_4sg_df, logfc_CRISPRoff_df),
                          axis_labels_list = list(expression("T.gonfio" ~ "(" * gamma * ")"),
                                                  expression("CRISPRoff" ~ "(" * gamma * ")")
                                                  )
                          ), common_args)
      do.call(ReplicateScatterPlot, args_list)


      args_list <- c(list(ScatterInputDf(logfc_4sg_df, logfc_original_df),
                          axis_labels_list = list(expression("T.gonfio" ~ "(" * gamma * ")"),
                                                  expression("Nu\u00f1ez et al." ~ "(" * gamma * ")")
                                                  )
                          ), common_args)
      do.call(ReplicateScatterPlot, args_list)


      args_list <- c(list(ScatterInputDf(logfc_CRISPRoff_df, logfc_original_df),
                          axis_labels_list = list(expression("CRISPRoff" ~ "(" * gamma * ")"),
                                                  expression("Nu\u00f1ez et al." ~ "(" * gamma * ")")
                                                  )
                          ), common_args)
      do.call(ReplicateScatterPlot, args_list)


      if (make_PDF) {
        dev.off()
      }
    }
  }
}



# Export scatter plots for the manuscript ---------------------------------

PDF_height <- 1.75
top_gap <- 0.1377392
left_gap <- 0.1085217
use_heights <- c(top_gap, 0.66, 0.2 + (0.14 - top_gap))
use_widths <- c(left_gap, 0.26, 0.045, 0.26, 0.045, 0.26, 0.02 + (0.11 - left_gap))
plot_height <- PDF_height * (use_heights[[2]] / sum(use_heights))
PDF_width <- (plot_height * 3) + ((sum(use_widths[c(1, 3, 5, 7)]) / use_widths[[2]]) * plot_height)

pdf(file = file.path(manuscript_dir, "Scatter plots - CRISPRoff vs. T.gonfio - both replicates.pdf"),
    width = PDF_width, height = PDF_height
    )
par(cex = 0.6, lwd = 0.8)
ThreeScatterPlots(logfc_4sg_df,
                  logfc_CRISPRoff_df,
                  "T.gonfio",
                  "CRISPRoff",
                  y_axis_label_line = 2.5,
                  layout_widths     = use_widths,
                  layout_heights    = use_heights,
                  embed_PNG         = TRUE
                  )
dev.off()



pdf(file = file.path(manuscript_dir, "Scatter plots - CRISPRoff vs. T.gonfio - replicate 2.pdf"),
    width = PDF_width, height = PDF_height
    )
par(cex = 0.6, lwd = 0.8)
ThreeScatterPlots(logfc_4sg_df,
                  logfc_CRISPRoff_df,
                  "T.gonfio",
                  "CRISPRoff",
                  choose_rep           = 2,
                  add_rep_to_label     = TRUE,
                  show_remaining_genes = TRUE,
                  y_axis_label_line    = 2.5,
                  layout_widths        = use_widths,
                  layout_heights       = use_heights,
                  embed_PNG            = TRUE
                  )
dev.off()



# Export scatter plots for the thesis -------------------------------------

emf_widths <- c(left_gap, 0.27, 0.03, 0.27, 0.03, 0.27, 0.02 + (0.11 - left_gap))
emf_width <- (plot_height * 3) + ((sum(emf_widths[c(1, 3, 5, 7)]) / emf_widths[[2]]) * plot_height)

df_names <- c(
  "T.gonfio"    = "logfc_4sg_df",
  "CRISPRoff"   = "logfc_CRISPRoff_df",
  "Re-analysis" = "logfc_original_df"
)
dataset_combos <- list(
  c("T.gonfio", "CRISPRoff"),
  c("T.gonfio", "Re-analysis"),
  c("CRISPRoff", "Re-analysis")
)

for (i in 1:6) {
  if (i %in% 1:3) {
    use_rep <- NULL
  } else {
    use_rep <- 2L
  }
  use_index <- rep(1:3, length.out = i)[[i]]
  file_name <- paste0("5A ", tolower(as.character(as.roman(i))), ") ",
                      "Scatter plots - CRISPRoff vs. T.gonfio - ",
                      if (is.null(use_rep)) "both replicates" else paste0("replicate ", use_rep)
                      )
  devEMF::emf(file = file.path(thesis_dir, paste0(file_name, ".emf")),
              width = PDF_width, height = PDF_height, emfPlus = FALSE, coordDPI = 1500
              )
  par(cex = 0.6, lwd = 0.7)
  data1_name <- dataset_combos[[use_index]][[1]]
  data2_name <- dataset_combos[[use_index]][[2]]
  ThreeScatterPlots(get(df_names[[data1_name]]),
                    get(df_names[[data2_name]]),
                    data1_name,
                    data2_name,
                    label_gene_sets      = use_index == 1,
                    choose_rep           = use_rep,
                    add_rep_to_label     = !(is.null(use_rep)),
                    show_remaining_genes = TRUE,
                    x_axis_label_line    = 1.7,
                    y_axis_label_line    = 2.1,
                    layout_widths        = use_widths,
                    layout_heights       = use_heights,
                    show_empty_ticks     = FALSE,
                    embed_PNG            = TRUE
                    )
  dev.off()
}




# Export combined ROC curves ----------------------------------------------

essential_entrezs     <- intersect(essentials_2020Q2_df[, "Entrez_ID"], non_NA_common_genes)
non_essential_entrezs <- intersect(non_essentials_2020Q2_df[, "Entrez_ID"], non_NA_common_genes)

ROC_original_input_df  <- ROCInputDf(logfc_original_df, essential_entrezs, non_essential_entrezs)
ROC_CRISPRoff_input_df <- ROCInputDf(logfc_CRISPRoff_df, essential_entrezs, non_essential_entrezs)
ROC_4sg_input_df       <- ROCInputDf(logfc_4sg_df, essential_entrezs, non_essential_entrezs)

common_ROC_original_df  <- ROCDfForColumn(ROC_original_input_df, "Mean_log2FC")
common_ROC_CRISPRoff_df <- ROCDfForColumn(ROC_CRISPRoff_input_df, "Mean_log2FC")
common_ROC_4sg_df       <- ROCDfForColumn(ROC_4sg_input_df, "Mean_log2FC")



PlotROCDf(ROC_original_df, show_AUC = FALSE,
          line_color = brewer.pal(9, "Greys")[[7]]
          )
PlotROCDf(ROC_CRISPRoff_df, add = TRUE,
          line_color = brewer.pal(9, "Blues")[[7]]
          )
PlotROCDf(ROC_4sg_df, add = TRUE,
          line_color = brewer.pal(9, "Reds")[[4]]
          )
box()




individual_ROC_df_list <- list(ROC_original_df, ROC_CRISPRoff_df, ROC_4sg_df)
ROC_df_list <- list(common_ROC_original_df, common_ROC_CRISPRoff_df, common_ROC_4sg_df)

ThreeLinesROC(individual_ROC_df_list)
ThreeLinesROC(ROC_df_list)


ThreeLinesROC(ROC_df_list, embed_PNG = TRUE, small_gap_size = 1.2, large_gap_multiplier = 1.3,
              use_lwd = 1.5, legend_lwd = 2.5,
              line_x_distance = -0.4, legend_order = c(3, 1, 2),
              axis_limits = c(-0.02, 1.02)
              )



devEMF::emf(file.path(thesis_dir, "3B) ROC curves.emf"),
            width = 2.72, height = 2, emfPlus = FALSE, coordDPI = 3000
            )
old_par <- par(mar = c(3, 4, 2, 7), cex = 0.6, lwd = 0.7)
ThreeLinesROC(ROC_df_list, embed_PNG = FALSE, small_gap_size = 1.25, large_gap_multiplier = 1.5,
              transparency = FALSE, use_lwd = 1.75, legend_lwd = 2.5,
              line_x_distance = -0.5, legend_order = c(3, 1, 2), legend_inside = FALSE,
              length_in_lines = 0.55, lines_x_start = 1
              )
par(old_par)
dev.off()




# Examine bidirectional promoters -----------------------------------------

common_logfc_original_df <- logfc_original_df[logfc_original_df[, "Entrez_ID"] %in% common_genes, ]
row.names(common_logfc_original_df) <- NULL

common_logfc_CRISPRoff_df <- logfc_CRISPRoff_df[logfc_CRISPRoff_df[, "Entrez_ID"] %in% common_genes, ]
row.names(common_logfc_CRISPRoff_df) <- NULL

common_logfc_4sg_df <- logfc_4sg_df[logfc_4sg_df[, "Entrez_ID"] %in% common_genes, ]
row.names(common_logfc_4sg_df) <- NULL


BidirectionalViolins(bidirectional_df, common_logfc_original_df, max_distance = 20000,
                     num_controls = 30L, compare_across = FALSE, point_cex = 0.7
                     )
BidirectionalViolins(bidirectional_df, common_logfc_CRISPRoff_df, max_distance = 20000,
                     num_controls = 30L, compare_across = FALSE, point_cex = 0.7
                     )
BidirectionalViolins(bidirectional_df, common_logfc_4sg_df, max_distance = 20000,
                     num_controls = 30L, compare_across = FALSE, point_cex = 0.7
                     )


datasets_vec <- c(
  "Data re-analysis"  = "common_logfc_original_df",
  "CRISPRoff library" = "common_logfc_CRISPRoff_df",
  "T.gonfio library"  = "common_logfc_4sg_df"
)


for (i in seq_along(datasets_vec)) {
  file_name <- paste0("4) Bidirectional violins - ",
                      tolower(as.character(as.roman(i))), ") ",
                      names(datasets_vec)[[i]]
                      )
  devEMF::emf(file.path(thesis_dir, paste0(file_name, ".emf")),
              width = 66, height = 42, emfPlus = FALSE
              )
  old_par <- par(mar = c(3, 4, 2, 1), cex = 10.6, lwd = 10.5)
  BidirectionalViolins(bidirectional_df, get(datasets_vec[[i]]), max_distance = 20000,
                       num_controls = 30L, compare_across = FALSE, point_cex = 0.8,
                       zero_lty = "solid", zero_color = "gray86",
                       quantiles_lty = c("dotted", "dashed", "dotted"), # compatibility with emf device
                       annotation_cex = 1, draw_groups_n = FALSE, swarm_method = "compactswarm",
                       y_limits = c(-0.52, 0.2), label_points = TRUE,
                       gap_ratio = 1.2, wex = 0.86,
                       side_gap = 0.4, right_gap = 0.3, use_spacing = 0.85,
                       y_start_adj = -0.6, y_line_adj = +0.1
                       )
  mtext(names(datasets_vec)[[i]], at = 0.7, adj = 0, line = -0.3,
        cex = par("cex"), font = 2
        )
  par(old_par)
  dev.off()
}






