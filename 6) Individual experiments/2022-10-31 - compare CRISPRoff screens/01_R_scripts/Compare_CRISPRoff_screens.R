### 2022-10-31


# Load packages and source code -------------------------------------------

library("devEMF")

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(R_functions_dir, "07_comparing_CRISPRoff_screens.R"))



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
  required_objects <- c("essentials_2020Q2_df", "non_essentials_2020Q2_df")
  stopifnot(all(required_objects %in% ls(envir = globalenv())))
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
  are_essential <- common_entrezs %in% essentials_2020Q2_df[, "Entrez_ID"]
  are_non_essential <- common_entrezs %in% non_essentials_2020Q2_df[, "Entrez_ID"]
  results_df <- data.frame(
    "Entrez_ID"      = common_entrezs,
    "Essentiality"   = ifelse(are_essential,
                              "essential",
                              ifelse(are_non_essential, "non-essential", "neither")
                              ),
    "Is_NT"          = FALSE,
    "Rep1_data"      = logfc_1_df[, logfc_column],
    "Rep2_data"      = logfc_2_df[, logfc_column],
    "Rep1_min_spec"  = logfc_1_df[, "Min_specificity"],
    "Rep2_min_spec"  = logfc_2_df[, "Min_specificity"],
    "Rep1_comb_spec" = logfc_1_df[, "Combined_specificity"],
    "Rep2_comb_spec" = logfc_2_df[, "Combined_specificity"],
    stringsAsFactors = FALSE
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
                              show_empty_ticks     = TRUE,
                              point_cex            = 0.5
                              ) {

  ## Prepare data
  scatter_df <- ScatterInputDf(logfc_1_df, logfc_2_df, choose_rep = choose_rep)
  xy_mat <- as.matrix(scatter_df[, c("Rep1_data", "Rep2_data")])
  if (show_phenotype_score) {
    xy_mat <- xy_mat / num_cell_divisions
  }

  ## Categorize genes
  are_essential <- scatter_df[, "Essentiality"] %in% "essential"
  are_non_essential <- scatter_df[, "Essentiality"] %in% "non-essential"

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
      current_device <- StartEmbedPNG(output_dir)
      MakeEmptyPlot(xy_lim, xy_lim)
    }

    ## Set up plot canvas
    abline(v = 0, h = 0,  col = "gray85", lend = "butt")
    abline(a = 0, b = 1,  col = "gray85", lend = "butt")
    abline(a = 0, b = -1, col = "gray85", lend = "butt")

    ## Draw points
    points(xy_mat[are_selected_mat[, i], ],
           pch = 16,
           cex = point_cex,
           col = point_colors[[i]],
           xpd = NA
           )

    if (embed_PNG) {
      StopEmbedPNG(current_device, output_dir, make_empty_plot = FALSE)
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
                              group_labels_y    = 1.7,
                              bracket_color     = "gray50"
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
             cex = 1 * point_size_factor,
             pch = lollipop_pch,
             xpd = NA
             )
    } else if (lollipop_mode == "bisected") {
      points(x   = bar_positions,
             y   = bars_vec,
             col = bar_color,
             bg  = stem_color,
             cex = 1.2 * point_size_factor,
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
             cex = 1.2 * point_size_factor,
             pch = 16,
             xpd = NA
             )
      points(x   = bar_positions,
             y   = bars_vec,
             col = stem_color,
             cex = 0.175 * point_size_factor,
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
             col = bracket_color,
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



# Prioritize plasmids included in the new CRISPRoff 2sg library -----------

intersect_plasmids <- intersect(logfc_CRISPRoff_df[, "sgID"], logfc_original_df[, "sgID"])
matches_vec <- match(logfc_original_df[, "sgID"], logfc_CRISPRoff_df[, "sgID"])
logfc_original_df <- logfc_original_df[order(matches_vec), ]
row.names(logfc_original_df) <- NULL



# Choose one plasmid for each Entrez gene ID ------------------------------

logfc_4sg_df       <- ChooseOnePerEntrez(logfc_4sg_df)
logfc_CRISPRoff_df <- ChooseOnePerEntrez(logfc_CRISPRoff_df)
logfc_original_df  <- ChooseOnePerEntrez(logfc_original_df)

new_common_genes <- intersect(logfc_CRISPRoff_df[, "Entrez_ID"], logfc_4sg_df[, "Entrez_ID"])
common_genes <- intersect(logfc_original_df[, "Entrez_ID"], new_common_genes)

present_CRISPRoff_genes <- logfc_CRISPRoff_df[, "Entrez_ID"][!(is.na(logfc_CRISPRoff_df[, "Mean_log2FC"]))]
present_4sg_genes <- logfc_4sg_df[, "Entrez_ID"][!(is.na(logfc_4sg_df[, "Mean_log2FC"]))]

non_NA_common_genes <- intersect(common_genes, present_CRISPRoff_genes)
non_NA_common_genes <- intersect(non_NA_common_genes, present_4sg_genes)


new_common_genes <- intersect(logfc_CRISPRoff_df[, "Entrez_ID"], logfc_4sg_df[, "Entrez_ID"])

new_non_NA_common_genes <- intersect(new_common_genes, present_CRISPRoff_genes)
new_non_NA_common_genes <- intersect(new_non_NA_common_genes, present_4sg_genes)




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
                  gap_ratio = 1.7, side_gap = 0.6, point_size_factor = 1.2,
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



pdf(file.path(manuscript_dir, "SSMD - bar chart.pdf"),
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


scaling_factor <- 20
devEMF::emf(file.path(thesis_dir, "3C) SSMD - bar chart.emf"),
            width = 2.3 * scaling_factor, height = 2 * scaling_factor,
            emfPlus = FALSE, coordDPI = 3000
            )
old_par <- par(mar = c(3, 4, 2, 1),
               cex = 0.6 * scaling_factor * 0.95,
               lwd = 0.7 * scaling_factor * 0.95
               )
this_bars_vec <- abs(c(separation_original_mat["Robust SSMD", ],
                       separation_CRISPRoff_mat["Robust SSMD", ],
                       separation_4sg_mat["Robust SSMD", ]
                       ))
bar_positions <- CompareScreenBars(
  this_bars_vec,
  short_labels = TRUE,
  use_tcl = 0.35, y_axis_label_line = 1.9, y_axis_mgp = 0.525,
  gini_index = FALSE, y_axis_label = "SSMD*",
  y_upper_limit = 1.5, draw_grid = TRUE,
  y_axis_n = 3, bar_color = "#407ab5", stem_color = "#e0ebf5",
  gap_ratio = 1.7, side_gap = 0.6,
  bar_width = 0.4, lollipop = TRUE, point_size_factor = 1.25,
  group_labels_y = 1.75, bracket_color = "gray60"
)
dev.off()



# Plot the separation between E and NE genes for the manuscript -----------

this_points_vec <- this_bars_vec[3:6]

pdf(file.path(manuscript_dir, "bioRxiv v2 - Figure 6H - SSMD.pdf"),
    width = 1.1, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
ComparePoints(this_points_vec, left_gap = 0.55, right_gap = 0.45)
par(old_par)
dev.off()


## Export source data
ssmd_source_data_df <- data.frame(
  "Screen" = rep(c("CRISPRoff", "T.gonfio"), each = 2),
  "Replicate" = rep(1:2, times = 2),
  "Robust_SSMD" = this_points_vec
)

write.table(ssmd_source_data_df,
            file.path(manuscript_dir, "bioRxiv v2 - Figure 6H - SSMD - source data.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE
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
              width = PDF_width, height = PDF_height, emfPlus = FALSE, coordDPI = 6000
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
                    embed_PNG            = TRUE,
                    point_cex            = 0.6
                    )
  dev.off()
}



# Export combined ROC curves ----------------------------------------------

essential_entrezs     <- intersect(essentials_2020Q2_df[, "Entrez_ID"], non_NA_common_genes)
non_essential_entrezs <- intersect(non_essentials_2020Q2_df[, "Entrez_ID"], non_NA_common_genes)

# essential_fraction <- essential_df[["CRISPR_num_essential"]] /
#                       essential_df[["CRISPR_num_cell_lines"]]
# are_inclusive_ess <- essential_fraction > 0.95
# are_inclusive_noness <- (essential_fraction == 0) &
#                         (essential_df[, "CRISPR_mean_probability"] <= 0.03)
# inclusive_essentials <- intersect(essential_df[, "Entrez_ID"][are_inclusive_ess],
#                                   non_NA_common_genes
#                                   )
# inclusive_non_essentials <- intersect(essential_df[, "Entrez_ID"][are_inclusive_noness],
#                                       non_NA_common_genes
#                                      )

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

MultiLinesROC(individual_ROC_df_list)
MultiLinesROC(ROC_df_list)


MultiLinesROC(ROC_df_list, embed_PNG = TRUE, small_gap_size = 1.2,
              large_gap_multiplier = 1.3, use_lwd = 1.5, legend_lwd = 2.5,
              line_x_distance = -0.4, legend_order = c(3, 1, 2),
              x_axis_limits = c(-0.02, 1.02),
              )



devEMF::emf(file.path(thesis_dir, "3B) ROC curves.emf"),
            width = 2.72, height = 2, emfPlus = FALSE, coordDPI = 3000
            )
old_par <- par(mar = c(3, 4, 2, 7), cex = 0.6, lwd = 0.7)
MultiLinesROC(ROC_df_list, embed_PNG = FALSE, small_gap_size = 1.25, large_gap_multiplier = 1.5,
              transparency = FALSE, use_lwd = 2.5, legend_lwd = 2.5,
              line_x_distance = -0.4, legend_order = c(3, 1, 2), legend_inside = FALSE,
              length_in_lines = 0.55, lines_x_start = 0.9
              )
par(old_par)
dev.off()



scaling_factor <- 20
scaled_width <- 2.72 * scaling_factor
scaled_height <- 2 * scaling_factor
scaled_cex <- 0.6 * scaling_factor * 0.95
scaled_lwd <- 0.7 * scaling_factor * 0.95

use_mar <- c(3, 4, 2, 7)

plot_height <- (2 * scaling_factor) - (sum(use_mar[c(1, 3)]) * scaled_cex * 0.2)
plot_width <- (2.72 * scaling_factor) - (sum(use_mar[c(2, 4)]) * scaled_cex * 0.2)
width_diff <- plot_width - plot_height
use_mar[[4]] <- use_mar[[4]] + (width_diff / scaled_cex * 5)



custom_colors <- c("black", "#0664ef", "#da0b0b")
devEMF::emf(file.path(thesis_dir, "3B) ROC curves - only annotation.emf"),
            width = scaled_width, height = scaled_height,
            emfPlus = FALSE, coordDPI = 3000
            )
old_par <- par(mar = use_mar, cex = scaled_cex, lwd = scaled_lwd)
MultiLinesROC(ROC_df_list, small_gap_size = 1.25, large_gap_multiplier = 1.5,
              transparency = TRUE, legend_lwd = 2.5,
              line_x_distance = -0.4, legend_order = c(3, 2, 1), legend_inside = FALSE,
              length_in_lines = 0.55, lines_x_start = 0.9,
              use_colors = custom_colors, black_alpha = 0.6, colors_alpha = 0.7,
              only_annotation = TRUE
              )
par(old_par)
dev.off()


svglite::svglite(file.path(thesis_dir, "3B) ROC curves.svg"),
            width = scaled_width, height = scaled_height, bg = "transparent"
            )
old_par <- par(mar = use_mar, cex = scaled_cex, lwd = scaled_lwd)
MultiLinesROC(ROC_df_list,
              transparency = TRUE, use_lwd = 2.5,
              use_colors = custom_colors, black_alpha = 0.6, colors_alpha = 0.7,
              only_annotation = FALSE
              )
par(old_par)
dev.off()



# Export combined ROC curves for the manuscript ---------------------------

new_essential_entrezs     <- intersect(essentials_2020Q2_df[, "Entrez_ID"], new_non_NA_common_genes)
new_non_essential_entrezs <- intersect(non_essentials_2020Q2_df[, "Entrez_ID"], new_non_NA_common_genes)

new_ROC_CRISPRoff_input_df <- ROCInputDf(logfc_CRISPRoff_df, new_essential_entrezs, new_non_essential_entrezs)
new_ROC_4sg_input_df       <- ROCInputDf(logfc_4sg_df, new_essential_entrezs, new_non_essential_entrezs)

new_common_ROC_CRISPRoff_df <- ROCDfForColumn(ROC_CRISPRoff_input_df, "Mean_log2FC")
new_common_ROC_4sg_df       <- ROCDfForColumn(ROC_4sg_input_df, "Mean_log2FC")


custom_colors <- c("black", "#0664ef")

pdf(file.path(manuscript_dir, "bioRxiv v2 - Figure 6I - ROC curves.pdf"),
    width = 2, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.3))
MultiLinesROC(list(new_common_ROC_CRISPRoff_df, new_common_ROC_4sg_df),
              transparency = TRUE, use_lwd = 2.2,
              use_colors = custom_colors, black_alpha = 0.6, colors_alpha = 0.7,
              y_label_line = 2.1, x_label_line = 1.7,
              middle_line = TRUE,
              legend_inside = TRUE, long_labels = FALSE,
              lines_x_start = -0.225, lines_y_start = 0.8,
              large_gap_multiplier = 1.2, small_gap_size = 0.88,
              text_cex = 0.9, AUC_num_digits = 3
              )
par(old_par)
dev.off()



# Export mean gamma violin plots for the manuscript -----------------------

rep_list <- c(split(new_common_ROC_CRISPRoff_df[, "Mean_log2FC"], !(new_common_ROC_CRISPRoff_df[, "Is_essential"])),
              split(new_common_ROC_4sg_df[, "Mean_log2FC"], !(new_common_ROC_4sg_df[, "Is_essential"]))
              )

pdf(file.path(manuscript_dir, "bioRxiv v2 - Figure 6G - violin plots.pdf"),
    width = 1.9, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
MeanSwarms(rep_list, show_truncation = FALSE)
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
              width = 66, height = 42, emfPlus = FALSE, coordDPI = 3000
              )
  par(cex = 10.6, lwd = 10.5)
  BidirectionalViolins(bidirect_df       = bidirectional_df,
                       logfc_df          = get(datasets_vec[[i]]),
                       max_distance      = 20000,
                       num_controls      = 30L,
                       compare_across    = FALSE,
                       point_cex         = 0.8,
                       draw_grid         = TRUE,
                       grid_color        = "gray88",
                       zero_lty          = "solid",
                       zero_color        = "gray68",
                       show_x_axis       = FALSE,
                       quantiles_lty     = c("dotted", "dashed", "dotted"), # compatibility with emf device
                       annotation_cex    = 1,
                       draw_groups_n     = FALSE,
                       swarm_method      = "compactswarm",
                       wex               = 0.86,
                       y_limits          = c(-0.6, 0.2),
                       lower_bound       = -0.6,
                       label_points      = TRUE,
                       gap_ratio         = 1.2,
                       side_gap          = 0.4,
                       right_gap         = 0.3,
                       use_spacing       = 0.85,
                       y_start_adj       = -0.8,
                       y_line_adj        = 0.1,
                       p_value_cex       = 0.8,
                       draw_group_labels = i == 3,
                       show_title        = names(datasets_vec)[[i]]
                       )
  dev.off()
}




# Try stuff ---------------------------------------------------------------

essential_fraction <- essential_df[["CRISPR_num_essential"]] /
                      essential_df[["CRISPR_num_cell_lines"]]

logfc_df <- common_logfc_CRISPRoff_df

are_essential <- logfc_df[, "Entrez_ID"] %in% essentials_2020Q2_df[, "Entrez_ID"]
are_non_essential <- logfc_df[, "Entrez_ID"] %in% non_essentials_2020Q2_df[, "Entrez_ID"]
are_unspecific <- logfc_df[, "Combined_specificity"] < 0.05

table(are_non_essential, are_unspecific)
beeswarm(split(logfc_df[are_non_essential, "Mean_log2FC"], are_unspecific[are_non_essential]), cex = 0.5, pch = 16)


pairs_df <- ScatterInputDf(common_logfc_CRISPRoff_df, common_logfc_4sg_df)
pairs_df <- pairs_df[pairs_df[, "Essentiality"] != "neither", ]
are_E <- pairs_df[, "Essentiality"] == "essential"
are_NE <- pairs_df[, "Essentiality"] == "non-essential"
pairs_df[, "Rep1_quantile"] <- rank(pairs_df[, "Rep1_data"]) / nrow(pairs_df)
pairs_df[, "Rep2_quantile"] <- rank(pairs_df[, "Rep2_data"]) / nrow(pairs_df)
pairs_df[, "Delta_quantile"] <- pairs_df[, "Rep1_quantile"] - pairs_df[, "Rep2_quantile"]
pairs_df[, "Delta_specificity"] <- pairs_df[, "Rep1_min_spec"] - pairs_df[, "Rep2_min_spec"]


MakeEmptyPlot(); abline(a = 0, b = 1); abline(v = 0.5, h = 0.5)
points(pairs_df[are_NE, "Rep1_quantile"], pairs_df[are_NE, "Rep2_quantile"])
MakeEmptyPlot(); abline(a = 0, b = 1); abline(v = 0.5, h = 0.5)
points(pairs_df[are_E, "Rep1_quantile"], pairs_df[are_E, "Rep2_quantile"])



MakeEmptyPlot()
abline(a = 0, b = 1)
abline(v = 0.5, h = 0.5, col = "gray50")
are_E <- pairs_df[, "Essentiality"] == "essential"
are_NE <- pairs_df[, "Essentiality"] == "non-essential"
points(pairs_df[are_E, "Rep1_quantile"], pairs_df[are_E, "Rep2_quantile"], pch = 16, cex = 0.5)
MakeEmptyPlot()
abline(a = 0, b = 1)
points(pairs_df[are_NE, "Rep1_quantile"], pairs_df[are_NE, "Rep2_quantile"], pch = 16, cex = 0.5)
abline(v = 0.5, h = 0.5, col = "gray50")

beeswarm(list(pairs_df[are_E, "Delta_quantile"], pairs_df[are_NE, "Delta_quantile"]),
         cex = 0.5, pch = 16
         )
boxplot(list(pairs_df[are_E, "Delta_quantile"], pairs_df[are_NE, "Delta_quantile"]), add = TRUE)

non_essential_entrezs <- essential_df[(essential_fraction == 0) &
                                      (essential_df[, "CRISPR_mean_probability"] <= 0.03), "Entrez_ID"]
are_non_essential <- logfc_df[, "Entrez_ID"] %in% non_essential_entrezs
table(are_non_essential, are_unspecific)
split_list <- split(logfc_df[are_non_essential, "Mean_log2FC"], are_unspecific[are_non_essential])
split_list <- lapply(split_list, function(x) x[is.finite(x)])
beeswarm(split_list, cex = 0.5, pch = 16)




# Analyze the effect of proximity to essential genes ----------------------

## Compute distances to the nearest essential gene
distances_original_data_df    <- DistanceToEssential(common_logfc_original_df, TSS_df,  use_data_for_essential = TRUE)
distances_CRISPRoff_data_df   <- DistanceToEssential(common_logfc_CRISPRoff_df, TSS_df, use_data_for_essential = TRUE)
distances_4sg_data_df         <- DistanceToEssential(common_logfc_4sg_df, TSS_df,       use_data_for_essential = TRUE)

distances_original_depmap_df  <- DistanceToEssential(common_logfc_original_df, TSS_df,  use_data_for_essential = FALSE)
distances_CRISPRoff_depmap_df <- DistanceToEssential(common_logfc_CRISPRoff_df, TSS_df, use_data_for_essential = FALSE)
distances_4sg_depmap_df       <- DistanceToEssential(common_logfc_4sg_df, TSS_df,       use_data_for_essential = FALSE)


## Compare the range of distances between datasets and
## different definitions of gene essentiality
range(abs(distances_original_data_df[, "Distance"]), na.rm = TRUE)
range(abs(distances_CRISPRoff_data_df[, "Distance"]), na.rm = TRUE)
range(abs(distances_4sg_data_df[, "Distance"]), na.rm = TRUE)
range(abs(distances_original_depmap_df[, "Distance"]), na.rm = TRUE)
range(abs(distances_CRISPRoff_depmap_df[, "Distance"]), na.rm = TRUE)
range(abs(distances_4sg_depmap_df[, "Distance"]), na.rm = TRUE)


## Visualize the phenotype score versus the distance to
## the nearest essential gene
DistanceScatter(distances_original_data_df,  x_limit = 7.625)
DistanceScatter(distances_CRISPRoff_data_df, x_limit = 7.625)
DistanceScatter(distances_4sg_data_df,       x_limit = 7.625)
DistanceScatter(distances_original_depmap_df)
DistanceScatter(distances_CRISPRoff_depmap_df)
DistanceScatter(distances_4sg_depmap_df)

DistanceScatter(distances_original_depmap_df, x_grid = TRUE)


## Perform visualizations for subsets of genes
DistanceScatter(distances_original_depmap_df[distances_original_depmap_df[, "DepMap_essentiality"] %in% "Essential", ],
                use_hue = "#6c4cbd"
                )
DistanceScatter(distances_original_depmap_df[distances_original_depmap_df[, "DepMap_essentiality"] %in% "Non-essential", ],
                use_hue = "#3c9f5d"
                )
DistanceScatter(distances_4sg_depmap_df[distances_4sg_depmap_df[, "DepMap_essentiality"] %in% "Essential", ],
                use_hue = "#6c4cbd"
                )
DistanceScatter(distances_4sg_depmap_df[distances_4sg_depmap_df[, "DepMap_essentiality"] %in% "Non-essential", ],
                use_hue = "#3c9f5d"
                )

are_low <- distances_original_data_df[, "Mean_log2FC"] < -2
DistanceScatter(distances_original_data_df[, "Distance"][are_low],
                distances_original_data_df[, "Nearest_mean_log2FC"][are_low],
                x_limit = 7.625
                )
DistanceScatter(distances_original_data_df[, "Distance"][are_low],
                distances_original_data_df[, "Mean_log2FC"][are_low],
                x_limit = 7.625
                )


## Export plots for a multi-figure panel

datasets_vec <- c(
  "Data re-analysis"  = "distances_original_depmap_df",
  "CRISPRoff library" = "distances_CRISPRoff_depmap_df",
  "T.gonfio library"  = "distances_4sg_depmap_df"
)

for (i in seq_along(datasets_vec)) {
  file_name <- paste0("4) Scatter plots of distance - ",
                      tolower(as.character(as.roman(i))), ") ",
                      names(datasets_vec)[[i]]
                      )
  devEMF::emf(file.path(thesis_dir, paste0(file_name, ".emf")),
              width = 66, height = 42, emfPlus = FALSE, coordDPI = 3000
              )
  par(mai = c(12.72, 8.48, 6.36, 3.18), cex = 10.6, lwd = 10.5)
  DistanceScatter(get(datasets_vec[[i]]),
                  point_cex = 0.6, #point_color = "#2a4d98", median_color = "black",
                  embed_PNG = TRUE, png_res = 50, png_padding = 0.7,
                  grid_lwd = 0.6, median_lwd = 1.4, x_grid = TRUE,
                  show_x_label = i == 3, grid_color = "gray88",
                  grid_highlight = "gray68"
                  )

  mtext(VerticalAdjust(names(datasets_vec)[[i]]),
        at = grconvertX(0.015, from = "npc", to = "user"),
        adj = 0, line = -0.05, cex = par("cex")
        )
  dev.off()
}

use_x_limit <- max(log10(abs(distances_4sg_depmap_df[, "Distance"])), na.rm = TRUE)
use_x_limit <- use_x_limit + ((use_x_limit - 1) * 0.015)

file_name <- paste0("4) Scatter plots of distance - ",
                    "iv) essential genes"
                    )
devEMF::emf(file.path(thesis_dir, paste0(file_name, ".emf")),
            width = 66, height = 42, emfPlus = FALSE, coordDPI = 3000
            )
par(mai = c(12.72, 8.48, 6.36, 3.18), cex = 10.6, lwd = 10.5)
DistanceScatter(distances_4sg_depmap_df[distances_4sg_depmap_df[, "DepMap_essentiality"] %in% "Essential", ],
                use_hue = "#6c4cbd",
                point_cex = 0.6, #point_color = "#2a4d98", median_color = "black",
                embed_PNG = TRUE, png_res = 50, png_padding = 0.7,
                grid_lwd = 0.6, median_lwd = 1.4, x_grid = TRUE,
                grid_color = "gray88",
                grid_highlight = "gray68",
                x_limit = use_x_limit
                )
mtext(VerticalAdjust("T.gonfio library: essential genes"),
      at = grconvertX(0.015, from = "npc", to = "user"),
      adj = 0, line = -0.05, cex = par("cex")
      )
dev.off()


file_name <- paste0("4) Scatter plots of distance - ",
                    "v) non-essential genes"
                    )
devEMF::emf(file.path(thesis_dir, paste0(file_name, ".emf")),
            width = 66, height = 42, emfPlus = FALSE, coordDPI = 3000
            )
par(mai = c(12.72, 8.48, 6.36, 3.18), cex = 10.6, lwd = 10.5)
DistanceScatter(distances_4sg_depmap_df[distances_4sg_depmap_df[, "DepMap_essentiality"] %in% "Non-essential", ],
                use_hue = "#3c9f5d",
                point_cex = 0.6, #point_color = "#2a4d98", median_color = "black",
                embed_PNG = TRUE, png_res = 50, png_padding = 0.7,
                grid_lwd = 0.6, median_lwd = 1.4, x_grid = TRUE,
                grid_color = "gray88",
                grid_highlight = "gray68",
                x_limit = use_x_limit
                )
mtext(VerticalAdjust("T.gonfio library: non-essential genes"),
      at = grconvertX(0.015, from = "npc", to = "user"),
      adj = 0, line = -0.05, cex = par("cex")
      )
dev.off()



