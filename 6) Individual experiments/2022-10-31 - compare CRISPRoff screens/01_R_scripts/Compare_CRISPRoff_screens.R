### 2022-10-31


# Load packages and source code -------------------------------------------

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
                              axis_cex_factor      = 1 / 0.7,
                              y_axis_label_line    = 2.2
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
  point_colors <- adjustcolor(point_colors, alpha.f = 0.4)
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
          bg       = "transparent"
          )
      par(lwd = current_par[["lwd"]])
      par(cex = current_par[["cex"]])
      par(mar = rep(0, 4))
      MakeEmptyPlot(xy_lim, xy_lim)
    }

    ## Set up plot canvas
    abline(v = 0, h = 0, col = "gray85", lend = "butt")
    abline(a = 0, b = 1, col = "gray85", lend = "butt")
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
    mtext(VerticalAdjust(top_labels[[i]]), line = 0.2, cex = par("cex"))
    mtext(VerticalAdjust(axis_labels_list[[1]]), side = 1, line = 1.8, cex = par("cex"))
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
      axis(j,
           at       = tick_locations,
           labels   = use_tick_labels,
           mgp      = c(3, if (j == 1) 0.35 else 0.5, 0),
           tcl      = -0.3,
           las      = 1,
           lwd      = par("lwd"),
           cex.axis = par("cex") * axis_cex_factor
           )
    }
    box()

    MakeEmptyPlot()
  }

  par(old_par)
  layout(1)

  return(invisible(NULL))
}



CompareGinis <- function(gini_vec,
                         use_title = "Plasmid count heterogeneity at baseline",
                         short_labels = FALSE,
                         use_tcl = 0.375,
                         y_axis_label_line = 2.25,
                         y_axis_mgp = 0.55
                         ) {
  stopifnot(length(gini_vec) == 4)

  ## Determine bar positions
  bar_positions <- RepositionByGroups(c(1, 1, 2, 2), gap_ratio = 1.6)
  bar_width <- 0.45
  num_bars <- length(bar_positions)
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - 0.5) - (num_bars * 0.04),
                    (max(bar_positions) + 0.5) + (num_bars * 0.04)
                     )

  ## Draw the bars
  MakeEmptyPlot(x_limits = group_limits, y_limits = c(0, 0.5))
  PlotBarplotMat(t(gini_vec),
                 colors_vec    = brewer.pal(9, "Blues")[[8]],
                 positions_vec = bar_positions,
                 bar_width     = bar_width
                 )


  ## Draw the y axis
  axis(2,
       las = 2,
       mgp = c(3, y_axis_mgp, 0),
       tcl = -(use_tcl),
       lwd = par("lwd")
       )
  mtext(VerticalAdjust("Gini index"),
        side = 2,
        line = y_axis_label_line,
        cex  = par("cex")
        )

  ## Annotate bars
  mtext(text = rep(paste0("R", 1:2), 2),
        at = bar_positions, side = 1, line = 0.3, cex = par("cex")
        )
  segments(x0  = bar_positions[c(1, 3)] - 0.2,
           x1  = bar_positions[c(2, 4)] + 0.2,
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.5), from = "lines", to = "user")),
           col = "gray50",
           xpd = NA
           )
  mtext(text = paste0(c("CRISPRoff", "T.gonfio"), if (short_labels) "" else " library"),
        at   = c(mean(bar_positions[1:2]), mean(bar_positions[3:4])),
        side = 1,
        line = 1.85,
        cex  = par("cex")
        )
  mtext(use_title, side = 3, line = 1, cex = par("cex"))
  box(bty = "l")
  return(invisible(NULL))
}




# Compare Gini indices at baseline ----------------------------------------

CompareGinis(c(gini_indices_CRISPRoff[c(1, 2)], gini_indices_4sg[c(1, 2)]))

pdf(file.path(manuscript_dir, "Gini index comparison.pdf"),
    width = 2, height = 1.75
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
CompareGinis(c(gini_indices_CRISPRoff[c(1, 2)], gini_indices_4sg[c(1, 2)]),
             short_labels = TRUE, use_title = "",
             use_tcl = 0.35, y_axis_label_line = 2.2, y_axis_mgp = 0.525
             )
title("Count heterogeneity", cex.main = 1, font.main = 1)
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

PDF_height <- 1.15 / sum((use_heights[[2]] / sum(use_heights)))
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
                  axis_cex_factor   = 1 / 0.45,
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
                  axis_cex_factor      = 1 / 0.45,
                  layout_widths        = use_widths,
                  layout_heights       = use_heights,
                  embed_PNG            = TRUE
                  )
dev.off()




# Examine bidirectional promoters -----------------------------------------

common_genes <- intersect(logfc_CRISPRoff_df[, "Entrez_ID"], logfc_4sg_df[, "Entrez_ID"])

common_logfc_CRISPRoff_df <- logfc_CRISPRoff_df[logfc_CRISPRoff_df[, "Entrez_ID"] %in% common_genes, ]
row.names(common_logfc_CRISPRoff_df) <- NULL

common_logfc_4sg_df <- logfc_4sg_df[logfc_4sg_df[, "Entrez_ID"] %in% common_genes, ]
row.names(common_logfc_4sg_df) <- NULL

BidirectionalViolins(bidirectional_df, common_logfc_CRISPRoff_df, max_distance = 20000,
                     num_controls = 30L, compare_across = FALSE
                     )
BidirectionalViolins(bidirectional_df, common_logfc_4sg_df, max_distance = 20000,
                     num_controls = 30L, compare_across = FALSE
                     )





