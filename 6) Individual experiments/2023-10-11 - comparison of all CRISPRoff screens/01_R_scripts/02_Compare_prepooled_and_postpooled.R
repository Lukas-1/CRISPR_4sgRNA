### 2023-10-17


# Load packages and source code -------------------------------------------

library("pBrackets")

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(R_functions_dir, "07_comparing_CRISPRoff_screens.R"))



# Define paths ------------------------------------------------------------

first_rdata_dir <- file.path(first_illumina_trial_dir, "03_R_objects")
pool_rdata_dir  <- file.path(experiments_directory,
                             "2023-09-28 - prepooled vs postpooled - Illumina",
                             "03_R_objects"
                             )
project_dir     <- file.path(experiments_directory, "2023-10-11 - comparison of all CRISPRoff screens")
output_dir      <- file.path(project_dir, "02_output")
figures_dir     <- output_dir



# Load data ---------------------------------------------------------------

load(file.path(pool_rdata_dir, "10_create_figures_from_count_data.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))



# Define functions --------------------------------------------------------

BracketsCoords <- function(x1, y1, x2, y2, h = NULL, ticks = 0.5, curvature = 0.5, type = 1) {
  ## copied from pBrackets::brackets
    if (!is.numeric(curvature))
        stop("curvature must be numeric")
    if (!is.numeric(type))
        stop("type must be a integer, 1 to 5")
    if (curvature < 0)
        curvature <- 0
    if (curvature > 1)
        curvature <- 1
    if (length(ticks) == 1)
        if (is.na(ticks))
            ticks <- NULL
    if (!is.numeric(ticks) & !is.null(ticks))
        stop("ticks must be numeric or NULL")
    if (length(ticks) > 1) {
        if (anyDuplicated(abs(ticks)))
            stop("duplicated ticks")
    }
    xm <- mean(c(x1, x2))
    ym <- mean(c(y1, y2))
    coords <- par("usr")
    xr <- abs(diff(coords[1:2]))
    yr <- abs(diff(coords[3:4]))
    if (is.null(h))
        h <- sqrt(xr * yr)/20
    dinps <- par("pin")
    xtox <- dinps[2]/dinps[1]
    rato <- xr/yr
    mysig <- sign(h)
    vx = (x2 - x1)/rato/xtox
    vy = (y2 - y1) * rato * xtox
    len = sqrt(vx^2 + vy^2)
    ux = vy/len
    uy = vx/len
    x3 = xm + h * ux
    y3 = ym - h * uy
    ret <- pBrackets:::a_get_relative(c(xm, x3), c(ym, y3))
    h <- mysig * dist(cbind(ret[, 1], ret[, 2]))
    rvalues <- pBrackets:::a_get_relative(c(x1, x2), c(y1, y2))
    x1 <- rvalues[1, 1]
    x2 <- rvalues[2, 1]
    y1 <- rvalues[1, 2]
    y2 <- rvalues[2, 2]
    brackets <- pBrackets:::a_cb_brackets(phi = curvature, ticks = ticks, type = type)
    x <- brackets[1, ] * dist(cbind(c(x1, x2), c(y1, y2)))
    y <- brackets[2, ] * h
    rout <- pBrackets:::a_rotate(cbind(x, y), atan2(y2 - y1, x2 - x1))
    xout <- rout[, 1] + x1
    yout <- rout[, 2] + y1
    native <- pBrackets:::a_get_native(xout, yout)
    return(native)
}



# Define gene selections --------------------------------------------------

ess_genes_blomen_hart    <- essentials_2020Q2_df[, "Entrez_ID"]
noness_genes_blomen_hart <- non_essentials_2020Q2_df[, "Entrez_ID"]

ess_genes_depmap <- essential_df[essential_df[, "Three_categories"] %in% "Essential", "Entrez_ID"]
noness_genes_depmap <- essential_df[essential_df[, "Three_categories"] %in% "Non-essential", "Entrez_ID"]

ess_genes_residual <- setdiff(ess_genes_depmap, ess_genes_blomen_hart)

common_plasmids_logfc_df_list <- CommonPlasmidsRocDfList(list(logfc_prepooled_df, logfc_postpooled_df))

blomen_hart_ROC_df_list <- LogFcDfListToRocDfList(
  common_plasmids_logfc_df_list,
  essential_entrezs = ess_genes_blomen_hart,
  non_essential_entrezs = noness_genes_blomen_hart
)

EssResidual_NEBlomenHart_ROC_df_list <- LogFcDfListToRocDfList(
  common_plasmids_logfc_df_list,
  essential_entrezs = ess_genes_residual,
  non_essential_entrezs = noness_genes_blomen_hart
)

EssDepMap_NEBlomenHart_ROC_df_list <- LogFcDfListToRocDfList(
  common_plasmids_logfc_df_list,
  essential_entrezs = ess_genes_depmap,
  non_essential_entrezs = noness_genes_blomen_hart
)

demap_ROC_df_list <- LogFcDfListToRocDfList(
  common_plasmids_logfc_df_list,
  essential_entrezs = ess_genes_depmap,
  non_essential_entrezs = noness_genes_depmap
)



# Export mean gamma violin plots for the manuscript -----------------------

rep_list <- c(split(blomen_hart_ROC_df_list[[1]][, "Mean_log2FC"], !(blomen_hart_ROC_df_list[[1]][, "Is_essential"])),
              split(blomen_hart_ROC_df_list[[2]][, "Mean_log2FC"], !(blomen_hart_ROC_df_list[[2]][, "Is_essential"]))
              )

pdf(file.path(output_dir, "Manuscript - Figure 6K - violin plots.pdf"),
          width = 1.9, height = 2
          )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
violin_positions <- MeanSwarms(rep_list, group_labels = c("prepool", "postpool"),
                               show_truncation = FALSE
                               )
par(old_par)
dev.off()



CustomBracketsInBackground <- function(group_positions, quantiles_list, use_lwd = 1, x_distance = 1.5) {

  x_pos <- group_positions[[1]] - diff(grconvertX(c(0, x_distance), from = "lines", to = "user"))
  coords_mat <- BracketsCoords(x1 = x_pos,
                               x2 = x_pos,
                               y1 = quantiles_list[[1]][["25%"]],
                               y2 = quantiles_list[[1]][["75%"]]
                               )
  polygon(x      = c(coords_mat[, 1], group_positions[[1]], group_positions[[1]]),
          y      = c(coords_mat[, 2], coords_mat[[length(coords_mat)]], coords_mat[, 2][[1]]),
          col    = "gray96",
          border = NA
          )
  segments(x0  = coords_mat[1, 1],
           x1  = group_positions[[1]],
           y0  = c(quantiles_list[[1]][["25%"]], quantiles_list[[1]][["75%"]]),
           lty = "21",
           lwd = par("lwd") * 0.7,
           col = "gray80"
           )
  lines(coords_mat, col = "black")

  x_pos <- group_positions[[2]] + diff(grconvertX(c(0, x_distance - 0.025), from = "lines", to = "user"))
  coords_mat <- BracketsCoords(x1 = x_pos,
                               x2 = x_pos,
                               y1 = quantiles_list[[2]][["75%"]],
                               y2 = quantiles_list[[2]][["25%"]]
                               )
  polygon(x      = c(coords_mat[, 1], group_positions[[2]], group_positions[[2]]),
          y      = c(coords_mat[, 2], coords_mat[[length(coords_mat)]], coords_mat[, 2][[1]]),
          col    = "gray96",
          border = NA
          )
  segments(x0  = coords_mat[1, 1],
           x1  = group_positions[[2]],
           y0  = c(quantiles_list[[2]][["25%"]], quantiles_list[[2]][["75%"]]),
           lty = "21",
           lwd = par("lwd") * 0.7,
           col = "gray80"
           )
  lines(coords_mat, col = "black")
}



IQR_vec <- vapply(blomen_hart_ROC_df_list, function(x) {
  IQR(x[, "Mean_log2FC"][x[, "Is_essential"]] / 10)
}, numeric(1))

two_quantiles_list <- lapply(blomen_hart_ROC_df_list, function(x) {
  quantile(x[, "Mean_log2FC"][x[, "Is_essential"]] / 10)
})


x_positions <- MeanSwarms(rep_list, group_labels = c("prepool", "postpool"),
                          show_truncation = FALSE, x_limits = c(-1.35, 4.6),
                          x_positions = violin_positions - c(rep(0.85, 2), rep(0, 2))
                          )

num_digits_IQR <- 3L


pdf(file.path(output_dir, "Violin plot variant - IQR indicated.pdf"),
    width = 2.5, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
x_positions <- MeanSwarms(rep_list, group_labels = c("prepool", "postpool"),
                          show_truncation = FALSE, x_limits = c(-1.45, 4.6),
                          x_positions = violin_positions - c(rep(0.85, 2), rep(0, 2)),
                          GridFunction = CustomBracketsInBackground(x_positions, two_quantiles_list)
                          )
old_lheight <- par("lheight" = 0.925)
text(x      = x_positions[[1]] - diff(grconvertX(c(0, 2.85), from = "lines", to = "user")),
     y      = mean(c(two_quantiles_list[[1]][c("25%", "75%")])),
     labels = paste0("IQR\n", format(round(IQR_vec[[1]], digits = num_digits_IQR), nsmall = num_digits_IQR)),
     cex    = 0.9,
     adj    = c(0.5, 0.5)
     )
text(x      = x_positions[[2]] + diff(grconvertX(c(0, 2.9), from = "lines", to = "user")),
     y      = mean(c(two_quantiles_list[[2]][c("25%", "75%")])),
     labels = paste0("IQR\n", format(round(IQR_vec[[2]], digits = num_digits_IQR), nsmall = num_digits_IQR)),
     cex    = 0.9,
     adj    = c(0.5, 0.5)
     )
par(old_lheight)
par(old_par)
dev.off()



rep_list <- c(split(EssResidual_NEBlomenHart_ROC_df_list[[1]][["Mean_log2FC"]], !(EssResidual_NEBlomenHart_ROC_df_list[[1]][, "Is_essential"])),
              split(EssResidual_NEBlomenHart_ROC_df_list[[2]][["Mean_log2FC"]], !(EssResidual_NEBlomenHart_ROC_df_list[[2]][, "Is_essential"]))
              )

pdf(file.path(output_dir, "Violin plot variant - residual essential genes.pdf"),
    width = 1.9, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
MeanSwarms(rep_list, group_labels = c("prepool", "postpool"),
           GridFunction = rect(xleft = 1, xright = 2, ybottom = -0.4, ytop = 0, border = NA, col = "red")
           )
par(old_par)
dev.off()
rm(rep_list)



# Create scatter plots for the manuscript ---------------------------------

scatter_df <- ScatterInputDf(logfc_prepooled_df, logfc_postpooled_df)
scatter_df[, "Rep1_data"] <- scatter_df[, "Rep1_data"] / 10
scatter_df[, "Rep2_data"] <- scatter_df[, "Rep2_data"] / 10


ReplicateScatterPlot(scatter_df,
                     highlight_NT = FALSE,
                     axis_labels_list = list(expression("Prepool" ~ "(" * gamma * ")"),
                                             expression("Postpool" ~ "(" * gamma * ")")
                                             ),
                     lower_bound = -0.6, upper_bound = 0.6
                     )

use_cex <- 0.6
use_mai <- c(0.7, 0.8, 0.38 / use_cex, 1.4)
base_height <- 1.2

pdf(file.path(output_dir, "Manuscript - Figure 6L - scatter plot.pdf"),
    width  = base_height + (sum(use_mai[c(2, 4)] * use_cex)),
    height = base_height + (sum(use_mai[c(1, 3)]) * use_cex)
    )
old_par <- par(cex = 0.6, lwd = 0.7)

ReplicateScatterPlot(scatter_df,
                     axis_labels_list     = list(expression("Prepool" ~ "(" * gamma * ")"),
                                                 expression("Postpool" ~ "(" * gamma * ")")
                                                 ),
                     highlight_NT         = FALSE,
                     lower_bound          = -0.6,
                     upper_bound          = 0.2,
                     show_axis_truncation = FALSE,
                     axis_ticks_pretty_n  = 5,
                     use_mar              = use_mai * 5,
                     embed_PNG            = TRUE,
                     x_axis_label_line    = 1.8,
                     y_axis_label_line    = 2.1,
                     use_tcl              = 0.3,
                     x_axis_mgp           = 0.35,
                     y_axis_mgp           = 0.5,
                     point_cex            = 0.45,
                     legend_lines_x_start = 0.65,
                     legend_point_x_start = 0.05,
                     capitalize_legend    = FALSE,
                     break_lines          = TRUE,
                     axis_line_color      = "gray80",
                     small_gap_size       = 1.15,
                     large_gap_multiplier = 1.5
                     )
dev.off()



# Create a custom scatter plot --------------------------------------------

use_vec_list <- lapply(1:2, function(x) {
  logfc_df <- common_plasmids_logfc_df_list[[x]]
  results_vec <- logfc_df[, "Mean_log2FC"][logfc_df[, "Entrez_ID"] %in% ess_genes_blomen_hart]
  results_vec <- results_vec / 10
  results_vec[results_vec > 0.2] <- 0.2
  results_vec[results_vec < -0.6] <- -0.6
  return(results_vec)
})


pdf(file.path(output_dir, "Scatter plot - only essential genes.pdf"),
    width  = 2.1,
    height = 2.1
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.3))
RegressionScatter(use_vec_list[[1]], use_vec_list[[2]], same_limits = TRUE,
                  points_color = "#4b357e",
                  band_color = brewer.pal(9, "Purples")[[3]],
                  line_color = "#3d2d67",
                  confint_level = 0.99, show_axes = FALSE
                  )

tick_locations <- seq(-0.6, 0, by = 0.2)
axis(1, at = tick_locations, tcl = -0.3, mgp = c(3, 0.35, 0), lwd = par("lwd"))
axis(2, at = tick_locations, tcl = -0.3, mgp = c(3, 0.5, 0), las = 1, lwd = par("lwd"))
mtext(VerticalAdjust(expression("Prepool" ~ "(" * gamma * ")")),
      side = 1, line = 1.8, cex = par("cex")
      )
mtext(VerticalAdjust(expression("Postpool" ~ "(" * gamma * ")")),
      side = 2, line = 2.1, cex = par("cex")
      )
dev.off()



# Draw a scatter plot with density plots on the sides ---------------------

pdf(file.path(output_dir, "Scatter plot - with side density plots.pdf"),
    width  = 2.15,
    height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.45))
DensityTrapezoidScatter(use_vec_list)
dev.off()



# Plot the separation between E and NE genes for the manuscript -----------

this_points_vec <- c(separation_prepooled_mat["Robust SSMD", ],
                     separation_postpooled_mat["Robust SSMD", ]
                     )
this_points_vec <- abs(this_points_vec)

pdf(file.path(output_dir, "Manuscript - Figure 6M - SSMD.pdf"),
    width = 1.1, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
ComparePoints(this_points_vec, left_gap = 0.55, right_gap = 0.45,
              y_upper_limit = 2.5, group_labels = c(NA, "prepool", "postpool")
              )
par(old_par)
dev.off()



# Export combined ROC curves for the manuscript ---------------------------

manuscript_ROC_args <- list(
  ROC_df_list = blomen_hart_ROC_df_list,
  use_colors = c("#0664ef", "#da0b0b"),
  black_alpha = 0.6, colors_alpha = 0.7,
  y_label_line = 2.1, x_label_line = 1.7,
  middle_line = TRUE,
  legend_inside = TRUE, long_labels = FALSE,
  lines_x_start = -0.225, lines_y_start = 0.8,
  large_gap_multiplier = 1.2, small_gap_size = 0.88,
  legend_vec = c(NA, "prepool", "postpool"), text_cex = 0.9,
  AUC_num_digits = 3
)



pdf(file.path(output_dir, "Manuscript - Figure 6N - ROC curves.pdf"),
    width = 2, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.3))
do.call(MultiLinesROC, manuscript_ROC_args)
par(old_par)
dev.off()



TallyEssNoness <- function(input_df) {
  num_essential <- sum(input_df[, "Is_essential"])
  num_nonessential <- sum(!(input_df[, "Is_essential"]))
  result_text <- paste0(num_essential, " essential and ", num_nonessential, " non-essential genes")
  return(result_text)
}


variant_list_list <- list(
  list(file_suffix = " 1) standard (Blomen & Hart genes)",
       df_list     = blomen_hart_ROC_df_list,
       sub_text    = paste0("Intersection of Blomen and Hart, as in the CRISPRoff study by\nNu\u00f1ez et al. (", TallyEssNoness(blomen_hart_ROC_df_list[[1]]), ")")
       ),
  list(file_suffix = " 2) residual essential genes",
       df_list     = EssResidual_NEBlomenHart_ROC_df_list,
       sub_text    = paste0("Remaining essential genes from DepMap that are not part\nof Blomen and Hart; non-essential genes from Blomen and Hart\n(", TallyEssNoness(EssResidual_NEBlomenHart_ROC_df_list[[1]]), ")")
       ),
  list(file_suffix = " 3) DepMap essential genes",
       df_list     = EssDepMap_NEBlomenHart_ROC_df_list,
       sub_text    = paste0("Essential genes from DepMap; essential genes from Blomen and\nHart (", TallyEssNoness(EssDepMap_NEBlomenHart_ROC_df_list[[1]]), ")")
       ),
  list(file_suffix = " 4) DepMap essential and non-essential genes",
       df_list     = demap_ROC_df_list,
       sub_text    = paste0("Essential and non-essential genes defined by DepMap\n(", TallyEssNoness(demap_ROC_df_list[[1]]), ")")
       )
)



QuickTitle <- function(title_text) {
  title(title_text, font.main = 1, cex.main = 1)
}



for (variant_list in variant_list_list) {
  DrawSubtext <- function() {
    old_lheight <- par("lheight" = 1.15)
    mtext(variant_list[["sub_text"]], side = 1, cex = par("cex") * 0.6,
          line = 3.25, adj = 0, padj = 1,
          at = par("usr")[[1]] - diff(grconvertX(c(0, 2), from = "lines", to = "user"))
          )
    par(old_lheight)
  }
  manuscript_ROC_args[["ROC_df_list"]] <- variant_list[["df_list"]]
  pdf(file.path(output_dir, paste0("ROC curve variants - ", variant_list[["file_suffix"]], ".pdf")),
      width = 2, height = 2.35
      )
  old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.77, 0.5, 0.38, 0.3))

  do.call(MultiLinesROC, manuscript_ROC_args)
  QuickTitle("ROC curve"); DrawSubtext()

  do.call(MultiLinesROC, c(manuscript_ROC_args[names(manuscript_ROC_args) != "y_label_line"],
                           list(x_axis_limits = c(0, 0.2), y_axis_limits = c(0.8, 1)), y_label_line = 2.5, draw_legend = FALSE))
  QuickTitle("ROC curve (zoomed in)")

  do.call(MultiLinesROC, c(manuscript_ROC_args, list(log_FPR = TRUE, draw_legend = FALSE)))
  QuickTitle("ROC curve (logarithmic)")

  do.call(MultiLinesROC, c(manuscript_ROC_args, list(precision_recall = TRUE)))
  QuickTitle("PR curve (non-essential genes)")

  do.call(MultiLinesROC, c(manuscript_ROC_args, list(precision_recall = TRUE, flip = FALSE)))
  QuickTitle("PR curve (essential genes)")

  par(old_par)
  dev.off()
}



# Tabulate essential and non-essential genes (not complete...) ------------

use_roc_df_list <- blomen_hart_ROC_df_list
use_cutoff <- -1

use_roc_df_list <- lapply(use_roc_df_list, function(x) {
  x <- x[order(x[, "Entrez_ID"]), ]
  row.names(x) <- NULL
  x
})

stopifnot(length(unique(lapply(use_roc_df_list, function(x) x[, "Plasmid_ID"]))) == 1)

mean_mean_log2fc <- rowMeans(do.call(cbind, lapply(use_roc_df_list, function(x) x[, "Mean_log2FC"])))
new_order <- order(!(use_roc_df_list[[1]][, "Is_essential"]),
                   mean_mean_log2fc
                   )
use_roc_df_list <- lapply(use_roc_df_list, function(x) {
  x <- x[new_order, ]
  row.names(x) <- NULL
  return(x)
})


