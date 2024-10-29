## 2024-06-10


# Load packages and source code -------------------------------------------

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
project_dir     <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
rdata_dir       <- file.path(project_dir, "03_R_objects")
output_dir      <- file.path(project_dir, "05_output")
figures_dir     <- file.path(output_dir, "Figures", "Comparison prepool-postpool")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "10_create_figures_from_count_data.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))



# Define gene selections --------------------------------------------------

ess_genes_blomen_hart    <- essentials_2020Q2_df[, "Entrez_ID"]
noness_genes_blomen_hart <- non_essentials_2020Q2_df[, "Entrez_ID"]

common_plasmids_logfc_df_list <- CommonPlasmidsRocDfList(list(logfc_prepooled_df, logfc_postpooled_df))

blomen_hart_ROC_df_list <- LogFcDfListToRocDfList(
  common_plasmids_logfc_df_list,
  essential_entrezs = ess_genes_blomen_hart,
  non_essential_entrezs = noness_genes_blomen_hart
)



# Export manuscript-style mean gamma violin plots -------------------------

rep_list <- c(split(blomen_hart_ROC_df_list[[1]][, "Mean_log2FC"], !(blomen_hart_ROC_df_list[[1]][, "Is_essential"])),
              split(blomen_hart_ROC_df_list[[2]][, "Mean_log2FC"], !(blomen_hart_ROC_df_list[[2]][, "Is_essential"]))
              )

pdf(file.path(figures_dir, "Manuscript-style violin plots.pdf"),
    width = 1.9, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
violin_positions <- MeanSwarms(rep_list, group_labels = c("prepool", "postpool"),
                               show_truncation = FALSE
                               )
par(old_par)
dev.off()



# Create scatter plots ----------------------------------------------------

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

pdf(file.path(figures_dir, "Manuscript-style scatter plot.pdf"),
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



# Draw a scatter plot with density plots on the sides ---------------------

use_vec_list <- lapply(1:2, function(x) {
  logfc_df <- common_plasmids_logfc_df_list[[x]]
  results_vec <- logfc_df[, "Mean_log2FC"][logfc_df[, "Entrez_ID"] %in% ess_genes_blomen_hart]
  results_vec <- results_vec / 10
  results_vec[results_vec > 0.2] <- 0.2
  results_vec[results_vec < -0.6] <- -0.6
  return(results_vec)
})

pdf(file.path(figures_dir, "Scatter plot - with side density plots.pdf"),
    width  = 2.15,
    height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.45))
DensityTrapezoidScatter(use_vec_list)
dev.off()



# Compare the separation between E and NE genes ---------------------------

this_points_vec <- c(separation_prepooled_mat["Robust SSMD", ],
                     separation_postpooled_mat["Robust SSMD", ]
                     )
this_points_vec <- abs(this_points_vec)

pdf(file.path(figures_dir, "Manuscript-style SSMD.pdf"),
    width = 1.1, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.1))
ComparePoints(this_points_vec, left_gap = 0.55, right_gap = 0.45,
              y_upper_limit = 2, group_labels = c(NA, "prepool", "postpool")
              )
par(old_par)
dev.off()



# Export manuscript-style combined ROC curves -----------------------------

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

pdf(file.path(figures_dir, "Manuscript-style ROC curves.pdf"),
    width = 2, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.7, mai = c(0.42, 0.5, 0.38, 0.3))
do.call(MultiLinesROC, manuscript_ROC_args)
par(old_par)
dev.off()


