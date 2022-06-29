# 2022-01-20


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data", "04_Histograms.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir, "4_output", "PrP")
manuscript_dir <- file.path(output_dir, "Figures", "Manuscript", "2) Component plots")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Modify labels (for PrPc screen) -----------------------------------------

AdjustLabels()
controls_labels[["Pos"]] <- c("Positive", "controls", expression("(" * italic("PRNP") * " gene)"))



# Adjust default plot settings --------------------------------------------

## Use the following code to modify the default version of the histogram
## to use more pronounced borders for the vertical bars:
#
# OriginalThreeHistograms <- ThreeHistograms
# ThreeHistograms <- function(...) {
#   all_args <- list(...)
#   if (!("gene_border_color" %in% names(all_args))) {
#     all_args[["gene_border_color"]] <- "gray30"
#   }
#   do.call(OriginalThreeHistograms, all_args)
# }



# Plot histograms ---------------------------------------------------------

ThreeHistograms(PrP_df, "Raw_rep1")
ThreeHistograms(PrP_df, "FoldNT_rep1")
ThreeHistograms(PrP_df, "DeltaNT_rep1")
ThreeHistograms(PrP_df, "PercActivation_rep1")

ThreeHistograms(PrP_df, "Log2FC_rep1")
ThreeHistograms(PrP_df, "Log2FC_rep1", gene_border_color = "gray30")

ThreeHistograms(PrP_df, "CellTiterGlo_raw")
ThreeHistograms(PrP_df, "CellTiterGlo_foldNT")
ThreeHistograms(PrP_df, "p_value_deltaNT")




# Export plots for the manuscript -----------------------------------------

old_controls_labels <- controls_labels

controls_labels <- list(
  "Gene" = c("TF sub-library", "(T.gonfio)"),
  "NT"   = c("Non-targeting", "controls"),
  "Pos"  = c("Positive controls", expression("(" * italic("PRNP") * " gene)"))
)

manuscript_width <- 2.7
manuscript_height <- 2.25
manuscript_mai <- c(0.4, 0.5, 0.3, 0.2)

pdf(file = file.path(manuscript_dir, "Figure 6D - histogram.pdf"),
    width = manuscript_width, height = manuscript_height
    )
par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
ThreeHistograms(PrP_df, "FoldNT_rep1", use_mai = manuscript_mai,
                legend_x_start = -7, legend_y_mid = 0.7,
                use_mgp = c(2.4, 0.5, 0), x_axis_mgp = c(1.725, 0.375, 0),
                x_axis_label = FormatPlotMath("PrPc levels (normalized by plate medians)"),
                small_gap_size = 1.05, large_gap_multiplier = 1.5
                )
dev.off()

controls_labels <- old_controls_labels




# Export histograms as PDF and PNG files ----------------------------------

plot_width <- 7
plot_height <- 5


pdf(file = file.path(output_dir, "Figures", "Histograms", "Histograms.pdf"),
    width = plot_width, height = plot_height
    )
for (use_column in names(column_file_names)) {
  ThreeHistograms(PrP_df, use_column)
}
dev.off()



for (i in seq_along(column_file_names)) {
  use_column <- names(column_file_names)[[i]]
  file_name <- paste0("Histogram - ", i,  ") ", column_file_names[[i]], ".png")
  png(filename = file.path(output_dir, "Figures", "Histograms", "PNGs", file_name),
      width = plot_width, height = plot_height, units = "in", res = 600
      )
  ThreeHistograms(PrP_df, use_column)
  dev.off()
}


