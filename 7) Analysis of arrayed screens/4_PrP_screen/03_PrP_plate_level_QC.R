# 2022-01-18


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "2_Analyzing_data",    "01_calculating_scores.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "01_Plate_level_QC.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir, "4_output", "PrP")
manuscript_dir <- file.path(output_dir, "Figures", "Manuscript", "2) Component plots")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))




# Calculate plate-wise quality metrics ------------------------------------

PlotZPrimes(PrP_df, filter_NT = TRUE)
PlotSSMDControls(PrP_df, filter_NT = TRUE)

PlotZPrimes(PrP_df, filter_NT = TRUE, reorder_plates = TRUE)
PlotSSMDControls(PrP_df, filter_NT = TRUE, reorder_plates = TRUE)



# Export plots as PDF or PNG files ----------------------------------------

plot_width <- 5.5
plot_height <- 3.8

pdf(file = file.path(output_dir, "Figures", "Quality metrics", "Quality metrics.pdf"),
    width = plot_width, height = plot_height
    )
PlotZPrimes(PrP_df, filter_NT = TRUE)
PlotSSMDControls(PrP_df, filter_NT = TRUE)
dev.off()


png(filename = file.path(output_dir, "Figures", "Quality metrics", "Z_prime.png"),
    width = plot_width, height = plot_height, units = "in", res = 600
    )
PlotZPrimes(PrP_df, filter_NT = TRUE)
dev.off()


png(filename = file.path(output_dir, "Figures", "Quality metrics", "SSMD.png"),
    width = plot_width, height = plot_height, units = "in", res = 600
    )
PlotSSMDControls(PrP_df, filter_NT = TRUE)
dev.off()



# Export plots for the manuscript -----------------------------------------

manuscript_width <- 2.2
manuscript_height <- 1.2
manuscript_mai <- c(0.4, 0.5, 0.05, 0.05)

pdf(file = file.path(manuscript_dir, "Figure 6B - Z-prime.pdf"),
    width = manuscript_width, height = manuscript_height
    )
par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
y_axis_vec <- seq(-0.25, 1, by = 0.25)
y_axis_labs <- ifelse(y_axis_vec %in% c(0, 0.5, 1), format(round(y_axis_vec, digits = 1)), NA)
y_axis_labs[y_axis_labs == "0.0"] <- "0"
plates_in_order <- PlotZPrimes(PrP_df, filter_NT = TRUE, use_mai = manuscript_mai,
                               y_label_line = 2.05, roman_plates = FALSE, reorder_plates = TRUE,
                               label_plates = FALSE,
                               plate_labels_line = 0.1, x_label_line = 1.2,
                               y_axis_ticks = y_axis_vec,
                               y_axis_labels = y_axis_labs,
                               point_cex = 0.6
                               )
dev.off()


pdf(file = file.path(manuscript_dir, "Figure 6C - controls SSMD.pdf"),
    width = manuscript_width, height = manuscript_height
    )
par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
PlotSSMDControls(PrP_df, filter_NT = TRUE, use_mai = manuscript_mai,
                 y_label_line = 2.1, y_limits_include = c(0, 12),
                 roman_plates = FALSE, reorder_plates = TRUE,
                 plates_in_order = plates_in_order,
                 y_axis_label = "SSMD (controls)",
                 plate_labels_line = 0.15, x_label_line = 1.4,
                 point_cex = 0.6, label_plates = FALSE
                 )
mtext(c(1, 24), side = 1, line = 0.15, at = c(0.75, 12.25), cex = par("cex"))
mtext(FormatPlotMath("Plate number"), side = 1, line = 1.15, cex = par("cex"))
dev.off()



