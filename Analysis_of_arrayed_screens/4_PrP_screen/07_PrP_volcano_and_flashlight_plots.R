# 2022-01-20


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "05_Volcano_and_flashlight_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir, "4_output", "PrP")
manuscript_dir <- file.path(output_dir, "Figures", "Manuscript", "2) Component plots")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))




# Export plots for the manuscript -----------------------------------------

old_controls_colors <- controls_colors
controls_colors <- controls_colors[c(3, 1, 2)]

controls_labels <- list(
  "Gene" = c("Genes in", "T.gonfio", "library"),
  "NT"   = c("Non-", "targeting", "controls"),
  "Pos"  = c("Positive", "controls", expression("(" * italic("PRNP")), "gene)")
)

manuscript_width <- 3.5
manuscript_height <- 2.6
manuscript_mai <- c(0.5, 0.5, 0.2, 1)

pdf(file = file.path(manuscript_dir, "Figure 6F - volcano plot.pdf"),
    width = manuscript_width, height = manuscript_height
    )
par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 label_points = FALSE, indicate_areas = TRUE,
                 indicate_lines = TRUE,
                 indicate_log2FCs = log2(2), indicate_p_values = 0.05,
                 use_mai = manuscript_mai,
                 use_mgp = c(1.92, 0.5, 0), x_axis_mgp = c(1.725, 0.375, 0),
                 small_gap_size = 1.2, large_gap_multiplier = 1.5,
                 point_x_start = 0.1, lines_x_start = 1
                 )
dev.off()

controls_colors <- old_controls_colors



# Differentiate between own and Tubingen NT controls ----------------------

mat_384 <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
are_own_NT <- PrP_df[, "Target_flag"] %in% "Own NT control"
are_Tubingen_NT <- PrP_df[, "Is_NT_ctrl"] & !(are_own_NT)

are_NT_ctrl <- PrP_df[, "Is_NT_ctrl"]

PrP_df[, "Is_NT_ctrl"] <- are_Tubingen_NT
PrP_df[, "Custom_color"] <- are_own_NT

controls_labels <- list(
  "o_NT" = c("Own", "non-targeting", "controls"),
  "T_NT" = c("Tubingen", "non-targeting", "controls"),
  "Pos"  = c("Positive", "controls", expression("(" * italic("PRNP") * " gene)")),
  "Gene" = c("Genes in ", "CRISPRa", "library")
)

controls_colors <- c(custom_color, NT_ctrl_color, pos_ctrl_color, "black")



# Plot data ---------------------------------------------------------------

VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_deltaNT",
                 show_title = "Volcano plot (p values from untransformed data)"
                 )

VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_act",
                 show_title = "Volcano plot (p values from % activation)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 show_title = "Volcano plot (p values from log2-transformed data)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_act_log2",
                 show_title = "Volcano plot (p values from % activation, log2 data)"
                 )

VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_deltaNT_Glo",
                 show_title = "Volcano plot (untransformed, CellTitreGlo-norm.)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_act_Glo",
                 show_title = "Volcano plot (% activation, CellTitreGlo-norm.)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_log2_Glo",
                 show_title = "Volcano plot (log2, CellTitreGlo-normalized)"
                 )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_act_log2_Glo",
                 show_title = "Volcano plot (% activation, log2, CellTitreGlo-norm.)"
                 )


VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "SSMD_deltaNT",
                 show_title = "Volcano Plot (SSMD from untransformed data)"
                 )

VolcanoFlashPlot(PrP_df, "PercActivation_rep1", "SSMD_deltaNT",
                 show_title = "Volcano Plot (SSMD from untransformed data)"
                 )



# Export individually customized plots ------------------------------------

base_width <- 5.5
base_height <- 5.1

selected_volcanoes_dir <- file.path(output_dir, "Figures", "Volcano plots", "Selected plots")


png(file.path(selected_volcanoes_dir, "1) Volcano plot - cutoffs shown.png"),
    width = base_width + 0.8, height = base_height, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 show_title = Embolden(FormatPlotMath("Volcano plot (p values, log2FC)")),
                 label_points = FALSE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2), indicate_p_values = 0.05
                 )
dev.off()


pdf(file.path(selected_volcanoes_dir, "2) Volcano plot - genes shown.pdf"),
    width = base_width + 0.8, height = base_height#, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_rep1", "p_value_log2",
                 show_title = Embolden(FormatPlotMath("Volcano plot (p values, log2FC)")),
                 label_points = TRUE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2), tiny_labels = TRUE
                 )
dev.off()


png(file.path(selected_volcanoes_dir, "3) Volcano plot - Glo-normalized - cutoffs shown.png"),
    width = base_width + 0.8, height = base_height, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_log2_Glo",
                 show_title = Embolden(FormatPlotMath("Volcano plot (normalized to CellTiter-Glo)")),
                 label_points = FALSE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2)
                 )
dev.off()



pdf(file.path(selected_volcanoes_dir, "4) Volcano plot - Glo-normalized - genes shown.pdf"),
    width = base_width + 0.8, height = base_height#, units = "in", res = 600
    )
VolcanoFlashPlot(PrP_df, "Log2FC_Glo_rep1", "p_value_log2_Glo",
                 show_title = Embolden(FormatPlotMath("Volcano plot (normalized to CellTiter-Glo)")),
                 label_points = TRUE, indicate_areas = TRUE, indicate_lines = TRUE,
                 indicate_log2FCs = log2(2), tiny_labels = TRUE
                 )
dev.off()




# Export all plots as PDF and PNG files -----------------------------------

ExportAllVolcanoAndFlashlightPlots(PrP_df)




