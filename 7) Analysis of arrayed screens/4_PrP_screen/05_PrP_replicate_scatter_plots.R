# 2022-01-20


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "03_Replicate_scatter_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir, "4_output", "PrP")
manuscript_dir <- file.path(output_dir, "Figures", "Manuscript", "2) Component plots")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Modify labels (for PrPc screen) -----------------------------------------

AdjustLabels()



# Examine the correlation between replicates ------------------------------

ReplicateScatter(PrP_df, "Raw_rep1")
ReplicateScatter(PrP_df, "PercActivation_log2_Glo_rep1")
ReplicateScatter(PrP_df, "Log2FC_rep1", same_scale = FALSE)




# Export plots for the manuscript -----------------------------------------

manuscript_width <- 2.25
manuscript_height <- 2.25
manuscript_mai <- c(0.4, 0.55, 0.3, 0.15)

input_df <- PrP_df
are_gene <- !(is.na(input_df[, "Entrez_ID"]))

x_vec <- PrP_df[, "Log2FC_rep1"][are_gene]
y_vec <- PrP_df[, "Log2FC_rep2"][are_gene]

axis_limits <- range(c(x_vec, y_vec))

pdf(file = file.path(manuscript_dir, "Figure 6E - scatter plot.pdf"),
    width = manuscript_width, height = manuscript_height
    )
par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
ScatterPlot(x_vec,
            y_vec,
            top_label = "",
            use_limits = axis_limits,
            use_mgp = c(1.65, 0.55, 0), x_axis_mgp = c(1.725, 0.375, 0),
            point_size = 0.8
            )
dev.off()




# Export plots as PDF and PNG files ---------------------------------------

ExportAllReplicateScatterPlots(PrP_df)



