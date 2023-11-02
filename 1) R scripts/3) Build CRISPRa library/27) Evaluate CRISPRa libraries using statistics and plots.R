### 29th January 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
output_plots_directory  <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa", "Plots")



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRa_RData_directory, "24) Find all TSSs targeted by each sgRNA - summary data.RData"))




# Add data on other (unintended) targeted TSSs ----------------------------

merged_replaced_CRISPRa_df <- AddOtherTargetBooleans(merged_replaced_CRISPRa_df,
                                                     TSS_targets_df,
                                                     TSS_protein_targets_df,
                                                     main_TSS_targets_df,
                                                     main_TSS_protein_targets_df
                                                     )



# Show some example plots -------------------------------------------------

BarPlot_Sources(merged_replaced_CRISPRa_df,
                "all22_SNP_AF_max_Kaviar",
                show_sublibraries = FALSE,
                filter_top4       = TRUE,
                use_cutoff        = 0.1
                )

BarPlot_Sources(merged_replaced_CRISPRa_df,
                "Expected_all22_SNP_AF_max_Kaviar",
                show_sublibraries = FALSE,
                filter_top4       = TRUE
                )

ViolinBox_Sources(merged_replaced_CRISPRa_df,
                  "Deviation_from_TSS_window",
                  show_sublibraries = FALSE
                  )

BarPlot_Sources(merged_replaced_CRISPRa_df,
                "Are_overlapping",
                filter_complete_genes = FALSE
                )

ViolinBox_Sources(merged_replaced_CRISPRa_df,
                  "GuideScan_specificity",
                  aggregate_scores = TRUE
                  )




# Tabulate data for the manuscript ----------------------------------------

CRISPRa_figure_list <- PrepareManuscriptPlots(merged_replaced_CRISPRa_df)




# Save data ---------------------------------------------------------------

save(list = "CRISPRa_figure_list",
     file = file.path(CRISPRa_RData_directory, "27) Evaluate CRISPRa libraries using statistics and plots.RData")
     )



# Export plots for the manuscript -----------------------------------------

DrawAllManuscriptPlots(CRISPRa_figure_list, rename_libraries = TRUE)
DrawAllManuscriptPlots(CRISPRa_figure_list, make_PNGs = TRUE)
DrawAllManuscriptPlots(CRISPRa_figure_list, make_EMFs = TRUE,
                       rename_libraries = TRUE, line_breaks = FALSE,
                       sina_plot = TRUE, draw_whiskers = TRUE
                       )


ManuscriptViolinBox(CRISPRa_figure_list[["df_list_filtered"]][["GuideScan_specificity"]],
                    axis_label           = "GuideScan specificity",
                    horizontal           = FALSE,
                    use_mai              = par("mai"),
                    use_width            = pdf_width,
                    use_height           = pdf_height,
                    abbreviate_libraries = TRUE,
                    modality_on_bottom   = FALSE,
                    use_cex              = 1,
                    use_lwd              = 1,
                    CRISPRa_colors       = manuscript_CRISPRa_colors,
                    CRISPRo_colors       = manuscript_CRISPRo_colors,
                    rename_libraries     = TRUE,
                    line_breaks          = TRUE,
                    sina_plot            = TRUE
                    )



## Export the TSS doughnut plot

TSS_colors <- carto_pal(3, "Emrld") # or pink: c("#fdd8eb", "#f990c6", "#d31279")
TSS_colors <- c(TSS_colors[[1]], colorRampPalette(TSS_colors)(100)[[40]], TSS_colors[[3]])
TSS_colors <- c(TSS_colors, rep(TSS_colors[[length(TSS_colors)]], 6 - length(TSS_colors)))

pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              "Doughnut plot - TSSs.pdf"
              ),
    width = 2.4, height = 1.7 # width = 3.4, height = 2.5
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
do.call(TSSDonutBar,
        c(manuscript_donut_args[!(names(manuscript_donut_args) %in% c("use_mai", "bar_label_line", "donut_radius", "donut_y_mid"))],
          list(CRISPR_df       = merged_replaced_CRISPRa_df,
               x_axis_label    = "Genes in T.gonfio library",
               use_labels      = 1:6,
               use_colors      = TSS_colors,
               use_mai         = c(0.2, 0.47, 0.4, 0.32),
               y_axis_label    = "Number of TSSs",
               text_dark_color = "black",
               space           = 0.3,
               bar_label_line  = 0.8,
               donut_radius    = 0.34,
               donut_y_mid     = 0.35
               )
          )
        )
dev.off()



devEMF::emf(file.path(output_plots_directory, "Manuscript", "Thesis",
                      "Doughnut plot - TSSs.emf"
                      ),
            width = 3.5, height = 2.1, emfPlus = FALSE, coordDPI = 1500  # width = 3.4, height = 2.5
            )
par(cex = manuscript_cex, lwd = 0.8)
do.call(TSSDonutBar,
        c(manuscript_donut_args[!(names(manuscript_donut_args) %in% c("use_mai", "bar_label_line", "donut_radius", "donut_x_mid", "donut_y_mid"))],
          list(CRISPR_df       = merged_replaced_CRISPRa_df,
               x_axis_label    = "Genes in T.gonfio library",
               use_labels      = 1:6,
               use_colors      = TSS_colors,
               use_mai         = c(0.2, 0.47, 0.4, 0.32),
               y_axis_label    = "Number of TSSs",
               text_dark_color = "black",
               space           = 0.3,
               bar_label_line  = 0.8,
               donut_radius    = 0.26,
               donut_x_mid     = 0.85,
               donut_y_mid     = 0.23,
               y_axis_label_line = 2.2,
               x_axis_label_line = 1.87
               )
          )
        )
dev.off()




## Export the plasmid targets doughnut plot

pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              paste0("Doughnut plots - CRISPRa plasmids.pdf")
              ),
    width = 2.6, height = 1.7
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
new_donut_args <- list(space             = 0.3,
                       use_mai           = c(0.05, 1, 0.4, 0.16),
                       donut_radius      = 0.34,
                       use_line_height   = 0.72,
                       donut_y_mid       = 0.35,
                       donut_label_y_adj = 0.05
                       )
do.call(SummaryDonutBar,
        c(manuscript_donut_args[!(names(manuscript_donut_args) %in% names(new_donut_args))],
          list(CRISPR_df    = merged_replaced_CRISPRa_df,
               targets_df   = TSS_targets_df,
               x_axis_label = "Plasmids in T.gonfio library",
               use_map_list = manuscript_map_list,
               percent_max  = 100
               ),
          new_donut_args
          )
        )
dev.off()



devEMF::emf(file.path(output_plots_directory, "Manuscript", "Thesis",
                      "Doughnut plots - CRISPRa plasmids.emf"
                      ),
            width = 3, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
par(cex = manuscript_cex, lwd = manuscript_lwd)
new_donut_args <- list(space             = 0.3,
                       use_mai           = c(0.05, 1, 0.4, 0.16),
                       donut_radius      = 0.34,
                       use_line_height   = 0.9,
                       donut_y_mid       = 0.35,
                       donut_label_y_adj = 0.05
                       )
do.call(SummaryDonutBar,
        c(manuscript_donut_args[!(names(manuscript_donut_args) %in% names(new_donut_args))],
          list(CRISPR_df    = merged_replaced_CRISPRa_df,
               targets_df   = TSS_targets_df,
               x_axis_label = "Genes in T.gonfio library",
               use_map_list = manuscript_map_list,
               percent_max  = 100
               ),
          new_donut_args
          )
        )
dev.off()




## Export the TSS histogram

TSS_distances_df <- TSSHistogramsForModality(merged_replaced_CRISPRa_df,
                                             "CRISPRoff",
                                             omit_outside_x_range = TRUE
                                             )
TSS_4sg_distances_vec <- FilterDistanceByGroup(TSS_distances_df, "4sg")


pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              "Histogram - CRISPRoff - window around TSS.pdf"
              ),
    width = 2.05, height = 2
    )
par(cex = 0.6, lwd = 0.7,
    mai = c(0.4, 0.6, 0.2, 0.2)
    )
TSSHistogram(distances_vec        = TSS_4sg_distances_vec,
             use_breaks           = 200,
             use_title            = "",
             modality_text        = "",
             highlight_range      = c(-500, 500),
             highlight_color      = colorRampPalette(brewer.pal(9, "Blues")[c(2, 3)])(3)[[2]],
             hist_color           = brewer.pal(9, "Blues")[[8]],
             omit_outside_x_range = TRUE,
             label_range          = FALSE,
             draw_grid            = FALSE,
             hardcoded_x_axis     = FALSE,
             x_label_line         = 1.8,
             x_axis_mgp           = 0.45,
             y_axis_mgp           = 0.55,
             y_label_line         = 2.7,
             use_tcl              = 0.4,
             x_axis_label         = "Position relative to the TSS"
             )
title("T.gonfio library", cex.main = 1, font.main = 1, line = 0.5)
dev.off()




devEMF::emf(file.path(output_plots_directory, "Manuscript", "Whole library",
                      "Histogram - CRISPRoff - window around TSS.emf"
                      ),
            width = 4.2, height = 2.2, emfPlus = FALSE
            )
par(cex = 0.7, lwd = 0.8,
    mai = c(0.4, 0.5, 0.2, 0.34)
    )
TSSHistogram(distances_vec        = TSS_4sg_distances_vec,
             use_breaks           = 200,
             use_title            = "",
             modality_text        = "",
             highlight_range      = c(-500, 500),
             highlight_color      = colorRampPalette(brewer.pal(9, "Blues")[c(2, 3)])(3)[[2]],
             omit_outside_x_range = TRUE,
             label_range          = FALSE,
             draw_grid            = FALSE,
             hardcoded_x_axis     = FALSE,
             x_label_line         = 1.7,
             x_axis_mgp           = 0.45,
             y_axis_mgp           = 0.55,
             y_label_line         = 2.8,
             use_tcl              = 0.4,
             x_axis_label         = "Position relative to the TSS"
             )
title("T.gonfio library", cex.main = 1, font.main = 1, line = 0.5)
dev.off()




# Draw plots describing the 4sg library as a whole ------------------------

PlotNumGenesInLibrary()

DrawAllDonutBars(merged_replaced_CRISPRa_df)

Plot4sgData(sgRNAs_overview_df, merged_replaced_CRISPRa_df)

PlotVennDiagrams(merged_replaced_CRISPRa_df)

CreateAllTSSHistograms(merged_replaced_CRISPRa_df)





# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_replaced_CRISPRa_df)




# Plot categorical data ---------------------------------------------------

UniqueSequencesBarPlots(merged_replaced_CRISPRa_df)

SourcesBarPlots(merged_replaced_CRISPRa_df)





# Plot numerical data -----------------------------------------------------

UniquePointsBoxPlots(merged_replaced_CRISPRa_df)

SourcesBoxPlots(merged_replaced_CRISPRa_df)




