### 29th January 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
output_plots_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko", "Plots")



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))
load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRko_RData_directory, "16) Find all genes targeted by each sgRNA - summary data.RData"))




# Add data on other (unintended) targeted genes ---------------------------

merged_CRISPRko_df <- AddOtherTargetBooleans(merged_CRISPRko_df,
                                             guides_CDS_df,
                                             guides_CDS_protein_df
                                             )



# Tabulate data for the manuscript ----------------------------------------

CRISPRko_figure_list <- PrepareManuscriptPlots(merged_CRISPRko_df)




# Save data ---------------------------------------------------------------

save(list = "CRISPRko_figure_list",
     file = file.path(CRISPRko_RData_directory, "19) Evaluate CRISPRko libraries using statistics and plots.RData")
     )



# Show some example plots -------------------------------------------------

BarPlot_Sources(merged_CRISPRko_df,
                "all22_SNP_AF_max_Kaviar",
                show_sublibraries = FALSE,
                filter_top4       = TRUE,
                use_cutoff        = 50
                )

ViolinBox_Sources(merged_CRISPRko_df, "CRISPOR_4MM_specificity",
                  aggregate_scores = TRUE
                  )
BarPlot_Sources(merged_CRISPRko_df, "GuideScan_specificity",
                filter_top4 = TRUE
                )
BarPlot_Sources(merged_CRISPRko_df, "Are_overlapping")

BarPlot_Sources(merged_CRISPRko_df, "Have_homologies", filter_complete_genes = FALSE)




# Export plots for the manuscript -----------------------------------------

DrawAllManuscriptPlots(CRISPRko_figure_list, rename_libraries = TRUE)
DrawAllManuscriptPlots(CRISPRko_figure_list, make_PNGs = TRUE)
DrawAllManuscriptPlots(CRISPRko_figure_list, make_EMFs = TRUE,
                       rename_libraries = TRUE, line_breaks = FALSE,
                       sina_plot = TRUE, draw_whiskers = TRUE
                       )



# emf("T.spiezzo scatter.emf",
#     width = 4.25, height = 4.25, emfPlus = FALSE
#     )
# ScatterPlot(merged_CRISPRko_df, only_top4 = TRUE,
#             "GuideScan_specificity", "CRISPOR_Doench_efficacy",
#             point_alpha = 0.2, point_cex = 0.4,
#             show_title = "T.spiezzo",
#             embed_PNG = TRUE
#             )
# dev.off()
#
#
# emf("T.spiezzo specificity scatter.emf",
#     width = 4.25, height = 4.25, emfPlus = FALSE
#     )
# ScatterPlot(merged_CRISPRko_df, only_top4 = TRUE,
#             "GuideScan_specificity", "CRISPOR_3MM_specificity",
#             point_alpha = 0.2, point_cex = 0.4,
#             show_title = "CRISPR knocko
#             ut",
#             embed_PNG = TRUE
#             )
# dev.off()




## Export the deletion histogram
pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              "Histogram - size of expected deletions.pdf"
              ),
    width = 2.4, height = 1.75
    )
par(cex = manuscript_cex,
    lwd = manuscript_lwd,
    mai = c(0.4, 0.47, 0.2, 0.32) # c(0.7, 1, 0.5, 0.15)
    )
hist_results <- DrawDeletionHistogram(sgRNAs_overview_df,
                                      use_title       = "",
                                      x_axis_label    = "Size of expected deletion (bp)",
                                      y_axis_label    = "Number of genes",
                                      fill_color      = "#3576c0", #colorRampPalette(CRISPRo_colors)(100)[[75]],
                                      x_label_line    = 1.8,
                                      y_label_line    = 2.53,
                                      use_y_max       = 600,
                                      abbreviate_1000 = TRUE
                                      )
mtext(line = 0.25, "T.spiezzo library", cex = par("cex"))
dev.off()




## Export the projected deletions histogram
devEMF::emf(file.path(output_plots_directory, "Manuscript", "Whole library",
                      "Histogram - size of expected deletions.emf"
                      ),
            width = 3.8, height = 2.2, emfPlus = FALSE, coordDPI = 1500
            )
par(cex = manuscript_cex,
    lwd = manuscript_lwd,
    mai = c(0.4, 0.47, 0.2, 0.32)
    )
hist_results <- DrawDeletionHistogram(sgRNAs_overview_df,
                                      use_title       = "",
                                      x_axis_label    = "Size of expected deletion (bp)",
                                      y_axis_label    = "Number of genes",
                                      fill_color      = "#3576c0", #colorRampPalette(CRISPRo_colors)(100)[[75]],
                                      x_label_line    = 1.8,
                                      y_label_line    = 2.53,
                                      use_y_max       = 600,
                                      abbreviate_1000 = TRUE
                                      )
mtext(line = 0.25, "T.spiezzo library", cex = par("cex"))
dev.off()




## Export the plasmid targets doughnut plot
pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              "Doughnut plots - CRISPRo plasmids.pdf"
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
          list(CRISPR_df    = merged_CRISPRko_df,
               targets_df   = guides_CDS_df,
               x_axis_label = "Genes in T.spiezzo library",
               use_map_list = manuscript_map_list
               ),
          new_donut_args
          )
        )
dev.off()



devEMF::emf(file.path(output_plots_directory, "Manuscript", "Whole library",
                      "Doughnut plots - CRISPRo plasmids.emf"
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
          list(CRISPR_df    = merged_CRISPRko_df,
               targets_df   = guides_CDS_df,
               x_axis_label = "Genes in T.spiezzo library",
               use_map_list = manuscript_map_list
               ),
          new_donut_args
          )
        )
dev.off()




# Draw plots describing the 4sg library as a whole ------------------------

PlotNumGenesInLibrary()

DrawAllDonutBars()

Plot4sgData(sgRNAs_overview_df, merged_CRISPRko_df)

PlotVennDiagrams(merged_CRISPRko_df)




# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_CRISPRko_df)




# Plot categorical data ---------------------------------------------------

UniqueSequencesBarPlots(merged_CRISPRko_df)

SourcesBarPlots(merged_CRISPRko_df)




# Plot numerical data -----------------------------------------------------

UniquePointsBoxPlots(merged_CRISPRko_df)

SourcesBoxPlots(merged_CRISPRko_df)




