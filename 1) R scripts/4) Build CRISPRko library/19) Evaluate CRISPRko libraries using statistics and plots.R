### 29th January 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
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



# Show some example plots -------------------------------------------------

ViolinBox_Sources(merged_CRISPRko_df, "CRISPOR_4MM_specificity",
                  aggregate_scores = TRUE
                  )
BarPlot_Sources(merged_CRISPRko_df, "GuideScan_specificity",
                filter_top4 = TRUE
                )
BarPlot_Sources(merged_CRISPRko_df, "Are_overlapping")





# Export plots for the manuscript -----------------------------------------

DrawAllManuscriptPlots(merged_CRISPRko_df)


## Draw the deletion histogram
pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              paste0("Histogram - size of expected deletions", ".pdf")),
    width = 3.2, height = 2.2
    )
par(cex = manuscript_cex, lwd = manuscript_lwd, mai = c(0.4, 0.47, 0.4, 0.2))
hist_results <- DrawDeletionHistogram(sgRNAs_overview_df,
                                      use_title    = "",
                                      x_axis_label = "Size of expected deletion (base pairs)",
                                      y_axis_label = "Number of genes",
                                      fill_color   = CRISPRo_colors[[2]],
                                      x_label_line = 1.8,
                                      y_label_line = 2.53
                                      )

text(x      = hist_results[["mids"]][which.max(hist_results[["counts"]])],
     y      = par("usr")[[4]] + diff(grconvertY(c(0, 1.2), from = "lines", to = "user")),
     labels = "CRISPRo library",
     adj    = 0.5,
     xpd    = NA
     )
dev.off()




## Draw the gene deletions doughnut plot
pdf(file.path(output_plots_directory, "Manuscript", "Whole library",
              paste0("Doughnut plots - CRISPRo plasmids", ".pdf")
              ),
    width = 3.3, height = 2.15
    )

par(cex = manuscript_cex, lwd = manuscript_lwd)
do.call(SummaryDonutBar,
        c(manuscript_donut_args,
          list(CRISPR_df    = merged_CRISPRko_df,
               targets_df   = guides_CDS_df,
               x_axis_label = "Genes in CRISPRo library",
               use_map_list = manuscript_map_list
               )
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







