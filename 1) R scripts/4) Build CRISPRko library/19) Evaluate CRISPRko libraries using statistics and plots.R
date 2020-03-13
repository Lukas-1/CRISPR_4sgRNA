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
load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRko_RData_directory, "13) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))






# # Try stuff ---------------------------------------------------------------


ViolinBox_Sources(merged_CRISPRko_df, "CRISPOR_4MM_specificity",
                  aggregate_scores = TRUE
                  )


BarPlot_Sources(merged_CRISPRko_df, "GuideScan_specificity", show_sublibraries = TRUE, filter_top4 = TRUE)
BarPlot_Sources(merged_CRISPRko_df, "Are_overlapping")
BarPlot_Sources(merged_CRISPRko_df, "Have_homologies")



# ViolinBox_Sources(merged_CRISPRko_df, "CRISPOR_4MM_specificity",
#                   aggregate_scores = TRUE,
#                   filter_for_complete = TRUE
#                   )
#
#
#
#
# use_width <- pdf_width * 0.9
# use_height <- pdf_height * 1.2
#
#
# png(file = file.path(output_plots_directory, paste0("4sg library overview - Box plots - C1) sources - rest vs. 4sg - filtered - test - GuideScan_specificity.png")),
#     width = use_width, height = use_height, units = "in", res = 300
#     )
# ViolinBox_Sources(merged_CRISPRko_df, "GuideScan_specificity",
#                   show_rest_v_4sg = TRUE, filter_top4 = TRUE
#                   )
# dev.off()
#
#
#
#
#
#
#
# SourcesBoxPlots(merged_CRISPRko_df)
#
#
# expanded_CRISPR_df <- ExpandedSubgroupsDf(FilterCRISPRDf(merged_CRISPRko_df))
#
#
#
#
#
#
# ViolinBox_Sources(merged_CRISPRko_df, "GuideScan_specificity", show_sublibraries = FALSE, filter_top4 = TRUE)
#
# ViolinBox_Sources(merged_CRISPRko_df, "GuideScan_specificity", show_rest_v_4sg = TRUE, filter_top4 = FALSE)
# ViolinBox_Sources(merged_CRISPRko_df, "GuideScan_specificity", show_rest_v_4sg = TRUE, filter_top4 = TRUE)
#
# ViolinBox_Sources(merged_CRISPRko_df, "GuideScan_specificity", show_sublibraries = TRUE)
#
#
#
#
#
# CRISPR_df <- merged_CRISPRko_df
# y_column <- "GuideScan_specificity"
#
# filtered_df <- FilterCRISPRDf(merged_CRISPRko_df)
# expanded_df <- OriginalSubgroupsDf(filtered_df, "GuideScan_specificity")
#
#
#
# CRISPR_df <- merged_CRISPRko_df
# png(file = file.path(output_plots_directory, paste0("4sg library overview - barplots - test - GuideScan_specificity.png")),
#     width = pdf_width, height = pdf_height, units = "in", res = 300
#     )
# layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
# ViolinBox_UniqueTwoGroups(CRISPR_df, "GuideScan_specificity", show_title = FALSE)
# ViolinBox_results <- ViolinBox_UniqueLibraries(CRISPR_df, "GuideScan_specificity", show_title = FALSE)
# OuterTitleForLayout(ViolinBox_results[["title_text"]])
# dev.off()
#
#
# CRISPR_df <- merged_CRISPRko_df
# layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
# TwoViolinBox(CRISPR_df, "GuideScan_specificity", show_title = FALSE)
# ViolinBox_results <- ViolinBox(CRISPR_df, "GuideScan_specificity", show_title = FALSE)
# OuterTitleForLayout(ViolinBox_results[["title_text"]])









# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_CRISPRko_df)





# Draw Venn diagrams ------------------------------------------------------

PlotVennDiagrams(merged_CRISPRko_df)






# Plot histograms describing the 4sg combination as a whole ---------------

Plot4sgData(sgRNAs_overview_df)

# head(sgRNAs_overview_df[order(sgRNAs_overview_df[["Longest_subsequence"]], decreasing = TRUE), ])





# Plot categorical data ---------------------------------------------------

UniqueSequencesBarPlots(merged_CRISPRko_df)

SourcesBarPlots(merged_CRISPRko_df)





# Plot numerical data -----------------------------------------------------

UniquePointsBoxPlots(merged_CRISPRko_df)

SourcesBoxPlots(merged_CRISPRko_df)







