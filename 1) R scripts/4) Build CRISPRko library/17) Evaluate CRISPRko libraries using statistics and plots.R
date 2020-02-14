### 29th January 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
output_plots_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko", "Plots")






# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))





# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_CRISPRko_df)





# Draw Venn diagrams ------------------------------------------------------

PlotVennDiagrams(merged_CRISPRko_df)





# Plot histograms describing the 4sg combination as a whole ---------------

Plot4sgData(sgRNAs_overview_df)

head(sgRNAs_overview_df[order(sgRNAs_overview_df[["Longest_subsequence"]], decreasing = TRUE), ])





# Plot categorical data ---------------------------------------------------

PlotCategoricalData(merged_CRISPRko_df)





# Plot numerical data -----------------------------------------------------

PlotNumericalData(merged_CRISPRko_df)


# CRISPR_df <- merged_CRISPRko_df
# png(file = file.path(output_plots_directory, paste0("4sg library overview - barplots - test - GuideScan_specificity.png")),
#     width = pdf_width, height = pdf_height, units = "in", res = 300
#     )
# layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
# TwoViolinBox(CRISPR_df, "GuideScan_specificity", show_title = FALSE)
# ViolinBox_results <- ViolinBox(CRISPR_df, "GuideScan_specificity", show_title = FALSE)
# OuterTitleForLayout(ViolinBox_results[["title_text"]])
# dev.off()






