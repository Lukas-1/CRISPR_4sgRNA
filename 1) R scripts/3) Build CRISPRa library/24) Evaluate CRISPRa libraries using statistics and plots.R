### 29th January 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
output_plots_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa", "Plots")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRa_RData_directory, "19) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))






# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_replaced_CRISPRa_df)





# Draw Venn diagrams ------------------------------------------------------

PlotVennDiagrams(merged_replaced_CRISPRa_df)





# Plot histograms describing the 4sg combination as a whole ---------------

Plot4sgData(sgRNAs_overview_df)

head(sgRNAs_overview_df[order(sgRNAs_overview_df[["Longest_subsequence"]], decreasing = TRUE), ])





# Plot categorical data ---------------------------------------------------

PlotCategoricalData(merged_replaced_CRISPRa_df)





# Plot numerical data -----------------------------------------------------

PlotNumericalData(merged_replaced_CRISPRa_df)






