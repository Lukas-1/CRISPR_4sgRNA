### 29th January 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
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





# Draw plots describing the 4sg library as a whole ------------------------

PlotNumGenesInLibrary()

DrawAllDonutBars(merged_replaced_CRISPRa_df)

Plot4sgData(sgRNAs_overview_df, merged_replaced_CRISPRa_df)

PlotVennDiagrams(merged_replaced_CRISPRa_df)





# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_replaced_CRISPRa_df)






# Plot categorical data ---------------------------------------------------

UniqueSequencesBarPlots(merged_replaced_CRISPRa_df)

SourcesBarPlots(merged_replaced_CRISPRa_df)






# Plot numerical data -----------------------------------------------------

UniquePointsBoxPlots(merged_replaced_CRISPRa_df)

SourcesBoxPlots(merged_replaced_CRISPRa_df)





