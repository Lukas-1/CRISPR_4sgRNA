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
load(file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))





# Try stuff ---------------------------------------------------------------

ViolinBox_Sources(merged_replaced_CRISPRa_df,
                  "Deviation_from_TSS_window",
                  "show_sublibraries"      = FALSE,
                  "filter_complete_genes"  = FALSE,
                  "filter_complete_scores" = FALSE
                  )


BarPlot_Sources(merged_replaced_CRISPRa_df,
                "Are_overlapping",
                "show_sublibraries"      = FALSE,
                "filter_complete_genes"  = FALSE,
                "filter_complete_scores" = TRUE
                )

BarPlot_Sources(merged_replaced_CRISPRa_df,
                "Have_homologies",
                "show_sublibraries"      = FALSE,
                "filter_complete_genes"  = FALSE,
                "filter_complete_scores" = TRUE
                )



BarPlot_Sources(merged_replaced_CRISPRa_df, "Are_overlapping")
BarPlot_Sources(merged_replaced_CRISPRa_df, "Have_homologies")
BarPlot_Sources(merged_replaced_CRISPRa_df, "Deviation_from_TSS_window")




ViolinBox_Sources(merged_replaced_CRISPRa_df,
                  "Deviation_from_TSS_window",
                  "show_sublibraries"      = FALSE,
                  "filter_complete_genes"  = FALSE,
                  "filter_complete_scores" = FALSE
                  )

ViolinBox_Sources(merged_replaced_CRISPRa_df,
                  "GuideScan_specificity",
                  "aggregate_scores"       = TRUE,
                  "show_sublibraries"      = FALSE,
                  "filter_complete_scores" = TRUE
                  )


BarPlot_UniqueLibraries(merged_replaced_CRISPRa_df,
                        "Expected_all22_SNP_AF_max_Kaviar"
                        )


# load(file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData"))
# for (unique_ID in unique(merged_replaced_CRISPRa_df[["Combined_ID"]][merged_replaced_CRISPRa_df[["Is_control"]] == "No"])) {
#   are_this_ID <- merged_replaced_CRISPRa_df[["Combined_ID"]] == unique_ID
#   best_TSSs <- unique(merged_replaced_CRISPRa_df[are_this_ID, "Best_TSS"])
#   best_TSSs <- best_TSSs[!(is.na(best_TSSs))]
#   if (length(best_TSSs) > 1) {
#     stop()
#   }
# }
# print("hurray!")
# stop()
#
# show_columns <- c("Entrez_ID", "Original_entrez", "Gene_symbol", "Original_symbol",
#                   "Source", "Chromosome", "Entrez_chromosome", "Best_TSS", "Cut_location",
#                   "Distance_from_TSS", "GuideScan_specificity"
#                   )
#
# merged_replaced_CRISPRa_df[are_this_ID, show_columns]
#
#
#
# CRISPR_df <- merged_replaced_CRISPRa_df
# CRISPR_df <- FilterCRISPRDf(CRISPR_df)
#
# all_sources_fac <- ReformatSourceToFactor(CRISPR_df[["Source"]])
#
# all_sources_vec <- as.character(all_sources_fac)
#
#
# are_from_this_source_list <- sapply(intersect(libraries_order, all_sources_vec),
#                                     function(x) grepl(x, all_sources_vec, fixed = TRUE),
#                                     simplify = FALSE
#                                     )
# sources_df <- data.frame(
#   do.call(cbind, are_from_this_source_list),
#   "Are_4sg" = factor(ifelse(CRISPR_df[["Rank"]] %in% 1:4,
#                             "4sg",
#                             "Rest"
#                             ),
#                      levels = c("Rest", "4sg")
#                      ),
#   check.names = FALSE
# )
#
# use_colors <- c(brewer.pal(9, "Blues")[[3]], brewer.pal(9, "Greens")[[3]], brewer.pal(9, "Reds")[[3]])
#
#
# # Plots using 'by' argument
# foo <- plot(euler(sources_df, by = list(Are_4sg)),
#             labels = list(cex = 0.5),
#             quantities = list(font = 2, cex = 0.2),
#             fills = list(fill = use_colors, alpha = 1)
#             )
# print(foo)





# # Try stuff ---------------------------------------------------------------
#
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity",
#                   aggregate_scores = TRUE
#                   )
#
#
# use_width <- pdf_width * 0.9
# use_height <- pdf_height * 1.2
#
#
# png(file = file.path(output_plots_directory, paste0("4sg library overview - Box plots - D2) sources - sublibraries - unfiltered - test - GuideScan_specificity.png")),
#     width = use_width, height = use_height, units = "in", res = 300
#     )
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity",
#                   show_sublibraries = TRUE, filter_top4 = FALSE
#                   )
# dev.off()
#
#
#
# SourcesBoxPlots(merged_replaced_CRISPRa_df)
#
#
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_sublibraries = FALSE, filter_top4 = FALSE)
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_sublibraries = FALSE, filter_top4 = TRUE)
#
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_rest_v_4sg = TRUE, filter_top4 = FALSE)
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_rest_v_4sg = TRUE, filter_top4 = TRUE)
#
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_sublibraries = TRUE, filter_top4 = FALSE)
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_sublibraries = TRUE, filter_top4 = FALSE, collapse_GPP = TRUE)
# ViolinBox_Sources(merged_replaced_CRISPRa_df, "GuideScan_specificity", show_sublibraries = TRUE, filter_top4 = TRUE)
#
#
# SourcesBoxPlots(merged_replaced_CRISPRa_df)
#
#
#
#
#
#
# CRISPR_df <- merged_replaced_CRISPRa_df
# y_column <- "GuideScan_specificity"
#
#
#
#
# filtered_df <- FilterCRISPRDf(merged_replaced_CRISPRa_df)
# expanded_df <- OriginalSubgroupsDf(filtered_df, "GuideScan_specificity")
#
# nrow(filtered_df); nrow(expanded_df)
#







# Draw scatter plots ------------------------------------------------------

MakeScatterPlots(merged_replaced_CRISPRa_df)





# Draw Venn diagrams ------------------------------------------------------

PlotVennDiagrams(merged_replaced_CRISPRa_df)






# Plot histograms describing the 4sg combination as a whole ---------------

Plot4sgData(sgRNAs_overview_df, merged_replaced_CRISPRa_df)

# head(sgRNAs_overview_df[order(sgRNAs_overview_df[["Longest_subsequence"]], decreasing = TRUE), ])






# Plot categorical data ---------------------------------------------------

UniqueSequencesBarPlots(merged_replaced_CRISPRa_df)

SourcesBarPlots(merged_replaced_CRISPRa_df)






# Plot numerical data -----------------------------------------------------

UniquePointsBoxPlots(merged_replaced_CRISPRa_df)

SourcesBoxPlots(merged_replaced_CRISPRa_df)





