### 5th November 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData"))






# Rank sgRNAs -------------------------------------------------------------

merged_replaced_CRISPRa_df <- RankCRISPRDf(merged_replaced_CRISPRa_df, ID_column = "AltTSS_ID")






# Find combinations of non-overlapping sgRNAs -----------------------------

# load(file.path(CRISPRa_RData_directory, "18) Compile a list of human transcription factors.RData"))
# replaced_TF_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]
#
#
# PrioritizeNonOverlapping(merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Gene_symbol"] %in% "AEBP1", ])
#
#




CRISPR_df <- merged_replaced_CRISPRa_df
combined_IDs <- GetUniqueIDs(CRISPR_df, "AltTSS_ID")
controls_df <- NonOverlappingDfForControls(CRISPR_df)


cl <- parallel::makeCluster(5)
parallel::clusterExport(cl, list("SortCombinations", "CreateCombinations", "MessageID",
                                  "ReorderSubDfByLocation", "NumHomologousPairs", "SplitIntoSubstrings",
                                  "CRISPR_df",
                                  "preferred_AF_max_column", "SNP_frequency_cutoff"
                                  )
                         )
reordered_df_list <- parallel::parLapply(cl,
                                         combined_IDs,
                                         function(x) SortCombinations(CRISPR_df[CRISPR_df[, "AltTSS_ID"] == x, , drop = FALSE],
                                                                      min_overlaps    = 0:10,
                                                                      min_spaces      = 50L,
                                                                      only_top_24_GPP = FALSE
                                                                      )
                                         )
parallel::stopCluster(cl)


results_df <- do.call(rbind.data.frame, c(reordered_df_list,
                                          list(controls_df),
                                          list(stringsAsFactors = FALSE, make.row.names = FALSE)
                                          )
                      )

merged_replaced_CRISPRa_df <- results_df








# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"
                      )
     )





# have_CRISPOR_combined_IDs <- unique(merged_replaced_CRISPRa_df[!(is.na(merged_replaced_CRISPRa_df[, "CRISPOR_CFD_specificity"])), "Combined_ID"])
# view_merged_replaced_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Combined_ID"] %in% have_CRISPOR_combined_IDs, ]
# view_CRISPOR_columns <- c(
#   "Entrez_ID",
#   "Gene_symbol",
#   "Rank",
#   "Original_rank",
#   "GuideScan_specificity",
#   "CRISPOR_CFD_specificity",
#   "Num_0MM",
#   "Num_5G_MM",
#   "Num_1MM",
#   preferred_AF_max_column,
#   preferred_rsID_column,
#   "PAM",
#   "CRISPOR_Graf_status"
#
# )
# View(view_merged_replaced_CRISPRa_df[, view_CRISPOR_columns])



