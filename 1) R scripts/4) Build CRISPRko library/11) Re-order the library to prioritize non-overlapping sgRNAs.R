### 5th November 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "10) Find overlaps between CRISPRko sgRNA sequences and genetic polymorphisms.RData"))






# Rank sgRNAs -------------------------------------------------------------

merged_CRISPRko_df <- RankCRISPRDf(merged_CRISPRko_df, ID_column = "Combined_ID")






# Find combinations of non-overlapping sgRNAs -----------------------------

# CRISPR_df <- merged_CRISPRko_df
# unique_IDs <- GetUniqueIDs(CRISPR_df, "Combined_ID")
# controls_df <- NonOverlappingDfForControls(CRISPR_df)
#
#
# reordered_df_list <- lapply(unique_IDs,
#                             function(x) SortCombinations(CRISPR_df[CRISPR_df[, "Combined_ID"] == x, , drop = FALSE],
#                                                          min_spaces = 50L, only_top_24_GPP = FALSE,
#                                                          min_overlaps = 0:10
#                                                          )
#                             )




CRISPR_df <- merged_CRISPRko_df
combined_IDs <- unique(CRISPR_df[CRISPR_df[, "Is_control"] %in% "No", "Combined_ID"])
controls_df <- NonOverlappingDfForControls(CRISPR_df)


cl <- parallel:::makeCluster(8)
parallel:::clusterExport(cl, list("SortCombinations", "CreateCombinations", "MessageID",
                                  "ReorderSubDfByLocation", "NumHomologousPairs", "SplitIntoSubstrings",
                                  "CRISPR_df",
                                  "preferred_AF_max_column", "SNP_frequency_cutoff"
                                  )
                         )
reordered_df_list <- parallel:::parLapply(cl,
                                          combined_IDs,
                                          function(x) SortCombinations(CRISPR_df[CRISPR_df[, "Combined_ID"] == x, , drop = FALSE],
                                                                       min_overlaps    = 0:10,
                                                                       min_spaces      = 50L,
                                                                       only_top_24_GPP = FALSE
                                                                       )
                                          )
parallel:::stopCluster(cl)


results_df <- do.call(rbind.data.frame, c(reordered_df_list,
                                          list(controls_df),
                                          list(stringsAsFactors = FALSE, make.row.names = FALSE)
                                          )
                      )

merged_CRISPRko_df <- results_df







# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory,
                      "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"
                      )
     )









