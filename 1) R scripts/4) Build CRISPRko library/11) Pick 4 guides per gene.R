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

# SortCombinations(merged_CRISPRko_df[merged_CRISPRko_df[, "Gene_symbol"] %in% "LOC102723382", ])

merged_CRISPRko_df <- PrioritizeNonOverlapping(merged_CRISPRko_df, ID_column = "Combined_ID", parallel_mode = TRUE)





# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory,
                      "11) Pick 4 guides per gene.RData"
                      )
     )








# selected_combined_IDs <- unique(merged_CRISPRko_df[["Combined_ID"]])[1:1000]
# library("microbenchmark")
# microbenchmark(
#   selected_IDs_old_df <- PrioritizeNonOverlapping(merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% selected_combined_IDs, ],
#                                                   ID_column = "Combined_ID", parallel_mode = FALSE
#                                                   ),
#   selected_IDs_new_df <- PrioritizeNonOverlapping(merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% selected_combined_IDs, ],
#                                                   ID_column = "Combined_ID", parallel_mode = FALSE, max_overlaps = Inf
#                                                   ),
#   times = 1
# )
#
#
# identical(selected_IDs_old_df[4305:4310, colnames(selected_IDs_old_df) != "Overlaps_tolerance"],
#           selected_IDs_new_df[4305:4310, colnames(selected_IDs_new_df) != "Overlaps_tolerance"]
#           )
#
#
#
# show_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "Symbol_overlapping_0MM",
#                   "sgRNA_sequence",
#                   "Entrez_chromosome", "Chromosome", "Strand",
#                   "CRISPOR_3MM_specificity", "CRISPOR_Num_0MM", "CRISPOR_Num_1MM",
#                   "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
#                   "Rank", "Best_combination_rank", "Spacing",# "Overlaps_tolerance",
#                   "Num_overlaps", "Original_rank"
#                   )
#
# identical(selected_IDs_old_df[, !(colnames(selected_IDs_old_df) %in% c("Overlaps_tolerance", "Best_combination_rank"))],
#           selected_IDs_new_df[, !(colnames(selected_IDs_old_df) %in% c("Overlaps_tolerance", "Best_combination_rank"))]
#           )
#
#
#
# cbind(selected_IDs_old_df[4305:4310, 159],
#   selected_IDs_new_df[4305:4310, 159]
# )
#
#
#
#
#
# microbenchmark(
#   selected_IDs_old_df <- SortCombinations(merged_CRISPRko_df[merged_CRISPRko_df[, "Gene_symbol"] %in% "A1BG", ]),
#   selected_IDs_new_df <- SortCombinations(merged_CRISPRko_df[merged_CRISPRko_df[, "Gene_symbol"] %in% "A1BG", ], max_overlaps = Inf),
#   times = 1
# )
#
#
# microbenchmark(
#   selected_IDs_old_df <- PrioritizeNonOverlapping(merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% selected_combined_IDs, ]),
#   selected_IDs_new_df <- PrioritizeNonOverlapping(merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% selected_combined_IDs, ], max_overlaps = 99),
#   times = 3
# )



