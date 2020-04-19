### 25th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "18) Pick 4 guides per TSS.RData"))








# Replace genes with overlaps >7bp ----------------------------------------

CRISPR_df <- merged_replaced_CRISPRa_df

## Prepare for identifying problematic genes
are_controls <- !(CRISPR_df[["Is_control"]] == "No")
combined_IDs_vec <- CRISPR_df[["Combined_ID"]][!(are_controls)]
combined_IDs_fac <- factor(combined_IDs_vec, levels = unique(combined_IDs_vec))
CheckThatFactorIsInOrder(combined_IDs_fac)


## Identify genes that are to be replaced
are_top4 <- CRISPR_df[["Rank"]] %in% 1:4
sgRNA_splits <- split(CRISPR_df[["sgRNA_sequence"]][are_top4],
                      factor(CRISPR_df[["AltTSS_ID"]][are_top4], levels = unique(CRISPR_df[["AltTSS_ID"]][are_top4]))
                      )
homologies_vec <- vapply(sgRNA_splits, LongestSharedSubsequence, integer(1))

problematic_TSS_IDs <- names(sgRNA_splits)[homologies_vec >= 8]
problematic_combined_IDs <- unique(CRISPR_df[["Combined_ID"]][CRISPR_df[["AltTSS_ID"]] %in% problematic_TSS_IDs])

message(paste0(length(problematic_combined_IDs), " problematic genes were identified ",
               "(without a valid 4sg combination)!"
               )
        )

## Subset CRISPR_df to include only the problematic genes
drop_columns <- c("Best_combination_rank", "Spacing", "Overlaps_tolerance",
                  "Num_overlaps", "Original_rank", "Rank"
                  )
CRISPR_sub_df <- CRISPR_df[CRISPR_df[["Combined_ID"]] %in% problematic_combined_IDs, !(colnames(CRISPR_df) %in% drop_columns)]
CRISPR_sub_df <- RankCRISPRDf(CRISPR_sub_df, ID_column = "Combined_ID")



## Pick 4 guides for problematic genes
picked_df_list <- lapply(problematic_combined_IDs,
                         function(x) SortCombinations(CRISPR_sub_df[CRISPR_sub_df[["Combined_ID"]] == x, , drop = FALSE],
                                                      tolerate_divergent_chromosomes  = TRUE
                                                      )
                         )
picked_df_list <- lapply(picked_df_list, function(x) {
  x[["TSS_number"]]    <- NA_integer_
  x[["Allocated_TSS"]] <- NA_character_
  x[["Num_TSSs"]]      <- 1L
  x[["TSS_ID"]]        <- NA_character_
  x[["AltTSS_ID"]]     <- x[["Combined_ID"]]
  return(x)
})



## Merge with the existing data frame

old_df_list <- split(CRISPR_df[!(are_controls), ], combined_IDs_fac)

are_to_replace <- vapply(old_df_list, function(x) unique(x[["Combined_ID"]] %in% problematic_combined_IDs), logical(1))
old_df_list[are_to_replace] <- picked_df_list

combined_df <- do.call(rbind.data.frame,
                       c(old_df_list,
                         list(CRISPR_df[are_controls, ]),
                         list(stringsAsFactors = FALSE, make.row.names = FALSE)
                         )
                       )









# Save data ---------------------------------------------------------------

save(list = "combined_df",
     file = file.path(CRISPRa_RData_directory,
                      "18.5) Replace problematic genes.RData"
                      )
     )







