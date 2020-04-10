### 9th April 2020 ###




# Define functions --------------------------------------------------------

ReplaceUnspacedGuides <- function(CRISPR_df) {

  ## Prepare for identifying problematic genes
  are_controls <- !(CRISPR_df[["Is_control"]] == "No")
  combined_IDs_vec <- CRISPR_df[["Combined_ID"]][!(are_controls)]
  combined_IDs_fac <- factor(combined_IDs_vec, levels = unique(combined_IDs_vec))
  CheckThatFactorIsInOrder(combined_IDs_fac)


  ## Identify genes that are to be replaced
  are_spaced <- tapply(CRISPR_df[["Spacing"]][!(are_controls)],
                       combined_IDs_fac,
                       function(x) any(x %in% c(12L, 50L))
                       )
  message(paste0(sum(!(are_spaced)), " problematic genes were identified (without a valid 4sg combination)!"))
  unspaced_gene_IDs <- names(which(!(are_spaced)))



  ## Subset CRISPR_df to include only the problematic genes
  drop_columns <- c("Best_combination_rank", "Spacing", "Overlaps_tolerance",
                    "Num_overlaps", "Original_rank", "Rank"
                    )
  CRISPR_sub_df <- CRISPR_df[CRISPR_df[["Combined_ID"]] %in% unspaced_gene_IDs, !(colnames(CRISPR_df) %in% drop_columns)]
  CRISPR_sub_df <- RankCRISPRDf(CRISPR_sub_df, ID_column = "Combined_ID")



  ## Pick 4 guides for problematic genes, without reference to the TSS
  picked_df_list <- lapply(unspaced_gene_IDs,
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


  ## Check for genes that are *STILL* problematic
  are_now_spaced <- vapply(picked_df_list, function(x) any(x[["Spacing"]] %in% c(12L, 50L)), logical(1))
  message(paste0(sum(!(are_now_spaced)), " genes were *STILL* problematic (and were not replaced)!"))
  are_to_replace <- !(are_spaced)
  are_to_replace[are_to_replace] <- are_now_spaced


  ## Merge with the existing data frame
  old_df_list <- split(CRISPR_df[!(are_controls), ], combined_IDs_fac)
  old_df_list[are_to_replace] <- picked_df_list[are_now_spaced]

  combined_df <- do.call(rbind.data.frame,
                         c(old_df_list,
                           list(CRISPR_df[are_controls, ]),
                           list(stringsAsFactors = FALSE, make.row.names = FALSE)
                           )
                         )
  return(combined_df)
}


