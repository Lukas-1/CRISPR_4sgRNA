### 1st March 2020 ###




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

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRa_RData_directory, "18) Pick 4 guides per TSS.RData"))
load(file.path(CRISPRa_RData_directory, "19) Pick 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"))





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




# Pick a new set of 4 guides for genes with invalid 4sgs ------------------

merged_replaced_CRISPRa_df <- ReplaceUnspacedGuides(merged_replaced_CRISPRa_df)
lax_CRISPRa_df <- ReplaceUnspacedGuides(lax_CRISPRa_df)






# # Examine an example gene -------------------------------------------------
#
# show_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "sgRNA_sequence",
#                   "Chromosome", "Strand", "Start", "End", "Cut_location",
#                   "TSS_number", "Allocated_TSS", "Num_TSSs",
#                   "TSS_ID", "AltTSS_ID", "Rank", "Best_combination_rank", "Spacing",
#                   "Overlaps_tolerance", "Num_overlaps", "Original_rank"
#                   )
#
# lax_CRISPRa_df[lax_CRISPRa_df[["Gene_symbol"]] %in% "TRIM74", show_columns][1:10, ]
# merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Gene_symbol"]] %in% "TRIM74", show_columns][1:10, ]
#
#
#
#
#
# # Prepare for isolating problematic genes ---------------------------------
#
# are_controls <- !(merged_replaced_CRISPRa_df[["Is_control"]] == "No")
# combined_IDs_vec <- merged_replaced_CRISPRa_df[["Combined_ID"]][!(are_controls)]
# combined_IDs_fac <- factor(combined_IDs_vec, levels = unique(combined_IDs_vec))
#
# CheckThatFactorIsInOrder(combined_IDs_fac)
#
#
#
#
# # Identify genes that are to be replaced ----------------------------------
#
# are_spaced <- tapply(merged_replaced_CRISPRa_df[["Spacing"]][!(are_controls)],
#                      combined_IDs_fac,
#                      function(x) any(x %in% c(12L, 50L))
#                      )
#
# unspaced_gene_IDs <- names(which(!(are_spaced)))
#
#
#
#
#
# # Examine and compare the identified genes --------------------------------
#
# length(unspaced_gene_IDs)
# length(intersect(unspaced_gene_IDs, collected_entrez_IDs))
#
# # untargetable_annotations <- c("Not protein-coding", "Only annotated on alternate loci", "Not in current annotation release")
# # are_targetable <- !(sgRNAs_overview_df[["Gene_annotation_status"]] %in% untargetable_annotations)
# # are_unspaced_overview <- sgRNAs_overview_df[["Spacing"]] %in% "None"
# #
# # unspaced_overview_entrezs <- sgRNAs_overview_df[["Entrez_ID"]][are_targetable & are_unspaced_overview]
# #
# # setdiff(intersect(unspaced_gene_IDs, collected_entrez_IDs), unspaced_overview_entrezs)
#
#
#
#
#
# # Replace problematic genes -----------------------------------------------
#
# are_to_replace <- merged_replaced_CRISPRa_df[["Combined_ID"]] %in% unspaced_gene_IDs
#
# stopifnot(identical(merged_replaced_CRISPRa_df[["Combined_ID"]], lax_CRISPRa_df[["Combined_ID"]]))
#
# merged_replaced_CRISPRa_df[["Relaxed_location"]] <- "No"
# lax_CRISPRa_df[["Relaxed_location"]] <- "Yes"
#
# strict_df_splits <- split(merged_replaced_CRISPRa_df[!(are_controls), ],
#                           combined_IDs_fac[!(are_controls)]
#                           )
# lax_df_splits <- split(lax_CRISPRa_df[!(are_controls), !(grepl("_strict$", colnames(lax_CRISPRa_df)))],
#                        combined_IDs_fac[!(are_controls)]
#                        )
#
# stopifnot(identical(names(strict_df_splits), names(are_spaced)))
#
# combined_df_splits <- lapply(seq_along(lax_df_splits), function(x) {
#   if (are_spaced[[x]]) {
#     strict_df_splits[[x]]
#   } else {
#     lax_df_splits[[x]]
#   }
# })
#
# combined_df <- do.call(rbind.data.frame,
#                        c(combined_df_splits,
#                          list(merged_replaced_CRISPRa_df[are_controls, ]),
#                          list(stringsAsFactors = FALSE, make.row.names = FALSE)
#                          )
#                        )
#
# stopifnot(identical(combined_df[["Combined_ID"]], merged_replaced_CRISPRa_df[["Combined_ID"]]))
#
# merged_replaced_CRISPRa_df <- combined_df







# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "20) For problematic genes, pick 4 guides without reference to the TSS - merged_replaced_CRISPRa_df.RData"
                      )
     )

save(list = "lax_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "20) For problematic genes, pick 4 guides without reference to the TSS - lax_CRISPRa_df.RData"
                      )
     )





