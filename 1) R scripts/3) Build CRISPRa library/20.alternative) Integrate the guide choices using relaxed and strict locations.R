### 1st March 2020 ###




# Import packages and source code -----------------------------------------







# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRa_RData_directory, "18.5.2) Pick the top 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"))
load(file.path(CRISPRa_RData_directory, "19) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))








# Examine an example gene -------------------------------------------------

setdiff(colnames(lax_CRISPRa_df), colnames(merged_replaced_CRISPRa_df))

show_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "sgRNA_sequence",
                  "Chromosome", "Strand", "Start", "End", "Cut_location",
                  "TSS_number", "Allocated_TSS", "Num_TSSs",
                  "TSS_ID", "AltTSS_ID", "Rank", "Best_combination_rank", "Spacing",
                  "Overlaps_tolerance", "Num_overlaps", "Original_rank"
                  )
lax_CRISPRa_df[lax_CRISPRa_df[["Gene_symbol"]] %in% "TRIM74", show_columns][1:10, ]

merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Gene_symbol"]] %in% "TRIM74", show_columns][1:10, ]





# Prepare for isolating problematic genes ---------------------------------

are_controls <- !(merged_replaced_CRISPRa_df[["Is_control"]] == "No")
combined_IDs_vec <- merged_replaced_CRISPRa_df[["Combined_ID"]][!(are_controls)]
combined_IDs_fac <- factor(combined_IDs_vec, levels = unique(combined_IDs_vec))

stopifnot(identical(length(unique(combined_IDs_fac)),
                    length(rle(as.integer(combined_IDs_fac))[["lengths"]])
                    )
          )





# Identify genes that are to be replaced ----------------------------------

are_spaced <- tapply(merged_replaced_CRISPRa_df[["Spacing"]][!(are_controls)],
                     combined_IDs_fac,
                     function(x) any(x %in% c(12L, 50L))
                     )

unspaced_gene_IDs <- names(which(!(are_spaced)))





# Examine and compare the identified genes --------------------------------

length(unspaced_gene_IDs)
length(intersect(unspaced_gene_IDs, collected_entrez_IDs))

untargetable_annotations <- c("Not protein-coding", "Only annotated on alternate loci", "Not in current annotation release")
are_targetable <- !(sgRNAs_overview_df[["Gene_annotation_status"]] %in% untargetable_annotations)
are_unspaced_overview <- sgRNAs_overview_df[["Spacing"]] %in% "None"

unspaced_overview_entrezs <- sgRNAs_overview_df[["Entrez_ID"]][are_targetable & are_unspaced_overview]

setdiff(intersect(unspaced_gene_IDs, collected_entrez_IDs), unspaced_overview_entrezs)




# Replace problematic genes -----------------------------------------------

merged_replaced_CRISPRa_df[["Relaxed_location"]] <- "No"
lax_CRISPRa_df[["Relaxed_location"]] <- "Yes"

strict_df_splits <- split(merged_replaced_CRISPRa_df[!(are_controls), ],
                          combined_IDs_fac[!(are_controls)]
                          )
lax_df_splits <- split(lax_CRISPRa_df, factor(lax_CRISPRa_df[["Combined_ID"]], levels = unique(lax_CRISPRa_df[["Combined_ID"]])))


combined_df_splits <- lapply(levels(combined_IDs_fac), function(x) {
  if (x %in% unspaced_gene_IDs) {
    lax_df_splits[[x]][, colnames(merged_replaced_CRISPRa_df)]
  } else {
    strict_df_splits[[x]]
  }
})

combined_df <- do.call(rbind.data.frame,
                       c(combined_df_splits,
                         list(merged_replaced_CRISPRa_df[are_controls, ]),
                         list(stringsAsFactors = FALSE, make.row.names = FALSE)
                         )
                       )

stopifnot(identical(combined_df[["Combined_ID"]], merged_replaced_CRISPRa_df[["Combined_ID"]]))

merged_replaced_CRISPRa_df <- combined_df





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "18.6.2) Integrate the guide choices using relaxed and strict locations.RData"
                      )
     )

