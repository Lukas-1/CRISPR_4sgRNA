### 20th July 2020 ###




# Define functions --------------------------------------------------------


PrepareTSSDf <- function(TSS_df,
                         distance_before             = 1001L,
                         distance_after              = 1001L,
                         only_consistent_chromosomes = TRUE,
                         only_protein_coding         = FALSE,
                         only_best_TSS               = FALSE,
                         check_entrezs_and_symbols   = TRUE
                         ) {

  all_chromosomes <- paste0("chr", c(1:23, "X", "Y", "M"))
  stopifnot(all(TSS_df[["Chromosome"]] %in% all_chromosomes))

  if (check_entrezs_and_symbols) {
    num_symbols_vec <- lengths(strsplit(TSS_df[["Gene_symbol"]], ", ", fixed = TRUE))
    num_entrezs_vec <- lengths(strsplit(TSS_df[["Entrez_ID"]], ", ", fixed = TRUE))
    are_discrepant <- num_symbols_vec != num_entrezs_vec
    if (any(are_discrepant)) {
      mapped_df <- MapToEntrezs(entrez_IDs_vec = TSS_df[["Entrez_ID"]][are_discrepant])
      stopifnot(identical(mapped_df[["Gene_symbol"]],
                          TSS_df[["Gene_symbol"]][are_discrepant]
                          )
                )
      stopifnot(identical(lengths(strsplit(mapped_df[["Entrez_ID"]], ", ", fixed = TRUE)),
                          num_symbols_vec[are_discrepant]
                          )
                )
      TSS_df[["Entrez_ID"]][are_discrepant] <- mapped_df[["Entrez_ID"]]
    }
  }

  are_selected <- rep(TRUE, nrow(TSS_df))

  if (only_consistent_chromosomes) {
    are_selected[are_selected] <- TSS_df[["Has_consistent_chromosome"]][are_selected]
  }
  if (only_best_TSS) {
    are_selected[are_selected] <- TSS_df[["Is_chosen_TSS"]][are_selected]
  }
  if (only_protein_coding) {
    are_selected[are_selected] <- grepl("protein-coding", TSS_df[["Gene_type"]][are_selected], fixed = TRUE)
  }

  names(TSS_df)[names(TSS_df) == "Entrez_ID"]       <- "Entrez_IDs"
  names(TSS_df)[names(TSS_df) == "Gene_symbol"]     <- "Gene_symbols"
  names(TSS_df)[names(TSS_df) == "Original_symbol"] <- "Original_symbols"

  TSS_df[["Number_of_Entrez_IDs"]] <- lengths(strsplit(TSS_df[["Entrez_IDs"]], ", ", fixed = TRUE))

  select_columns <- c("Entrez_IDs", "Number_of_Entrez_IDs",
                      "Gene_symbols", "Original_symbols", "Gene_type",
                      "Source", "Score", "Is_main_TSS", "TSS",
                      "Entrez_chromosome", "Chromosome", "Strand"
                      )
  if (only_consistent_chromosomes) {
    select_columns <- setdiff(select_columns, "Entrez_chromosome")
  }

  results_df <- data.frame(
    TSS_df[are_selected, select_columns],
    "Start" = TSS_df[["TSS"]][are_selected] - distance_before,
    "End"   = TSS_df[["TSS"]][are_selected] + distance_after,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



FindNearbyTSSs <- function(ranges_df,
                           distance_before     = 1001L,
                           distance_after      = 1001L,
                           only_protein_coding = FALSE,
                           only_best_TSS       = FALSE
                           ) {



}


FindOverlappingGenes <- function(ranges_df,
                                 genes_df,
                                 only_protein_coding = FALSE,
                                 exclude_pseudogenes = FALSE
                                 ) {

}