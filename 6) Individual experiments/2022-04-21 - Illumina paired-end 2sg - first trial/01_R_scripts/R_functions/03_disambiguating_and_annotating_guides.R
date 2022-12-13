## 2022-05-03


# Define functions --------------------------------------------------------

FindAffectedGenes <- function(CRISPR_df, sg_number) {
  use_columns <- c("Entrez_IDs", "Gene_symbols", paste0(c("Locations_0MM", "Num_0MM"), "_sg", sg_number))
  input_df <- CRISPR_df[, use_columns]
  names(input_df) <- c("Entrez_ID", "Gene_symbol", "Locations_0MM", "Num_0MM")
  for (location_column in c("Chromosome", "Strand", "Start", "End")) {
    input_df[, location_column] <- NA
  }
  input_df[, "Entrez_ID"] <- gsub(", ", "/", input_df[, "Entrez_ID"], fixed = TRUE)
  nearby_list <- AlignSummaryDf(FindNearbyTSSs, input_df, all_TSS_df)
  affected_df <- nearby_list[["summary_df"]][, c("Affected_Entrez_IDs", "Affected_gene_symbols", "Affects_intended_main_TSS")]
  names(affected_df) <- paste0(names(affected_df), paste0("_sg", sg_number))
  return(affected_df)
}



Disambiguate_IDs <- function(input_df, from_symbols_column, sg1_column, sg2_column) {
  from_symbols_splits <- strsplit(input_df[, from_symbols_column], ", ", fixed = TRUE)
  affected_sg1_splits <- strsplit(input_df[, sg1_column], "[,;] ")
  affected_sg2_splits <- strsplit(input_df[, sg2_column], "[,;] ")
  results_list <- mapply(function(x, y, z) {
    results_vec <- intersect(x, c(y, z))
    if (length(results_vec) == 0) {
      results_vec <- x
    }
    return(results_vec)
  }, from_symbols_splits, affected_sg1_splits, affected_sg2_splits, SIMPLIFY = FALSE)
  results_vec <- sapply(results_list, "[[", 1)
  return(results_vec)
}



GetGCcontent <- function(char_vec) {
  char_vec <- toupper(char_vec)
  if (all(substr(char_vec, 1, 1) == "G")) {
    char_vec <- substr(char_vec, 2, nchar(char_vec))
  }
  char_mat <- do.call(rbind, strsplit(char_vec, "", fixed = TRUE))
  are_GC_mat <- (char_mat == "G") | (char_mat == "C")
  num_GC_vec <- as.integer(rowSums(are_GC_mat))
  return(num_GC_vec)
}



ChooseBetweenMultiplePlasmids <- function(input_df) {
  num_plasmids <- table(input_df[, "Entrez_ID"])[as.character(input_df[, "Entrez_ID"])]
  are_first_vec <- ifelse(num_plasmids %in% 1, TRUE, NA)
  are_preferred_vec <- are_first_vec
  TSS_columns <- paste0("Affects_intended_main_TSS_sg", 1:2)
  target_main_TSS <- rowSums(as.matrix(input_df[, TSS_columns]))
  duplicated_entrezs <- unique(input_df[, "Entrez_ID"][(num_plasmids > 1) %in% TRUE])
  for (this_entrez in duplicated_entrezs) {
    are_this_entrez <- input_df[, "Entrez_ID"] %in% this_entrez
    this_seq <- seq_len(sum(are_this_entrez))
    are_first_vec[are_this_entrez] <- this_seq == 1
    are_preferred_vec[are_this_entrez] <- this_seq == which.max(target_main_TSS[are_this_entrez])
  }
  input_df[, "Num_plasmids_for_Entrez"] <- num_plasmids
  input_df[, "Is_first_plasmid"] <- are_first_vec
  input_df[, "Is_preferred_plasmid"] <- are_preferred_vec
  return(input_df)
}



