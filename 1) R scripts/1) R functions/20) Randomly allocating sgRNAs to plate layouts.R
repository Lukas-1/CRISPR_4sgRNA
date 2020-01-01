### 1st January 2020 ###



# Define functions --------------------------------------------------------


AreCompleteTranscripts <- function(CRISPR_df) {
  are_incomplete <- rep.int(NA, nrow(CRISPR_df))
  unique_TSS_IDs <- unique(CRISPR_df[CRISPR_df[, "Is_control"] == "No", "AltTSS_ID"])
  for (unique_TSS_ID in unique_TSS_IDs) {
    are_this_TSS <- CRISPR_df[, "AltTSS_ID"] %in% unique_TSS_ID
    are_unmapped_TSS <- all(CRISPR_df[are_this_TSS, "Num_TSSs"] >= 2) && all(is.na(CRISPR_df[are_this_TSS, "Start"]))
    if (!(are_unmapped_TSS)) { # AltTSS_IDs corresponding to unmapped sgRNAs are left to be 'NA'
      is_complete <- all(1:4 %in% CRISPR_df[are_this_TSS, "Rank"])
      are_incomplete[are_this_TSS] <- is_complete
    }
  }
  return(are_incomplete)
}



RandomizeAllIndices <- function(n_total = NULL, n_per_plate_vec = NULL, n_per_plate = 384L) {

  if (is.null(n_total) && is.null(n_per_plate_vec)) {
    stop("Either n_total or n_per_plate_vec must be specified!")
  }

  if (is.null(n_per_plate_vec)) {
    num_full_plates <- floor(n_total / n_per_plate)
    last_plate_n <- n_total - (num_full_plates * n_per_plate)
    n_per_plate_vec <- rep(n_per_plate, num_full_plates)
    if (last_plate_n > 0) {
      n_per_plate_vec <- c(n_per_plate_vec, last_plate_n)
    }
  } else {
    n_total <- sum(n_per_plate_vec)
  }

  indices_pool <- seq_len(n_total)
  num_plates <- length(n_per_plate_vec)
  results_list <- vector(mode = "list", length = num_plates)

  for (plate_number in seq_len(num_plates - 1L)) {
    this_sample <- sample(indices_pool, n_per_plate_vec[[plate_number]])
    results_list[[plate_number]] <- this_sample
    indices_pool <- indices_pool[!(indices_pool %in% this_sample)]
  }
  results_list[[num_plates]] <- sample(indices_pool, n_per_plate_vec[[num_plates]])
  return(results_list)
}


RetrieveIndices <- function(CRISPR_df, indices_vec, ID_column) {
  unique_combined_IDs <- unique(CRISPR_df[, ID_column])
  combined_ID_matches <- match(CRISPR_df[, ID_column], unique_combined_IDs)
  results_df <- CRISPR_df[combined_ID_matches %in% indices_vec, ]
  return(results_df)
}







