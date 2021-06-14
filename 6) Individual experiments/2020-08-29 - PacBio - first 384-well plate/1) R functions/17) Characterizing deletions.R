### 2nd June 2021 ###



# Define functions --------------------------------------------------------

AnnotateDeletions <- function(deletions_df) {

  stopifnot("features_list" %in% ls(envir = globalenv()))

  ### Prepare the features data frame
  use_features_df <- FeaturesListToDf(features_list)

  are_promoter <- grepl("promoter", use_features_df[["Feature"]])
  promoter_splits <- strsplit(use_features_df[["Feature"]][are_promoter], "_", fixed = TRUE)
  use_features_df[["Feature"]][are_promoter] <- sapply(promoter_splits, "[[", 1)

  are_sg_cr <- use_features_df[["Feature"]] %in% paste0("sg", 1:4, "_cr", 1:4)
  use_features_df[["Feature"]][are_sg_cr] <- paste0("sg_cr_", 1:4)


  ### Find spanning deletions
  deletions_df[["Span_tracrRNAs"]] <- FindSpanningDeletions(deletions_df, use_features_df, "tracrRNA")
  deletions_df[["Span_promoters"]] <- FindSpanningDeletions(deletions_df, use_features_df, "promoter")
  deletions_df[["Span_sg_cr"]]     <- FindSpanningDeletions(deletions_df, use_features_df, "sg_cr_")
  deletions_df[["Span_sgRNAs"]]    <- FindSpanningDeletions(deletions_df, use_features_df, "sg")
  return(deletions_df)
}


CheckThatFactorIsInOrder <- function(my_factor) {
  stopifnot(identical(length(unique(my_factor)),
                      length(rle(as.integer(my_factor))[["lengths"]])
                      )
            )
}



FindLargestDeletion <- function(deletions_df, set_seed = TRUE) {
  zmw_fac <- factor(deletions_df[["ZMW"]])
  CheckThatFactorIsInOrder(zmw_fac)
  if (set_seed) {
    set.seed(1)
  }
  are_largest_random <- tapply(deletions_df[["Deletion_size"]], zmw_fac, function(x) {
    are_max <- x == max(x)
    if (sum(are_max) > 1) {
      candidates_vec <- which(are_max)
      use_index <- sample(candidates_vec, 1)
      are_max <- seq_along(x) == use_index
    }
    return(are_max)
    stopifnot(sum(are_max) == 1)
    return(are_max)
  }, simplify = FALSE)
  deletions_df[["Is_largest_selected"]] <- unlist(are_largest_random, use.names = FALSE)
  are_largest <- tapply(deletions_df[["Deletion_size"]], zmw_fac, function(x) {
    x == max(x)
  }, simplify = FALSE)
  deletions_df[["Is_largest"]] <- unlist(are_largest, use.names = FALSE)
  return(deletions_df)
}


PrioritizeDeletions <- function(deletions_df, set_seed = TRUE) {
  zmw_fac <- factor(deletions_df[["ZMW"]])
  wells_fac <- factor(deletions_df[["Combined_ID"]])
  CheckThatFactorIsInOrder(zmw_fac)
  CheckThatFactorIsInOrder(wells_fac)
  if (set_seed) {
    set.seed(1)
  }
  random_del_ranks <- tapply(deletions_df[["Deletion_size"]], zmw_fac, function(x) {
    sample(seq_along(x))
  }, simplify = FALSE)
  random_well_ranks <- tapply(deletions_df[["Is_largest_selected"]], wells_fac, function(x) {
    results_vec <- rep(NA, length(x))
    results_vec[x] <- sample(seq_len(sum(x)))
    return(results_vec)
  }, simplify = FALSE)
  deletions_df[["ZMW_random_rank"]]  <- unlist(random_del_ranks,  use.names = FALSE)
  deletions_df[["Well_random_rank"]] <- unlist(random_well_ranks, use.names = FALSE)
  return(deletions_df)
}



FindWithin <- function(positions_vec, features_df, feature) {
  are_within_mat <- do.call(cbind, lapply(1:4, function(x) {
    this_index <- which(features_df[["Feature"]] == paste0(feature, x))
    are_within <- (positions_vec >= features_df[["Start"]][[this_index]]) &
                  (positions_vec <= features_df[["End"]][[this_index]])
    return(are_within)
  }))
  return(are_within_mat)
}




FindSpanningDeletions <- function(deletions_df, features_df, use_feature) {

  start_in_feature_mat <- FindWithin(deletions_df[["Deletion_start"]],
                                     features_df,
                                     use_feature
                                     )
  end_in_feature_mat   <- FindWithin(deletions_df[["Deletion_end"]],
                                     features_df,
                                     use_feature
                                     )

  stopifnot(all(rowSums(start_in_feature_mat) %in% c(0, 1)))
  stopifnot(all(rowSums(end_in_feature_mat) %in% c(0, 1)))

  span_feature <- (start_in_feature_mat[, 1] &
                   (rowSums(end_in_feature_mat[, 2:4]) == 1)) |
                  (start_in_feature_mat[, 2] &
                   (rowSums(end_in_feature_mat[, 3:4]) == 1)) |
                  (start_in_feature_mat[, 3] &
                   end_in_feature_mat[, 4])

  return(span_feature)
}






CompileDeletions <- function(alignments_df,
                             ID_column = "Well_number",
                             unique_IDs = seq_len(384)
                             ) {

  deletion_size <- 20L
  deletion_string <- paste0(rep("-", deletion_size), collapse = "")
  have_deletion <- grepl(deletion_string, alignments_df[["Aligned_read"]], fixed = TRUE)

  well_df_list <- lapply(seq_along(unique_IDs), function(x) {

    current_ID <- unique_IDs[[x]]

    message(paste0("Processing the well: ", current_ID, "..."))

    are_this_ID <- alignments_df[[ID_column]] == current_ID
    these_have_deletion <- have_deletion[are_this_ID]
    if (!(any(these_have_deletion))) {
      return(NULL)
    }

    aligned_plasmid_vec <- alignments_df[["Aligned_plasmid"]][are_this_ID]
    aligned_read_vec <- alignments_df[["Aligned_read"]][are_this_ID]

    aligned_plasmid_char_list <- strsplit(aligned_plasmid_vec, "")
    aligned_read_char_list <- strsplit(aligned_read_vec, "")

    this_ID_zmws <- alignments_df[["ZMW"]][are_this_ID]
    num_reads <- length(this_ID_zmws)

    if (is.character(current_ID)) {
      use_index <- x
    } else {
      use_index <- current_ID
    }

    deletions_mat_list <- lapply(which(these_have_deletion), function(x) {

      plasmid_are_gap <- aligned_plasmid_char_list[[x]] == "-"
      plasmid_char_numbers <- cumsum(!(plasmid_are_gap))

      read_are_gap <- aligned_read_char_list[[x]] == "-"
      read_char_numbers <- cumsum(!(read_are_gap))

      gap_rle <- rle(read_are_gap)

      ends <- cumsum(gap_rle[["lengths"]])
      starts <- c(1L, ends[seq_len(length(ends) - 1)] + 1L)

      are_gaps <- (gap_rle[["lengths"]] >= deletion_size) &
                  (gap_rle[["values"]] == TRUE)

      deletion_starts <- plasmid_char_numbers[starts[are_gaps]]
      deletion_ends <- plasmid_char_numbers[ends[are_gaps]]

      deletions_mat <- cbind(
        "ZMW"            = this_ID_zmws[[x]],
        "Deletion_start" = deletion_starts,
        "Deletion_end"   = deletion_ends
      )
      return(deletions_mat)
    })
    well_mat <- do.call(rbind, deletions_mat_list)

    well_df <- data.frame(
      "Combined_ID" = current_ID,
      well_mat,
      stringsAsFactors = FALSE
    )
    return(well_df)
  })
  message("Collating the final data frame...")
  results_df <- do.call(rbind.data.frame, c(well_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  results_df[["Deletion_size"]] <- results_df[["Deletion_end"]] - results_df[["Deletion_start"]]
  message("Finding the largest deletion for each read...")
  results_df <- FindLargestDeletion(results_df)
  message("Annotating deletions...")
  results_df <- AnnotateDeletions(results_df)
  return(results_df)
}





