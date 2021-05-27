### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("ShortRead") # For processing the quality scores






# Functions for subsampling reads (for simulations) -----------------------

ProcessWithSubsampling <- function(ccs_df,
                                   barcodes_df,
                                   extracted_df,
                                   wells_vec = seq_len(384),
                                   use_fractions = c(1, 0.75, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.025, 0.01),
                                   num_repetitions = 10L
                                   ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  set.seed(1)

  results_list <- lapply(use_fractions, function(use_fraction) {

    num_reads <- round(nrow(ccs_df) * use_fraction)
    if (use_fraction == 1) {
      reps_vec <- 1
    } else {
      reps_vec <- seq_len(num_repetitions)
    }

    rep_list <- lapply(reps_vec, function(i) {

      sample_indices <- sample(seq_len(nrow(ccs_df)), size = num_reads)
      subsampled_ccs_df <- ccs_df[sample_indices, ]
      row.names(subsampled_ccs_df) <- NULL

      if ("Plate_number" %in% names(ccs_df)) {
        use_unique_IDs <- sg_sequences_df[["Combined_ID"]]
        use_ID_column <- "Combined_ID"
        subsampled_ccs_df[["Passed_filters"]] <- subsampled_ccs_df[["Plate_passed_filters"]] &
                                                (subsampled_ccs_df[["Well_passed_filters"]] %in% TRUE)
        ccs3_lima_zmws <- GetCCS3_ZMWs(subsampled_ccs_df)
        analysis_list <- AnalyzePlates(subsampled_ccs_df,
                                       sg_sequences_df,
                                       barcodes_df,
                                       extracted_df,
                                       set_seed = FALSE
                                       )
      } else {
        use_unique_IDs <- wells_vec
        use_ID_column <- "Well_number"
        ccs3_lima_zmws <- NULL
        analysis_list <- AnalyzeWells(subsampled_ccs_df,
                                      sg_sequences_df,
                                      barcodes_df,
                                      extracted_df,
                                      set_seed = FALSE
                                      )
      }
      ccs5_lima_zmws <- GetCCS5_ZMWs(subsampled_ccs_df, wells_vec = wells_vec)
      ccs7_lima_zmws <- GetCCS7_ZMWs(subsampled_ccs_df, wells_vec = wells_vec)

      ccs3_df_list <- SummarizeWells(analysis_list,
                                     use_zmws = ccs3_lima_zmws,
                                     unique_IDs = use_unique_IDs,
                                     ID_column = use_ID_column
                                     )
      ccs5_df_list <- SummarizeWells(analysis_list,
                                     use_zmws = ccs5_lima_zmws,
                                     unique_IDs = use_unique_IDs,
                                     ID_column = use_ID_column
                                     )
      ccs7_df_list <- SummarizeWells(analysis_list,
                                     use_zmws = ccs7_lima_zmws,
                                     unique_IDs = use_unique_IDs,
                                     ID_column = use_ID_column
                                     )
      df_list_list <- list(
        "ccs3" = ccs3_df_list,
        "ccs5" = ccs5_df_list,
        "ccs7" = ccs7_df_list
      )
      return(df_list_list)
    })
    names(rep_list) <- paste0("rep", reps_vec)
    return(rep_list)
  })
  names(results_list) <- paste0(use_fractions * 100, "% sampled")
  return(results_list)
}






# Functions for calculating distances on 384-well plates ------------------

MakeDistanceList <- function(manhattan_distance = FALSE, matrix_format = FALSE) {

  plate_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)

  current_well <- 1
  distance_list <- lapply(seq_len(384), function(current_well) {
    distance_vec <- vapply(seq_len(384), function(target_well) {

      are_current_well <- plate_mat == current_well
      are_target_well <- plate_mat == target_well

      current_row_index <- which(rowSums(are_current_well) == 1)
      target_row_index  <- which(rowSums(are_target_well) == 1)
      current_col_index <- which(colSums(are_current_well) == 1)
      target_col_index  <- which(colSums(are_target_well) == 1)

      col_distance <- abs(current_col_index - target_col_index)
      row_distance <- abs(current_row_index - target_row_index)

      if (manhattan_distance) {
        distance <- col_distance + row_distance
      } else {
        distance <- sqrt(col_distance^2 + row_distance^2) # Compute the Euclidian distance
      }
      return(distance)
    }, if (manhattan_distance) integer(1) else numeric(1))
    if (matrix_format) {
      return(matrix(distance_vec, nrow = 16, ncol = 24, byrow = TRUE))
    } else {
      return(distance_vec)
    }
  })
  return(distance_list)
}






# Functions for processing reads ------------------------------------------

ContainSequences <- function(query_seq, target_seq) {

  assign("delete_query_seq", query_seq, envir = globalenv())
  assign("delete_target_seq", target_seq, envir = globalenv())

  num_queries <- length(query_seq)
  num_targets <- length(target_seq)

  ## The functions are from the Biostrings package

  if (!("DNAStringSet" %in% class(query_seq))) {
    query_seq <- DNAStringSet(query_seq)
  }
  if (!("DNAStringSet" %in% class(target_seq))) {
    target_seq <- DNAStringSet(target_seq)
  }

  rev_search_seq <- reverseComplement(query_seq)

  contain_fwd_mat <- do.call(cbind, lapply(seq_len(num_queries), function(x) {
    vcountPattern(query_seq[[x]], target_seq)
  }))
  contain_rev_mat <- do.call(cbind, lapply(seq_len(num_queries), function(x) {
    vcountPattern(rev_search_seq[[x]], target_seq)
  }))

  are_fwd <- rowSums(contain_fwd_mat) > 0
  are_rev <- rowSums(contain_rev_mat) > 0

  num_duplicated_fwd <- sum(contain_fwd_mat > 1)
  num_duplicated_rev <- sum(contain_rev_mat > 1)

  sent_message <- FALSE

  if (num_duplicated_fwd != 0) {
    message(paste0("The same sequence was found more than once for ",
                   num_duplicated_fwd, " read",
                   if (num_duplicated_fwd == 1) "" else "s",
                   " in the forward direction!"
                   ))
    sent_message <- TRUE
  }
  if (num_duplicated_rev != 0) {
    message(paste0("The same sequence was found more than once for ",
                   num_duplicated_rev, " read",
                   if (num_duplicated_rev == 1) "" else "s",
                   " in the reverse direction!"
                   ))
    sent_message <- TRUE
  }
  # stopifnot(!(any(contain_fwd_mat == 2)))
  # stopifnot(!(any(contain_rev_mat == 2)))

  contain_mat <- (contain_fwd_mat > 0) + (contain_rev_mat > 0)

  same_guide_twice_mat <- contain_mat > 1
  have_duplication <- rowSums(same_guide_twice_mat) > 0

  if (any(have_duplication)) {
    message(paste0("For ", sum(same_guide_twice_mat), " out of 4 * ",
                   num_targets, " searches, the same query matched in both the",
                   " forward and reverse direction (an inverted duplication)..."
                   ))
    fwd_helper_mat <- contain_fwd_mat
    fwd_helper_mat[same_guide_twice_mat] <- 0L
    have_inversion <- (rowSums(contain_fwd_mat) > 0) &
                      (rowSums(contain_rev_mat) > 0)
    sent_message <- TRUE
  } else {
    have_inversion <- are_fwd & are_rev
  }
  if (any(have_inversion)) {
    message(paste0("There was a likely inversion for ", sum(have_inversion),
                   " out of ", num_targets, " reads in this well ",
                   "(some queries matched in the forward and others in ",
                   "the reverse direction)..."
                   ))
    sent_message <- TRUE
  }
  if (sent_message) {
    message("")
  }

  contain_mat[same_guide_twice_mat] <- 1L

  colnames(contain_mat) <- paste0("sequence_", seq_len(num_queries))

  are_fwd <- ifelse(are_fwd,
                    ifelse(are_rev, NA, TRUE),
                    ifelse(are_rev, FALSE, NA)
                    )

  num_correct_vec <- rowSums(contain_mat)
  contain_num_mat <- do.call(cbind,
                             lapply(seq_len(num_queries),
                                    function(x) num_correct_vec >= x
                                    )
                             )
  colnames(contain_num_mat) <- paste0("at_least_", seq_len(num_queries))

  results_mat <- cbind(contain_mat,
                       contain_num_mat,
                       "orientation_fwd_sg" = are_fwd,
                       "inconsistent_orientation" = have_duplication | have_inversion
                       )
  return(results_mat)
}




StatsForWell <- function(well_index, ccs_df, sg_df) {

  assign("delete_well_index", well_index, envir = globalenv())
  assign("delete_ccs_df", ccs_df, envir = globalenv())
  assign("delete_sg_df", sg_df, envir = globalenv())

  well_number <- sg_df[well_index, "Well_number"]

  are_this_well <- (ccs_df[["Well_number"]] %in% well_number) & ccs_df[["Passed_filters"]]
  num_reads <- sum(are_this_well)

  well_sequences <- DNAStringSet(substr(ccs_df[["Sequence"]][are_this_well],
                                        ccs_df[["Clip_start"]][are_this_well],
                                        ccs_df[["Clip_end"]][are_this_well]
                                        ))

  sg_sequences <- vapply(1:4, function(x) sg_df[well_index, paste0("sg_cr_", x)], "")
  sg_with_promoter_seq <- vapply(1:4, function(x) sg_df[well_index, paste0("sg_cr_pr_", x)], "")

  contain_sg_cr_mat   <- ContainSequences(sg_sequences, well_sequences)
  contain_sg_prom_mat <- ContainSequences(sg_with_promoter_seq, well_sequences)
  contain_plasmid_mat <- ContainSequences(sg_df[well_index, "Whole_plasmid"], well_sequences)

  colnames(contain_sg_cr_mat)[1:4] <- paste0("sg", 1:4, "_cr", 1:4)
  colnames(contain_sg_cr_mat)[colnames(contain_sg_cr_mat) == "at_least_4"] <- "all_4"
  results_mat <- cbind(
    contain_sg_cr_mat[, 1:8, drop = FALSE],
    "all_4_promoters" = contain_sg_prom_mat[, "at_least_4"],
    "whole_plasmid"   = contain_plasmid_mat[, 1],
    contain_sg_cr_mat[, "orientation_fwd_sg", drop = FALSE]
  )
  rownames(results_mat) <- NULL
  stopifnot(!(any(results_mat[, "all_4_promoters"] < results_mat[, "whole_plasmid"])))

  return(results_mat)
}




GetMode <- function(vec) {
   uniqv <- unique(vec)
   uniqv[which.max(tabulate(match(vec, uniqv)))]
}



OrientationContaminationMat <- function(query_seq, reads_seq) {

  are_usual_length <- nchar(query_seq) == GetMode(nchar(query_seq))

  PDict_object <- PDict(query_seq[are_usual_length], max.mismatch = 0)

  counts_mat <- vcountPDict(PDict_object, reads_seq)

  row_list <- vector(mode = "list", length = length(query_seq))
  row_list[are_usual_length] <- lapply(seq_len(nrow(counts_mat)),
                                       function(x) counts_mat[x, ]
                                       )

  unusual_seq <- DNAStringSet(query_seq[!(are_usual_length)])

  unusual_list <- lapply(seq_len(sum(!(are_usual_length))),
                         function(x) vcountPattern(unusual_seq[[x]], reads_seq)
                         )

  row_list[!(are_usual_length)] <- unusual_list

  all_counts_mat <- do.call(rbind, row_list)
  stopifnot(all(unique(as.vector(all_counts_mat))) %in% 0:1)
  mode(all_counts_mat) <- "integer"
  return(all_counts_mat)
}




GetContaminationMat <- function(query_seq, ccs_df, sg_df) {

  assign("delete_query_seq_1", query_seq, envir = globalenv())
  assign("delete_ccs_df", ccs_df, envir = globalenv())
  assign("delete_sg_df", sg_df, envir = globalenv())

  stopifnot(nrow(sg_df) == length(query_seq))

  have_well <- !(is.na(ccs_df[["Well_number"]])) & ccs_df[["Passed_filters"]]

  read_sequences <- DNAStringSet(substr(ccs_df[["Sequence"]][have_well],
                                        ccs_df[["Clip_start"]][have_well],
                                        ccs_df[["Clip_end"]][have_well]
                                        ))


  rev_query_seq <- as.character(reverseComplement(DNAStringSet(query_seq)))

  fwd_counts_mat <- OrientationContaminationMat(query_seq, read_sequences)
  rev_counts_mat <- OrientationContaminationMat(rev_query_seq, read_sequences)

  counts_mat <- fwd_counts_mat + rev_counts_mat

  contamination_mat <- counts_mat
  contamination_mat[contamination_mat == 2L] <- 1L

  ## If the sequence of two wells is exactly the same,
  ## perhaps it shouldn't count as a contamination...
  sg_cr <- sg_df[, paste0("sg_cr_", 1:4)]
  same_seq_list <- lapply(seq_len(nrow(sg_df)),
                          function(x) which(query_seq %in% vapply(sg_cr, function(y) y[[x]], ""))
                          )
  same_seq_wells_list <- lapply(same_seq_list, function(x) sg_df[["Well_number"]][x])
  are_this_well_list <- lapply(same_seq_wells_list,
                               function(x) ccs_df[["Well_number"]][have_well] %in% x
                               )

  for (i in seq_len(nrow(sg_df))) {
    contamination_mat[i, are_this_well_list[[i]]] <- 0L
  }
  return(contamination_mat)
}



ContaminationsSummaryDf <- function(reads_df, wells_vec, close_well_range = 3L) {

  assign("delete_reads_df", reads_df, envir = globalenv())
  assign("delete_wells_vec", wells_vec, envir = globalenv())

  wells_fac <- factor(reads_df[["Well_number"]], levels = wells_vec)

  are_contaminated <- reads_df[, "Num_contaminating_guides"] >= 1
  are_close_rand <- (reads_df[, "Random_distance"] <= close_well_range) %in% TRUE

  num_other_genes_vec <- tapply(reads_df[, "Contaminating_well"],
                                wells_fac,
                                function(x) {
                                  this_list <- strsplit(x[!(is.na(x))], ", ", fixed = TRUE)
                                  this_vec <- unlist(this_list, use.names = FALSE)
                                  length(unique(this_vec))
                                })

  num_contaminants <- tapply(are_contaminated, wells_fac, sum)
  num_close_contaminants_rand <- tapply(are_contaminated & are_close_rand, wells_fac, sum)

  plate_num_close_vec <- vapply(manhattan_dist_list,
                                function(x) sum(x <= close_well_range),
                                integer(1)
                                ) - 1L

  num_expected_close_wells <- plate_num_close_vec[wells_vec] * (num_contaminants / (length(wells_vec) - 1))

  per_well_distances_vec <- tapply(reads_df[, "Mean_distance"],
                                   wells_fac,
                                   mean,
                                   na.rm = TRUE
                                   )

  expected_distances_vec <- vapply(manhattan_dist_list,
                                   function(x) mean(x[x != 0]),
                                   numeric(1)
                                   )

  distances_list <- split(reads_df[, "Mean_distance"], wells_fac)

  p_val_vec <- vapply(seq_along(distances_list),
                      function(x) {
                        dist_vec <- distances_list[[x]]
                        if (all(is.na(dist_vec))) {
                          return(NA_real_)
                        } else {
                          mean_expected <- expected_distances_vec[[x]]
                          return(suppressWarnings(wilcox.test(x = dist_vec,
                                                              mu = mean_expected,
                                                              alternative = "less"
                                                              )[["p.value"]]
                                                  ))
                        }
                      },
                      numeric(1)
                      )

  results_df <- data.frame(
    "Num_contaminating_genes"   = num_other_genes_vec,
    "Num_contaminated_reads"    = num_contaminants,
    "Num_from_close_wells"      = num_close_contaminants_rand,
    "Expected_from_close_wells" = num_expected_close_wells,
    "Mean_distance"             = per_well_distances_vec,
    "Expected_distance"         = expected_distances_vec[wells_vec],
    "Distance_p_value"          = p_val_vec,
    stringsAsFactors = FALSE
  )
  return(results_df)
}




CreateSummaryDf <- function(reads_df,
                            filter_reads = FALSE,
                            filter_subsequences = FALSE,
                            close_well_range = 3L,
                            ID_column = "Well_number",
                            unique_IDs = seq_len(384)
                            ) {

  assign("delete_reads_df", reads_df, envir = globalenv())
  assign("delete_ID_column", ID_column, envir = globalenv())
  assign("delete_unique_IDs", unique_IDs, envir = globalenv())

  if (filter_reads) {
    reads_df <- reads_df[reads_df[["Passes_filters"]] == 1, ]
  }
  if (filter_subsequences) {
    reads_df <- reads_df[reads_df[["Passes_sg_quality"]] == 1, ]
  }

  wells_fac <- factor(reads_df[[ID_column]], levels = unique_IDs)


  ## Process 100% correct read counts and percentages

  binary_columns <- c(paste0("sg", 1:4, "_cr", 1:4),
                      paste0("at_least_", 1:3),
                      "all_4", "all_4_promoters", "whole_plasmid"
                      )

  well_stats_mat <- as.matrix(reads_df[, binary_columns])

  stats_vec_list <- tapply(seq_len(nrow(reads_df)),
                           wells_fac,
                           function(x) colSums(well_stats_mat[x, , drop = FALSE])
                           )
  are_empty <- vapply(stats_vec_list, is.null, logical(1))
  if (any(are_empty)) {
    stats_vec_list[are_empty] <- list(integer(ncol(well_stats_mat)))
  }

  assign("delete_stats_vec_list", stats_vec_list, envir = globalenv())


  summary_counts_mat <- do.call(rbind, stats_vec_list)
  mode(summary_counts_mat) <- "integer"

  total_vec <- as.integer(table(wells_fac))
  summary_perc_mat <- apply(summary_counts_mat,
                            2,
                            function(x) x / total_vec * 100
                            )

  colnames(summary_counts_mat) <- paste0("Count_", colnames(summary_counts_mat))
  colnames(summary_perc_mat) <- paste0("Perc_", colnames(summary_perc_mat))


  ## Process alteration categories ##

  has_alterations <- "sg1_category" %in% names(reads_df)
  if (has_alterations) {
    alterations_mat <- AlterationCategoriesToIntegerMat(reads_df)
    alterations_vec_list <- tapply(seq_len(nrow(reads_df)),
                                   wells_fac,
                                   function(x) colSums(alterations_mat[x, , drop = FALSE])
                                   )
    are_empty <- vapply(alterations_vec_list, is.null, logical(1))
    if (any(are_empty)) {
      alterations_vec_list[are_empty] <- list(integer(ncol(alterations_mat)))
    }
    alteration_counts_mat <- do.call(rbind, alterations_vec_list)
  }


  ## Process contaminations / wrong barcodes

  if ("Plate_number" %in% names(reads_df)) {
    all_plates <- setdiff(reads_df[, "Plate_number"], NA)
    contam_list <- lapply(all_plates, function(x) {
      wells_sub_vec <- sg_sequences_df[["Well_number"]][sg_sequences_df[["Plate_number"]] %in% x]
      sub_df <- reads_df[reads_df[, "Plate_number"] %in% x, ]
      row.names(sub_df) <- NULL
      ContaminationsSummaryDf(sub_df,
                              wells_vec = wells_sub_vec,
                              close_well_range = close_well_range
                              )
    })
    contam_df <- do.call(rbind.data.frame,
                         c(contam_list,
                           stringsAsFactors = FALSE,
                           make.row.names = FALSE
                         ))
  } else {
    contam_df <- ContaminationsSummaryDf(reads_df,
                                         wells_vec = unique_IDs,
                                         close_well_range = close_well_range
                                         )
  }


  num_under_2kb <- tapply(reads_df[, "Clipped_read_length"],
                          wells_fac,
                          function(x) sum(x < 2000)
                          )

  results_df <- data.frame(
    "Combined_ID"            = unique_IDs,
    "Plate_number"           = NA,
    "Well_number"            = NA,
    "Count_total"            = total_vec,
    "Num_under_2kb"          = num_under_2kb,
    "Num_low_barcode_scores" = tapply(reads_df[["Passes_barcode_filters"]],
                                      wells_fac, function(x) sum(x == 0)
                                      ),
    "Num_low_quality_scores" = tapply(reads_df[["Passes_read_quality"]],
                                      wells_fac, function(x) sum(x == 0)
                                      ),
    "Num_low_bc_or_qual"     = tapply(reads_df[["Passes_filters"]],
                                      wells_fac, function(x) sum(x == 0)
                                      ),
    "Mean_read_quality"      = tapply(reads_df[["Mean_quality"]],
                                      wells_fac, mean
                                      ),
    summary_counts_mat,
    summary_perc_mat,
    contam_df,
    stringsAsFactors         = FALSE,
    row.names                = NULL
  )

  if (ID_column == "Well_number") {
    results_df[["Well_number"]] <- results_df[["Combined_ID"]]
    results_df <- results_df[, !(names(results_df) %in% c("Combined_ID", "Plate_number"))]
  } else {
    results_df[["Plate_number"]] <- as.integer(substr(results_df[["Combined_ID"]], 6, 7))
    results_df[["Well_number"]]  <- as.integer(substr(results_df[["Combined_ID"]], 13, 15))
  }

  if (has_alterations) {
    results_df <- data.frame(results_df,
                             alteration_counts_mat,
                             stringsAsFactors = FALSE,
                             row.names = NULL
                             )
  }
  return(results_df)
}



CreateCrossContamMat <- function(contamin_mat_list, well_numbers_vec, wells_vec = seq_len(384)) {

  assign("delete_contamin_mat_list", contamin_mat_list, envir = globalenv())
  assign("delete_well_numbers_vec", well_numbers_vec, envir = globalenv())
  assign("delete_wells_vec", wells_vec, envir = globalenv())

  well_numbers_fac <- factor(well_numbers_vec, levels = wells_vec)

  all_comb_mat <- t(combn(wells_vec, 2))
  colnames(all_comb_mat) <- paste0("well", 1:2, "_ID")

  counts_wells_mat_list <- lapply(1:4, function(sg_number) {
    results_list <- tapply(seq_along(well_numbers_fac),
                           well_numbers_fac,
                           function(x) rowSums(contamin_mat_list[[sg_number]][, x, drop = FALSE]),
                           simplify = FALSE
                           )
    are_null <- vapply(results_list, is.null, logical(1))
    results_list[are_null] <- list(rep(NA_integer_, nlevels(well_numbers_fac)))
    results_mat <- do.call(rbind, results_list)
    return(results_mat)
  })

  cross_contam_mat <- do.call(cbind, lapply(1:4, function(sg_number) {
    use_mat <- counts_wells_mat_list[[sg_number]]
    sub_mat <- t(mapply(function(x, y) {
      assign("delete_x", x, envir = globalenv())
      assign("delete_y", y, envir = globalenv())
      assign("delete_use_mat", use_mat, envir = globalenv())
      x_index <- wells_vec == x
      y_index <- wells_vec == y
      c(use_mat[x_index, y_index], use_mat[y_index, x_index])
    },
    all_comb_mat[, 1],
    all_comb_mat[, 2]
    ))
    colnames(sub_mat) <- c(paste0("w2_sg", sg_number, "_in_well1"),
                           paste0("w1_sg", sg_number, "_in_well2")
                           )
    return(sub_mat)
  }))

  well_1_columns <- paste0("w2_sg", 1:4, "_in_well1")
  well_2_columns <- paste0("w1_sg", 1:4, "_in_well2")
  cross_contam_mat <- cbind("w2_in_well1" = rowSums(cross_contam_mat[, well_1_columns], na.rm = TRUE),
                            "w1_in_well2" = rowSums(cross_contam_mat[, well_2_columns], na.rm = TRUE),
                            cross_contam_mat
                            )
  mode(cross_contam_mat) <- "integer"
  cross_contam_mat <- cbind(all_comb_mat, cross_contam_mat)

  totals_vec <- rowSums(cross_contam_mat[, c("w2_in_well1", "w1_in_well2")])
  use_order <- order(totals_vec, decreasing = TRUE)
  are_zero <- totals_vec == 0

  cross_contam_mat <- cross_contam_mat[use_order[!(are_zero[use_order])], ]
  return(cross_contam_mat)
}




CheckThatFactorIsInOrder <- function(my_factor) {
  stopifnot(identical(length(unique(my_factor)),
                      length(rle(as.integer(my_factor))[["lengths"]])
                      )
            )
}



GetFeaturesData <- function(extracted_df, use_ZMWs = NULL, get_quality = FALSE) {

  CheckThatFactorIsInOrder(extracted_df[["ZMW"]])
  stopifnot(length(unique(table(extracted_df[["ZMW"]]))) == 1)
  stopifnot(length(unique(table(extracted_df[["Feature"]]))) == 1)

  if (!(is.null(use_ZMWs))) {
    unique_zmws <- unique(extracted_df[["ZMW"]])
    zmw_matches <- match(use_ZMWs, unique_zmws)
    stopifnot(!(anyNA(zmw_matches)))
  }

  use_features <- c(paste0("sg", 1:4),
                    paste0("sg", 1:4, "_cr", 1:4)
                    )

  if (get_quality) {
    extract_column <- "Mean_quality"
    column_suffix <- "quality"
  } else {
    use_features <- c(use_features, "TpR_DHFR")
    extract_column <- "Category"
    column_suffix <- "category"
  }

  categories_vec_list <- sapply(use_features, function(x) {
    are_this_feature <- extracted_df[["Feature"]] == x
    categories_sub_vec <- extracted_df[[extract_column]][are_this_feature]
    if (!(is.null(use_ZMWs))) {
      categories_sub_vec <- categories_sub_vec[zmw_matches]
    }
    return(categories_sub_vec)
  }, simplify = FALSE)

  results_mat <- do.call(cbind, categories_vec_list)
  colnames(results_mat) <- paste0(colnames(results_mat), "_", column_suffix)
  return(results_mat)
}





AlterationCategoriesToIntegerMat <- function(input_df) {
  suffix_regex <- "_category$"
  categories_columns <- grep(suffix_regex, names(input_df), value = TRUE)
  all_features <- sub(suffix_regex, "", categories_columns)
  features_mat <- as.matrix(input_df[, categories_columns])
  features_mat[features_mat == "Flanking insertion"] <- "Correct"

  basic_categories <- c("Correct", "Mutation", "Deletion")
  contam_categories <- c(basic_categories, "Contamination")
  contamination_categories <- c(paste0("sg", 1:4),
                                paste0("sg", 1:4, "_cr", 1:4)
                                )
  results_mat_list <- lapply(seq_along(all_features), function(x) {
    if (all_features[[x]] %in% contamination_categories) {
      use_categories <- contam_categories
    } else {
      use_categories <- basic_categories
    }
    categories_list <- lapply(use_categories, function(y) {
      features_mat[, categories_columns[[x]]] == y
    })
    categories_mat <- do.call(cbind, categories_list)
    colnames(categories_mat) <- paste0(use_categories, "_", all_features[[x]])
    return(categories_mat)
  })
  results_mat <- do.call(cbind, results_mat_list)
  mode(results_mat) <- "integer"
  return(results_mat)
}




AnalyzeWells <- function(ccs_df,
                         sg_df,
                         barcodes_df,
                         extracted_df,
                         bc_min_comb_score   = 60,
                         bc_min_score_lead   = 30,
                         bc_length_cutoff    = 8L,
                         min_mean_quality    = 85 / 93 * 100,
                         set_seed            = TRUE,
                         use_guides_ref_list = guides_ref_list
                         ) {

  stopifnot("manhattan_dist_list" %in% ls(envir = globalenv()))

  ccs_df[["Mean_quality"]] <- GetMeanQuality(substr(ccs_df[["Quality"]],
                                                    ccs_df[["Clip_start"]],
                                                    ccs_df[["Clip_end"]]
                                                    )
                                             )

  ccs_df[["Pass_CCS5"]] <- as.integer(ccs_df[["Pass_CCS5"]])


  ## Add data on contaminations / wrong barcodes

  contamin_mat_list <- lapply(sg_df[, paste0("sg_cr_", 1:4)],
                              function(x) GetContaminationMat(x, ccs_df, sg_df)
                              )

  contamin_mat <- Reduce(`+`, contamin_mat_list)
  assign("delete_contamin_mat_list", contamin_mat_list, envir = globalenv())
  assign("delete_contamin_mat", contamin_mat, envir = globalenv())

  have_contamin_mat <- contamin_mat >= 1
  have_well <- !(is.na(ccs_df[["Well_number"]])) &
               ccs_df[["Passed_filters"]]
  have_valid_well <- have_well & ccs_df[["Well_number"]] %in% sg_df[, "Well_number"]
  well_numbers_vec <- ifelse(have_well, ccs_df[["Well_number"]], NA_integer_)
  assign("delete_well_numbers_vec", well_numbers_vec, envir = globalenv())
  distance_mat <- DistanceForContam(contamin_mat,
                                    well_numbers_vec,
                                    wells_vec = sg_df[, "Well_number"]
                                    )

  ccs_df[["Num_contaminating_guides"]] <- NA_integer_
  ccs_df[["Num_contaminating_genes"]]  <- NA_integer_
  ccs_df[["Contaminating_well"]]       <- NA_integer_
  ccs_df[["Min_distance"]]             <- NA_integer_
  ccs_df[["Random_distance"]]          <- NA_integer_
  ccs_df[["Mean_distance"]]            <- NA_real_
  ccs_df[["Distances"]]                <- NA_character_

  ccs_df[["Num_contaminating_guides"]][have_well] <- colSums(contamin_mat)
  ccs_df[["Num_contaminating_genes"]][have_well] <- colSums(have_contamin_mat)

  min_distances_vec <- colMins(distance_mat, na.rm = TRUE)
  min_distances_vec[is.infinite(min_distances_vec)] <- NA_integer_
  ccs_df[["Min_distance"]][have_well] <- min_distances_vec

  mean_distances_vec <- colMeans(distance_mat, na.rm = TRUE)
  ccs_df[["Mean_distance"]][have_well] <- mean_distances_vec

  have_contam <- !(is.na(min_distances_vec))

  if (set_seed) {
    set.seed(1)
  }
  random_distances_vec <- vapply(which(have_contam),
                                 function(x) {
                                   dist_vec <- distance_mat[, x]
                                   dist_vec <- dist_vec[!(is.na(dist_vec))]
                                   if (length(dist_vec) == 1) {
                                     return(dist_vec)
                                   } else {
                                     return(sample(dist_vec, 1))
                                   }},
                                 integer(1)
                                 )
  ccs_df[["Random_distance"]][have_well][have_contam] <- random_distances_vec


  all_distances_vec <- vapply(which(have_contam),
                              function(x) {
                                dist_vec <- distance_mat[, x]
                                dist_vec <- dist_vec[!(is.na(dist_vec))]
                                return(paste0(dist_vec, collapse = ", "))
                              }, "")
  ccs_df[["Distances"]][have_well][have_contam] <- all_distances_vec


  contams_vec <- vapply(which(have_contam),
                        function(x) paste0(sg_df[["Well_number"]][which(have_contamin_mat[, x])], collapse = ", "),
                        ""
                        )
  ccs_df[["Contaminating_well"]][have_well][have_contam] <- contams_vec


  ## Check for the presence of sequences

  well_stats_mat_list <- lapply(seq_len(nrow(sg_df)),
                                function(x) StatsForWell(x, ccs_df, sg_df)
                                )


  ## Create a data frame of individual ZMWs

  ZMW_mat <- do.call(rbind, well_stats_mat_list)

  sg_features <- c(paste0("sg", 1:4), paste0("sg", 1:4, "_cr", 1:4))

  only_wells <- "Barcode_combined_score" %in% names(ccs_df)
  if (only_wells) {
    barcode_columns <- c("Barcode_combined_score", "Barcode_score_lead")
    ID_columns <- "Well_number"
  } else {
    barcode_columns <- c("Plate_barcode_combined_score", "Plate_barcode_score_lead",
                         "Well_barcode_combined_score", "Well_barcode_score_lead"
                         )
    ID_columns <- c("Combined_ID", "Plate_number", "Well_number")
  }

  all_columns <- c(
    "ZMW", ID_columns,
    "Clipped_read_length",

    "Read_quality", "Num_full_passes", "Pass_CCS5",

    "Correct_barcodes", "Correct_row", "Correct_column",
    "Correct_flanks", "Correct_row_flank", "Correct_column_flank",

    "Starts_with_row_barcode", "Ends_with_column_barcode",

    "Row_bc_length", "Column_bc_length",

    "Row_mean_quality",  "Column_mean_quality", "Row_barcode", "Column_barcode",
    "Row_quality", "Column_quality",

    paste0(c("TpR_DHFR", sg_features), "_category"),

    paste0(sg_features, "_quality"),

    barcode_columns,

    "Mean_quality",
    "Passes_filters", "Passes_barcode_filters", "Passes_read_quality",
    "Passes_sg_quality",
    "Num_contaminating_guides", "Num_contaminating_genes", "Contaminating_well",
    "Min_distance", "Random_distance", "Mean_distance", "Distances",
    "sg1_cr1", "sg2_cr2", "sg3_cr3", "sg4_cr4", "at_least_1", "at_least_2",
    "at_least_3", "all_4", "all_4_promoters", "whole_plasmid",
    "Orientation_fwd"
  )

  wells_ccs_df <- do.call(rbind.data.frame,
                          c(split(ccs_df[, names(ccs_df) %in% all_columns],
                                  well_numbers_vec
                                  ),
                          list(stringsAsFactors = FALSE,
                               make.row.names = FALSE
                               )
                          ))

  are_real_wells <- wells_ccs_df[["Well_number"]] %in% sg_df[["Well_number"]]

  individual_ZMWs_df <- data.frame(
    wells_ccs_df[are_real_wells, ],
    ZMW_mat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  contamin_mat_list <- lapply(contamin_mat_list, function(x) x[, are_real_wells])

  if (is.null(barcodes_df) && is.null(extracted_df)) {
    results_list <- list(
      "individual_reads_df" = individual_ZMWs_df[, intersect(all_columns, names(individual_ZMWs_df))],
      "contamin_mat_list"   = contamin_mat_list
    )
    return(results_list)
  }

  matches_vec <- match(individual_ZMWs_df[["ZMW"]], barcodes_df[["ZMW"]])
  stopifnot(!(anyNA(matches_vec)))
  barcodes_df <- barcodes_df[matches_vec, !(names(barcodes_df) %in% names(individual_ZMWs_df))]
  individual_ZMWs_df <- data.frame(individual_ZMWs_df,
                                   barcodes_df,
                                   stringsAsFactors = FALSE,
                                   row.names = NULL
                                   )

  have_correct_row <- individual_ZMWs_df[, "Starts_with_row_barcode"] |
                      ((individual_ZMWs_df[, "Row_bc_length"] >= bc_length_cutoff) &
                        individual_ZMWs_df[, "Correct_row_flank"])

  have_correct_column <- individual_ZMWs_df[, "Ends_with_column_barcode"] |
                         ((individual_ZMWs_df[, "Column_bc_length"] >= bc_length_cutoff) &
                           individual_ZMWs_df[, "Correct_column_flank"])


  if ("Barcode_combined_score" %in% names(individual_ZMWs_df)) {
    pass_bc <- have_correct_row & have_correct_column &
              (individual_ZMWs_df[, "Barcode_combined_score"] >= bc_min_comb_score) &
              (individual_ZMWs_df[, "Barcode_score_lead"] >= bc_min_score_lead)
  } else {
    pass_bc <- have_correct_row & have_correct_column &
              (individual_ZMWs_df[, "Plate_barcode_combined_score"] >= bc_min_comb_score) &
              (individual_ZMWs_df[, "Well_barcode_combined_score"] >= bc_min_comb_score) &
              (individual_ZMWs_df[, "Plate_barcode_score_lead"] >= bc_min_score_lead) &
              (individual_ZMWs_df[, "Well_barcode_score_lead"] >= bc_min_score_lead)
  }

  pass_rq <- (individual_ZMWs_df[["Mean_quality"]]) >= min_mean_quality
  pass_filters <- pass_bc & pass_rq

  individual_ZMWs_df[["Passes_filters"]]         <- as.integer(pass_filters)
  individual_ZMWs_df[["Passes_barcode_filters"]] <- as.integer(pass_bc)
  individual_ZMWs_df[["Passes_read_quality"]]    <- as.integer(pass_rq)

  if (is.null(extracted_df)) {
    results_list <- list(
      "individual_reads_df" = individual_ZMWs_df[, intersect(all_columns, names(individual_ZMWs_df))],
      "contamin_mat_list"   = contamin_mat_list
    )
    return(results_list)
  }

  alterations_mat <- GetFeaturesData(extracted_df, individual_ZMWs_df[["ZMW"]])
  subqualities_mat <- GetFeaturesData(extracted_df, individual_ZMWs_df[["ZMW"]], get_quality = TRUE)
  subqualities_mat <- subqualities_mat / 93 * 100

  are_poor_quality_mat <- subqualities_mat < min_mean_quality
  are_poor_quality_mat[is.na(are_poor_quality_mat)] <- FALSE

  individual_ZMWs_df[["Passes_sg_quality"]] <- as.integer(rowSums(are_poor_quality_mat) == 0)

  individual_ZMWs_df <- data.frame(
    individual_ZMWs_df,
    alterations_mat,
    subqualities_mat,
    stringsAsFactors = FALSE
  )

  individual_ZMWs_df <- individual_ZMWs_df[, all_columns]

  results_list <- list(
    "individual_reads_df" = individual_ZMWs_df,
    "contamin_mat_list"   = contamin_mat_list
  )
  return(results_list)
}




SummarizeWells <- function(analysis_list,
                           use_zmws = NULL,
                           close_well_range = 3L,
                           ID_column = "Well_number",
                           unique_IDs = seq_len(384)
                           ) {

  reads_df <- analysis_list[["individual_reads_df"]]

  has_contam <- "contamin_mat_list" %in% names(analysis_list)
  if (has_contam) {
    contamin_mat_list <- analysis_list[["contamin_mat_list"]]
  }

  if (!(is.null(use_zmws))) {
    zmw_matches <- match(use_zmws, reads_df[["ZMW"]])
    stopifnot(!(anyNA(zmw_matches)))
    reads_df <- reads_df[zmw_matches, ]
    row.names(reads_df) <- NULL
    if (has_contam) {
      contamin_mat_list <- lapply(contamin_mat_list, function(x) x[, zmw_matches])
    }
  }

  unfiltered_summary_df <- CreateSummaryDf(reads_df,
                                           filter_reads = FALSE,
                                           close_well_range = close_well_range,
                                           ID_column = ID_column,
                                           unique_IDs = unique_IDs
                                           )
  filtered_summary_df   <- CreateSummaryDf(reads_df,
                                           filter_reads = TRUE,
                                           close_well_range = close_well_range,
                                           ID_column = ID_column,
                                           unique_IDs = unique_IDs
                                           )

  results_list <- list(
    "original_summary_df" = unfiltered_summary_df,
    "filtered_summary_df" = filtered_summary_df
  )
  if ("Passes_sg_quality" %in% names(reads_df)) {
    filtered_gRNAs_df <- CreateSummaryDf(reads_df,
                                         filter_reads = TRUE,
                                         filter_subsequences = TRUE,
                                         close_well_range = close_well_range,
                                         ID_column = ID_column,
                                         unique_IDs = unique_IDs
                                         )
    results_list <- c(results_list, list("filtered_gRNAs_df" = filtered_gRNAs_df))
  }
  results_list <- c(results_list, list("individual_reads_df" = reads_df))

  if (has_contam) {
    stopifnot(ID_column == "Well_number")
    stopifnot(identical(unique(vapply(contamin_mat_list, ncol, integer(1))),
                        nrow(reads_df)
                        ))
    cross_contam_mat <- CreateCrossContamMat(contamin_mat_list,
                                             reads_df[["Well_number"]],
                                             wells_vec = unique_IDs
                                             )
    results_list <- c(results_list, list("contaminations_mat"  = cross_contam_mat))
  }
  return(results_list)
}




DistanceForContam <- function(contamination_mat, long_wells_vec, wells_vec = seq_len(384)) {

  stopifnot("manhattan_dist_list" %in% ls(envir = globalenv()))
  have_wells <- !(is.na(long_wells_vec))

  stopifnot(ncol(contamination_mat) == sum(have_wells))

  distance_mat <- matrix(nrow = nrow(contamination_mat),
                         ncol = ncol(contamination_mat)
                         )
  for (i in wells_vec) {
    are_this_well <- long_wells_vec[have_wells] %in% i
    are_contam <- contamination_mat[, are_this_well, drop = FALSE] >= 1
    sub_mat <- distance_mat[, are_this_well, drop = FALSE]
    contam_rows <- which(rowSums(are_contam) >= 1)
    for (row_i in contam_rows) {
      sub_mat[row_i, are_contam[row_i, ]] <- manhattan_dist_list[[i]][[row_i]]
    }
    distance_mat[, are_this_well] <- sub_mat
  }
  return(distance_mat)
}





# Functions for processing multiple plates --------------------------------

AnalyzePlates <- function(use_ccs_df,
                          use_sg_df,
                          use_barcodes_df,
                          use_extracted_df,
                          set_seed = TRUE
                          ) {

  are_eligible <- (use_ccs_df[["Well_exists"]] %in% TRUE) &
                  (use_ccs_df[["Read_quality"]] > 0)

  use_ccs_df[["Passed_filters"]] <- use_ccs_df[["Plate_passed_filters"]] &
                                    (use_ccs_df[["Well_passed_filters"]] %in% TRUE)

  use_ccs_df[["Pass_CCS5"]] <- NA_integer_
  names(use_ccs_df)[names(use_ccs_df) == "Well_clip_start"]          <- "Clip_start"
  names(use_ccs_df)[names(use_ccs_df) == "Well_clip_end"]            <- "Clip_end"
  names(use_ccs_df)[names(use_ccs_df) == "Well_clipped_read_length"] <- "Clipped_read_length"

  all_plates <- sort(setdiff(use_ccs_df[["Plate_number"]], NA))

  sub_list_list <- lapply(all_plates, function(x) {
    message(paste0("Analyzing plate #", x, "...\n"))
    ccs_sub_df <- use_ccs_df[are_eligible & (use_ccs_df[["Plate_number"]] %in% x), ]
    are_this_plate <- sg_sequences_df[["Plate_number"]] == x
    sg_sub_df <- use_sg_df[are_this_plate, ]
    row.names(sg_sub_df) <- NULL
    sub_list <- AnalyzeWells(ccs_sub_df,
                             sg_sub_df,
                             use_barcodes_df,
                             use_extracted_df,
                             set_seed = set_seed
                             )
    message("\n\n")
    return(sub_list)
  })

  combined_df <- do.call(rbind.data.frame,
                         c(lapply(sub_list_list, function(x) x[["individual_reads_df"]]),
                           stringsAsFactors = FALSE, make.row.names = FALSE
                           ))

  contamin_mat_list_list <- lapply(sub_list_list, function(x) x[["contamin_mat_list"]])
  names(contamin_mat_list_list) <- paste0("Plate", all_plates)

  results_list <- list(
    "individual_reads_df" = combined_df,
    "contamin_mat_list_list" = contamin_mat_list_list
  )
  return(results_list)
}




# Functions for exporting data --------------------------------------------

ExportTable <- function(export_df,
                        file_name_only,
                        file_directory = tables_output_directory
                        ) {

  have_NA_columns <- vapply(export_df, anyNA, logical(1))
  for (i in which(have_NA_columns)) {
    export_df[[i]] <- ifelse(is.na(export_df[[i]]),
                             "",
                             export_df[[i]]
                             )
  }
  print(file_directory)
  write.table(export_df,
              file      = file.path(file_directory, paste0(file_name_only, ".tsv")),
              sep       = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE
              )
  return(invisible(NULL))
}



ExportIndivTable <- function(indiv_reads_df, ...) {
  round_columns <- grep( "_quality$", names(indiv_reads_df), value = TRUE)
  round_columns <- setdiff(round_columns,
                           c("Read_quality", "Passes_read_quality",
                             "Row_quality", "Column_quality"
                             )
                           )
  for (quality_column in round_columns) {
    indiv_reads_df[[quality_column]] <- round(indiv_reads_df[[quality_column]], digits = 1)
  }
  exclude_columns <- c("Min_distance", "Random_distance", "Mean_distance",
                       "Passes_barcode_filters", "Passes_read_quality"
                       )
  if (all(is.na(indiv_reads_df[["Pass_CCS5"]]))) {
    exclude_columns <- c(exclude_columns, "Pass_CCS5")
  }
  indiv_reads_df <- indiv_reads_df[, !(names(indiv_reads_df) %in% exclude_columns)]
  ExportTable(indiv_reads_df, ...)
}



ExportSummaryTable <- function(summary_df, ...) {
  if ("Distance_p_value" %in% names(summary_df)) {
    p_val_vec <- signif(summary_df[["Distance_p_value"]], 1)
    p_val_vec <- vapply(signif(p_val_vec, 1), format, scientific = 9L, "")
    summary_df[["Distance_p_value"]] <- p_val_vec
  }
  are_integer <- vapply(summary_df, is.integer, logical(1))
  are_numeric <- vapply(summary_df, is.numeric, logical(1))
  for (i in which(!(are_integer) & are_numeric)) {
    summary_df[[i]] <- round(summary_df[[i]], digits = 1)
  }
  ExportTable(summary_df, ...)
}




