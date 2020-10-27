### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("ShortRead") # For processing the quality scores



# Define functions --------------------------------------------------------

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



AddBarCodeCombos <- function(lima_report_df) {
  combo_IDs_vec <- paste0(lima_report_df[["IdxLowestNamed"]], "--",
                          lima_report_df[["IdxHighestNamed"]]
                          )
  lima_report_df[["Barcode_combo"]] <- combo_IDs_vec
  lima_report_df[["Well_number"]] <- barcodes_to_wells_map[combo_IDs_vec]
  return(lima_report_df)
}



ContainSequences <- function(query_seq, target_seq) {

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

  contain_fwd_mat <- vapply(seq_len(num_queries), function(x) {
    vcountPattern(query_seq[[x]], target_seq)
  }, integer(num_targets))
  contain_rev_mat <- vapply(seq_len(num_queries), function(x) {
    vcountPattern(rev_search_seq[[x]], target_seq)
  }, integer(num_targets))

  are_fwd <- rowSums(contain_fwd_mat) > 0
  are_rev <- rowSums(contain_rev_mat) > 0


  num_duplicated_fwd <- sum(contain_fwd_mat > 1)
  num_duplicated_rev <- sum(contain_rev_mat > 1)

  sent_message <- FALSE

  if (num_duplicated_fwd != 0) {
    message(paste0("The same sequence was found more than once for ",
                   num_duplicated_fwd, " read",
                   if (num_duplicated_fwd == 1) "" else "",
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
  contain_num_mat <- vapply(seq_len(num_queries),
                            function(x) num_correct_vec >= x,
                            integer(num_targets)
                            )
  colnames(contain_num_mat) <- paste0("at_least_", seq_len(num_queries))

  results_mat <- cbind(contain_mat,
                       contain_num_mat,
                       "orientation_fwd_sg" = are_fwd,
                       "inconsistent_orientation" = have_duplication | have_inversion
                       )
  return(results_mat)
}




StatsForWell <- function(well_number, reads_list, bam_wells_vec, wells_vec = seq_len(384)) {

  well_index <- which(wells_vec == well_number)

  are_this_well <- bam_wells_vec %in% well_number
  num_reads <- sum(are_this_well)
  well_sequences <- reads_list[["seq"]][are_this_well, ]

  sg_sequences <- vapply(1:4, function(x) guides_ref_list[[x]][[well_index]], "")
  sg_with_promoter_seq <- vapply(1:4, function(x) guides_with_promoters_list[[x]][[well_index]], "")

  contain_sg_cr_mat   <- ContainSequences(sg_sequences, well_sequences)
  contain_sg_prom_mat <- ContainSequences(sg_with_promoter_seq, well_sequences)
  contain_plasmid_mat <- ContainSequences(plasmids_vec[[well_index]], well_sequences)

  colnames(contain_sg_cr_mat)[1:4] <- paste0("sg", 1:4, "_cr", 1:4)
  colnames(contain_sg_cr_mat)[colnames(contain_sg_cr_mat) == "at_least_4"] <- "all_4"
  results_mat <- cbind(
    contain_sg_cr_mat[, 1:8],
    "all_4_promoters" = contain_sg_prom_mat[, "at_least_4"],
    "whole_plasmid"   = contain_plasmid_mat[, 1],
    contain_sg_cr_mat[, "orientation_fwd_sg", drop = FALSE]
  )
  stopifnot(!(any(results_mat[, "all_4_promoters"] < results_mat[, "whole_plasmid"])))

  return(results_mat)
}






GetZMWs <- function(char_vec) {
   as.integer(sapply(strsplit(char_vec, "/", fixed = TRUE), "[[", 2))
}




MakeDfForReads <- function(reads_list,
                           report_df,
                           counts_df = NULL,
                           counts_384_df = NULL
                           ) {

  ## Create barcode IDs for the CCS report files

  report_df <- AddBarCodeCombos(report_df)


  ## Check for equivalence

  if (!(is.null(counts_df))) {
    counts_IDs_vec <- paste0(counts_df[["IdxFirstNamed"]], "--",
                             counts_df[["IdxCombinedNamed"]]
                             )

    pass_filters <- report_df[["PassedFilters"]] == 1

    IDs_to_counts <- table(factor(report_df[["Barcode_combo"]][pass_filters],
                                  levels = counts_IDs_vec
                                  )
                           )

    stopifnot(identical(as.integer(IDs_to_counts), counts_df[["Counts"]]))
  }




  ## Assign the sequences to wells

  profile_ZMWs <- GetZMWs(report_df[["ZMW"]])
  bam_ZMWs <- GetZMWs(reads_list[["qname"]])

  matches_vec <- match(bam_ZMWs, profile_ZMWs)

  results_df <- data.frame(
    "ZMW"               = bam_ZMWs,
    "Well_number"       = report_df[["Well_number"]][matches_vec],
    "Length"            = width(reads_list[["seq"]]),
    "BC_combined_score" = report_df[["ScoreCombined"]][matches_vec],
    "BC_score_lead"     = report_df[["ScoreLead"]][matches_vec],
    stringsAsFactors    = FALSE
  )


  ## Check for equivalence

  if (!(is.null(counts_384_df))) {
    stopifnot(identical(as.integer(table(results_df[["Well_number"]])),
                        counts_384_df[["Counts"]]
                        )
              )
  }
  return(results_df)
}





GetMode <- function(vec) {
   uniqv <- unique(vec)
   uniqv[which.max(tabulate(match(vec, uniqv)))]
}




GetContaminationMat <- function(query_seq, reads_list, bam_wells_vec, wells_vec = seq_len(384)) {

  stopifnot("guides_ref_list" %in% ls(envir = globalenv()))

  are_usual_length <- nchar(query_seq) == GetMode(nchar(query_seq))

  PDict_object <- PDict(query_seq[are_usual_length], max.mismatch = 0)

  have_well <- !(is.na(bam_wells_vec))

  counts_mat <- vcountPDict(PDict_object, reads_list[["seq"]][have_well, ])

  row_list <- vector(mode = "list", length = length(wells_vec))
  row_list[are_usual_length] <- lapply(seq_len(nrow(counts_mat)),
                                       function(x) counts_mat[x, ]
                                       )

  unusual_seq <- DNAStringSet(query_seq[!(are_usual_length)])

  unusual_list <- lapply(seq_len(sum(!(are_usual_length))),
                         function(x) vcountPattern(unusual_seq[[x]],
                                                   reads_list[["seq"]][have_well, ]
                                                   )
                         )

  row_list[!(are_usual_length)] <- unusual_list

  all_counts_mat <- do.call(rbind, row_list)
  stopifnot(all(unique(as.vector(all_counts_mat))) %in% 0:1)
  mode(all_counts_mat) <- "integer"

  ## If the sequence of two wells is exactly the same,
  ## perhaps it shouldn't count as a contamination...
  same_seq_list <- lapply(seq_along(wells_vec),
                          function(x) which(query_seq %in% vapply(guides_ref_list, function(y) y[[x]], ""))
                          )
  are_this_well_list <- lapply(same_seq_list,
                               function(x) bam_wells_vec[have_well] %in% x
                               )

  contamination_mat <- all_counts_mat
  for (i in seq_along(wells_vec)) {
    contamination_mat[i, are_this_well_list[[i]]] <- 0L
  }
  return(contamination_mat)
}





CreateSummaryDf <- function(reads_df,
                            filter_reads = FALSE,
                            close_well_range = 3L,
                            wells_vec = seq_len(384)
                            ) {

  assign("delete_reads_df", reads_df, envir = globalenv())
  assign("delete_wells_vec", wells_vec, envir = globalenv())

  if (filter_reads) {
    reads_df <- reads_df[reads_df[["Passes_filters"]] == 1, ]
  }

  wells_fac <- factor(reads_df[["Well_number"]])


  ## Process 100% correct read counts and percentages

  binary_columns <- c(paste0("sg", 1:4, "_cr", 1:4),
                      paste0("at_least_", 1:3),
                      "all_4", "all_4_promoters", "whole_plasmid"
                      )

  well_stats_mat <- as.matrix(reads_df[, binary_columns])

  summary_counts_mat <- do.call(rbind,
                                tapply(seq_len(nrow(reads_df)),
                                       wells_fac,
                                       function(x) colSums(well_stats_mat[x, ])
                                       )
                                )
  mode(summary_counts_mat) <- "integer"

  total_vec <- as.integer(table(wells_fac))
  summary_perc_mat <- apply(summary_counts_mat,
                            2,
                            function(x) x / total_vec * 100
                            )

  colnames(summary_counts_mat) <- paste0("Count_", colnames(summary_counts_mat))
  colnames(summary_perc_mat) <- paste0("Perc_", colnames(summary_perc_mat))


  ## Process alteration categories ##

  alterations_mat <- AlterationCategoriesToIntegerMat(reads_df)
  alteration_counts_mat <- do.call(rbind,
                                   tapply(seq_len(nrow(reads_df)),
                                          wells_fac,
                                          function(x) colSums(alterations_mat[x, ])
                                          )
                                   )

  ## Process contaminations / wrong barcodes

  are_contaminated <- reads_df[["Contam_guides"]] >= 1
  are_close_rand <- (reads_df[["Random_distance"]] <= close_well_range) %in% TRUE

  num_other_genes_vec <- tapply(reads_df[["Contam_well"]],
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

  per_well_distances_vec <- tapply(reads_df[["Mean_distance"]],
                                   wells_fac,
                                   mean,
                                   na.rm = TRUE
                                   )

  expected_distances_vec <- vapply(manhattan_dist_list,
                                   function(x) mean(x[x != 0]),
                                   numeric(1)
                                   )

  distances_list <- split(reads_df[["Mean_distance"]], wells_fac)

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

  num_under_2kb <- tapply(reads_df[["Length"]],
                          wells_fac,
                          function(x) sum(x < 2000)
                          )

  results_df <- data.frame(
    "Well_number"            = wells_vec,
    "Count_total"            = total_vec,
    "Num_under_2kb"          = num_under_2kb,
    "Num_low_barcode_scores" = tapply(reads_df[["Passes_barcode_filters"]],
                                      wells_fac, function(x) sum(x == 0)
                                      ),
    "Num_low_quality_scores" = tapply(reads_df[["Passes_read_filters"]],
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
    "Num_contam_genes"       = num_other_genes_vec,
    "Num_contam_wells"       = num_contaminants,
    "Num_from_close_wells"   = num_close_contaminants_rand,
    "Expected_close_wells"   = num_expected_close_wells,
    "Mean_distance"          = per_well_distances_vec,
    "Expected_distance"      = expected_distances_vec[wells_vec],
    "Distance_p_value"       = p_val_vec,
    alteration_counts_mat,
    stringsAsFactors         = FALSE,
    row.names                = NULL
  )
  return(results_df)
}



CreateCrossContamMat <- function(contamin_mat_list, well_numbers_vec, wells_vec = seq_len(384)) {

  assign("delete_contamin_mat_list", contamin_mat_list, envir = globalenv())
  assign("delete_well_numbers_vec", well_numbers_vec, envir = globalenv())
  assign("delete_wells_vec", wells_vec, envir = globalenv())

  all_comb_mat <- t(combn(wells_vec, 2))
  colnames(all_comb_mat) <- paste0("well", 1:2, "_ID")

  counts_wells_mat_list <- lapply(1:4, function(sg_number) {
    do.call(rbind,
            tapply(seq_along(well_numbers_vec),
                   well_numbers_vec,
                   function(x) rowSums(contamin_mat_list[[sg_number]][, x])
                   )
            )
  })

  cross_contam_mat <- do.call(cbind, lapply(1:4, function(sg_number) {
    use_mat <- counts_wells_mat_list[[sg_number]]
    sub_mat <- t(mapply(function(x, y) {
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
  cross_contam_mat <- cbind("w2_in_well1" = rowSums(cross_contam_mat[, well_1_columns]),
                            "w1_in_well2" = rowSums(cross_contam_mat[, well_2_columns]),
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



AddAlterationCategories <- function(extracted_df, use_ZMWs = NULL) {

  CheckThatFactorIsInOrder(extracted_df[["ZMW"]])
  stopifnot(length(unique(table(extracted_df[["ZMW"]]))) == 1)
  stopifnot(length(unique(table(extracted_df[["Feature"]]))) == 1)

  if (!(is.null(use_ZMWs))) {
    unique_zmws <- unique(extracted_df[["ZMW"]])
    zmw_matches <- match(use_ZMWs, unique_zmws)
    stopifnot(!(anyNA(zmw_matches)))
  }

  use_features <- c(paste0("sg", 1:4),
                    paste0("sg", 1:4, "_cr", 1:4),
                    "TpR_DHFR"
                    )

  categories_vec_list <- sapply(use_features, function(x) {
    are_this_feature <- extracted_df[["Feature"]] == x
    categories_sub_vec <- extracted_df[["Category"]][are_this_feature]
    if (!(is.null(use_ZMWs))) {
      categories_sub_vec <- categories_sub_vec[zmw_matches]
    }
    return(categories_sub_vec)
  }, simplify = FALSE)

  results_mat <- do.call(cbind, categories_vec_list)
  colnames(results_mat) <- paste0(colnames(results_mat), "_category")
  return(results_mat)
}



AlterationCategoriesToIntegerMat <- function(input_df) {
  suffix_regex <- "_category$"
  categories_columns <- grep(suffix_regex, colnames(input_df), value = TRUE)
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




AnalyzeWells <- function(reads_list,
                         report_df,
                         barcodes_df,
                         extracted_df,
                         counts_df         = NULL,
                         counts_384_df     = NULL,
                         bc_min_comb_score = 60,
                         bc_min_score_lead = 30,
                         bc_length_cutoff  = 8L,
                         min_mean_quality  = 85,
                         wells_vec         = seq_len(384)
                         ) {

  stopifnot("manhattan_dist_list" %in% ls(envir = globalenv()))

  reads_report_df <- MakeDfForReads(reads_list, report_df)
  reads_report_df[["Mean_quality"]] <- GetMeanQuality(reads_list[["qual"]])


  ## Add data on contaminations / wrong barcodes

  contamin_mat_list <- lapply(guides_ref_list,
                              function(x) GetContaminationMat(x,
                                                              reads_list,
                                                              reads_report_df[["Well_number"]],
                                                              wells_vec = wells_vec
                                                              )
                              )


  contamin_mat <- Reduce(`+`, contamin_mat_list)
  have_contamin_mat <- contamin_mat >= 1
  distance_mat <- DistanceForContam(contamin_mat,
                                    reads_report_df[["Well_number"]],
                                    wells_vec = wells_vec
                                    )

  reads_report_df[["Contam_guides"]]   <- NA_integer_
  reads_report_df[["Contam_genes"]]    <- NA_integer_
  reads_report_df[["Contam_well"]]     <- NA_integer_
  reads_report_df[["Min_distance"]]    <- NA_integer_
  reads_report_df[["Random_distance"]] <- NA_integer_
  reads_report_df[["Mean_distance"]]   <- NA_real_
  reads_report_df[["Distances"]]       <- NA_character_

  have_well <- !(is.na(reads_report_df[["Well_number"]]))

  reads_report_df[["Contam_guides"]][have_well] <- colSums(contamin_mat)
  reads_report_df[["Contam_genes"]][have_well] <- colSums(have_contamin_mat)

  min_distances_vec <- colMins(distance_mat, na.rm = TRUE)
  min_distances_vec[is.infinite(min_distances_vec)] <- NA_integer_
  reads_report_df[["Min_distance"]][have_well] <- min_distances_vec

  mean_distances_vec <- colMeans(distance_mat, na.rm = TRUE)
  reads_report_df[["Mean_distance"]][have_well] <- mean_distances_vec

  have_contam <- !(is.na(min_distances_vec))

  set.seed(1)
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
  reads_report_df[["Random_distance"]][have_well][have_contam] <- random_distances_vec


  all_distances_vec <- vapply(which(have_contam),
                              function(x) {
                                dist_vec <- distance_mat[, x]
                                dist_vec <- dist_vec[!(is.na(dist_vec))]
                                return(paste0(dist_vec, collapse = ", "))
                              }, "")
  reads_report_df[["Distances"]][have_well][have_contam] <- all_distances_vec


  contams_vec <- vapply(which(have_contam),
                        function(x) paste0(which(have_contamin_mat[, x]), collapse = ", "),
                        ""
                        )
  reads_report_df[["Contam_well"]][have_well][have_contam] <- contams_vec



  ## Check for the presence of sequences

  well_stats_mat_list <- lapply(wells_vec,
                                function(x) StatsForWell(x,
                                                         reads_list,
                                                         reads_report_df[["Well_number"]],
                                                         wells_vec = wells_vec
                                                         )
                                )



  ## Create a data frame of individual ZMWs

  ZMW_mat <- do.call(rbind, well_stats_mat_list)


  report_for_wells_df <- do.call(rbind.data.frame,
                                 c(split(reads_report_df, reads_report_df[["Well_number"]]),
                                   list(stringsAsFactors = FALSE,
                                        make.row.names = FALSE
                                        )
                                 ))

  are_real_wells <- report_for_wells_df[["Well_number"]] %in% wells_vec

  individual_ZMWs_df <- data.frame(
    report_for_wells_df[are_real_wells, ],
    ZMW_mat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  stopifnot(identical(individual_ZMWs_df[["ZMW"]], barcodes_df[["ZMW"]]))

  barcodes_df <- barcodes_df[, !(names(barcodes_df) %in% names(individual_ZMWs_df))]

  # assign("delete_df", data.frame(individual_ZMWs_df,
  #                                  barcodes_df[are_real_wells, ],
  #                                  stringsAsFactors = FALSE
  #                                  ),
  #        envir = globalenv())

  individual_ZMWs_df <- data.frame(individual_ZMWs_df,
                                   barcodes_df,
                                   stringsAsFactors = FALSE
                                   )



  have_correct_row <- individual_ZMWs_df[["Starts_with_row_barcode"]] |
                      ((individual_ZMWs_df[["Row_bc_length"]] >= bc_length_cutoff) &
                        individual_ZMWs_df[["Correct_row_flank"]])

  have_correct_column <- individual_ZMWs_df[["Ends_with_column_barcode"]] |
                         ((individual_ZMWs_df[["Column_bc_length"]] >= bc_length_cutoff) &
                           individual_ZMWs_df[["Correct_column_flank"]])


  pass_bc <- have_correct_row & have_correct_column &
             (individual_ZMWs_df[["BC_combined_score"]] >= bc_min_comb_score) &
             (individual_ZMWs_df[["BC_score_lead"]] >= bc_min_score_lead)

  pass_rq <- individual_ZMWs_df[["Mean_quality"]] >= min_mean_quality
  pass_filters <- pass_bc & pass_rq

  individual_ZMWs_df[["Passes_filters"]]         <- as.integer(pass_filters)
  individual_ZMWs_df[["Passes_barcode_filters"]] <- as.integer(pass_bc)
  individual_ZMWs_df[["Passes_read_filters"]]    <- as.integer(pass_rq)

  alterations_mat <- AddAlterationCategories(extracted_df, individual_ZMWs_df[["ZMW"]])
  individual_ZMWs_df <- data.frame(
    individual_ZMWs_df,
    alterations_mat,
    stringsAsFactors = FALSE
  )

  all_columns <- c(
    "ZMW", "Well_number", "Length",

    "Correct_barcodes", "Correct_row", "Correct_column",
    "Correct_flanks", "Correct_row_flank", "Correct_column_flank",

    "Starts_with_row_barcode", "Ends_with_column_barcode",

    "Row_bc_length", "Column_bc_length",

    "Row_mean_quality",  "Column_mean_quality", "Row_barcode", "Column_barcode",
    "Row_quality", "Column_quality",

    paste0(c("TpR_DHFR", paste0("sg", 1:4), paste0("sg", 1:4, "_cr", 1:4)),
           "_category"
           ),

    "BC_combined_score", "BC_score_lead", "Mean_quality",
    "Passes_filters", "Passes_barcode_filters", "Passes_read_filters",
    "Contam_guides", "Contam_genes", "Contam_well",
    "Min_distance", "Random_distance", "Mean_distance", "Distances",
    "sg1_cr1", "sg2_cr2", "sg3_cr3", "sg4_cr4", "at_least_1", "at_least_2",
    "at_least_3", "all_4", "all_4_promoters", "whole_plasmid",
    "Orientation_fwd"
  )
  assign("delete_individual_ZMWs_df", individual_ZMWs_df, envir = globalenv())
  individual_ZMWs_df <- individual_ZMWs_df[, all_columns]


  contamin_mat_list <- lapply(contamin_mat_list, function(x) x[, are_real_wells])

  results_list <- list(
    "individual_reads_df" = individual_ZMWs_df,
    "contamin_mat_list"   = contamin_mat_list
  )
  return(results_list)
}




SummarizeWells <- function(analysis_list,
                           use_zmws = NULL,
                           close_well_range = 3L,
                           wells_vec = seq_len(384)
                           ) {

  stopifnot(identical(unique(vapply(analysis_list[["contamin_mat_list"]], ncol, integer(1))),
                      nrow(analysis_list[["individual_reads_df"]])
                      ))

  reads_df <- analysis_list[["individual_reads_df"]]
  contamin_mat_list <- analysis_list[["contamin_mat_list"]]
  if (!(is.null(use_zmws))) {
    zmw_matches <- match(use_zmws, reads_df[["ZMW"]])
    stopifnot(!(anyNA(zmw_matches)))
    reads_df <- reads_df[zmw_matches, ]
    row.names(reads_df) <- NULL
    contamin_mat_list <- lapply(contamin_mat_list, function(x) x[, zmw_matches])
  }

  unfiltered_summary_df <- CreateSummaryDf(reads_df,
                                           filter_reads = FALSE,
                                           close_well_range = close_well_range,
                                           wells_vec = wells_vec
                                           )
  filtered_summary_df   <- CreateSummaryDf(reads_df,
                                           filter_reads = TRUE,
                                           close_well_range = close_well_range,
                                           wells_vec = wells_vec
                                           )

  cross_contam_mat <- CreateCrossContamMat(contamin_mat_list,
                                           reads_df[["Well_number"]],
                                           wells_vec = wells_vec
                                           )

  results_list <- list(
    "original_summary_df" = unfiltered_summary_df,
    "filtered_summary_df" = filtered_summary_df,
    "individual_reads_df" = reads_df,
    "contaminations_mat"  = cross_contam_mat
  )
  return(results_list)
}







DistanceForContam <- function(contamination_mat, bam_wells_vec, wells_vec = seq_len(384)) {

  have_wells <- !(is.na(bam_wells_vec))
  stopifnot("manhattan_dist_list" %in% ls(envir = globalenv()))
  distance_mat <- matrix(nrow = nrow(contamination_mat),
                         ncol = ncol(contamination_mat)
                         )
  for (i in wells_vec) {
    are_this_well <- bam_wells_vec[have_wells] %in% i
    are_contam <- contamination_mat[, are_this_well] >= 1
    sub_mat <- distance_mat[, are_this_well]
    contam_rows <- which(rowSums(are_contam) >= 1)
    for (row_i in contam_rows) {
      sub_mat[row_i, are_contam[row_i, ]] <- manhattan_dist_list[[i]][[row_i]]
    }
    distance_mat[, are_this_well] <- sub_mat
  }
  return(distance_mat)
}





# Define functions for exporting data -------------------------------------

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
  indiv_reads_df[["Mean_quality"]] <- round(indiv_reads_df[["Mean_quality"]], digits = 1)
  exclude_columns <- c("Min_distance", "Random_distance", "Mean_distance",
                       "Passes_barcode_filters", "Passes_read_filters"
                       )
  indiv_reads_df <- indiv_reads_df[, !(colnames(indiv_reads_df) %in% exclude_columns)]
  ExportTable(indiv_reads_df, ...)
}



ExportSummaryTable <- function(summary_df, ...) {
  p_val_vec <- signif(summary_df[["Distance_p_value"]], 1)
  p_val_vec <- vapply(signif(p_val_vec, 1), format, scientific = 9L, "")
  summary_df[["Distance_p_value"]] <- p_val_vec
  are_integer <- vapply(summary_df, is.integer, logical(1))
  are_numeric <- vapply(summary_df, is.numeric, logical(1))
  for (i in which(!(are_integer) & are_numeric)) {
    summary_df[[i]] <- round(summary_df[[i]], digits = 1)
  }
  ExportTable(summary_df, ...)
}




