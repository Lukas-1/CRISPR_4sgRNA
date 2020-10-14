### 8th September 2020 ###



# Import packages and source code -----------------------------------------

library("ShortRead") # For processing the quality scores




# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
file_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory      <- file.path(file_directory, "2) Input")
R_objects_directory       <- file.path(file_directory, "3) R objects")

file_output_directory     <- file.path(file_directory, "5) Output")
tables_output_directory   <- file.path(file_output_directory, "Tables")

raw_data_directory        <- file.path(file_input_directory, "Raw data")
reanalysis_directory      <- file.path(raw_data_directory, "LukasAnalysis")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "3) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads.RData"))
load(file.path(R_objects_directory, "7) Extract barcode sequences and quality scores.RData"))





# Read in extra data (for checks) -----------------------------------------

ReadTable <- function(file_path, header = TRUE, sep = "\t") {
  read.table(file_path, sep = sep, quote = "", stringsAsFactors = FALSE,
             header = header, row.names = NULL, check.names = FALSE
             )
}

sl7_CCS3_lima_counts_file <- file.path(reanalysis_directory, "SmrtLink7_CCS3/lima", "lukaslimaCCS3SL7output.lima.counts")
sl7_CCS5_lima_counts_file <- file.path(reanalysis_directory, "SmrtLink7_CCS5/lima", "lukaslimaCCS5SL7output.lima.counts")

sl7_ccs3_counts_df <- ReadTable(sl7_CCS3_lima_counts_file)
sl7_ccs5_counts_df <- ReadTable(sl7_CCS5_lima_counts_file)

sl7_ccs3_counts_384_file <- file.path(raw_data_directory, "FGCZ/ccs_3_99/lima_26/lima.lima.384.counts.txt")
sl7_ccs5_counts_384_file <- file.path(raw_data_directory, "FGCZ/ccs_5_999/lima_26/lima.lima.384.counts.txt")

sl7_ccs3_384_counts_df <- ReadTable(sl7_ccs3_counts_384_file)
sl7_ccs5_384_counts_df <- ReadTable(sl7_ccs5_counts_384_file)

ccs3_fgcz_accuracies_df <- ReadTable(file.path(raw_data_directory, "FGCZ/ccs_3_99/lima_26/blastn/identical.perc.txt"))
ccs5_fgcz_accuracies_df <- ReadTable(file.path(raw_data_directory, "FGCZ/ccs_5_999/lima_26/blastn/identical.perc.txt"))







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

  stopifnot(!(any(contain_fwd_mat == 2)))
  stopifnot(!(any(contain_rev_mat == 2)))

  contain_mat <- contain_fwd_mat + contain_rev_mat

  same_guide_twice_mat <- contain_mat > 1
  have_duplication <- rowSums(same_guide_twice_mat) > 0

  sent_message <- FALSE
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
                       "orientation_fwd" = are_fwd,
                       "inconsistent_orientation" = have_duplication | have_inversion
                       )
  return(results_mat)
}




StatsForWell <- function(well_number, use_reads, bam_wells_vec) {

  are_this_well <- bam_wells_vec %in% well_number
  num_reads <- sum(are_this_well)
  well_sequences <- use_reads[["seq"]][are_this_well, ]

  sg_sequences <- vapply(1:4, function(x) guides_ref_list[[x]][[well_number]], "")
  sg_with_promoter_seq <- vapply(1:4, function(x) guides_with_promoters_list[[x]][[well_number]], "")

  contain_sg_cr_mat   <- ContainSequences(sg_sequences, well_sequences)
  contain_sg_prom_mat <- ContainSequences(sg_with_promoter_seq, well_sequences)
  contain_plasmid_mat <- ContainSequences(plasmids_vec[[well_number]], well_sequences)

  colnames(contain_sg_cr_mat)[1:4] <- paste0("sg", 1:4, "_cr", 1:4)
  colnames(contain_sg_cr_mat)[colnames(contain_sg_cr_mat) == "at_least_4"] <- "all_4"
  results_mat <- cbind(
    contain_sg_cr_mat[, 1:8],
    "all_4_promoters" = contain_sg_prom_mat[, "at_least_4"],
    "whole_plasmid"   = contain_plasmid_mat[, 1],
    contain_sg_cr_mat[, "orientation_fwd", drop = FALSE]
  )
  stopifnot(!(any(results_mat[, "all_4_promoters"] < results_mat[, "whole_plasmid"])))

  return(results_mat)
}






GetZMWs <- function(char_vec) {
   as.integer(sapply(strsplit(char_vec, "/", fixed = TRUE), "[[", 2))
}




MakeDfForReads <- function(use_reads,
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
  bam_ZMWs <- GetZMWs(use_reads[["qname"]])

  matches_vec <- match(bam_ZMWs, profile_ZMWs)

  results_df <- data.frame(
    "ZMW"               = bam_ZMWs,
    "Well_number"       = report_df[["Well_number"]][matches_vec],
    "Length"            = width(use_reads[["seq"]]),
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



GetMeanQualities <- function(use_reads) {
  ShortRead::alphabetScore(use_reads[["qual"]]) /
  width(use_reads[["qual"]])
}



## DELETE THIS!!!!
# GetOriginalMeanQualities <- function(lima_reads, ccs_reads) {
#   matches_vec <- match(lima_reads[["qname"]], ccs_reads[["qname"]])
#   stopifnot(!(anyNA(matches_vec)))
#   mean_qualities <- ShortRead::alphabetScore(ccs_reads[["qual"]]) /
#                     width(ccs_reads[["qual"]])
#   return(mean_qualities[matches_vec])
# }


ProcessBarcodesDf <- function(barcodes_df) {
  # matches_vec <- match(lima_reads[["qname"]], ccs_reads[["qname"]])
  # stopifnot(!(anyNA(matches_vec)))
  # mean_qualities <- ShortRead::alphabetScore(ccs_reads[["qual"]]) /
  #                   width(ccs_reads[["qual"]])
  # return(mean_qualities[matches_vec])
}




GetMode <- function(vec) {
   uniqv <- unique(vec)
   uniqv[which.max(tabulate(match(vec, uniqv)))]
}




GetContaminationMat <- function(query_seq, use_reads, bam_wells_vec) {

  stopifnot("guides_ref_list" %in% ls(envir = globalenv()))

  are_usual_length <- nchar(query_seq) == GetMode(nchar(query_seq))

  PDict_object <- PDict(query_seq[are_usual_length], max.mismatch = 0)

  have_well <- !(is.na(bam_wells_vec))

  counts_mat <- vcountPDict(PDict_object, use_reads[["seq"]][have_well, ])

  row_list <- vector(mode = "list", length = 384)
  row_list[are_usual_length] <- lapply(seq_len(nrow(counts_mat)),
                                       function(x) counts_mat[x, ]
                                       )

  unusual_seq <- DNAStringSet(query_seq[!(are_usual_length)])

  unusual_list <- lapply(seq_len(sum(!(are_usual_length))),
                         function(x) vcountPattern(unusual_seq[[x]],
                                                   use_reads[["seq"]][have_well, ]
                                                   )
                         )

  row_list[!(are_usual_length)] <- unusual_list

  all_counts_mat <- do.call(rbind, row_list)
  stopifnot(all(unique(as.vector(all_counts_mat))) %in% 0:1)
  mode(all_counts_mat) <- "integer"

  ## If the sequence of two wells is exactly the same,
  ## perhaps it shouldn't count as a contamination...
  same_seq_list <- lapply(seq_len(384),
                          function(x) which(query_seq %in% vapply(guides_ref_list, function(y) y[[x]], ""))
                          )
  are_this_well_list <- lapply(same_seq_list,
                               function(x) bam_wells_vec[have_well] %in% x
                               )

  contamination_mat <- all_counts_mat
  for (i in seq_len(384)) {
    contamination_mat[i, are_this_well_list[[i]]] <- 0L
  }
  return(contamination_mat)
}





CreateSummaryDf <- function(reads_df,
                            filter_reads = FALSE,
                            close_well_range = 3L
                            ) {

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


  ## Process contaminations / wrong barcodes

  are_contaminated <- reads_df[["Contam_guides"]] >= 1
  are_close_rand <- (reads_df[["Random_distance"]] <= close_well_range) %in% TRUE

  num_other_genes_vec <- tapply(reads_df[["Contam_well"]],
                                wells_fac,
                                function(x) {
                                  wells_list <- strsplit(x[!(is.na(x))], ", ", fixed = TRUE)
                                  wells_vec <- unlist(wells_list, use.names = FALSE)
                                  length(unique(wells_vec))
                                })

  num_contaminants <- tapply(are_contaminated, wells_fac, sum)
  num_close_contaminants_rand <- tapply(are_contaminated & are_close_rand, wells_fac, sum)

  plate_num_close_vec <- vapply(manhattan_dist_list,
                                function(x) sum(x <= close_well_range),
                                integer(1)
                                ) - 1L

  num_expected_close_wells <- plate_num_close_vec * (num_contaminants / 383)

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
    "Well_number"            = seq_len(384),
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
    "Expected_distance"      = expected_distances_vec,
    "Distance_p_value"       = p_val_vec,
    stringsAsFactors         = FALSE,
    row.names                = NULL
  )
  return(results_df)
}



CreateCrossContamMat <- function(contamin_mat_list, reads_report_df) {

  have_well <- !(is.na(reads_report_df[["Well_number"]]))

  all_comb_mat <- t(combn(seq_len(384), 2))
  colnames(all_comb_mat) <- paste0("well", 1:2, "_ID")

  counts_wells_mat_list <- lapply(1:4, function(sg_number) {
    do.call(rbind,
            tapply(seq_len(sum(have_well)),
                   reads_report_df[["Well_number"]][have_well],
                   function(x) rowSums(contamin_mat_list[[sg_number]][, x])
                   )
            )
  })

  cross_contam_mat <- do.call(cbind, lapply(1:4, function(sg_number) {
    use_mat <- counts_wells_mat_list[[sg_number]]
    sub_mat <- t(mapply(function(x, y) c(use_mat[x, y], use_mat[y, x]),
                        all_comb_mat[, 1],
                        all_comb_mat[, 2]
                        )
                 )
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



AnalyzeWells <- function(use_reads,
                         report_df,
                         ccs_reads         = NULL,
                         counts_df         = NULL,
                         counts_384_df     = NULL,
                         close_well_range  = 3L,
                         bc_min_comb_score = 80,
                         bc_min_score_lead = 40,
                         min_mean_quality  = 85
                         ) {

  stopifnot("manhattan_dist_list" %in% ls(envir = globalenv()))

  reads_report_df <- MakeDfForReads(use_reads, report_df)
  reads_report_df[["Mean_quality"]] <- GetMeanQualities(use_reads)

  pass_bc <- (reads_report_df[["BC_combined_score"]] >= bc_min_comb_score) &
             (reads_report_df[["BC_score_lead"]] >= bc_min_score_lead)
  pass_rq <- reads_report_df[["Mean_quality"]] >= min_mean_quality
  pass_filters <- pass_bc & pass_rq
  reads_report_df[["Passes_filters"]]         <- as.integer(pass_filters)
  reads_report_df[["Passes_barcode_filters"]] <- as.integer(pass_bc)
  reads_report_df[["Passes_read_filters"]]    <- as.integer(pass_rq)


  ## Add data on contaminations / wrong barcodes

  contamin_mat_list <- lapply(guides_ref_list,
                              function(x) GetContaminationMat(x,
                                                              use_reads,
                                                              reads_report_df[["Well_number"]]
                                                              )
                              )


  contamin_mat <- Reduce(`+`, contamin_mat_list)
  have_contamin_mat <- contamin_mat >= 1
  distance_mat <- DistanceForContam(contamin_mat, reads_report_df[["Well_number"]])

  cross_contam_mat <- CreateCrossContamMat(contamin_mat_list, reads_report_df)

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

  well_stats_mat_list <- lapply(seq_len(384),
                                function(x) StatsForWell(x,
                                                         use_reads,
                                                         reads_report_df[["Well_number"]]
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

  individual_ZMWs_df <- data.frame(
    report_for_wells_df,
    ZMW_mat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  ## Create summary data frames

  unfiltered_summary_df <- CreateSummaryDf(individual_ZMWs_df,
                                           filter_reads = FALSE,
                                           close_well_range = close_well_range
                                           )
  filtered_summary_df   <- CreateSummaryDf(individual_ZMWs_df,
                                           filter_reads = TRUE,
                                           close_well_range = close_well_range
                                           )

  results_list <- list(
    "original_summary_df" = unfiltered_summary_df,
    "filtered_summary_df" = filtered_summary_df,
    "individual_reads_df" = individual_ZMWs_df,
    "contaminations_mat"  = cross_contam_mat
  )
  return(results_list)
}



DistanceForContam <- function(contamination_mat, bam_wells_vec) {

  have_wells <- !(is.na(bam_wells_vec))
  stopifnot("manhattan_dist_list" %in% ls(envir = globalenv()))
  distance_mat <- matrix(nrow = nrow(contamination_mat),
                         ncol = ncol(contamination_mat)
                         )
  for (i in seq_len(384)) {
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








# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)







# Create the summary data frames ------------------------------------------

sl7_ccs3_df_list <- AnalyzeWells(sl7_ccs3_lima,
                                 sl7_ccs3_report_df,
                                 sl7_ccs3_ccs,
                                 sl7_ccs3_counts_df,
                                 sl7_ccs3_384_counts_df
                                 )
sl7_ccs5_df_list <- AnalyzeWells(sl7_ccs5_lima,
                                 sl7_ccs5_report_df,
                                 sl7_ccs5_ccs,
                                 sl7_ccs5_counts_df,
                                 sl7_ccs5_384_counts_df
                                 )

sl9_ccs3_df_list <- AnalyzeWells(sl9_ccs3_lima, sl9_ccs3_report_df, sl9_ccs3_ccs)
sl9_ccs5_df_list <- AnalyzeWells(sl9_ccs5_lima, sl9_ccs5_report_df, sl9_ccs3_ccs)


## Check for equivalence

CheckCountsAreEqual <- function(reanalysis_df, fgcz_df) {
  stopifnot(identical(reanalysis_df[["Count_total"]], fgcz_df[["tot"]]))
  for (i in 1:4) {
    reanalysis_col <- paste0("Count_sg", i, "_cr", i)
    fgcz_col <- paste0("sg", i, "_cr", i, "_identical")
    stopifnot(identical(reanalysis_df[[reanalysis_col]], fgcz_df[[fgcz_col]]))
  }
  return(invisible(NULL))
}


CheckCountsAreEqual(sl7_ccs3_df_list[["original_summary_df"]],
                    ccs3_fgcz_accuracies_df
                    )
CheckCountsAreEqual(sl7_ccs5_df_list[["original_summary_df"]],
                    ccs5_fgcz_accuracies_df
                    )




# Export tables -----------------------------------------------------------

ExportSummaryTable(sl7_ccs3_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_original"
                   )
ExportSummaryTable(sl7_ccs3_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS3_99_summary_per_well_filtered"
                   )
ExportIndivTable(sl7_ccs3_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS3_99_individual_reads"
                 )
ExportTable(sl7_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS3_99_contaminations"
            )


ExportSummaryTable(sl7_ccs5_df_list[["original_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_original"
                   )
ExportSummaryTable(sl7_ccs5_df_list[["filtered_summary_df"]],
                   "SmrtLink7/SmrtLink7_CCS5_999_summary_per_well_filtered"
                   )
ExportIndivTable(sl7_ccs5_df_list[["individual_reads_df"]],
                 "SmrtLink7/SmrtLink7_CCS5_999_individual_reads"
                 )
ExportTable(sl7_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink7/SmrtLink7_CCS5_999_contaminations"
            )



ExportSummaryTable(sl9_ccs3_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_original"
                   )
ExportSummaryTable(sl9_ccs3_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS3_99_summary_per_well_filtered"
                   )
ExportIndivTable(sl9_ccs3_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS3_99_individual_reads"
                 )
ExportTable(sl9_ccs3_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS3_99_contaminations"
            )

ExportSummaryTable(sl9_ccs5_df_list[["original_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_original"
                   )
ExportSummaryTable(sl9_ccs5_df_list[["filtered_summary_df"]],
                   "SmrtLink9/SmrtLink9_CCS5_999_summary_per_well_filtered"
                   )
ExportIndivTable(sl9_ccs5_df_list[["individual_reads_df"]],
                 "SmrtLink9/SmrtLink9_CCS5_999_individual_reads"
                 )
ExportTable(sl9_ccs5_df_list[["contaminations_mat"]],
            "SmrtLink9/SmrtLink9_CCS5_999_contaminations"
            )









# Save data ---------------------------------------------------------------

save(list = c("sl7_ccs3_df_list", "sl7_ccs5_df_list",
              "sl9_ccs3_df_list", "sl9_ccs5_df_list"
              ),
     file = file.path(R_objects_directory, "8) Process demultiplexed PacBio reads.RData")
     )





# Check for associations with contaminations ------------------------------

use_reads_df <- sl7_ccs3_df_list[["individual_reads_df"]]

bc_comb_cutoff <- 80
bc_lead_cutoff <- 40
close_well_range <- 3L
are_poor_barcodes <- (use_reads_df[["BC_combined_score"]] < bc_comb_cutoff) &
                     (use_reads_df[["BC_score_lead"]] < bc_lead_cutoff)

are_close_contams <- use_reads_df[["Random_distance"]] <= close_well_range

are_standard_lengths <- use_reads_df[["Length"]] %in% 2223:2229

are_contaminated <- use_reads_df[["Contam_guides"]] >= 1


fisher.test(table("are_contam" = are_contaminated,
                  "bad_bc"     = are_poor_barcodes
                  ))

fisher.test(table("are_close" = are_close_contams,
                  "bad_bc"    = are_poor_barcodes
                  ))
wilcox.test(use_reads_df[["Random_distance"]] ~
              are_poor_barcodes
            )



fisher.test(table("are_contam"      = are_contaminated,
                  "standard_length" = are_standard_lengths
                  ))

fisher.test(table("standard_length" = are_standard_lengths,
                  "bad_bc"          = are_poor_barcodes
                  ))




# sl7_ccs3_df_list[[2]][["BC_score_lead"]][sl7_ccs3_df_list[[2]][["Contam_guides"]] >= 1]
#
#
# counts_wells_mat <- do.call(rbind,
#                             tapply(seq_len(ncol(all_counts_mat)),
#                                    reads_report_df[["Well_number"]][have_well],
#                                    function(x) rowSums(all_counts_mat[, x])
#                                    )
#                             )
#
#
# counts_wells_copy_mat <- counts_wells_mat
# diag(counts_wells_copy_mat) <- 0
#
#
#   MakeEmptyPlot()
#
#
#   ## Draw the heatmap
#
#   numeric_mat <- counts_wells_copy_mat / max(counts_wells_copy_mat)
#
#   ColorFunction <- colorRampPalette(brewer.pal(9, "Blues"))
#
#   my_breaks <- c(0, 0.8, 0.9, 1.01)
#
#   my_cmap <- makecmap(numeric_mat, colFn = ColorFunction)
#   my_color_mat <- cmap(numeric_mat, my_cmap)
#   Do_cimage(my_color_mat)
#   x_range <- par("usr")[[2]] - par("usr")[[1]]
#
#










