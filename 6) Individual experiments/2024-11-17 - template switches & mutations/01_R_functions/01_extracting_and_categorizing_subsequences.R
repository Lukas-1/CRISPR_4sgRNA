## 2024-11-17

## This script contains helper functions which rely on functions from the following scripts:
# 2020-08-29 - PacBio - first 384-well plate/1) R functions/07) Categorizing subsequences of reads aligned to the reference.R
# 2022-01-05 - first nanopore sequencing run/01_R_scripts/1_R_functions/03_extracting_aligned_sgRNAs.R


# Define maps -------------------------------------------------------------

features_list <- list(
  "promoter1_hU6"  = c(177, 426),
  "sg1"            = c(427, 446),
  "tracrRNA1"      = c(447, 532),
  "polyT_1"        = c(533, 539),

  "EM7_promoter"   = c(540, 587),
  "pre_TpR"        = c(588, 605),
  "TpR_DHFR"       = c(606, 842),
  "polyT_TpR"      = c(843, 849),

  "promoter2_mU6"  = c(850, 1165),
  "sg2"            = c(1166, 1185),
  "tracrRNA2"      = c(1186, 1273),
  "polyT_2"        = c(1274, 1280),

  "promoter3_hH1"  = c(1281, 1504),
  "sg3"            = c(1505, 1524),
  "tracrRNA3"      = c(1525, 1612),
  "polyT_3"        = c(1613, 1619),

  "promoter4_h7SK" = c(1620, 1863),
  "sg4"            = c(1864, 1883),
  "tracrRNA4"      = c(1884, 1969),
  "polyT_4"        = c(1970, 1976)
)



# General utility functions -----------------------------------------------

CheckThatIntegerVectorIsInOrder <- function(my_factor) {
  stopifnot(identical(length(unique(my_factor)),
                      length(rle(my_factor)[["lengths"]])
                      )
            )
}



# Preparatory functions ---------------------------------------------------

TweakFeaturesDf <- function(input_df) {
  exclude_features <- c("column_barcode", "column_primer", "row_primer", "row_barcode",
                      paste0("sg", 1:4, "_cr", 1:4)
                      )
  features_df <- input_df[!(input_df[, "Feature"] %in% exclude_features), ]
  row.names(features_df) <- NULL

  for (column in c("Start", "End")) {
    features_df[[column]] <- features_df[[column]] - 10L
  }
  features_df[, "Length"] <- features_df[, "End"] - features_df[, "Start"] + 1L
  if ("amplicon_ref" %in% ls(envir = globalenv())) {
    features_df[, "Template"] <- mapply(function(x, y) substr(amplicon_ref, x, y),
                                        features_df[, "Start"],
                                        features_df[, "End"]
                                        )
  }
  return(features_df)
}



# Functions that use output from the pairwiseAlignment R function ---------

FilterAlignmentsDf <- function(use_align_df, mapped_df) {
  are_included <- (mapped_df[, "Num_matched_sgRNAs"] >= 2) &
                  (mapped_df[, "Num_full_passes"] >= 7) &
                  (mapped_df[, "Read_quality"] >= 0.9999)
  included_indices <- mapped_df[, "Read_number"][are_included]

  use_align_df <- data.frame(
    "Read_number" = seq_len(nrow(use_align_df)),
    use_align_df
  )
  use_align_df <- use_align_df[included_indices, ]
  row.names(use_align_df) <- NULL
  return(use_align_df)
}


ExtractAllFeatures <- function(use_align_df) {

  use_align_df[, "Mean_quality"] <- GetMeanQuality(use_align_df[, "Read_sequence"])

  num_reads <- nrow(use_align_df)
  reads_per_chunk <- 10000
  num_chunks <- ceiling(num_reads / reads_per_chunk)
  chunks_vec <- rep(seq_len(num_chunks), each = reads_per_chunk)[seq_len(num_reads)]
  chunks_list <- vector(mode = "list", length = num_chunks)
  first_vec <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[1]]))
  last_vec  <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[length(x)]]))
  chunk_numbers <- format(seq_len(num_chunks))
  for (i in seq_len(num_chunks)) {
    are_this_chunk <- chunks_vec == i
    message("Processing chunk #", chunk_numbers[[i]], " of ",
            chunk_numbers[[length(chunk_numbers)]],  " (extracting reads ",
            first_vec[[i]], " to ", last_vec[[i]], ")..."
            )
    sub_df <- ExtractAlignedSequences(use_align_df[are_this_chunk, ])
    chunks_list[[i]] <- sub_df
  }
  extracted_df <- do.call(rbind.data.frame,
                          c(chunks_list,
                            stringsAsFactors = FALSE,
                            make.row.names = FALSE
                          ))
  extracted_df <- data.frame(
    "Read_number" = rep(use_align_df[, "Read_number"], each = nrow(features_df)),
    extracted_df
  )
  return(extracted_df)
}



# Helper functions for categorizing aligned sequences ---------------------

AddThreeCategories <- function(extracted_df, features_df) {
  matches_vec <- match(extracted_df[, "Feature"], features_df[, "Feature"])
  extracted_df[, "Template"] <- features_df[, "Template"][matches_vec]
  categorized_df <- ThreeBasicCategories(extracted_df, verbose = TRUE)
  categorized_df[, "Template"] <- NULL
  return(categorized_df)
}


AddGuideData <- function(categor_df, mapped_df) {
  CheckThatIntegerVectorIsInOrder(categor_df[, "Read_number"])
  read_numbers <- unique(categor_df[, "Read_number"])
  matches_vec <- match(read_numbers, mapped_df[, "Read_number"])
  stopifnot(!(anyNA(matches_vec)))
  for (i in 1:4) {
    are_correct <- !(is.na(mapped_df[matches_vec, paste0("Plasmid_sg", i)]))
    are_this_sg <- categor_df[, "Feature"] %in% paste0("sg", i)
    categor_df[, "Category"][are_this_sg][are_correct] <- "Correct"
    categor_df[, "Is_correct"][are_this_sg][are_correct] <- TRUE
    categor_df[, "Num_incorrect"][are_this_sg][are_correct] <- 0L
    categor_df[, "Over_5_percent_incorrect"][are_this_sg][are_correct] <- FALSE
  }
  return(categor_df)
}



# Functions for categorizing subsequences based on single-base data -------

UnMatrix <- function(input_mat) {
  dim(input_mat) <- NULL
  input_mat
}


CategorizeFeaturesFromBases <- function(correct_mat, del_mat) {

  stopifnot("features_df" %in% ls(envir = globalenv()))

  indices_list <- Map(function(x, y) x:y, features_df[, "Start"], features_df[, "End"])
  names(indices_list) <- features_df[, "Feature"]

  all_correct_mat <- do.call(cbind, lapply(indices_list, function(x) {
    rowSums(correct_mat[, x]) == length(x)
  }))
  at_least_95_percent_correct_mat <- do.call(cbind, lapply(indices_list, function(x) {
    rowSums(correct_mat[, x]) >= ceiling(length(x) * 0.95)
  }))
  num_incorrect_mat <- do.call(cbind, lapply(indices_list, function(x) {
    rowSums(!(correct_mat[, x]))
  }))
  num_missing_mat <- do.call(cbind, lapply(indices_list, function(x) {
    rowSums(del_mat[, x])
  }))
  mostly_deleted_mat <- do.call(cbind, lapply(indices_list, function(x) {
    rowSums(del_mat[, x]) > floor(length(x) * 0.5)
  }))
  category_mat <- matrix(nrow = nrow(all_correct_mat),
                         ncol = ncol(all_correct_mat)
                         )
  colnames(category_mat) <- features_df[, "Feature"]
  category_mat[all_correct_mat] <- "Correct"
  category_mat[!(all_correct_mat)] <- "Mutation"
  category_mat[mostly_deleted_mat] <- "Deletion"

  results_df <- data.frame(
    "Feature"                  = rep(features_df[, "Feature"], times = nrow(correct_mat)),
    "Is_correct"               = UnMatrix(t(all_correct_mat)),
    "Num_incorrect"            = UnMatrix(t(num_incorrect_mat)),
    "Over_5_percent_incorrect" = UnMatrix(t(!(at_least_95_percent_correct_mat))),
    "Num_missing"              = UnMatrix(t(num_missing_mat)),
    "Mostly_deleted"           = UnMatrix(t(mostly_deleted_mat)),
    "Category"                 = UnMatrix(t(category_mat))
  )
  return(results_df)
}





