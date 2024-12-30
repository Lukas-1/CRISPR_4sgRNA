## 2024-11-17

## This script contains helper functions which rely on functions from the following scripts:
# 2020-08-29 - PacBio - first 384-well plate/1) R functions/07) Categorizing subsequences of reads aligned to the reference.R
# 2022-01-05 - first nanopore sequencing run/01_R_scripts/1_R_functions/03_extracting_aligned_sgRNAs.R


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
  features_df[, "Template"] <- mapply(function(x, y) substr(amplicon_ref, x, y),
                                      features_df[, "Start"],
                                      features_df[, "End"]
                                      )
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


