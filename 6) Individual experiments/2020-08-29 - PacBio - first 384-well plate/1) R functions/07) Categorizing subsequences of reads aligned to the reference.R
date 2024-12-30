### 20th October 2020 ###



# Import packages and source code -----------------------------------------

library("Biostrings")
library("data.table")



# Define feature positions ------------------------------------------------

features_list <- list(
  "column_barcode" = c(1, 10),
  "column_primer"  = c(11, 30),

  "promoter1_hU6"  = c(177, 426),
  "sg1"            = c(427, 446),
  "sg1_cr1"        = c(427, 532),
  "tracrRNA1"      = c(447, 532),

  "TpR_DHFR"       = c(606, 842),

  "promoter2_mU6"  = c(850, 1165),
  "sg2"            = c(1166, 1185),
  "sg2_cr2"        = c(1166, 1273),
  "tracrRNA2"      = c(1186, 1273),

  "promoter3_hH1"  = c(1281, 1504),
  "sg3"            = c(1505, 1524),
  "sg3_cr3"        = c(1505, 1612),
  "tracrRNA3"      = c(1525, 1612),

  "promoter4_h7SK" = c(1620, 1863),
  "sg4"            = c(1864, 1883),
  "sg4_cr4"        = c(1864, 1969),
  "tracrRNA4"      = c(1884, 1969),

  "row_primer"     = c(2216, 2235),
  "row_barcode"    = c(2236, 2245)
)




# Define functions --------------------------------------------------------

AdjustForsgRNALength <- function(features_df, sg_vec) {
  stopifnot(length(sg_vec) == 4)
  sg_lengths <- nchar(sg_vec)
  results_mat <- as.matrix(features_df[, c("Start", "End")])
  if (!(all(sg_lengths == 20))) {
    for (i in 1:4) {
      if (sg_lengths[[i]] != 20) {
        sg_row <- match(paste0("sg", i), features_df[, "Feature"])
        addend <- -20L + sg_lengths[[i]]
        results_mat[sg_row, "End"] <- results_mat[sg_row, "End"] + addend
        results_mat[(sg_row + 1):nrow(results_mat), ] <- results_mat[(sg_row + 1):nrow(results_mat), ] + addend
      }
    }
  }
  return(results_mat)
}


FeaturesListToDf <- function(features_list) {
  features_mat <- do.call(rbind, features_list)
  mode(features_mat) <- "integer"
  rownames(features_mat) <- NULL
  colnames(features_mat) <- c("Start", "End")
  features_df <- data.frame("Feature" = names(features_list),
                            features_mat,
                            stringsAsFactors = FALSE
                            )
  return(features_df)
}



ExtractAlignedSequences <- function(ccs_df,
                                    alignments_df,
                                    ID_column = "Well_number",
                                    unique_IDs = seq_len(384),
                                    process_results = TRUE
                                    ) {

  stopifnot(all(c("features_df", "features_templates_list", "features_indices_list") %in% ls(envir = globalenv())))

  num_features <- nrow(features_df)
  features_vec <- features_df[["Feature"]]

  well_df_list <- lapply(seq_along(unique_IDs), function(x) {

    current_ID <- unique_IDs[[x]]

    message(paste0("Processing the well: ", current_ID, "..."))

    are_this_ID <- alignments_df[[ID_column]] == current_ID

    if (!(any(are_this_ID))) {
      message("No reads are available for this well ID! It was skipped.")
      return(NULL)
    }

    aligned_plasmid_vec <- alignments_df[["Aligned_plasmid"]][are_this_ID]
    aligned_read_vec <- alignments_df[["Aligned_read"]][are_this_ID]

    aligned_plasmid_char_list <- strsplit(aligned_plasmid_vec, "")
    aligned_read_char_list <- strsplit(aligned_read_vec, "")

    this_ID_zmws <- alignments_df[["ZMW"]][are_this_ID]
    zmw_matches <- match(this_ID_zmws, ccs_df[["ZMW"]])

    unaligned_qual_vec <- ccs_df[["Quality"]][zmw_matches]
    this_ID_fwd <- alignments_df[["Orientation_fwd"]][are_this_ID]

    unaligned_qual_vec[!(this_ID_fwd)] <- reverse(unaligned_qual_vec[!(this_ID_fwd)])

    unaligned_qual_char_list <- strsplit(unaligned_qual_vec, "")

    num_reads <- length(this_ID_zmws)

    if (is.character(current_ID)) {
      use_index <- x
    } else {
      use_index <- current_ID
    }

    extracted_mat_list <- lapply(seq_len(num_reads), function(x) {

      plasmid_char_numbers <- cumsum(aligned_plasmid_char_list[[x]] != "-")

      are_gaps <- aligned_read_char_list[[x]] == "-"
      read_char_numbers <- cumsum(!(are_gaps))
      aligned_qual_char_vec <- unaligned_qual_char_list[[x]][read_char_numbers]
      aligned_qual_char_vec[are_gaps] <- " "

      original_read_length <- length(unaligned_qual_char_list[[x]])
      stopifnot(max(read_char_numbers) == original_read_length)
      stopifnot(sum(!(are_gaps)) == original_read_length)

      feature_indices <- lapply(features_indices_list[[use_index]],
                                function(x) which(plasmid_char_numbers %in% x)
                                )
      aligned_templates <- vapply(feature_indices,
                                  function(y) paste0(aligned_plasmid_char_list[[x]][y], collapse = ""),
                                  ""
                                  )
      extracted_sequences <- vapply(feature_indices,
                                    function(y) paste0(aligned_read_char_list[[x]][y], collapse = ""),
                                    ""
                                    )
      extracted_qualities <- vapply(feature_indices,
                                    function(x) paste0(aligned_qual_char_vec[x], collapse = ""),
                                    ""
                                    )
      extracted_mat <- cbind(
        "Aligned_template" = aligned_templates,
        "Aligned_read"     = extracted_sequences,
        "Quality"          = extracted_qualities
      )
      return(extracted_mat)
    })
    well_mat <- do.call(rbind, extracted_mat_list)

    well_df <- data.frame(
      "Combined_ID" = current_ID,
      "ZMW"         = rep(this_ID_zmws, each = num_features),
      "Feature"     = rep(features_vec, times = num_reads),
      "Template"    = rep(features_templates_list[[use_index]], times = num_reads),
      well_mat,
      stringsAsFactors = FALSE
    )
    return(well_df)
  })
  message("Collating the final data frame...")
  results_df <- data.table::rbindlist(well_df_list)
  data.table::setDF(results_df)

  if (ID_column == "Combined_ID") {
    ID_matches <- match(results_df[, "Combined_ID"], ccs_df[, "Combined_ID"])
    results_df[["Plate_number"]] <- ccs_df[, "Plate_number"][ID_matches]
    results_df[["Well_number"]] <- ccs_df[, "Well_number"][ID_matches]
    names(results_df)[names(results_df) == "Combined_ID"] <- "Combined_ID"
  } else if (ID_column == "Well_number") {
    names(results_df)[names(results_df) == "Combined_ID"] <- "Well_number"
  }
  results_df[["Template"]] <- toupper(results_df[["Template"]])

  if ("Original_ZMW" %in% names(ccs_df)) {
    zmw_matches <- match(results_df[["ZMW"]], ccs_df[["ZMW"]])
    results_df[["Original_ZMW"]] <- ccs_df[zmw_matches, "Original_ZMW"]
    results_df[["SmrtCell"]] <- ccs_df[zmw_matches, "SmrtCell"]
  }

  first_columns <- c("Combined_ID", "Plate_number", "Well_number",
                     "ZMW", "Original_ZMW", "SmrtCell"
                     )
  first_columns <- intersect(first_columns, names(results_df))
  all_columns <- c(first_columns, setdiff(names(results_df), first_columns))
  results_df <- results_df[, all_columns]

  if (process_results) {
    message("Processing the final data frame...")
    results_df <- ProcessExtractedDf(results_df, unique_IDs)
  }
  return(results_df)
}



GetMeanQualityInChunks <- function(phred_vec) {
  chunk_size <- 500 * 1000
  num_chunks <- ceiling(length(phred_vec) / chunk_size)
  chunks_vec <- rep(seq_len(num_chunks), each = chunk_size)[seq_along(phred_vec)]
  mean_quality_list <- tapply(phred_vec,
                              chunks_vec,
                              GetMeanQualityNoGaps,
                              simplify = FALSE
                              )
  mean_quality_vec <- unlist(mean_quality_list, use.names = FALSE)
  return(mean_quality_vec)
}



GetMeanQualityNoGaps <- function(phred_vec) {
  quality_int_list <- as(PhredQuality(phred_vec), "IntegerList")
  quality_int_list <- lapply(quality_int_list, function(x) x[x != -1L])
  mean_quality_vec <- vapply(quality_int_list,
                             function(x) if (length(x) == 0) NA_real_ else mean(x),
                             numeric(1)
                             )
  return(mean_quality_vec)
}




ProcessByPlate <- function(input_df, UseFunction, is_df = TRUE, message_prefix = NULL) {
  if (is.null(message_prefix)) {
    message_prefix <- "Processing "
  }
  all_plates <- setdiff(input_df[, "Plate_number"], NA)
  results_list <- lapply(all_plates, function(x) {
    message(paste0(message_prefix, "plate #", x, "..."))
    sub_df <- input_df[input_df[, "Plate_number"] %in% x, ]
    row.names(sub_df) <- NULL
    UseFunction(sub_df)
  })
  if (is_df) {
    result_object <- do.call(rbind.data.frame,
                             c(results_list,
                               stringsAsFactors = FALSE,
                               make.row.names = FALSE
                               )
                             )
  } else {
    result_object <- unlist(results_list, use.names = FALSE)
  }
  return(result_object)
}



ThreeBasicCategories <- function(extracted_df, verbose = FALSE) {

  ShowMessage <- function(x) if (verbose) message(x)

  ShowMessage("Computing mean qualities...")
  mean_quality_vec <- GetMeanQualityNoGaps(extracted_df[, "Quality"])

  ShowMessage("Splitting sequences into single-character strings...")
  aligned_read_chars <- strsplit(extracted_df[, "Aligned_read"], "", fixed = TRUE)
  aligned_template_chars <- strsplit(extracted_df[, "Aligned_template"], "", fixed = TRUE)

  ShowMessage("Counting the number of incorrect bases...")
  num_bases_incorrect <- mapply(function(x, y) sum(x != y),
                                aligned_read_chars,
                                aligned_template_chars
                                )

  ShowMessage("Counting the number of deleted bases...")
  num_bases_lost <- vapply(aligned_read_chars, function(x) sum(x == "-"), integer(1))

  ShowMessage("Combining data...")
  sequence_lengths <- nchar(extracted_df[["Aligned_read"]])
  results_df <- data.frame(
    extracted_df,
    "Mean_quality"             = mean_quality_vec,
    "Is_correct"               = extracted_df[, "Template"] == extracted_df[, "Aligned_read"],
    "Alignment_length"         = sequence_lengths,
    "Num_incorrect"            = num_bases_incorrect,
    "Over_5_percent_incorrect" = (num_bases_incorrect / sequence_lengths) > 0.05,
    "Num_missing"              = num_bases_lost,
    "Mostly_deleted"           = num_bases_lost >= (sequence_lengths / 2)
  )

  ShowMessage("Assigning 3 categories...")
  category_vec <- ifelse(results_df[["Is_correct"]],
                         "Correct",
                         ifelse(results_df[["Mostly_deleted"]],
                                "Deletion",
                                "Mutation"
                                )
                         )
  results_df[["Category"]] <- category_vec
  return(results_df)
}



GetFlankingInsertions <- function(extracted_df) {
  wells_vec <- unique(extracted_df[, "Well_number"])
  insertion_possible <- !(extracted_df[, "Is_correct"] |
                          extracted_df[, "Mostly_deleted"] |
                          extracted_df[, "Is_contamination"]
                          )
  have_flanking_insertion <- rep(FALSE, nrow(extracted_df))
  for (include_tracrRNA in c(FALSE, TRUE)) {
    for (i in 1:4) {
      if (include_tracrRNA) {
        this_feature <- paste0("sg", i, "_cr", i)
      } else {
        this_feature <- paste0("sg", i)
      }
      message(paste0("Checking for flanking insertions for ", this_feature, "..."))
      are_this_feature <- extracted_df[, "Feature"] == this_feature
      for (well_number in wells_vec) {
        are_this_well <- (extracted_df[, "Well_number"] == well_number) & are_this_feature
        well_template <- unique(extracted_df[are_this_well, "Template"])
        are_selected <- are_this_well & insertion_possible
        have_flanking_insertion[are_selected] <- grepl(well_template,
                                                       extracted_df[are_selected, "Aligned_read"],
                                                       fixed = TRUE
                                                       )
      }
    }
  }
  results_vec <- extracted_df[, "Category"]
  results_vec[have_flanking_insertion] <- "Flanking insertion"
  return(results_vec)
}




ProcessExtractedDf <- function(extracted_df, unique_IDs = seq_len(384)) {

  stopifnot(all("sg_sequences_df" %in% ls(envir = globalenv())))

  if ("Plate_number" %in% names(extracted_df)) {
    results_df <- ProcessByPlate(extracted_df,
                                 ThreeBasicCategories,
                                 message_prefix  = "Dividing features into three basic categories for "
                                 )
  } else {
    results_df <- ThreeBasicCategories(extracted_df)
  }

  assign("delete_1_results_df", results_df, envir = globalenv())

  ## Check for contaminations
  contamination_possible <- !(results_df[["Is_correct"]] | results_df[["Mostly_deleted"]])
  sg_list <- lapply(1:4, function(x) sg_sequences_df[[paste0("Sequence_sg", x)]])
  all_guides_vec <- toupper(unlist(sg_list))

  tracrRNA_matches <- match(paste0("tracrRNA", 1:4), features_df[["Feature"]])
  tracrRNA_sequences <- features_templates_list[[1]][tracrRNA_matches]
  all_tracRNAs_vec <- unlist(lapply(tracrRNA_sequences, function(x) paste0(all_guides_vec, x)))

  are_contamination_sg <- rep(FALSE, nrow(results_df))
  for (i in 1:4) {
    are_this_sg <- extracted_df[["Feature"]] == paste0("sg", i)
    are_eligible <- are_this_sg & contamination_possible
    are_contamination_sg[are_eligible] <- extracted_df[["Aligned_read"]][are_eligible] %in% all_guides_vec
  }
  are_contamination_sg_cr <- rep(FALSE, nrow(results_df))
  for (i in 1:4) {
    are_this_sg_cr <- extracted_df[["Feature"]] == paste0("sg", i, "_cr", i)
    are_eligible <- are_this_sg_cr & contamination_possible
    are_contamination_sg_cr[are_eligible] <- extracted_df[["Aligned_read"]][are_eligible] %in% all_tracRNAs_vec
  }
  stopifnot(!(any(are_contamination_sg & are_contamination_sg_cr)))
  are_contamination <- are_contamination_sg | are_contamination_sg_cr
  stopifnot(all(results_df[["Category"]][are_contamination] == "Mutation"))
  results_df[["Is_contamination"]] <- are_contamination
  results_df[["Category"]][are_contamination] <- "Contamination"

  ## Check for small insertions that may or may not be considered to fall within
  ## the sgRNA + tracrRNA region
  if ("Plate_number" %in% names(results_df)) {
    category_vec <- ProcessByPlate(results_df,
                                   GetFlankingInsertions,
                                   is_df = FALSE,
                                   message_prefix = "Checking for flanking insertions for "
                                   )
  } else {
    category_vec <- GetFlankingInsertions(results_df)
  }
  results_df[["Category"]] <- category_vec
  return(results_df)
}



