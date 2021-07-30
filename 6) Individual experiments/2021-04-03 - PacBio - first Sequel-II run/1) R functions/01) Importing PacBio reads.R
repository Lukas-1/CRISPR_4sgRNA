### 11th April 2021 ###




# Define functions --------------------------------------------------------

IntegrateSingleReportDf <- function(sam_df, tidied_report_df, check_zmws = TRUE) {

  ZMW_matches <- match(sam_df[["ZMW"]], tidied_report_df[["ZMW"]])

  matched_report_df <- tidied_report_df[ZMW_matches, ]

  if (check_zmws) {
    stopifnot(identical(sam_df[["ZMW"]], matched_report_df[["ZMW"]]))
  }

  sequence_lengths <- nchar(sam_df[["Sequence"]])
  stopifnot(identical(sequence_lengths, nchar(sam_df[["Quality"]])))
  if ("Pre_clip_start" %in% names(matched_report_df)) {
    start_vec <- matched_report_df[["Pre_clip_start"]]
  } else {
    start_vec <- 1
  }
  if ("Pre_clip_end" %in% names(tidied_report_df)) {
    end_vec <- matched_report_df[["Pre_clip_end"]]
  } else {
    end_vec <- sequence_lengths
  }

  for (column_name in c("Sequence", "Quality")) {
    vec_1 <- substr(sam_df[[column_name]], start_vec, matched_report_df[["Clip_start"]] - 1)
    vec_2 <- substr(sam_df[[column_name]], matched_report_df[["Clip_end"]] + 1, end_vec)
    new_column_names <- paste0("Barcode_", 1:2, "_", tolower(column_name))
    sam_df[[new_column_names[[1]]]] <- vec_1
    sam_df[[new_column_names[[2]]]] <- vec_2
  }
  sam_df[["Barcode_1_length"]] <- nchar(vec_1)
  sam_df[["Barcode_2_length"]] <- nchar(vec_2)

  matched_report_df <- matched_report_df[, names(matched_report_df) != "ZMW"]
  results_df <- data.frame(
    sam_df,
    matched_report_df,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(results_df)
}



IntegrateReportDfs <- function(sam_df, plates_report_df, wells_report_df, use_plates_df) {

  plate_combos <- paste0(use_plates_df[["Barcode_ID"]], "--", use_plates_df[["Barcode_ID"]])
  barcodes_to_plates_map <- use_plates_df[["Plate_number"]]
  names(barcodes_to_plates_map) <- plate_combos
  tidied_plates_df <- TidyReportDf(plates_report_df, barcodes_to_plates_map)
  tidied_wells_df <- TidyReportDf(wells_report_df, barcodes_to_wells_map)

  zmw_matches <- match(tidied_wells_df[["ZMW"]], tidied_plates_df[["ZMW"]])

  add_vec <- tidied_plates_df[["Clip_start"]][zmw_matches] - 1L
  tidied_wells_df[["Clip_start"]]     <- tidied_wells_df[["Clip_start"]] + add_vec
  tidied_wells_df[["Clip_end"]]       <- tidied_wells_df[["Clip_end"]] + add_vec
  tidied_wells_df[["Pre_clip_start"]] <- tidied_plates_df[["Clip_start"]][zmw_matches]
  tidied_wells_df[["Pre_clip_end"]]   <- tidied_plates_df[["Clip_end"]][zmw_matches]

  plate_df <- IntegrateSingleReportDf(sam_df, tidied_plates_df)
  well_df <- IntegrateSingleReportDf(sam_df, tidied_wells_df, check_zmws = FALSE)

  names(plate_df)[names(plate_df) == "Combo_number"] <- "number"
  names(well_df)[names(well_df) == "Combo_number"] <- "number"
  plate_df <- plate_df[, names(plate_df) != "Barcodes_orientation_fwd"]

  AddPrefix <- function(use_df, prefix) {
    are_new <- !(names(use_df) %in% names(sam_df))
    names(use_df)[are_new] <- paste0(prefix, tolower(names(use_df)[are_new]))
    return(use_df)
  }

  plate_df <- AddPrefix(plate_df, "Plate_")
  well_df <- AddPrefix(well_df, "Well_")

  results_df <- data.frame(
    sam_df,
    plate_df[, !(names(plate_df) %in% names(sam_df))],
    well_df[, !(names(well_df) %in% names(sam_df))],
    stringsAsFactors = FALSE
  )

  plate_number_width <- max(nchar(as.character(use_plates_df[["Plate_number"]])))
  combined_IDs <- paste0("Plate", formatC(results_df[["Plate_number"]], width = plate_number_width, flag = "0"),
                         "_Well",  formatC(results_df[["Well_number"]],  width = 3, flag = "0")
                         )
  combined_IDs <- ifelse(is.na(results_df[["Plate_number"]]) |
                         is.na(results_df[["Well_number"]]),
                         NA_character_,
                         combined_IDs
                         )
  results_df[["Combined_ID"]] <- combined_IDs

  common_columns <- c(
    "Clipped_read_length",
    "Passed_filters",
    "Barcode_combined_score", "Barcode_score_lead",
    "Clip_start", "Clip_end", "Barcodes_orientation_fwd",
    "Barcode_1_length",   "Barcode_2_length",
    "Barcode_1_sequence", "Barcode_2_sequence",
    "Barcode_1_quality",  "Barcode_2_quality"
  )
  all_columns <- c(
    "ZMW", "Combined_ID", "Plate_number", "Well_number",
    "Num_full_passes", "Read_quality",
    "Sequence", "Quality",
    paste0("Plate_", tolower(common_columns)),
    paste0("Well_", tolower(common_columns))
  )
  all_columns <- setdiff(all_columns, "Plate_barcodes_orientation_fwd")

  results_df <- results_df[, all_columns]

  are_valid <- is.na(results_df[["Well_clip_end"]]) |
               (results_df[["Well_clip_end"]] <= results_df[["Plate_clip_end"]])
  stopifnot(all(are_valid))

  new_order <- order(results_df[["Plate_number"]],
                     results_df[["Well_number"]]
                     )
  results_df <- results_df[new_order, ]
  row.names(results_df) <- NULL

  stopifnot(!(any(duplicated(results_df[["ZMW"]]))))
  return(results_df)
}




AddWellExistsColumn <- function(use_ccs_df, use_library_df) {
  use_ccs_df[["Well_exists"]] <- ifelse(is.na(use_ccs_df[["Combined_ID"]]),
                                        NA,
                                        use_ccs_df[["Combined_ID"]] %in% use_library_df[["Combined_ID"]]
                                        )

  are_preceding <- seq_len(ncol(use_ccs_df)) < which(names(use_ccs_df) == "Well_number")
  new_columns <- c(names(ccs_df)[are_preceding],
                   "Well_exists",
                   setdiff(names(use_ccs_df)[!(are_preceding)], "Well_exists")
                   )
  use_ccs_df <- use_ccs_df[, new_columns]
  return(use_ccs_df)
}



