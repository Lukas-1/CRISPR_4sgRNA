### 14th June 2021 ###



# Define functions --------------------------------------------------------

CharacterizeContaminations <- function(extract_df, sg_df) {

  ## Tidy guide RNA sequences
  guides_mat <- as.matrix(sg_df[, paste0("Sequence_sg", 1:4)])
  all_guides <- as.character(guides_mat)


  ## Create helper variables
  are_sg <- extract_df[["Feature"]] %in% paste0("sg", 1:4)
  are_sg_contam <- extract_df[["Is_contamination"]] & are_sg


  ## Check for cross-plate contaminations
  multi_plate <- "Plate_number" %in% names(sg_df)
  if (multi_plate) {
    are_cross_plate <- rep(NA, nrow(extract_df))
    for (plate_number in unique(sg_df[["Plate_number"]])) {
      are_eligible <- (extract_df[["Plate_number"]] == plate_number) & are_sg_contam
      this_plate_guides <- as.character(guides_mat[sg_df[["Plate_number"]] == plate_number, ])
      other_plate_guides <- setdiff(all_guides, this_plate_guides)
      sub_cross_plate <- extract_df[["Aligned_read"]][are_eligible] %in% other_plate_guides
      are_cross_plate[are_eligible] <- sub_cross_plate
    }
  }

  ## Create a helper data frame
  if (multi_plate) {
    ID_columns <- c("Combined_ID", "Plate_number", "Well_number")
  } else {
    ID_columns <- c("Well_number")
  }
  long_df_list <- lapply(1:4,
                         function(x) data.frame(sg_df[, ID_columns, drop = FALSE],
                                                "Sg_number"       = x,
                                                "Sequence"        = guides_mat[, x],
                                                "Num_occurrences" = sg_df[, paste0("Num_occurrences_sg", x)],
                                                stringsAsFactors  = FALSE,
                                                row.names         = NULL
                                                )
                         )
  long_df <- do.call(rbind.data.frame, c(long_df_list, stringsAsFactors = FALSE))
  if (multi_plate) {
    new_order <- order(long_df[, "Plate_number"], long_df[, "Well_number"])
  } else {
    new_order <- order(long_df[["Well_number"]])
  }
  long_df <- long_df[new_order, ]
  row.names(long_df) <- NULL


  ## Integrate data on the contaminating well
  contam_df_list <- lapply(unique(sg_df[["Plate_number"]]), function(x) {
    sub_order <- order(long_df[["Plate_number"]] == x,
                       decreasing = TRUE
                       )
    reordered_df <- long_df[sub_order, ]
    are_this_plate <- extract_df[["Plate_number"]][are_sg_contam] == x
    this_plate_vec <- extract_df[["Aligned_read"]][are_sg_contam][are_this_plate]
    matches_vec <- match(this_plate_vec, toupper(reordered_df[["Sequence"]]))
    results_df <- reordered_df[matches_vec, ]
    return(results_df)
  })
  contam_df <- do.call(rbind.data.frame,
                       c(contam_df_list,
                         stringsAsFactors = FALSE,
                         make.row.names = FALSE
                         )
                       )
  names(contam_df) <- paste0("Contaminating_",
                             tolower(substr(names(contam_df), 1, 1)),
                             substr(names(contam_df), 2, nchar(names(contam_df)))
                             )

  redundant_columns <- c("Num_missing", "Mostly_deleted", "Is_correct",
                         "Category", "Is_contamination", "Aligned_template",
                         "Alignment_length", "Over_5_percent_incorrect"
                         )
  contam_df <- data.frame(extract_df[are_sg_contam, !(names(extract_df) %in% redundant_columns)],
                          "Are_cross_plate" = if (multi_plate) are_cross_plate[are_sg_contam] else NA,
                          contam_df[, names(contam_df) != "Contaminating_sequence"],
                          stringsAsFactors = FALSE,
                          row.names = NULL
                          )


  ## Tidy the data frame on contaminations
  contam_df[["Reference_sg_number"]] <- as.integer(substr(contam_df[["Feature"]], 3, 3))
  contam_df <- contam_df[, names(contam_df) != "Feature"]

  rename_vec <- c(
                                      "ZMW",
                                      "Original_ZMW",
                                      "SmrtCell",
                                      "Are_cross_plate",
    "Combined_ID"                   = "Reference_ID",
    "Contaminating_combined_ID"     = "Contaminating_ID",
    "Plate_number"                  = "Reference_plate_number",
                                      "Contaminating_plate_number",
    "Well_number"                   = "Reference_well_number",
                                      "Contaminating_well_number",
                                      "Reference_sg_number",
                                      "Contaminating_sg_number",
    "Template"                      = "Reference_sequence",
    "Aligned_read"                  = "Contaminating_sequence",
    "Num_incorrect"                 = "Num_mismatched_bases",
                                      "Quality",
                                      "Mean_quality",
    "Contaminating_num_occurrences" = "Contaminating_sequence_occurrences_in_library"
  )
  if (!(multi_plate)) {
    multi_plate_columns <- c(
      "Reference_ID", "Contaminating_ID",
      "Reference_plate_number", "Contaminating_plate_number"
    )
    rename_vec <- rename_vec[!(rename_vec %in% multi_plate_columns)]
  }
  if (!("Original_ZMW" %in% names(contam_df))) {
    rename_vec <- rename_vec[!(rename_vec %in% c("Original_ZMW", "SmrtCell"))]
  }

  for (column_name in setdiff(names(rename_vec), "")) {
    names(contam_df)[names(contam_df) == column_name] <- rename_vec[[column_name]]
  }
  contam_df <- contam_df[, order(match(names(contam_df), rename_vec))]
  contam_df[["Mean_quality"]] <- contam_df[["Mean_quality"]] / 93

  if (!(multi_plate)) {
    contam_df <- contam_df[, names(contam_df) != "Are_cross_plate"]
  }

  return(contam_df)
}








