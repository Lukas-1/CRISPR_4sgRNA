### 20th October 2020 ###





# Import packages and source code -----------------------------------------

library("Biostrings")



# Define functions --------------------------------------------------------

GetBarcodes <- function(use_sl7 = TRUE, wells_vec = seq_len(384)) {

  barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

  if (use_sl7) {
    ccs_list <- sl7_ccs3_ccs
    lima_list <- sl7_ccs3_lima
    report_df <- sl7_ccs3_report_df
    alignments_df <- sl7_alignments_df
  } else {
    ccs_list <- sl9_ccs3_ccs
    lima_list <- sl9_ccs3_lima
    report_df <- sl9_ccs3_report_df
    alignments_df <- sl9_alignments_df
  }

  lima_zmws <- as.integer(substr(lima_list[["qname"]], 22, nchar(lima_list[["qname"]]) - 4))
  ccs_zmws  <- as.integer(substr(ccs_list[["qname"]],  22, nchar(ccs_list[["qname"]])  - 4))

  ccs_well_numbers <- GetWellNumbers(report_df)
  lima_well_numbers <- ccs_well_numbers[match(lima_zmws, ccs_zmws)]

  lima_seq_vec <- as.character(lima_list[["seq"]])
  lima_qual_vec <- as.character(lima_list[["qual"]])

  ccs_seq_vec <- as.character(ccs_list[["seq"]])
  ccs_qual_vec <- as.character(ccs_list[["qual"]])

  stopifnot(identical(ccs_zmws, as.integer(substr(report_df[[1]], 22, nchar(report_df[[1]])))))

  alignment_matches <- match(ccs_zmws, alignments_df[["ZMW"]])
  are_fwd_reordered <- alignments_df[["Orientation_fwd"]][alignment_matches]

  well_barcodes_df_list <- lapply(wells_vec, function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- lima_well_numbers %in% well_number
    this_well_zmws <- lima_zmws[are_this_well]

    lima_seq <- lima_seq_vec[are_this_well]
    lima_qual <- lima_qual_vec[are_this_well]

    ccs_matches <- match(this_well_zmws, ccs_zmws)
    stopifnot(!(anyNA(ccs_matches)))
    ccs_seq <- ccs_seq_vec[ccs_matches]
    ccs_qual <- ccs_qual_vec[ccs_matches]

    are_forward_vec <- are_fwd_reordered[ccs_matches]
    stopifnot(!(anyNA(are_forward_vec)))

    row_template <- row_bc_vec[[well_number]]
    column_template <- column_bc_vec[[well_number]]

    barcode_df_list <- lapply(seq_along(this_well_zmws), function(x) {

      match_pos <- as.integer(regexpr(lima_seq[[x]], ccs_seq[[x]], fixed = TRUE))
      lima_length <- nchar(lima_seq[[x]])
      ccs_length <- nchar(ccs_seq[[x]])
      stopifnot(match_pos != -1)

      first_bc_seq <- substr(ccs_seq[[x]], 1, match_pos - 1)
      second_bc_seq <- substr(ccs_seq[[x]], match_pos + lima_length, ccs_length)

      first_bc_qual <- substr(ccs_qual[[x]], 1, match_pos - 1)
      second_bc_qual <- substr(ccs_qual[[x]], match_pos + lima_length, ccs_length)

      stopifnot(identical(as.character(lima_qual[[x]]),
                          substr(ccs_qual[[x]], match_pos, match_pos + lima_length - 1)
                          )
                ) ## DELETE THIS LATER

      if (are_forward_vec[[x]]) {
        column_barcode <- first_bc_seq
        column_quality <- first_bc_qual
        row_barcode <- second_bc_seq
        row_quality <- second_bc_qual
      } else {
        column_barcode <- as.character(reverseComplement(DNAString(second_bc_seq)))
        column_quality <- reverse(second_bc_qual)
        row_barcode <- as.character(reverseComplement(DNAString(first_bc_seq)))
        row_quality <- reverse(first_bc_qual)
      }

      match_template <- (row_barcode == row_template) && (column_barcode == column_template)

      result_list <- list(
        "ZMW"            = this_well_zmws[[x]],
        "Match_template" = match_template,
        "Is_forward"     = are_forward_vec[[x]],
        "Row_barcode"    = row_barcode,
        "Row_quality"    = row_quality,
        "Column_barcode" = column_barcode,
        "Column_quality" = column_quality
      )
      return(result_list)
    })
    barcode_well_df <- do.call(rbind.data.frame, c(barcode_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
    barcode_well_df <- data.frame("Well" = well_number, barcode_well_df, stringsAsFactors = FALSE)
    return(barcode_well_df)
  })

  message("Collating the final data frame...")
  results_df <- do.call(rbind.data.frame, c(well_barcodes_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  results_df <- ProcessBarcodesDf(results_df)
  return(results_df)
}



ProcessBarcodesDf <- function(barcodes_df, wells_vec = seq_len(384)) {

  contains_mat <- matrix(nrow = nrow(barcodes_df), ncol = 6)
  colnames(contains_mat) <- c("Contains_row_barcode",
                             "Starts_with_row_barcode",
                             "Ends_with_row_barcode",
                             "Contains_column_barcode",
                             "Starts_with_column_barcode",
                             "Ends_with_column_barcode"
                             )
  for (i in wells_vec) {

    are_this_well <- barcodes_df[["Well"]] %in% i
    row_bc <- row_bc_vec[[i]]
    column_bc <- column_bc_vec[[i]]

    contains_mat[are_this_well, "Contains_row_barcode"]       <- grepl(row_bc, barcodes_df[["Row_barcode"]][are_this_well], fixed = TRUE)
    contains_mat[are_this_well, "Starts_with_row_barcode"]    <- grepl(paste0("^", row_bc), barcodes_df[["Row_barcode"]][are_this_well])
    contains_mat[are_this_well, "Ends_with_row_barcode"]      <- grepl(paste0(row_bc, "$"), barcodes_df[["Row_barcode"]][are_this_well])
    contains_mat[are_this_well, "Contains_column_barcode"]    <- grepl(column_bc, barcodes_df[["Column_barcode"]][are_this_well], fixed = TRUE)
    contains_mat[are_this_well, "Starts_with_column_barcode"] <- grepl(paste0("^", column_bc), barcodes_df[["Column_barcode"]][are_this_well])
    contains_mat[are_this_well, "Ends_with_column_barcode"]   <- grepl(paste0(column_bc, "$"), barcodes_df[["Column_barcode"]][are_this_well])
  }

  unexpected_col_bc <- (contains_mat[, "Starts_with_column_barcode"] |
                        contains_mat[, "Contains_column_barcode"]) &
                       !(contains_mat[, "Ends_with_column_barcode"])

  unexpected_row_bc <- (contains_mat[, "Ends_with_row_barcode"] |
                        contains_mat[, "Contains_row_barcode"]) &
                       !(contains_mat[, "Starts_with_row_barcode"])

  if (any(unexpected_col_bc)) {
    stop("Unexpected column barcode!")
  }
  if (any(unexpected_row_bc)) {
    stop("Unexpected row barcode!")
  }
  redundant_columns <- c("Starts_with_column_barcode", "Ends_with_row_barcode",
                         "Contains_column_barcode", "Contains_row_barcode"
                         )
  contains_mat <- contains_mat[, !(colnames(contains_mat) %in% redundant_columns)]

  mode(contains_mat) <- "integer"

  results_df <- data.frame(
    barcodes_df[, c("Well", "ZMW")],
    "Correct_barcodes"    = as.integer(barcodes_df[["Match_template"]]),
    "Row_bc_length"       = nchar(barcodes_df[["Row_barcode"]]),
    "Column_bc_length"    = nchar(barcodes_df[["Column_barcode"]]),
    "Row_mean_quality"    = GetMeanQuality(barcodes_df[["Row_quality"]]),
    "Column_mean_quality" = GetMeanQuality(barcodes_df[["Column_quality"]]),
    barcodes_df[, c("Row_barcode", "Column_barcode", "Row_quality", "Column_quality")],
    "Orientation_fwd"     = as.integer(barcodes_df[["Is_forward"]]),
    contains_mat
  )
  return(results_df)
}



