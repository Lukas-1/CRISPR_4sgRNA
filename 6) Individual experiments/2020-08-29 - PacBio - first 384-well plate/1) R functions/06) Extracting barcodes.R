### 20th October 2020 ###




# Import packages and source code -----------------------------------------

library("Biostrings")



# Define functions --------------------------------------------------------

GetBarcodes <- function(ccs_df, alignments_df, wells_vec = seq_len(384)) {

  stopifnot(all(c("row_bc_vec", "column_bc_vec", "row_constant_region", "column_constant_region") %in% ls(envir = globalenv())))

  matches_vec <- match(alignments_df[["ZMW"]], ccs_df[["ZMW"]])

  stopifnot(!(anyNA(matches_vec)))

  are_forward <- alignments_df[["Orientation_fwd"]]
  ccs_seq <- ccs_df[["Sequence"]][matches_vec]

  bc_seq_1 <- ccs_df[["Barcode_1_sequence"]][matches_vec]
  bc_seq_2 <- ccs_df[["Barcode_2_sequence"]][matches_vec]

  bc_qual_1 <- ccs_df[["Barcode_1_quality"]][matches_vec]
  bc_qual_2 <- ccs_df[["Barcode_2_quality"]][matches_vec]

  column_seq <- ifelse(are_forward, bc_seq_1, NA_character_)
  column_seq[!(are_forward)] <- as.character(reverseComplement(DNAStringSet(bc_seq_2[!(are_forward)])))
  row_seq <- ifelse(are_forward, bc_seq_2, NA_character_)
  row_seq[!(are_forward)] <- as.character(reverseComplement(DNAStringSet(bc_seq_1[!(are_forward)])))

  column_qual <- ifelse(are_forward, bc_qual_1, NA_character_)
  column_qual[!(are_forward)] <- reverse(bc_qual_2[!(are_forward)])
  row_qual <- ifelse(are_forward, bc_qual_2, NA_character_)
  row_qual[!(are_forward)] <- reverse(bc_qual_1[!(are_forward)])

  flanking_1 <- substr(ccs_seq,
                       ccs_df[["Barcode_1_length"]][matches_vec] + 1,
                       ccs_df[["Barcode_1_length"]][matches_vec] + 20
                       )
  flanking_2 <- substr(ccs_seq,
                       nchar(ccs_seq) - ccs_df[["Barcode_2_length"]][matches_vec] - 19,
                       nchar(ccs_seq) - ccs_df[["Barcode_2_length"]][matches_vec]
                       )

  column_is_flanked <- (flanking_1 == column_constant_region) |
                       (flanking_2 == as.character(reverseComplement(DNAString(column_constant_region))))

  row_is_flanked <- (flanking_1 == row_constant_region) |
                    (flanking_2 == as.character(reverseComplement(DNAString(row_constant_region))))

  row_templates_vec <- row_bc_vec[ccs_df[["Well_number"]][matches_vec]]
  column_templates_vec <- column_bc_vec[ccs_df[["Well_number"]][matches_vec]]

  results_df <- data.frame(
    alignments_df[, c("ZMW", "Well_number")],
    "Match_template_row"    = row_seq == row_templates_vec,
    "Match_template_column" = column_seq == column_templates_vec,
    "Match_flanking_row"    = row_is_flanked,
    "Match_flanking_column" = column_is_flanked,
    "Is_forward"            = are_forward,
    "Row_barcode"           = row_seq,
    "Row_quality"           = row_qual,
    "Column_barcode"        = column_seq,
    "Column_quality"        = column_qual,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  results_df <- ProcessBarcodesDf(results_df, wells_vec = wells_vec)
  return(results_df)
}





ProcessBarcodesDf <- function(barcodes_df, wells_vec = seq_len(384)) {

  stopifnot(all(c("row_bc_vec", "column_bc_vec") %in% ls(envir = globalenv())))

  contains_mat <- matrix(nrow = nrow(barcodes_df), ncol = 6)
  colnames(contains_mat) <- c("Contains_row_barcode",
                              "Starts_with_row_barcode",
                              "Ends_with_row_barcode",
                              "Contains_column_barcode",
                              "Starts_with_column_barcode",
                              "Ends_with_column_barcode"
                              )
  for (i in wells_vec) {
    are_this_well <- barcodes_df[["Well_number"]] %in% i
    row_bc <- row_bc_vec[[i]]
    column_bc <- column_bc_vec[[i]]

    contains_mat[are_this_well, "Contains_row_barcode"]       <- grepl(row_bc, barcodes_df[["Row_barcode"]][are_this_well], fixed = TRUE)
    contains_mat[are_this_well, "Starts_with_row_barcode"]    <- startsWith(barcodes_df[["Row_barcode"]][are_this_well], row_bc)
    contains_mat[are_this_well, "Ends_with_row_barcode"]      <- endsWith(barcodes_df[["Row_barcode"]][are_this_well], row_bc)
    contains_mat[are_this_well, "Contains_column_barcode"]    <- grepl(column_bc, barcodes_df[["Column_barcode"]][are_this_well], fixed = TRUE)
    contains_mat[are_this_well, "Starts_with_column_barcode"] <- startsWith(barcodes_df[["Column_barcode"]][are_this_well], column_bc)
    contains_mat[are_this_well, "Ends_with_column_barcode"]   <- endsWith(barcodes_df[["Column_barcode"]][are_this_well], column_bc)
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
    barcodes_df[, c("Well_number", "ZMW")],
    "Correct_barcodes"     = as.integer(barcodes_df[["Match_template_row"]] & barcodes_df[["Match_template_column"]]),
    "Correct_row"          = as.integer(barcodes_df[["Match_template_row"]]),
    "Correct_column"       = as.integer(barcodes_df[["Match_template_column"]]),
    "Correct_flanks"       = as.integer(barcodes_df[["Match_flanking_row"]] & barcodes_df[["Match_flanking_column"]]),
    "Correct_row_flank"    = as.integer(barcodes_df[["Match_flanking_row"]]),
    "Correct_column_flank" = as.integer(barcodes_df[["Match_flanking_column"]]),
    "Row_bc_length"        = nchar(barcodes_df[["Row_barcode"]]),
    "Column_bc_length"     = nchar(barcodes_df[["Column_barcode"]]),
    "Row_mean_quality"     = GetMeanQuality(barcodes_df[["Row_quality"]]),
    "Column_mean_quality"  = GetMeanQuality(barcodes_df[["Column_quality"]]),
    barcodes_df[, c("Row_barcode", "Column_barcode", "Row_quality", "Column_quality")],
    "Orientation_fwd"      = as.integer(barcodes_df[["Is_forward"]]),
    contains_mat,
    stringsAsFactors = FALSE
  )
  return(results_df)
}



