### 15th October 2020 ###



# Import packages and source code -----------------------------------------

library("Rsamtools")

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "2) Functions for analyzing reads.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed - ccs3.RData"))





# Define functions --------------------------------------------------------

GetBarcodes <- function(use_ccs3 = TRUE, use_sl7 = TRUE) {

  barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

  if (use_ccs3) {
    if (use_sl7) {
      ccs_list <- sl7_ccs3_ccs
      lima_list <- sl7_ccs3_lima
      report_df <- sl7_ccs3_report_df
    } else {
      ccs_list <- sl9_ccs3_ccs
      lima_list <- sl9_ccs3_lima
      report_df <- sl9_ccs3_report_df
    }
  } else {
    if (use_sl7) {
      ccs_list <- sl7_ccs5_ccs
      lima_list <- sl7_ccs5_lima
      report_df <- sl7_ccs5_report_df
    } else {
      ccs_list <- sl9_ccs5_ccs
      lima_list <- sl9_ccs5_lima
      report_df <- sl9_ccs5_report_df
    }
  }

  lima_zmws <- as.integer(substr(lima_list[["qname"]], 22, nchar(lima_list[["qname"]]) - 4))
  ccs_zmws  <- as.integer(substr(ccs_list[["qname"]],  22, nchar(ccs_list[["qname"]])  - 4))

  ccs_well_numbers <- GetWellNumbers(report_df)
  lima_well_numbers <- ccs_well_numbers[match(lima_zmws, ccs_zmws)]

  stopifnot(identical(ccs_zmws, as.integer(substr(report_df[[1]], 22, nchar(report_df[[1]])))))

  well_barcodes_df_list <- lapply(seq_len(384), function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- lima_well_numbers %in% well_number
    this_well_zmws <- lima_zmws[are_this_well]

    lima_seq <- lima_list[["seq"]][are_this_well]
    lima_qual <- lima_list[["qual"]][are_this_well]

    ccs_matches <- match(this_well_zmws, ccs_zmws)
    stopifnot(!(anyNA(ccs_matches)))
    ccs_seq <- ccs_list[["seq"]][ccs_matches]
    ccs_qual <- ccs_list[["qual"]][ccs_matches]

    plasmid <- DNAStringSet(barcoded_plasmids[[well_number]])
    row_template <- row_bc_vec[[well_number]]
    column_template <- column_bc_vec[[well_number]]

    fwd_alignments <- pairwiseAlignment(ccs_seq, plasmid, type = "global")
    rev_alignments <- pairwiseAlignment(reverseComplement(ccs_seq), plasmid, type = "global")
    are_forward_vec <- score(fwd_alignments) > score(rev_alignments)

    barcode_df_list <- lapply(seq_along(this_well_zmws), function(x) {

      this_lima <- as.character(lima_seq[[x]])
      this_ccs <- as.character(ccs_seq[[x]])
      this_ccs_qual <- as.character(ccs_qual[[x]])

      match_pos <- as.integer(regexpr(this_lima, this_ccs, fixed = TRUE))
      lima_length <- nchar(this_lima)
      ccs_length <- nchar(this_ccs)
      stopifnot(match_pos != -1)

      first_bc_seq <- substr(this_ccs, 1, match_pos - 1)
      second_bc_seq <- substr(this_ccs, match_pos + lima_length, ccs_length)

      first_bc_qual <- substr(this_ccs_qual, 1, match_pos - 1)
      second_bc_qual <- substr(this_ccs_qual, match_pos + lima_length, ccs_length)

      stopifnot(identical(as.character(lima_qual[[x]]),
                          substr(this_ccs_qual, match_pos, match_pos + lima_length - 1)
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

  results_df <- do.call(rbind.data.frame, c(well_barcodes_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}



ProcessBarcodesDf <- function(barcodes_df, include_original_columns = T) {

  contains_mat <- matrix(nrow = nrow(barcodes_df), ncol = 6)
  colnames(contains_mat) <- c("Contains_row_barcode",
                             "Starts_with_row_barcode",
                             "Ends_with_row_barcode",
                             "Contains_column_barcode",
                             "Starts_with_column_barcode",
                             "Ends_with_column_barcode"
                             )
  for (i in seq_len(384)) {

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







# Extract barcodes --------------------------------------------------------

sl7_ccs3_barcodes_df <- GetBarcodes(use_sl7 = TRUE, use_ccs3 = TRUE)
sl7_ccs5_barcodes_df <- GetBarcodes(use_sl7 = TRUE, use_ccs3 = FALSE)
sl9_ccs3_barcodes_df <- GetBarcodes(use_sl7 = FALSE, use_ccs3 = TRUE)
sl9_ccs5_barcodes_df <- GetBarcodes(use_sl7 = FALSE, use_ccs3 = FALSE)

sl7_ccs3_barcodes_df <- ProcessBarcodesDf(sl7_ccs3_barcodes_df)
sl7_ccs5_barcodes_df <- ProcessBarcodesDf(sl7_ccs5_barcodes_df)
sl9_ccs3_barcodes_df <- ProcessBarcodesDf(sl9_ccs3_barcodes_df)
sl9_ccs5_barcodes_df <- ProcessBarcodesDf(sl9_ccs5_barcodes_df)






# Save data ---------------------------------------------------------------

save(list = paste0(c("sl7_ccs3", "sl7_ccs5", "sl9_ccs3", "sl9_ccs5"), "_barcodes_df"),
     file = file.path(R_objects_directory, "8) Extract barcode sequences and quality scores.RData")
     )








