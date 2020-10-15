### 11th October 2020 ###



# Import packages and source code -----------------------------------------

library("Rsamtools")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
file_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory      <- file.path(file_directory, "2) Input")
R_objects_directory       <- file.path(file_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed.RData"))
load(file.path(R_objects_directory, "7) Process demultiplexed PacBio reads.RData"))





# Define functions --------------------------------------------------------

GetWellNumbers <- function(lima_report_df) {
  combo_IDs_vec <- paste0(lima_report_df[["IdxLowestNamed"]], "--",
                          lima_report_df[["IdxHighestNamed"]]
                          )
  barcodes_to_wells_map[combo_IDs_vec]
}



GetBarcodes <- function(use_ccs3 = TRUE, use_sl7 = TRUE) {

  barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

  if (use_ccs3) {
    if (use_sl7) {
      ccs_list <- sl7_ccs3_ccs
      lima_list <- sl7_ccs3_lima
      report_df <- sl7_ccs3_report_df
      output_folder <- "SmrtLink7_CCS3"
    } else {
      ccs_list <- sl9_ccs3_ccs
      lima_list <- sl9_ccs3_lima
      report_df <- sl9_ccs3_report_df
      output_folder <- "SmrtLink7_CCS5"
    }
  } else {
    if (use_sl7) {
      ccs_list <- sl7_ccs5_ccs
      lima_list <- sl7_ccs5_lima
      report_df <- sl7_ccs5_report_df
      output_folder <- "SmrtLink9_CCS3"
    } else {
      ccs_list <- sl9_ccs5_ccs
      lima_list <- sl9_ccs5_lima
      report_df <- sl9_ccs5_report_df
      output_folder <- "SmrtLink9_CCS3"
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

      match_object <- matchPattern(lima_seq[[x]], ccs_seq[[x]])
      stopifnot(length(match_object) == 1)

      first_bc_seq <- substr(as.character(ccs_seq[[x]]), 1, start(match_object) - 1)
      second_bc_seq <- substr(as.character(ccs_seq[[x]]), end(match_object) + 1, nchar(ccs_seq[[x]]))

      first_bc_qual <- substr(as.character(ccs_qual[[x]]), 1, start(match_object) - 1)
      second_bc_qual <- substr(as.character(ccs_qual[[x]]), end(match_object) + 1, nchar(ccs_qual[[x]]))

      stopifnot(identical(as.character(lima_qual[[x]]),
                          substr(as.character(ccs_qual[[x]]), start(match_object), end(match_object))
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





# Extract barcodes --------------------------------------------------------

sl7_ccs3_barcodes_df <- GetBarcodes(use_sl7 = TRUE, use_ccs3 = TRUE)
sl7_ccs5_barcodes_df <- GetBarcodes(use_sl7 = TRUE, use_ccs3 = FALSE)
sl9_ccs3_barcodes_df <- GetBarcodes(use_sl7 = FALSE, use_ccs3 = TRUE)
sl9_ccs5_barcodes_df <- GetBarcodes(use_sl7 = FALSE, use_ccs3 = FALSE)





# Save data ---------------------------------------------------------------

save(list = paste0(c("sl7_ccs3", "sl7_ccs5", "sl9_ccs3", "sl9_ccs5"), "_barcodes_df"),
     file = file.path(R_objects_directory, "7) Extract barcode sequences and quality scores.RData")
     )











