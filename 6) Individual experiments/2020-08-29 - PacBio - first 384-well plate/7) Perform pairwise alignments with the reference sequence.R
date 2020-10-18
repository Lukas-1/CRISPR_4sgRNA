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

PairWiseAlignments <- function(use_sl7 = TRUE) {

  barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

  if (use_sl7) {
    ccs_list <- sl7_ccs3_ccs
    report_df <- sl7_ccs3_report_df
  } else {
    ccs_list <- sl9_ccs3_ccs
    report_df <- sl9_ccs3_report_df
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





ExtractAlignedSequences <- function(use_ccs3 = TRUE, use_sl7 = TRUE) {

  barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

  if (use_sl7) {
    ccs_list <- sl7_ccs3_ccs
    report_df <- sl7_ccs3_report_df
  } else {
    ccs_list <- sl9_ccs3_ccs
    report_df <- sl9_ccs3_report_df
  }

  ccs_zmws <- as.integer(substr(ccs_list[["qname"]],  22, nchar(ccs_list[["qname"]])  - 4))
  ccs_well_numbers <- GetWellNumbers(report_df)

  stopifnot(identical(ccs_zmws, as.integer(substr(report_df[[1]], 22, nchar(report_df[[1]])))))

  alignments_list <- lapply(seq_len(384), function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- ccs_well_numbers %in% well_number
    ccs_seq <- ccs_list[["seq"]][are_this_well]
    ccs_qual <- ccs_list[["qual"]][are_this_well]

    plasmid <- DNAStringSet(barcoded_plasmids[[well_number]])
    row_template <- row_bc_vec[[well_number]]
    column_template <- column_bc_vec[[well_number]]

    fwd_alignments <- pairwiseAlignment(ccs_seq, plasmid, type = "global")
    rev_alignments <- pairwiseAlignment(reverseComplement(ccs_seq), plasmid, type = "global")
    are_forward_vec <- score(fwd_alignments) > score(rev_alignments)

    meta_df <- data.frame("ZMW"             = ccs_zmws[are_this_well],
                          "Orientation_fwd" = are_forward_vec,
                          "Score_fwd"       = score(fwd_alignments),
                          "Score_rev"       = score(rev_alignments),
                          stringsAsFactors = FALSE
                          )

    assign("delete_fwd_alignments", fwd_alignments, envir = globalenv())
    assign("delete_rev_alignments", rev_alignments, envir = globalenv())
    assign("delete_are_forward_vec", are_forward_vec, envir = globalenv())
    combined_alignments <- do.call(c, lapply(seq_along(are_forward_vec),
                                             function(x) {
                                               if (are_forward_vec[[x]]) {
                                                 fwd_alignments[x]
                                               } else {
                                                 rev_alignments[x]
                                               }
                                             })
                                   )
    results_list <- list(
      "meta_df" = meta_df,
      "alignments" = combined_alignments
    )
    return(results_list)
  })
  return(alignments_list)
}




# Extract barcodes --------------------------------------------------------

sl7_alignments_list <- ExtractAlignedSequences(use_sl7 = TRUE)
sl9_alignments_list <- ExtractAlignedSequences(use_sl7 = FALSE)




# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_alignments_list"),
     file = file.path(R_objects_directory, "7) Perform pairwise alignments with the reference sequence.RData")
     )








