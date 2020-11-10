### 20th October 2020 ###



# Import packages and source code -----------------------------------------

library("Biostrings")



# Define functions --------------------------------------------------------

ExtractAlignedSequences <- function(use_sl7 = TRUE, wells_vec = seq_len(384)) {

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

  alignments_df_list <- lapply(seq_along(wells_vec), function(x) {

    well_number <- wells_vec[[x]]

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- ccs_well_numbers %in% well_number
    ccs_seq <- ccs_list[["seq"]][are_this_well]

    plasmid <- DNAStringSet(barcoded_plasmids[[x]])
    row_template <- row_bc_vec[[well_number]]
    column_template <- column_bc_vec[[well_number]]

    fwd_alignments <- pairwiseAlignment(ccs_seq, plasmid, type = "global")
    rev_alignments <- pairwiseAlignment(reverseComplement(ccs_seq), plasmid, type = "global")
    are_forward_vec <- score(fwd_alignments) > score(rev_alignments)

    new_order <- order(c(which(are_forward_vec), which(!(are_forward_vec))))

    aligned_plasmid_vec <- c(as.character(alignedSubject(fwd_alignments)[are_forward_vec]),
                             as.character(alignedSubject(rev_alignments)[!(are_forward_vec)])
                             )[new_order]
    aligned_read_vec <- c(as.character(alignedPattern(fwd_alignments)[are_forward_vec]),
                          as.character(alignedPattern(rev_alignments)[!(are_forward_vec)])
                          )[new_order]

    alignments_df <- data.frame("Well_number"     = well_number,
                                "ZMW"             = ccs_zmws[are_this_well],
                                "Orientation_fwd" = are_forward_vec,
                                "Score_fwd"       = score(fwd_alignments),
                                "Score_rev"       = score(rev_alignments),
                                "Aligned_plasmid" = aligned_plasmid_vec,
                                "Aligned_read"    = aligned_read_vec,
                                stringsAsFactors = FALSE
                                )
    return(alignments_df)
  })
  results_df <- do.call(rbind.data.frame, c(alignments_df_list,
                                            list(stringsAsFactors = FALSE,
                                                 make.row.names = FALSE
                                                 )
                                            )
                        )
  return(results_df)
}

