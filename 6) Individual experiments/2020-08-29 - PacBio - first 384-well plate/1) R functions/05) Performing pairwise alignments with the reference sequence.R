### 20th October 2020 ###



# Import packages and source code -----------------------------------------

library("Biostrings")
library("data.table")


# Define functions --------------------------------------------------------

ProcessWell <- function(ccs_df,
                        sg_df,
                        ID_index,
                        ID_column = "Well_number",
                        opening_penalty = 30
                        ) {

  current_ID <- sg_df[ID_index, ID_column]

  message(paste0("Processing the well: ", current_ID, "..."))

  are_this_ID <- ccs_df[[ID_column]] %in% current_ID

  if (!(any(are_this_ID))) {
    message("No reads are available for this well ID! It was skipped.")
    return(NULL)
  }

  ccs_seq <- DNAStringSet(ccs_df[["Sequence"]][are_this_ID])
  plasmid <- DNAStringSet(sg_df[ID_index, "Barcoded_plasmid"])

  fwd_alignments <- pairwiseAlignment(ccs_seq, plasmid,type = "global",
                                      gapOpening = opening_penalty
                                      )
  rev_alignments <- pairwiseAlignment(reverseComplement(ccs_seq), plasmid,
                                      type = "global",
                                      gapOpening = opening_penalty
                                      )
  are_forward_vec <- score(fwd_alignments) > score(rev_alignments)

  new_order <- order(c(which(are_forward_vec), which(!(are_forward_vec))))

  aligned_plasmid_vec <- c(as.character(alignedSubject(fwd_alignments)[are_forward_vec]),
                           as.character(alignedSubject(rev_alignments)[!(are_forward_vec)])
                           )[new_order]
  aligned_read_vec <- c(as.character(alignedPattern(fwd_alignments)[are_forward_vec]),
                        as.character(alignedPattern(rev_alignments)[!(are_forward_vec)])
                        )[new_order]

  alignments_df <- data.frame("Source_ID"       = current_ID,
                              "ZMW"             = ccs_df[["ZMW"]][are_this_ID],
                              "Orientation_fwd" = are_forward_vec,
                              "Score_fwd"       = score(fwd_alignments),
                              "Score_rev"       = score(rev_alignments),
                              "Aligned_plasmid" = aligned_plasmid_vec,
                              "Aligned_read"    = aligned_read_vec,
                              stringsAsFactors  = FALSE
                              )
  return(alignments_df)
}




ExtractAlignedSequences <- function(ccs_df,
                                    sg_df,
                                    ID_column = "Well_number",
                                    parallel_mode = FALSE,
                                    num_cores = NULL
                                    ) {

  if (parallel_mode) {
    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores() - 2
    }
    ccs_df <- ccs_df[, c("ZMW", "Sequence", ID_column)]

    ccs_df <- ccs_df[!(is.na(ccs_df[[ID_column]])), ]
    row.names(ccs_df) <- NULL
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterExport(cl,
                            varlist = c("DNAStringSet", "pairwiseAlignment",
                                        "reverseComplement", "score",
                                        "alignedSubject", "alignedPattern",
                                        "ProcessWell",
                                        "sg_df", "ccs_df", "ID_column"
                                        ),
                            envir = environment()
                            )
    alignments_df_list <- parallel::parLapply(cl,
                                              seq_len(nrow(sg_df)),
                                              function(x) ProcessWell(ccs_df,
                                                                      sg_df,
                                                                      x,
                                                                      ID_column
                                                                      )
                                              )
    parallel::stopCluster(cl)
  } else {
    alignments_df_list <- lapply(seq_len(nrow(sg_df)), function(x) {
      ProcessWell(ccs_df, sg_df, x, ID_column)
    })
  }

  message("Combining results into a data frame...")
  results_df <- data.table::rbindlist(alignments_df_list)
  data.table::setDF(results_df)

  names(results_df)[names(results_df) == "Source_ID"] <- ID_column
  return(results_df)
}

