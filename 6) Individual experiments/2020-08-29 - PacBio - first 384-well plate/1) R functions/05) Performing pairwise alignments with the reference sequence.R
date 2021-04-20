### 20th October 2020 ###



# Import packages and source code -----------------------------------------

library("Biostrings")



# Define functions --------------------------------------------------------


ProcessWell <- function(ccs_df,
                        unique_IDs,
                        ID_index,
                        ID_column = "Well_number"
                        ) {

  current_ID <- unique_IDs[[ID_index]]

  message(paste0("Processing the well: ", current_ID, "..."))

  are_this_ID <- ccs_df[[ID_column]] %in% current_ID
  ccs_seq <- DNAStringSet(ccs_df[["Sequence"]][are_this_ID])

  if (is.null(names(barcoded_plasmids))) {
    use_index <- ID_index
  } else {
    use_index <- current_ID
  }
  plasmid <- DNAStringSet(barcoded_plasmids[[use_index]])

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
                                    ID_column = "Well_number",
                                    unique_IDs = seq_len(384),
                                    parallel_mode = FALSE,
                                    num_cores = NULL
                                    ) {

  stopifnot("barcoded_plasmids" %in% ls(envir = globalenv()))

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
                                        "ProcessWell", "barcoded_plasmids",
                                        "ccs_df", "unique_IDs", "ID_column"
                                        ),
                            envir = environment()
                            )
    alignments_df_list <- parallel::parLapply(cl,
                                              seq_along(unique_IDs),
                                              function(x) ProcessWell(ccs_df,
                                                                      unique_IDs,
                                                                      x,
                                                                      ID_column
                                                                      )
                                              )
    parallel::stopCluster(cl)
  } else {
    alignments_df_list <- lapply(seq_along(unique_IDs), function(x) {
      ProcessWell(ccs_df, unique_IDs, x, ID_column)
    })
  }

  message("Combining results into a data frame...")
  results_df <- data.table::rbindlist(alignments_df_list)
  data.table::setDF(results_df)
  return(results_df)
}

