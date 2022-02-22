## 2022-01-04


# Load packages and source code -------------------------------------------

library("ShortRead")
library("Biostrings")
library("parallel")
library("data.table")



# Define paths ------------------------------------------------------------

input_dir <- file.path(project_dir, "02_input_data")
reads_dir <- file.path(input_dir, "raw_reads")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Read in data ------------------------------------------------------------

amplicon_ref <- read.table(file.path(input_dir, "amplicon_4sg.txt"),
                           quote = "", stringsAsFactors = FALSE
                           )[, 1]



# Detect number of cores --------------------------------------------------

number_of_cores <- parallel::detectCores() - 3L



# Define functions --------------------------------------------------------

AlignReads <- function(use_reference, use_sequences, opening_penalty = 30) {

  assign("delete_use_reference", use_reference, envir = globalenv())
  assign("delete_use_sequences", use_sequences, envir = globalenv())
  assign("delete_opening_penalty", opening_penalty, envir = globalenv())

  use_reference <- DNAStringSet(use_reference)
  if (!("DNAStringSet" %in% class(use_sequences))) {
    use_sequences <- DNAStringSet(use_sequences)
  }

  message("Computing forward alignments...")
  fwd_alignments <- pairwiseAlignment(use_sequences, use_reference,
                                      type = "global",
                                      gapOpening = opening_penalty
                                      )

  message("Computing reverse alignments...")
  rev_alignments <- pairwiseAlignment(reverseComplement(use_sequences), use_reference,
                                      type = "global",
                                      gapOpening = opening_penalty
                                      )


  message("Compiling results...")
  are_forward_vec <- score(fwd_alignments) > score(rev_alignments)

  new_order <- order(c(which(are_forward_vec), which(!(are_forward_vec))))

  aligned_plasmid_vec <- c(as.character(alignedSubject(fwd_alignments)[are_forward_vec]),
                           as.character(alignedSubject(rev_alignments)[!(are_forward_vec)])
                           )[new_order]
  aligned_read_vec <- c(as.character(alignedPattern(fwd_alignments)[are_forward_vec]),
                        as.character(alignedPattern(rev_alignments)[!(are_forward_vec)])
                        )[new_order]

  alignments_df <- data.frame("Read_number"     = seq_along(use_sequences),
                              "Orientation_fwd" = are_forward_vec,
                              "Score_fwd"       = score(fwd_alignments),
                              "Score_rev"       = score(rev_alignments),
                              "Aligned_ref"     = aligned_plasmid_vec,
                              "Aligned_read"    = aligned_read_vec,
                              stringsAsFactors  = FALSE
                              )
  old_order <- order(alignments_df[, "Read_number"])
  alignments_df <- alignments_df[old_order, names(alignments_df) != "Read_number"]
  row.names(alignments_df) <- NULL
  return(alignments_df)
}



ParallelAlign <- function(use_reference, use_sequences, opening_penalty = 30, num_cores = NULL) {

  ### See https://github.com/NathanSkene/EWCE/issues/5#issuecomment-497616095
  if (is.null(num_cores)) {
    num_cores <- parallel::detectCores() - 3
  }

  num_reads <- length(use_sequences)
  reads_per_chunk <- 500
  num_chunks <- ceiling(num_reads / reads_per_chunk)
  chunks_vec <- rep(seq_len(num_chunks), each = reads_per_chunk)[seq_len(num_reads)]

  use_seq_list <- lapply(seq_len(num_chunks), function(x) {
    are_this_chunk <- x == chunks_vec
    use_sequences[are_this_chunk]
  })

  cl <- parallel::makeCluster(num_cores, setup_strategy = "sequential")
  parallel::clusterExport(cl,
                          varlist = c("DNAStringSet", "pairwiseAlignment",
                                      "reverseComplement", "score",
                                      "alignedSubject", "alignedPattern",
                                      "AlignReads", "use_reference", "opening_penalty"
                          ),
                          envir = environment()
                          )
  alignments_df_list <- parallel::parLapply(cl,
                                            use_seq_list,
                                            function(x) AlignReads(use_reference,
                                                                   x,
                                                                   opening_penalty = opening_penalty
                                                                   )
                                            )
  parallel::stopCluster(cl)

  message("Combining results into a data frame...")
  results_df <- data.table::rbindlist(alignments_df_list)
  data.table::setDF(results_df)
  return(results_df)
}



ParallelAlignInChunks <- function(all_reads) {
  ## 'all_reads' is a ShortReadQ object
  num_reads <- length(all_reads)
  reads_per_chunk <- 250000
  num_chunks <- ceiling(num_reads / reads_per_chunk)
  chunks_vec <- rep(seq_len(num_chunks), each = reads_per_chunk)[seq_len(num_reads)]
  chunks_list <- vector(mode = "list", length = num_chunks)
  first_vec <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[1]]))
  last_vec  <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[length(x)]]))
  chunk_numbers <- format(seq_len(num_chunks))

  for (i in seq_len(num_chunks)) {
    are_this_chunk <- chunks_vec == i
    message("Processing chunk #", chunk_numbers[[i]], " of ",
            chunk_numbers[[length(chunk_numbers)]],  " (aligning reads ",
            first_vec[[i]], " to ", last_vec[[i]], ")..."
            )
    sub_df <- ParallelAlign(amplicon_ref,
                            ShortRead::sread(all_reads)[are_this_chunk],
                            num_cores = number_of_cores # See https://github.com/NathanSkene/EWCE/issues/5#issuecomment-497616095
                            )
    chunks_list[[i]] <- sub_df
  }

  alignments_df <- data.table::rbindlist(chunks_list)
  data.table::setDF(alignments_df)

  alignments_df <- data.frame(
    "Read_sequence" = as.character(ShortRead::sread(all_reads)),
    "Read_quality" = as.character(as(Biostrings::quality(all_reads), "PhredQuality")),
    alignments_df,
    stringsAsFactors = FALSE
  )
  return(alignments_df)
}



