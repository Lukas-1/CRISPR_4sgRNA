## 2023-10-05



# Define functions --------------------------------------------------------

DemultiplexReads <- function(fastq_yield) {

  required_objects <- c(
    "fwd_replicate_barcodes", "fwd_condition_barcodes", "condition_names",
    "include_flanking"
  )
  stopifnot(all(required_objects %in% ls(envir = globalenv())))
  rev_replicate_barcodes <- as.character(reverseComplement(DNAStringSet(fwd_replicate_barcodes)))
  rev_condition_barcodes <- as.character(reverseComplement(DNAStringSet(fwd_condition_barcodes)))

  use_reads <- ShortRead::sread(fastq_yield)
  use_quality <- Biostrings::quality(fastq_yield)

  message("Searching for barcodes in the forward direction...")
  replicate_fwd_mat <- t(vcountPDict(DNAStringSet(fwd_replicate_barcodes), use_reads, max.mismatch = 1))
  condition_fwd_mat <- t(vcountPDict(DNAStringSet(fwd_condition_barcodes), use_reads, max.mismatch = 1))

  message("Searching for barcodes in the reverse direction...")
  replicate_rev_mat <- t(vcountPDict(DNAStringSet(rev_replicate_barcodes), use_reads, max.mismatch = 1))
  condition_rev_mat <- t(vcountPDict(DNAStringSet(rev_condition_barcodes), use_reads, max.mismatch = 1))


  message("Identifying read orientation...")
  replicate_fwd_mat <- replicate_fwd_mat == 1
  condition_fwd_mat <- condition_fwd_mat == 1
  replicate_rev_mat <- replicate_rev_mat == 1
  condition_rev_mat <- condition_rev_mat == 1

  are_fwd <- (rowSums(replicate_fwd_mat) == 1) & (rowSums(condition_fwd_mat) == 1)
  are_rev <- (rowSums(replicate_rev_mat) == 1) & (rowSums(condition_rev_mat) == 1)

  are_fwd <- are_fwd & (rowSums(replicate_rev_mat) == 0) & (rowSums(condition_rev_mat) == 0)
  are_rev <- are_rev & (rowSums(replicate_fwd_mat) == 0) & (rowSums(condition_fwd_mat) == 0)

  are_matched <- are_fwd | are_rev

  num_matched <- sum(are_matched)
  orientation_vec <- rep(NA, num_matched)
  orientation_vec <- are_fwd[are_matched]


  message("Assigning reads to samples...")
  replicate_vec <- rep(NA, num_matched)
  condition_vec <- rep(NA, num_matched)

  replicate_vec[orientation_vec] <- apply(replicate_fwd_mat[are_matched, ][orientation_vec, ], 1, which)
  replicate_vec[!(orientation_vec)] <- apply(replicate_rev_mat[are_matched, ][!(orientation_vec), ], 1, which)

  condition_vec[orientation_vec] <- apply(condition_fwd_mat[are_matched, ][orientation_vec, ], 1, which)
  condition_vec[!(orientation_vec)] <- apply(condition_rev_mat[are_matched, ][!(orientation_vec), ], 1, which)


  message("Standardize read orientation...")
  valid_seq <- use_reads[are_matched]
  valid_seq[!(orientation_vec)] <- reverseComplement(valid_seq[!(orientation_vec)])

  valid_qual <- as(Biostrings::quality(use_quality[are_matched]), "PhredQuality")
  valid_qual[!(orientation_vec)] <- reverse(valid_qual[!(orientation_vec)])


  trimmed_seq <- rep(NA_character_, num_matched)
  trimmed_qual <- trimmed_seq

  for (i in seq_along(fwd_condition_barcodes)) {
    message("Trimming barcodes for the condition '", condition_names[[i]], "'...")
    are_this_condition <- condition_vec == i
    matches_object <- vmatchPattern(DNAString(fwd_condition_barcodes[[i]]),
                                    valid_seq[are_this_condition],
                                    max.mismatch = 1
                                    )
    starts_vec <- unlist(start(matches_object))
    stopifnot(length(starts_vec) == sum(are_this_condition))

    this_seq <- substr(valid_seq[are_this_condition], 1, starts_vec + (include_flanking - 1))
    trimmed_seq[are_this_condition] <- this_seq
    this_qual <- substr(valid_qual[are_this_condition], 1, starts_vec + (include_flanking - 1))
    trimmed_qual[are_this_condition] <- this_qual
  }


  for (i in seq_along(fwd_replicate_barcodes)) {
    message("Trimming barcodes for replicate ", i, "...")
    are_this_replicate <- replicate_vec == i
    matches_object <- vmatchPattern(DNAString(fwd_replicate_barcodes[[i]]),
                                    trimmed_seq[are_this_replicate],
                                    max.mismatch = 1
                                    )
    starts_list <- start(matches_object)
    are_valid <- lengths(starts_list) == 1 # Some reads violate the assumption that the replicate barcode comes before the condition barcode
    starts_vec <- unlist(starts_list[are_valid])

    trimmed_seq[are_this_replicate][!(are_valid)] <- NA_character_
    trimmed_qual[are_this_replicate][!(are_valid)] <- NA_character_

    string_widths <- nchar(trimmed_seq[are_this_replicate][are_valid])
    this_seq <- substr(trimmed_seq[are_this_replicate][are_valid], starts_vec + 10, string_widths)
    trimmed_seq[are_this_replicate][are_valid] <- this_seq
    this_qual <- substr(trimmed_qual[are_this_replicate][are_valid], starts_vec + 10, string_widths)
    trimmed_qual[are_this_replicate][are_valid] <- this_qual
  }

  message("Constructing the final data frame...")

  lengths_vec <- nchar(trimmed_seq)
  are_valid <- (!(is.na(trimmed_seq))) & (lengths_vec <= 3000)

  results_df <- data.frame(
    "Sequence"     = trimmed_seq[are_valid],
    "Quality"      = trimmed_qual[are_valid],
    "Mean_quality" = GetMeanQuality(trimmed_qual[are_valid]),
    "Is_fwd"       = orientation_vec[are_valid],
    "Replicate"    = replicate_vec[are_valid],
    "Condition"    = condition_vec[are_valid]
  )
  return(results_df)
}


RangesForChunks <- function(num_reads, chunk_size) {
  if (missing(chunk_size)) {
    stop("Please provide the 'chunk_size' argument!")
  }
  chunk_from <- seq(from = 1, to = num_reads, by = chunk_size)
  chunk_to <- c(chunk_from[seq_len(length(chunk_from) - 1) + 1] - 1L,
                num_reads
                )
  chunk_mat <- cbind(
    "chunk_number" = seq_along(chunk_from),
    "from"         = chunk_from,
    "to"           = chunk_to
  )
  return(chunk_mat)
}



DemultiplexFastqFiles <- function(fastq_files,
                                  file_numbers    = NULL,
                                  chunk_size      = 7 * 10^6,
                                  max_num_chunks  = NULL,
                                  start_at_chunk  = 1,
                                  fastq_num_reads = NULL,
                                  readerBlockSize = 1e8
                                  ) {
  if (is.null(file_numbers)) {
    file_numbers <- seq_along(fastq_files)
  } else {
    stopifnot(length(fastq_files) == length(file_numbers))
  }
  for (i in file_numbers) {

    fastq_file <- fastq_files[[which(file_numbers == i)]]

    use_n <- chunk_size
    if (!(is.null(max_num_chunks))) {
      if (is.null(fastq_num_reads)) {
        fastq_num_reads <- countFastq(fastq_file)
      }
      if (fastq_num_reads > (chunk_size * max_num_chunks)) {
        ranges_mat <- RangesForChunks(fastq_num_reads, chunk_size)
        use_chunks <- seq(from = start_at_chunk,
                          to = min(start_at_chunk + max_num_chunks - 1,
                                   nrow(ranges_mat)
                                   ),
                          by = 1
                          )
        use_n <- IRanges(start = ranges_mat[, "from"][use_chunks],
                         end = ranges_mat[, "to"][use_chunks]
                         )
      }
    }

    fastq_streamer <- ShortRead::FastqStreamer(fastq_file, n = use_n, readerBlockSize = readerBlockSize)

    chunk_i <- start_at_chunk - 1L

    message("Reading in the FASTQ reads for file ", i, " and chunk ", chunk_i + 1L, "...")

    while (length(fastq_yield_chunk <- ShortRead::yield(fastq_streamer))) {

      chunk_i <- chunk_i + 1L

      message("Demultiplexing ", length(fastq_yield_chunk), " reads...")
      chunk_df <- DemultiplexReads(fastq_yield_chunk)

      message("Saving demultiplexed reads for file ", i, " and chunk ", chunk_i, "...")
      chunk_df_name <- paste0("file_", i, "_demuxed_", chunk_i, "_df")
      assign(chunk_df_name, chunk_df, envir = globalenv())
      save(list = chunk_df_name,
           file = file.path(rdata_dir, paste0("02_demultiplex_reads__file_", i, "_chunk_", chunk_i, ".RData"))
           )
      rm(list = c("chunk_df_name", "chunk_df"))
      gc()

      if ((length(fastq_yield_chunk) == chunk_size) &&
          (is.null(max_num_chunks) || (chunk_i < (start_at_chunk + max_num_chunks - 1)))
          ) {
        message("\nReading in the FASTQ reads for file ", i, " and chunk ", chunk_i + 1, "...")
      }
    }
    message("No more reads found! Demultiplexing complete!")
  }
}


