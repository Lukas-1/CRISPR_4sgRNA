## 2022-02-15


# Load packages and source code -------------------------------------------

library("Biostrings")
library("ShortRead")
library("stringdist")



# Define paths ------------------------------------------------------------

project_dir <- "~/CRISPR/6) Individual experiments/2022-01-05 - first nanopore sequencing run"
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "05_extract_aligned_sgRNAs.RData"))
load(file.path(rdata_dir, "06_reformat_CRISPRa_library.RData"))



# Define functions --------------------------------------------------------

AllowOneMismatch <- function(query_sequence, template_vec) {

  distance_mat <- stringdistmatrix(query_sequence, template_vec, method = "lv")

  have_match_mat <- distance_mat <= 1
  stopifnot(!(any(distance_mat == 0)))
  have_one_match <- rowSums(have_match_mat) == 1
  have_multiple_matches <- rowSums(have_match_mat) > 1
  sub_vec <- rep(NA, length(query_sequence))
  sub_vec[have_one_match] <- template_vec[apply(have_match_mat[have_one_match, ], 1, which)]
  sub_vec[have_multiple_matches] <- vapply(which(have_multiple_matches), function(x) {
    sgRNAs_vec <- template_vec[have_match_mat[x, ]]
    sgRNAs_vec <- unique(sgRNAs_vec)
    if (length(sgRNAs_vec) == 1) {
      return(sgRNAs_vec)
    } else {
      return(NA_character_)
    }
  }, "")

  long_multiple <- rep(NA, length(sub_vec))
  long_multiple[have_multiple_matches] <- is.na(sub_vec[have_multiple_matches])

  results_df <- data.frame(
    "Sequence"             = sub_vec,
    "Has_multiple_matches" = long_multiple,
    stringsAsFactors = FALSE
  )
  return(results_df)
}


MakeChunks <- function(n_total, n_per_chunk = 1000) {
  num_chunks <- ceiling(n_total / n_per_chunk)
  chunks_vec <- rep(seq_len(num_chunks), each = n_per_chunk)[seq_len(n_total)]
  return(chunks_vec)
}


CheckGuides <- function(extract_df, sg_number = 1) {

  aligned_vec <- extracted_df[, paste0("Aligned_read_sg", sg_number)]
  template_vec <- toupper(sg_sequences_df[, paste0("Sequence_sg", sg_number)])

  ## Check for sequences that match the corresponding sgRNA perfectly
  are_identical <- aligned_vec %in% template_vec
  message("Processing sg", sg_number, "...")
  message(paste0(sum(are_identical), " out of ", nrow(extract_df), " sequences",
                 " (", formatC(sum(are_identical) / nrow(extract_df) * 100,
                               digits = 1, format = "f"
                               ),
                 "%) had perfect matches in extract_df."
                 ))

  ## Check for sequences with two or more gaps in the alignment
  num_gaps <- nchar(aligned_vec) - nchar(gsub("-", "", aligned_vec, fixed = TRUE))
  too_many_gaps <- num_gaps >= 2
  message(paste0(sum(too_many_gaps), " out of ", nrow(extract_df), " sequences",
                 " (", formatC(sum(too_many_gaps) / nrow(extract_df) * 100,
                               digits = 1, format = "f"
                               ),
                 "%) had two or more gaps (missing bases)."
                 ))

  ## Check for aligned sequences that are two or more base pairs too long
  are_too_long <- nchar(aligned_vec) > 21
  message(paste0(sum(are_too_long), " out of ", nrow(extract_df), " sequences",
                 " (", formatC(sum(are_too_long) / nrow(extract_df) * 100,
                               digits = 1, format = "f"
                               ),
                 "%) were too long (12 base pairs or longer)."
                 ))

  ## Choose the remaining sequences (to look for hits with one mismatch)
  are_remaining <- !(are_identical | too_many_gaps | are_too_long)
  remaining_vec <- aligned_vec[are_remaining]

  ## Initialize vectors
  num_remaining <- length(remaining_vec)
  sequences_vec <- rep(NA, num_remaining)
  multiple_matches <- rep(NA, num_remaining)

  ## Prepare 'large' chunks (for reporting on the progress...)
  reads_per_chunk <- 200000
  chunks_vec <- MakeChunks(num_remaining, reads_per_chunk)
  first_vec <- format(tapply(seq_len(num_remaining), chunks_vec, function(x) x[[1]]))
  last_vec  <- format(tapply(seq_len(num_remaining), chunks_vec, function(x) x[[length(x)]]))
  num_chunks <- ceiling(num_remaining / reads_per_chunk)
  chunk_numbers <- format(seq_len(num_chunks))

  ## Prepare for parallelization
  num_cores <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(num_cores)
  parallel::clusterExport(cl,
                          varlist = c("template_vec", "remaining_vec",
                                      "AllowOneMismatch", "stringdistmatrix"
                                      ),
                          envir = environment()
                          )

  ## Run a parallelized loop (for which data is split into smaller chunks...)
  for (i in seq_len(num_chunks)) {
    message("Processing chunk #", chunk_numbers[[i]], " of ",
            chunk_numbers[[length(chunk_numbers)]],  " (checking gRNAs ",
            first_vec[[i]], " to ", last_vec[[i]], ")..."
            )

    are_this_chunk <- chunks_vec == i
    small_chunks_vec <- MakeChunks(sum(are_this_chunk), 10000)
    small_chunks_list <- split(seq_len(num_remaining)[are_this_chunk],
                               small_chunks_vec
                               )

    sg_df_list <- parallel::parLapply(cl, small_chunks_list, function(x) {
      AllowOneMismatch(remaining_vec[x], template_vec)
    })
    sg_df <- do.call(rbind.data.frame, c(sg_df_list, stringsAsFactors = FALSE, make.row.names = FALSE))
    multiple_matches[are_this_chunk] <- sg_df[, "Has_multiple_matches"]
    sequences_vec[are_this_chunk] <- sg_df[, "Sequence"]
  }

  ## Initialize vectors for containing the results
  NA_vec <- rep(NA, nrow(extract_df))
  long_sequences_vec <- NA_vec
  long_multiple_vec <- NA_vec
  numMM_vec <- NA_vec

  ## Compile results
  long_sequences_vec[are_remaining] <- sequences_vec
  long_multiple_vec[are_remaining] <- multiple_matches
  numMM_vec[are_identical] <- 0L
  numMM_vec[are_remaining] <- ifelse(is.na(sequences_vec), NA, 1L)
  results_df <- data.frame(
    "Num_MM"        = numMM_vec,
    "Correct_sgRNA" = long_sequences_vec,
    "Multiple_1MM"  = long_multiple_vec,
    stringsAsFactors = FALSE
  )
  names(results_df) <- paste0(names(results_df), "_sg", sg_number)
  return(results_df)
}



# Look up aligned sgRNAs --------------------------------------------------

extract_df_list <- lapply(1:4, function(x) {
  CheckGuides(extracted_df, x)
})

matched_df <- do.call(data.frame, c(list(extracted_df), extract_df_list))



# Save data ---------------------------------------------------------------

save(matched_df,
     file = file.path(rdata_dir, "07_look_up_aligned_sgRNAs.RData")
     )






