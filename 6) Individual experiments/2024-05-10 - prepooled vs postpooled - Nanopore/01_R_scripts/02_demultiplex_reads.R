## 2024-05-30


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
pacbio_seq_functions_directory  <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate", "1) R functions")
source(file.path(pacbio_seq_functions_directory, "02) Analyzing reads.R")) # For GetMeanQuality

first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
project_dir <- first_nanopore_dir
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))
rm(project_dir)



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
input_dir <- file.path(project_dir, "02_input")
rdata_dir <- file.path(project_dir, "03_R_objects")

fastq_files <- list.files(file.path(input_dir, "Raw_reads"), full.names = TRUE)
stopifnot(all(grepl(".pass.fastq.gz", fastq_files, fixed = TRUE)))



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_prepare_barcodes.RData"))



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



# Demultiplex reads -------------------------------------------------------

chunk_size <- 7 * 10^6

for (i in seq_along(fastq_files)) {

  fastq_file <- fastq_files[[i]]

  fastq_streamer <- ShortRead::FastqStreamer(fastq_file, n = chunk_size)

  chunk_i <- 0L

  message("Reading in the FASTQ reads for file ", i, " and chunk 1...")

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

    if (length(fastq_yield_chunk) == chunk_size) {
      message("\nReading in the FASTQ reads for file ", i, " and chunk ", chunk_i + 1, "...")
    }
  }
  message("No more reads found! Demultiplexing complete!")
}



