## 2022-10-25


# Define functions --------------------------------------------------------

ExtendSamplesDf <- function(samples_df) {
  sample_splits <- strsplit(samples_df[, "Sample_name"], "_", fixed = TRUE)
  replicate_numbers <- as.integer(sub("^R", "", sapply(sample_splits, "[[", 2)))
  timepoint_numbers <- match(sapply(sample_splits, "[[", 1),
                             c("Tbefore", "T0", "T12")
                             )
  read_numbers <- rep(1:2, times = nrow(samples_df) / 2)
  results_df <- data.frame(
    samples_df,
    "Timepoint"        = sapply(sample_splits, "[[", 1),
    "Timepoint_number" = timepoint_numbers,
    "Replicate"        = replicate_numbers,
    "Read"             = read_numbers,
    "Rank"             = order(order(timepoint_numbers, replicate_numbers, read_numbers)),
    stringsAsFactors = FALSE
  )
  results_df[, "Long_name"] <- paste0(results_df[, "Timepoint"],
                                      "_rep", results_df[, "Replicate"],
                                      "_read", results_df[, "Read"]
                                      )
  return(results_df)
}



MeanBaseQuality <- function(all_reads_list) {
  samples_df <- ExtendSamplesDf(ProcessSamples(all_reads_list))
  base_vec_list <- lapply(all_reads_list, function(x) {
    qualities_mat <- as(Biostrings::quality(x), "matrix")
    colMeans(qualities_mat)
  })
  results_mat <- do.call(cbind, base_vec_list)
  colnames(results_mat) <- samples_df[, "Long_name"]
  results_mat <- results_mat[, order(samples_df[, "Rank"])]
  return(results_mat)
}



QualityDensities <- function(all_reads_list, use_reads = 1:2) {
  samples_df <- ExtendSamplesDf(ProcessSamples(all_reads_list))
  all_reads_list <- all_reads_list[order(samples_df[, "Rank"])]
  samples_df <- samples_df[order(samples_df[, "Rank"]), ]
  samples_df[, "Sample_number"] <- match(samples_df[, "Sample_number"], unique(samples_df[, "Sample_number"]))
  density_output_list <- lapply(unique(samples_df[, "Sample_number"]), function(x) {
    use_indices <- which((samples_df[, "Sample_number"] == x) &
                         (samples_df[, "Read"] %in% use_reads)
                         )
    qualities_vec <- unlist(lapply(use_indices, function(y) {
      qualities <- as(Biostrings::quality(all_reads_list[[y]]), "PhredQuality")
      ShortRead::alphabetScore(qualities) / Biostrings::width(qualities)
    }))
    assign("delete_use_reads", use_reads, envir = globalenv())
    assign("delete_use_indices", use_indices, envir = globalenv())
    assign("delete_qualities_vec", qualities_vec, envir = globalenv())
    density(qualities_vec, adj = 5)
  })
  names(density_output_list) <- unique(paste0(samples_df[, "Timepoint"], "_rep", samples_df[, "Replicate"]))
  return(density_output_list)
}



GC_Densities <- function(all_reads_list, use_reads = 1:2) {
  samples_df <- ExtendSamplesDf(ProcessSamples(all_reads_list))
  all_reads_list <- all_reads_list[order(samples_df[, "Rank"])]
  samples_df <- samples_df[order(samples_df[, "Rank"]), ]
  samples_df[, "Sample_number"] <- match(samples_df[, "Sample_number"], unique(samples_df[, "Sample_number"]))
  density_output_list <- lapply(unique(samples_df[, "Sample_number"]), function(x) {
    use_indices <- which((samples_df[, "Sample_number"] == x) &
                         (samples_df[, "Read"] %in% use_reads)
                         )
    GC_vec <- unlist(lapply(use_indices, function(y) {
      reads <- ShortRead::sread(all_reads_list[[y]])
      num_N <- letterFrequency(reads, "N")
      num_GC <- letterFrequency(reads, "GC")
      num_GC / (20 - num_N)
    }))
    density(GC_vec, adj = 10, from = 0, to = 1)
  })
  names(density_output_list) <- unique(paste0(samples_df[, "Timepoint"], "_rep", samples_df[, "Replicate"]))
  return(density_output_list)
}




