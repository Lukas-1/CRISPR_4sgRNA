### 6th June 2021 ###




# Define functions --------------------------------------------------------

MakeRandomDeletionDf <- function(range_spans, deletion_sizes, barcode_length = 27L) {

  stopifnot(identical(length(range_spans), length(deletion_sizes)))

  random_ranges <- vapply(range_spans, function(x) sample(seq_len(x), 1), integer(1))
  random_starts <- random_ranges + barcode_length

  random_df <- data.frame(
    "Deletion_start" = random_starts,
    "Deletion_end"   = random_starts + deletion_sizes
  )

  random_df <- AnnotateDeletions(random_df)
  return(random_df)
}




GetSimulatedMat <- function(deletion_sizes,
                            plasmid_length  = 2225L,
                            barcode_length  = 27L,
                            set_seed        = TRUE,
                            num_simulations = 1000L
                            ) {

  range_spans <- plasmid_length - deletion_sizes
  if (set_seed) {
    set.seed(1)
  }

  random_del_df_list <- lapply(seq_len(num_simulations), function(x) {
    if (x %in% c(1, seq(20, num_simulations, by = 20))) {
      message(paste0("Working on iteration #", x, "..."))
    }
    MakeRandomDeletionDf(range_spans, deletion_sizes, barcode_length = barcode_length)
  })

  use_columns <- c("Span_tracrRNAs", "Span_sgRNAs", "Span_sg_cr", "Span_promoters")
  num_spanning_mat <- do.call(rbind, lapply(random_del_df_list, function(x) {
    colSums(as.matrix(x[, use_columns]))
  }))
  percent_spanning_mat <- num_spanning_mat / length(deletion_sizes)
  return(percent_spanning_mat)
}




FilterDeletionsDf <- function(del_df, use_df_list) {

  pass_filters <- use_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
  use_zmws <- use_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]
  are_selected <- del_df[["ZMW"]] %in% use_zmws

  results_df <- del_df[are_selected, ]
  row.names(results_df) <- NULL
  results_df <- PrioritizeDeletions(results_df)
  return(results_df)
}




