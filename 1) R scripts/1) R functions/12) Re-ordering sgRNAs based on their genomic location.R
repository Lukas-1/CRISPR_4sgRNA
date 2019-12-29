### 2nd December 2019 ###



# Define functions --------------------------------------------------------

MessageID <- function(CRISPR_sub_df) {
  current_ID <- CRISPR_sub_df[1, "Combined_ID"]
  message(paste0("Processing sgRNAs with the combined ID: ", current_ID))
  return(invisible(NULL))
}



ReorderSubDfByLocation <- function(CRISPR_sub_df) {

  were_mapped <- !(apply(do.call(cbind,
                                 sapply(c("Chromosome", "Strand", "Start", "End"),
                                        function(x) is.na(CRISPR_sub_df[, x]),
                                        simplify = FALSE
                                        )
                                 ),
                         1, any
                         )
                   )
  if (length(unique(CRISPR_sub_df[were_mapped, "Chromosome"])) > 1) {
    chromosome_table <- sort(table(CRISPR_sub_df[were_mapped, "Chromosome"]))
    if ((chromosome_table[[1]] == 1) && (chromosome_table[[2]] >= 3) && (length(chromosome_table) == 2)) {
      minority_chromosome <- names(chromosome_table)[[1]]
      were_mapped <- ifelse(CRISPR_sub_df[, "Chromosome"] %in% minority_chromosome, FALSE, were_mapped)
      message(paste0("For the combined ID ", CRISPR_sub_df[1, "Combined_ID"], ",
                     sgRNAs mapped to multiple chromosomes, and the minority chromosome was ignored!"
                     )
              )
    } else {
      stop("The sgRNAs mapped to more than one chromosome, and this inconsistency could not be resolved!")
    }
  }

  new_order <- order(CRISPR_sub_df[, "Start"])
  sub_df_reordered <- CRISPR_sub_df[new_order, ]
  were_mapped <- were_mapped[new_order]

  results_list <- list(
    "reordered_df"    = sub_df_reordered,
    "were_mapped_vec" = were_mapped
  )
  return(results_list)
}

