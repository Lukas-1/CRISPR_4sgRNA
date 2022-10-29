### 2022-04-09



# Define functions --------------------------------------------------------

GetCounts2sg <- function(mapped_df,
                         only_0MM = FALSE,
                         no_template_switch = FALSE,
                         no_template_switch_requires_both_sgRNAs = TRUE,
                         choose_sample = NULL
                         ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  are_selected <- rep(TRUE, nrow(mapped_df))

  if (only_0MM) {
    have_0MM <- (mapped_df[, "Num_MM_sg1"] %in% 0) &
                (mapped_df[, "Num_MM_sg2"] %in% 0)
    are_selected <- have_0MM
  }
  if (no_template_switch) {
    have_no_switch <- (!(mapped_df[, "Has_template_switch"]))
    if (no_template_switch_requires_both_sgRNAs) {
      have_no_switch <- have_no_switch & (mapped_df[, "Num_matched_sgRNAs"] == 2L)
    }
    are_selected <- are_selected & have_no_switch
  }
  if (!(is.null(choose_sample))) {
    are_this_sample <- mapped_df[, "Sample_number"] == choose_sample
    are_selected <- are_selected & are_this_sample
  }

  sel_mapped_df <- mapped_df[are_selected, c("Plasmid_sg1", "Plasmid_sg2")]
  row.names(sel_mapped_df) <- NULL
  counts_vec <- GetCounts(sg_sequences_df[, "Plasmid_ID"], sel_mapped_df,
                          sg_numbers = 1:2
                          )
  return(counts_vec)
}



AllSamplesCounts <- function(mapped_df, only_0MM = FALSE, no_template_switch = FALSE) {

  sample_numbers <- unique(mapped_df[, "Sample_number"])
  sample_names <- mapped_df[, "Sample_name"][match(sample_numbers, mapped_df[, "Sample_number"])]
  short_names <- sub("before_off", "before", sample_names, fixed = TRUE)
  short_names <- sub("2sg_", "", short_names, fixed = TRUE)

  sample_counts_list <- lapply(seq_along(sample_numbers), function(x) {
    message(paste0("Calculating counts for the sample: ", short_names[[x]], "..."))
    GetCounts2sg(mapped_df, choose_sample = sample_numbers[[x]],
                 only_0MM = only_0MM, no_template_switch = no_template_switch
                 )
  })
  results_mat <- do.call(cbind, sample_counts_list)
  colnames(results_mat) <- short_names
  rownames(results_mat) <- NULL
  return(results_mat)
}



NumMappedReads <- function(use_lumi_df, sg_numbers = 1:2) {
  samples_fac <- factor(use_lumi_df[, "Sample_name"],
                        levels = unique(use_lumi_df[, "Sample_name"])
                        )
  results_df <- data.frame(
    "Sample_name"             = levels(samples_fac),
    "Num_either_read_mapped"  = tabulate(samples_fac),
    "Num_both_reads_mapped"   = tapply(use_lumi_df[, "Num_matched_sgRNAs"] == 2, samples_fac, sum),
    "Num_unmapped_read1_only" = tapply(is.na(use_lumi_df[, paste0("Plasmid_sg", sg_numbers[[1]])]), samples_fac, sum),
    "Num_unmapped_read2_only" = tapply(is.na(use_lumi_df[, paste0("Plasmid_sg", sg_numbers[[2]])]), samples_fac, sum),
    "Num_template_switch"     = tapply(use_lumi_df[, "Has_template_switch"], samples_fac, sum),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



