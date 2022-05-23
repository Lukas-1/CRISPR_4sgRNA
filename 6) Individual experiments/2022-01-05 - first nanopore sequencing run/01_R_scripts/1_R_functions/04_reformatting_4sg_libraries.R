## 2022-02-15



# Define functions --------------------------------------------------------

ReformatLibrary <- function(input_df) {

  ### Create a new data frame with one row per CRISPR library plasmid

  input_df[, "Plate_ID"] <- sapply(strsplit(input_df[, "Plate_string"], "_"), "[[", 2)
  plasmids_vec <- paste0(input_df[, "Entrez_ID"],
                         "__", input_df[, "Plate_ID"],
                         "__", input_df[, "Well_number"]
                         )
  input_df[, "Plasmid_ID"] <- plasmids_vec
  plasmids_fac <- factor(plasmids_vec, levels = unique(plasmids_vec))
  are_plasmid_specific <- vapply(names(input_df), function(x) {
    message("Checking column '", x, "' whether it relates to the plasmid (rather than to sgRNAs)... ")
    all(tapply(input_df[, x], plasmids_fac, function(y) length(unique(y)) == 1))
  }, logical(1))

  use_columns <- c("sgRNA_sequence", names(input_df)[are_plasmid_specific])
  plasmids_df_list <- split(input_df[, use_columns], plasmids_fac)

  unique(split(input_df[, "Rank"], plasmids_fac))

  plasmids_df_list <- lapply(plasmids_df_list, function(x) {
    results_list <- as.list(x[, "sgRNA_sequence"])
    names(results_list) <- paste0("Sequence_sg", 1:4)
    results_list <- c(as.list(x[1, names(x) != "sgRNA_sequence"]), results_list)
    return(results_list)
  })

  sg_sequences_df <- do.call(rbind.data.frame,
                             c(plasmids_df_list,
                               stringsAsFactors = FALSE,
                               make.row.names = FALSE
                             ))

  return(sg_sequences_df)
}






