## 2022-02-17


# Load packages and source code -------------------------------------------

root_directory        <- "~/CRISPR/6) Individual experiments"
R_functions_directory <- file.path(root_directory, "2020-08-29 - PacBio - first 384-well plate/1) R functions")
source(file.path(R_functions_directory, "02) Analyzing reads.R"))



# Define functions --------------------------------------------------------

Assign_gRNAs <- function(sg_df, match_df, sg_numbers = 1:4, include_columns = NULL) {

  message("Preparing the CRISPR library for identifying non-unique plasmids...")
  upper_library_mat <- toupper(as.matrix(sg_df[, paste0("Sequence_sg", sg_numbers)]))
  plasmids_list_list <- lapply(seq_along(sg_numbers), function(x) {
    split(sg_df[, "Plasmid_ID"],
          factor(upper_library_mat[, x], levels = unique(upper_library_mat[, x]))
          )
  })
  are_unique_library_mat <- do.call(cbind, lapply(seq_along(sg_numbers), function(x) {
    seq_vec <- upper_library_mat[, x]
    !(duplicated(seq_vec) | duplicated(seq_vec, fromLast = TRUE))
  }))

  sg_combos_vec <- do.call(paste, c(as.list(data.frame(upper_library_mat)), sep = "_"))
  combo_plasmids_list <- split(sg_df[, "Plasmid_ID"], sg_combos_vec)
  combo_plasmids_vec <- vapply(combo_plasmids_list, paste0, collapse = ", ", "")

  message("Identifying sequences that map to more than one plasmid...")
  sequences_mat <- do.call(cbind, lapply(sg_numbers, function(x) {
    ifelse(match_df[, paste0("Num_MM_sg", x)] %in% 0,
           match_df[, paste0("Aligned_read_sg", x)],
           match_df[, paste0("Correct_sgRNA_sg", x)]
           )
  }))
  colnames(sequences_mat) <- paste0("Sequence_sg", sg_numbers)

  have_guide <- rowSums(is.na(sequences_mat)) != length(sg_numbers)
  have_guide_mat <- sequences_mat[have_guide, ]

  all_unique_mat <- do.call(cbind, lapply(seq_along(sg_numbers), function(x) {
    matches_vec <- match(have_guide_mat[, x], upper_library_mat[, x])
    are_unique_library_mat[matches_vec, x]
  }))
  are_unique <- rowSums(is.na(all_unique_mat) | all_unique_mat) == length(sg_numbers)


  message("Mapping reads to plasmids...")
  possible_plasmids_list_list <- lapply(seq_along(sg_numbers), function(x) {
    lapply(have_guide_mat[!(are_unique), x], function(y) {
      if (is.na(y)) {
        NA
      } else {
        plasmids_list_list[[x]][[y]]
      }
    })
  })

  short_NA_mat <- do.call(cbind, lapply(seq_along(sg_numbers), function(x) {
    is.na(possible_plasmids_list_list[[x]])
  }))

  IntersectLists <- function(index_1, index_2) {
    both_non_NA <- !(short_NA_mat[, index_1] | short_NA_mat[, index_2])
    intersect_list <- mapply(intersect,
                             possible_plasmids_list_list[[index_1]][both_non_NA],
                             possible_plasmids_list_list[[index_2]][both_non_NA]
                             )
    have_intersect <- lengths(intersect_list) > 0
    possible_plasmids_list_list[[index_1]][both_non_NA][have_intersect] <- intersect_list[have_intersect]
    possible_plasmids_list_list[[index_2]][both_non_NA][have_intersect] <- intersect_list[have_intersect]
    return(possible_plasmids_list_list)
  }


  if (length(sg_numbers) == 4) {
    old_possible_plasmids_list_list <- NULL
    for (h in 1:10) {
      if (identical(old_possible_plasmids_list_list, possible_plasmids_list_list)) {
        break
      }
      message("Iteration #", h, "...")
      old_possible_plasmids_list_list <- possible_plasmids_list_list
      for (i in 1:3) {
        possible_plasmids_list_list <- IntersectLists(i, i + 1)
      }
      for (i in 1:2) {
        possible_plasmids_list_list <- IntersectLists(i, i + 2)
      }
      possible_plasmids_list_list <- IntersectLists(1, 4)
    }
  } else if (length(sg_numbers) == 2) {
    possible_plasmids_list_list <- IntersectLists(1, 2)
  } else {
    stop("Unexpected number of sgRNAs!")
  }


  message("Collapsing lists of plasmids into comma-separated vectors...")
  plasmids_mat <- do.call(cbind, lapply(seq_along(sg_numbers), function(x) {
    plasmids_vec <- rep(NA, nrow(have_guide_mat))
    matches_vec <- match(have_guide_mat[are_unique, x], upper_library_mat[, x])
    plasmids_vec[are_unique] <- sg_df[, "Plasmid_ID"][matches_vec]
    collapsed_vec <- vapply(possible_plasmids_list_list[[x]],
                            function(y) if (all(is.na(y))) NA_character_ else paste0(y, collapse = ", "),
                            ""
                            )
    plasmids_vec[!(are_unique)] <- collapsed_vec
    return(plasmids_vec)
  }))
  colnames(plasmids_mat) <- paste0("Plasmid_sg", sg_numbers)


  message("Mapping plasmids to genes...")
  MakeGenesMat <- function(genes_column) {
    results_mat <- do.call(cbind, lapply(seq_along(sg_numbers), function(x) {
      plasmid_splits <- strsplit(plasmids_mat[, x], ", ", fixed = TRUE)
      long_vec <- unlist(plasmid_splits, use.names = FALSE)
      indices_vec <- rep(seq_along(plasmid_splits), lengths(plasmid_splits))
      matches_vec <- match(long_vec, sg_df[, "Plasmid_ID"])
      lookup_vec <- sg_df[, genes_column]
      if (genes_column == "Entrez_ID") {
        lookup_vec <- ifelse(is.na(lookup_vec), sg_df[, "Gene_symbol"], lookup_vec)
      }
      genes_vec <- lookup_vec[matches_vec]
      results_vec <- rep(NA, length(plasmid_splits))
      genes_list <- split(genes_vec, indices_vec)
      are_NA <- is.na(genes_list)
      results_vec[!(are_NA)] <- vapply(genes_list[!(are_NA)],
                                       function(x) paste0(unique(x), collapse = ", "),
                                       ""
                                       )
      return(results_vec)
    }))
    return(results_mat)
  }
  symbols_mat <- MakeGenesMat("Gene_symbol")
  colnames(symbols_mat) <- paste0("Symbol_sg", sg_numbers)
  entrezs_mat <- MakeGenesMat("Entrez_ID")
  colnames(entrezs_mat) <- paste0("Entrez_sg", sg_numbers)


  message("Identifying template switches...")
  multiple_list <- as.list(data.frame(t(plasmids_mat), stringsAsFactors = FALSE))
  unique_list <- lapply(multiple_list, function(x) unique(x[!(is.na(x))]))

  rle_list <- unique_list
  have_switch <- lengths(unique_list) > 1
  rle_list[have_switch] <- lapply(multiple_list[have_switch],
                                  function(x) rle(x[!(is.na(x))])[["values"]]
                                  )
  switch_back <- lengths(rle_list) != lengths(unique_list)
  switch_twice <- lengths(rle_list) > (lengths(unique_list) + 1)


  message("Preparing additional relevant data...")
  num_MM_mat <- as.matrix(match_df[, paste0("Num_MM_sg", sg_numbers)])[have_guide, ]
  rownames(num_MM_mat) <- NULL

  quality_columns <- paste0("Quality_sg", sg_numbers)
  if (all(quality_columns %in% names(match_df))) {
    qual_string_mat <- as.matrix(match_df[, quality_columns])[have_guide, ]
    qual_string_mat <- gsub(" ", "", qual_string_mat, fixed = TRUE)
    NA_mat <- is.na(have_guide_mat)
    qualities_vec <- GetMeanQuality(qual_string_mat[!(NA_mat)])
    qualities_mat <- matrix(nrow = nrow(NA_mat), ncol = length(sg_numbers))
    qualities_mat[!(NA_mat)] <- qualities_vec
    colnames(qualities_mat) <- paste0("Quality_sg", sg_numbers)
  } else {
    qualities_mat <- matrix(nrow = sum(have_guide), ncol = 0)
  }

  message("Creating a data frame combining all relevant data...")
  if (is.null(include_columns)) {
    column_index <- which(names(match_df) == "Aligned_template_sg1")
    include_columns <- seq_len(column_index - 1)
  }
  if ("Read_number" %in% names(match_df)) {
    read_numbers <- match_df[, "Read_number"][have_guide]
  } else {
    read_numbers <- which(have_guide)
  }
  results_df <- data.frame(
    "Read_number" = read_numbers,
    match_df[have_guide, include_columns],
    have_guide_mat, qualities_mat, num_MM_mat,
    plasmids_mat, symbols_mat, entrezs_mat,
    "Num_matched_sgRNAs"    = as.integer(rowSums(!(is.na(have_guide_mat)))),
    "Num_targeted_plasmids" = lengths(unique_list),
    "Num_template_switches" = lengths(rle_list) - 1L,
    stringsAsFactors = FALSE
  )
  if (length(sg_numbers) == 4) {
    results_df <- data.frame(
      results_df,
      "Switch_sg1_to_sg2" = plasmids_mat[, 1] != plasmids_mat[, 2],
      "Switch_sg2_to_sg3" = plasmids_mat[, 2] != plasmids_mat[, 3],
      "Switch_sg3_to_sg4" = plasmids_mat[, 3] != plasmids_mat[, 4],
      "Num_switch_backs"  = ifelse(switch_twice, 2L,
                                   ifelse(switch_back == 1, 1L, 0L)
                                   ),
      stringsAsFactors = FALSE
    )
  }
  return(results_df)
}



GetCounts <- function(library_plasmids, reads_df, sg_numbers = 1:4) {

  plasmids_mat <- as.matrix(reads_df[, paste0("Plasmid_sg", sg_numbers)])

  multiple_list <- as.list(data.frame(t(plasmids_mat), stringsAsFactors = FALSE))
  unique_list <- lapply(multiple_list, function(x) unique(x[!(is.na(x))]))

  non_unique_vec <- table(unlist(unique_list, use.names = FALSE))
  unique_plasmids_list <- strsplit(names(non_unique_vec), ", ", fixed = TRUE)
  plasmids_vec <- unlist(unique_plasmids_list)
  lengths_vec <- lengths(unique_plasmids_list)
  counts_vec <- rep(non_unique_vec, lengths_vec)
  weights_vec <- 1 / rep(lengths_vec, lengths_vec)
  counts_vec <- weights_vec * counts_vec

  matched_counts_vec <- sapply(library_plasmids, function(x) {
    are_this_plasmid <- plasmids_vec == x
    if (any(are_this_plasmid)) {
      sum(counts_vec[are_this_plasmid])
    } else {
      0L
    }
  })

  return(matched_counts_vec)
}



MakeCountsDf <- function(sg_df, mapped_df) {

  have_no_switch    <- mapped_df[, "Num_template_switches"] == 0
  no_switch_all_4   <- have_no_switch & (mapped_df[, "Num_matched_sgRNAs"] == 4)
  no_switch_2ormore <- have_no_switch & (mapped_df[, "Num_matched_sgRNAs"] >= 2)
  sg2_matches_sg3   <- (mapped_df[, "Plasmid_sg2"] == mapped_df[, "Plasmid_sg3"]) %in% TRUE
  sg3_matches_sg4   <- (mapped_df[, "Plasmid_sg3"] == mapped_df[, "Plasmid_sg4"]) %in% TRUE

  counts_all4_vec       <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                     mapped_df[no_switch_all_4, ],
                                     )
  counts_2ormore_vec    <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                     mapped_df[no_switch_2ormore, ]
                                     )
  counts_noswitch_vec   <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                     mapped_df[have_no_switch, ]
                                     )
  counts_unselected_vec <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                     mapped_df
                                     )
  counts_sg2_sg3_vec    <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                     mapped_df[sg2_matches_sg3, ]
                                     )
  counts_sg3_sg4_vec    <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                     mapped_df[sg3_matches_sg4, ]
                                     )

  use_columns <- c("Plasmid_ID", "Gene_symbol", "Entrez_ID",
                   "TSS_ID", "Is_main_TSS", "Plate_ID",
                   "Well_number", "Is_obsolete"
                   )
  counts_df <- data.frame(
    "Count_perfect"       = counts_all4_vec,
    "Count_2_filtered"    = counts_2ormore_vec,
    "Count_filtered"      = counts_noswitch_vec,
    "Count_unfiltered"    = counts_unselected_vec,
    "Count_sg2_match_sg3" = counts_sg2_sg3_vec,
    "Count_sg3_match_sg4" = counts_sg3_sg4_vec,
    sg_sequences_df[, use_columns],
    stringsAsFactors = FALSE
  )
  return(counts_df)
}


