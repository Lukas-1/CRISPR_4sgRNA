### 9th September 2019 ####




# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R"))




# Define constants --------------------------------------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End", "PAM")




# Functions for integrating with a genome search and GuideScan ------------

ExtendWithGenomeSearch <- function(CRISPR_df, search_df, allow_5pG = FALSE) {
  genome_matches <- match(toupper(CRISPR_df[["sgRNA_sequence"]]), toupper(search_df[["Sequence"]]))
  stopifnot(!(anyNA(genome_matches)))
  results_df <- data.frame(
    CRISPR_df,
    search_df[genome_matches, !(names(search_df) %in% c("Sequence", "Found"))],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  stopifnot(!(anyNA(results_df[["Num_0MM"]])))
  have_hit <- results_df[["Num_0MM"]] > 0
  if (allow_5pG) {
    have_hit <- ifelse(grepl("hCRISPR", CRISPR_df[["Source"]], fixed = TRUE) & (CRISPR_df[["Is_control"]] == "No"),
                       have_hit | ((results_df[["Num_0MM"]] == 0) & (results_df[["Num_5G_MM"]] > 0)),
                       have_hit
                       )
  }
  for (column_name in c("Chromosome", "Strand", "Start", "End")) {
    results_df[[column_name]] <- ifelse(have_hit, results_df[[column_name]], NA)
    names(results_df)[names(results_df) == column_name] <- paste0("Hits_", tolower(column_name))
  }
  return(results_df)
}



EntrezIDsToChromosomes <- function(entrez_vec) {
  unique_entrez_ID_strings <- unique(entrez_vec)
  expanded_entrez_IDs_df <- ExpandList(strsplit(unique_entrez_ID_strings, ", ", fixed = TRUE))
  expanded_entrez_IDs_df[["Chromosome"]] <- entrez_to_symbol_df[["Chromosome"]][match(expanded_entrez_IDs_df[["Value"]], entrez_to_symbol_df[["Entrez_ID"]])]
  chromosomes_list <- split(expanded_entrez_IDs_df[["Chromosome"]], expanded_entrez_IDs_df[["List_index"]])
  chromosomes_list <- lapply(chromosomes_list, unique)
  chromosomes_vec <- vapply(chromosomes_list, function(x) if (length(x) > 1) NA_character_ else x, "")
  results_vec <- chromosomes_vec[match(entrez_vec, unique_entrez_ID_strings)]
  return(results_vec)
}



GetCases <- function(string_vec) {
  is_lower_vec <- (toupper(string_vec) != string_vec) %in% TRUE
  is_upper_vec <- (tolower(string_vec) != string_vec) %in% TRUE
  results_vec <- ifelse(is_lower_vec,
                        "lower",
                        ifelse(is_upper_vec, "upper", "neither")
                        )
  return(results_vec)
}


MergeTSSandGuideScan <- function(CRISPR_df, guidescan_df, combined_TSS_df) {

  CRISPR_df[["Entrez_chromosome"]] <- EntrezIDsToChromosomes(CRISPR_df[["Entrez_ID"]])

  combined_TSS_df[["GuideScan_input"]] <- TSSStringForGuideScan(combined_TSS_df)
  tidy_guidescan_df <- TidyGuideScanColumns(guidescan_df)

  TSS_columns <- c("Strand", "Best_TSS", "First_TSS", "Last_TSS")

  IDs_per_region_list <- sapply(unique(tidy_guidescan_df[["Region"]]), function(x) {
    unique(combined_TSS_df[combined_TSS_df[["GuideScan_input"]] == x, c("Combined_ID", "Entrez_ID", "Gene_symbol", TSS_columns)])
  }, simplify = FALSE)
  guidescan_per_region_list <- split(tidy_guidescan_df, factor(tidy_guidescan_df[["Region"]], levels = unique(tidy_guidescan_df[["Region"]])))
  guidescan_per_region_extended_list <- sapply(names(guidescan_per_region_list), function(x) {
    num_reps_df <- IDs_per_region_list[[x]]
    gs_df <- guidescan_per_region_list[[x]]
    sub_list <- lapply(seq_len(nrow(num_reps_df)), function(x) data.frame(num_reps_df[rep.int(x, nrow(gs_df)), ], gs_df))
    do.call(rbind.data.frame, c(sub_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  }, simplify = FALSE)
  extended_guidescan_df <- do.call(rbind.data.frame, c(guidescan_per_region_extended_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  vec_for_matching_CRISPR               <- paste0(CRISPR_df[["Combined_ID"]], "__", toupper(CRISPR_df[["sgRNA_sequence"]]))
  vec_for_matching_GuideScan            <- paste0(extended_guidescan_df[["Combined_ID"]], "__", toupper(extended_guidescan_df[["gRNA"]]))
  vec_for_matching_CRISPR_chromosome    <- paste0(CRISPR_df[["Entrez_chromosome"]], "__", vec_for_matching_CRISPR)
  vec_for_matching_GuideScan_chromosome <- paste0(extended_guidescan_df[["GuideScan_chromosome"]], "__", vec_for_matching_GuideScan)

  matches_vec <- ifelse(is.na(CRISPR_df[["Entrez_chromosome"]]),
                        match(vec_for_matching_CRISPR, vec_for_matching_GuideScan),
                        match(vec_for_matching_CRISPR_chromosome, vec_for_matching_GuideScan_chromosome)
                        )

  names(extended_guidescan_df)[names(extended_guidescan_df) == "Entrez_ID"]   <- "GuideScan_entrez_ID"
  names(extended_guidescan_df)[names(extended_guidescan_df) == "Gene_symbol"] <- "GuideScan_symbol"
  names(extended_guidescan_df)[names(extended_guidescan_df) == "Strand"]      <- "Strand_of_TSS"
  names(extended_guidescan_df)[names(extended_guidescan_df) == "Region"]      <- "GuideScan_region"

  combined_df <- CRISPR_df
  for (column_name in names(extended_guidescan_df)[names(extended_guidescan_df) != "Combined_ID"]) {
    combined_df[[column_name]] <- extended_guidescan_df[[column_name]][matches_vec]
  }
  return(combined_df)
}





# Functions for tidying and summarizing the data from GuideScan -----------

GetOffTargetCategory <- function(merged_CRISPR_df) {
  if ("Num_5G_MM" %in% names(merged_CRISPR_df)) {
    num_perfect_hits_vec <- rowSums(merged_CRISPR_df[, c("Num_0MM", "Num_5G_MM")])
  } else {
    num_perfect_hits_vec <- merged_CRISPR_df[["Num_0MM"]]
  }
  is_unspecific_vec <- (num_perfect_hits_vec > 1) | (merged_CRISPR_df[["Num_1MM"]] > 0)
  offtarget_levels <- c(
    "5 or fewer",
    "10 or fewer",
    "20 or fewer",
    "50 or fewer",
    "51 or more",
    "Unspecific",
    "Unknown",
    "Not scanned"
  )
  results_vec <- vapply(seq_len(nrow(merged_CRISPR_df)), function(x) {
    num_2or3MM <- merged_CRISPR_df[["GuideScan_Num_2or3MM"]][[x]]
    if (is.na(num_2or3MM) && is_unspecific_vec[[x]]) {
      return(offtarget_levels[[6]])
    } else if (is.na(num_2or3MM)) {
      return(offtarget_levels[[7]])
    } else if (num_2or3MM <= 5) {
      return(offtarget_levels[[1]])
    } else if (num_2or3MM <= 10) {
      return(offtarget_levels[[2]])
    } else if (num_2or3MM <= 20) {
      return(offtarget_levels[[3]])
    } else if (num_2or3MM <= 50) {
      return(offtarget_levels[[4]])
    } else {
      return(offtarget_levels[[5]])
    }
  }, "")

  results_fac <- factor(results_vec, levels = offtarget_levels, ordered = TRUE)
  results_fac[is.na(results_fac)] <- "Unknown"

  if ("TSS_searched_by_GuideScan" %in% names(merged_CRISPR_df)) {
    were_not_scanned <- (results_fac %in% "Unknown") & (merged_CRISPR_df[["TSS_searched_by_GuideScan"]] %in% c("No", "Not this gene"))
    results_fac[were_not_scanned] <- "Not scanned"
  }
  return(results_fac)
}





# Functions for resolving ambiguous gene IDs or sgRNA locations -----------

AssignToGeneByNearbyTSS <- function(CRISPR_df, combined_TSS_df, prefix = "") {

  GRanges_object_TSSs <- GRanges(
    seqnames = sub("chr", "", combined_TSS_df[["Chromosome"]], fixed = TRUE),
    ranges   = IRanges(start = combined_TSS_df[["First_TSS"]] - 1500L, end = combined_TSS_df[["Last_TSS"]] + 1500L),
    strand   = combined_TSS_df[["Strand"]]
  )
  combined_IDs <- unique(CRISPR_df[["Combined_ID"]][CRISPR_df[["Is_control"]] == "No"])

  new_entrezs_vec <- rep.int(NA_character_, nrow(CRISPR_df))
  assignment_vec <- rep.int(NA_character_, nrow(CRISPR_df))

  for (combined_ID in combined_IDs) {
    are_this_ID <- CRISPR_df[["Combined_ID"]] == combined_ID
    sub_df <- CRISPR_df[are_this_ID, ]
    mapped_entrezs <- unique(sub_df[["GuideScan_entrez_ID"]])
    mapped_entrezs <- mapped_entrezs[!(is.na(mapped_entrezs))]
    if (length(mapped_entrezs) == 1) {
      assignment <- "mapped using GuideScan"
    } else if (length(mapped_entrezs) == 0) {
      overlapping_entrezs_list <- lapply(which(!(is.na(sub_df[["Start"]]))), function(x) {
        this_GRange_object <- GRanges(
          seqnames = sub("chr", "", sub_df[["Chromosome"]][[x]], fixed = TRUE),
          ranges   = IRanges(start = sub_df[["Start"]][[x]], end = sub_df[["End"]][[x]]),
          strand   = sub_df[["Strand"]][[x]]
        )
        overlap_matches_df <- as.data.frame(findOverlaps(this_GRange_object, GRanges_object_TSSs, ignore.strand = TRUE))
        matched_entrezs <- unique(combined_TSS_df[["Entrez_ID"]][overlap_matches_df[[2]]])
        return(matched_entrezs)
      })
      mapped_entrezs <- unique(unlist(overlapping_entrezs_list))
      mapped_entrezs <- mapped_entrezs[!(is.na(mapped_entrezs))]
      if (length(mapped_entrezs) == 1) {
        assignment <- "mapped by proximity to TSS"
      }
    }
    if (length(mapped_entrezs) == 1) {
      new_entrezs_vec[are_this_ID] <- mapped_entrezs
    } else if (length(mapped_entrezs) == 0) {
      assignment <- "could not be mapped to any gene"
    } else {
      assignment <- "could not be mapped unambiguously"
    }
    assignment_vec[are_this_ID] <- paste0(prefix, assignment)
  }
  results_df <- MapToEntrezs(new_entrezs_vec)
  results_df <- data.frame(results_df, "Assignment" = assignment_vec, stringsAsFactors = FALSE, row.names = NULL)
  return(results_df)
}



GetCutLocations <- function(ranges_df) {
  ifelse(ranges_df[["Strand"]] == "+", ranges_df[["End"]] - 2L, ranges_df[["Start"]] + 3L)
}


LocationStringToDf <- function(location_char_vec) {
  chromosome_splits <- strsplit(location_char_vec, "(", fixed = TRUE)
  strand_splits <- strsplit(sapply(chromosome_splits, "[[", 2), ")", fixed = TRUE)
  location_vec <- sapply(strand_splits, "[[", 2)
  location_vec <- substr(location_vec, 2, nchar(location_vec))
  location_splits <- strsplit(location_vec, "-", fixed = TRUE)
  results_df <- data.frame(
    "Chromosome" = sapply(chromosome_splits, "[[", 1),
    "Strand"     = sapply(strand_splits, "[[", 1),
    "Start"      = as.integer(sapply(location_splits, "[[", 1)),
    "End"        = as.integer(sapply(location_splits, "[[", 2)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



NoMatchForChromosome <- function(force_stop) {
  if (force_stop) {
    stop("No hit with a matching chromosome  was found!")
  }
  results_df <- data.frame("Chromosome" = NA_character_,
                           "Strand"     = NA_character_,
                           "Start"      = NA_integer_,
                           "End"        = NA_integer_,
                           "PAM"        = NA_character_,
                           stringsAsFactors = FALSE
                           )
  return(results_df)
}



FindSingleLocationFromTSS <- function(location_df, chromosome, TSS_location = NA, force_stop = TRUE) {
  are_valid <- substr(location_df[["PAM"]], 2, 3) == "GG"
  if (!(is.na(chromosome))) {
    are_valid <- are_valid & (location_df[["Chromosome"]] == chromosome)
  }
  if (!(any(are_valid))) {
    return(NoMatchForChromosome(force_stop))
  } else {
    location_df <- location_df[are_valid, ]
  }
  if (is.na(TSS_location)) {
    location_df <- location_df[1, ]
  } else {
    cut_location_vec <- GetCutLocations(location_df)
    TSS_distances <- abs(TSS_location - cut_location_vec)
    location_df <- location_df[which.min(TSS_distances), ]
  }
  return(location_df)
}




FindSingleLocationFromGene <- function(location_df, chromosome, gene_GRanges_object, force_stop = TRUE) {
  are_valid <- substr(location_df[["PAM"]], 2, 3) == "GG"
  if (!(is.na(chromosome))) {
    are_valid <- are_valid & (location_df[["Chromosome"]] == chromosome)
  }
  if (!(any(are_valid))) {
    return(NoMatchForChromosome(force_stop))
  } else {
    location_df <- location_df[are_valid, ]
  }
  if (is.null(gene_GRanges_object)) {
    location_df <- location_df[1, ]
  } else {
    sgRNAs_GRanges_object <- RangesDfToGRangesObject(location_df)
    distance_vec <- mcols(distanceToNearest(sgRNAs_GRanges_object, gene_GRanges_object, ignore.strand = TRUE, select = "all"))[, 1]
    if (length(distance_vec) == 0) {
      location_df <- location_df[1, ] # i.e. gene_GRanges_object is on a different chromosome
    } else {
      stopifnot(length(distance_vec) == nrow(location_df))
      location_df <- location_df[which.min(distance_vec), ]
    }
  }
  return(location_df)
}






Choose0MMLocation <- function(CRISPR_df, use_TSS) {

  have_multiple_0MM <- CRISPR_df[["Num_0MM"]] > 1

  print(unique(CRISPR_df[, c("Entrez_ID", "Gene_symbol")]))
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  assign("delete_have_multiple_0MM", have_multiple_0MM, envir = globalenv())
  assign("delete_use_TSS", use_TSS, envir = globalenv())

  location_splits <- strsplit(CRISPR_df[["Locations_0MM"]][have_multiple_0MM], "; ", fixed = TRUE)
  PAM_splits <- strsplit(CRISPR_df[["PAM_0MM"]][have_multiple_0MM], "; ", fixed = TRUE)

  location_df <- LocationStringToDf(unlist(location_splits, use.names = FALSE))
  location_df[["PAM"]] <- unlist(PAM_splits, use.names = FALSE)
  location_df_splits <- split(location_df, rep(seq_along(location_splits), lengths(location_splits)))

  assign("delete_location_df_splits",    location_df_splits,    envir = globalenv())

  entrez_chromosome <- unique(CRISPR_df[["Entrez_chromosome"]][have_multiple_0MM])
  if (entrez_chromosome %in% "chrX, chrY") {
    entrez_chromosome <- "chrX"
  } else if (is.na(entrez_chromosome)) {
    chromosome_table <- table(location_df[["Chromosome"]])
    num_occurrences <- as.integer(chromosome_table)
    entrez_chromosome <- names(chromosome_table)[num_occurrences == max(num_occurrences)][[1]]
  }

  if (use_TSS) {
    best_TSS <- unique(GetBestTSSPositions(CRISPR_df[have_multiple_0MM, ]))
    assign("delete_have_multiple_0MM", have_multiple_0MM,    envir = globalenv())
    assign("delete_CRISPR_df",         CRISPR_df,            envir = globalenv())
    assign("delete_best_TSS",          best_TSS,             envir = globalenv())
    assign("delete_entrez_chromosome", entrez_chromosome,    envir = globalenv())
    assign("delete_splits_list",       location_df_splits,   envir = globalenv())
    locations_df_list <- lapply(location_df_splits, function(x) FindSingleLocationFromTSS(x, entrez_chromosome, best_TSS, force_stop = FALSE))
  } else {
    entrez_ID <- unique(CRISPR_df[["Entrez_ID"]][have_multiple_0MM])
    if (is.na(entrez_ID)) {
      gene_GRanges_object <- NULL
    } else {
      are_this_entrez <- (mcols(human_genes_GRanges)[, "gene_id"] == entrez_ID) &
                         (as.character(seqnames(human_genes_GRanges)) == entrez_chromosome)
      if (any(are_this_entrez)) {
        gene_GRanges_object <- human_genes_GRanges[are_this_entrez]
      } else {
        gene_GRanges_object <- NULL
      }
    }
    locations_df_list <- lapply(location_df_splits, function(x) FindSingleLocationFromGene(x, entrez_chromosome, gene_GRanges_object, force_stop = FALSE))
  }
  replaced_results_df <- do.call(rbind.data.frame, c(locations_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  results_df <- CRISPR_df


  #### Delete this!!!! ###
  stopifnot("PAM" %in% location_columns)





  for (column_name in location_columns) {
    results_df[[column_name]][have_multiple_0MM] <- replaced_results_df[[column_name]]
  }
  stopifnot(identical(nrow(results_df), nrow(CRISPR_df)))
  assign("delete_results_sub_df", results_df, envir = globalenv())
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  return(results_df)
}








# Functions for refining the genomic locations of sgRNAs ------------------

ReassignTSS <- function(merged_CRISPR_df, combined_TSS_df) {
  # Assigns a TSS to all sgRNAs based on the gene ID they are associated with (also to those sgRNAs which were not found on GuideScan!)

  chromosome_vec <- merged_CRISPR_df[["Entrez_chromosome"]]
  chromosome_vec <- ifelse(chromosome_vec == "chrX, chrY", "chrX", chromosome_vec)
  chromosome_vec <- ifelse(is.na(chromosome_vec), merged_CRISPR_df[["Chromosome"]], chromosome_vec)

  vec_for_matching_CRISPR <- paste0(merged_CRISPR_df[["Combined_ID"]], "__",chromosome_vec)
  vec_for_matching_TSS <- paste0(combined_TSS_df[["Combined_ID"]], "__", combined_TSS_df[["Chromosome"]])
  matches_vec <- match(vec_for_matching_CRISPR, vec_for_matching_TSS)
  for (column_name in c("Strand_of_TSS", "Best_TSS", "First_TSS", "Last_TSS")) {
    TSS_column_names <- c("Strand_of_TSS" = "Strand")
    if (column_name %in% names(TSS_column_names)) {
      TSS_column_name <- TSS_column_names[[column_name]]
    } else {
      TSS_column_name <- column_name
    }
    merged_CRISPR_df[[column_name]] <- ifelse(is.na(matches_vec),
                                              merged_CRISPR_df[[column_name]],
                                              combined_TSS_df[[TSS_column_name]][matches_vec]
                                              )
  }
  return(merged_CRISPR_df)
}




GetBestTSSPositions <- function(CRISPR_df) {
  ifelse(is.na(CRISPR_df[["Best_TSS"]]),
         ifelse(CRISPR_df[["Strand_of_TSS"]] == "+", CRISPR_df[["First_TSS"]], CRISPR_df[["Last_TSS"]]),
         CRISPR_df[["Best_TSS"]]
         )
}



GetDistanceFromTSS <- function(CRISPR_df) {
  TSS_position_vec <- GetBestTSSPositions(CRISPR_df)
  have_neg_strand_TSS <- CRISPR_df[["Strand_of_TSS"]] == "-"
  results_vec <- (CRISPR_df[["Cut_location"]] - TSS_position_vec) * ifelse(have_neg_strand_TSS, -1L, 1L) +
                 ifelse(have_neg_strand_TSS, +2L, -1L)
  return(results_vec)
}




MergeLocations <- function(merged_CRISPR_df, combined_TSS_df) {

  ### Check for discrepancies between various methods of finding the location of an sgRNA
  discordant_mat <- cbind(
    "Chromosome_entrez_hits"      = !(is.na(merged_CRISPR_df[["Entrez_chromosome"]])) & !(is.na(merged_CRISPR_df[["Hits_chromosome"]]))      & (merged_CRISPR_df[["Entrez_chromosome"]] != merged_CRISPR_df[["Hits_chromosome"]]),
    "Chromosome_entrez_GuideScan" = !(is.na(merged_CRISPR_df[["Entrez_chromosome"]])) & !(is.na(merged_CRISPR_df[["GuideScan_chromosome"]])) & (merged_CRISPR_df[["Entrez_chromosome"]] != merged_CRISPR_df[["GuideScan_chromosome"]]),
    "Chromosome_hits_GuideScan"   = !(is.na(merged_CRISPR_df[["Hits_chromosome"]]))   & !(is.na(merged_CRISPR_df[["GuideScan_chromosome"]])) & (merged_CRISPR_df[["Hits_chromosome"]]   != merged_CRISPR_df[["GuideScan_chromosome"]]),
    "Strand"                      = !(is.na(merged_CRISPR_df[["Hits_strand"]]))       & !(is.na(merged_CRISPR_df[["GuideScan_strand"]]))     & (merged_CRISPR_df[["Hits_strand"]]       != merged_CRISPR_df[["GuideScan_strand"]]),
    "Start"                       = !(is.na(merged_CRISPR_df[["Hits_start"]]))        & !(is.na(merged_CRISPR_df[["GuideScan_start"]]))      & (merged_CRISPR_df[["Hits_start"]] != merged_CRISPR_df[["GuideScan_start"]]),
    "End"                         = !(is.na(merged_CRISPR_df[["Hits_end"]]))          & !(is.na(merged_CRISPR_df[["GuideScan_end"]]))        & (merged_CRISPR_df[["Hits_end"]]   != merged_CRISPR_df[["GuideScan_end"]])
  )
  merged_CRISPR_df[["Discordant_locations"]] <- rowSums(discordant_mat) >= 1
  is_discordant_location <- rowSums(discordant_mat[, c("Chromosome_hits_GuideScan", "Strand", "Strand", "Start", "End")]) >= 1
  stopifnot(all(!(is_discordant_location), na.rm = TRUE))

  ### Assign each sgRNA to a location
  were_hits <- !(is.na(merged_CRISPR_df[["Hits_start"]]))

  merged_CRISPR_df[["Chromosome"]] <- ifelse(were_hits,
                                             merged_CRISPR_df[["Hits_chromosome"]],
                                             merged_CRISPR_df[["Entrez_chromosome"]]
                                             )
  merged_CRISPR_df[["Strand"]] <- ifelse(were_hits, merged_CRISPR_df[["Hits_strand"]], merged_CRISPR_df[["GuideScan_strand"]])
  merged_CRISPR_df[["Start"]]  <- ifelse(were_hits, merged_CRISPR_df[["Hits_start"]],  merged_CRISPR_df[["GuideScan_start"]])
  merged_CRISPR_df[["End"]]    <- ifelse(were_hits, merged_CRISPR_df[["Hits_end"]],    merged_CRISPR_df[["GuideScan_end"]])

  merged_CRISPR_df <- ReassignTSS(merged_CRISPR_df, combined_TSS_df)

  return(merged_CRISPR_df)
}





AdjustPositionColumns <- function(merged_CRISPR_df, guidescan_df, combined_TSS_df, reorder_by_rank = TRUE, allow_5pG_MM = TRUE, minimal_version = FALSE) {

  # Assign sgRNAs to their genomic locations
  merged_CRISPR_df <- MergeLocations(merged_CRISPR_df, combined_TSS_df)

  # Resolve the issue of sgRNAs that could not be assigned to a single Entrez ID based on gene annotation alone (i.e. the gene symbol or Entrez ID provided)
  # (These will be assigned to a gene whose TSS is located within 1500 bp of the sgRNA)
  are_ambiguous <- grepl(",", merged_CRISPR_df[["Entrez_ID"]], fixed = TRUE)
  are_NA <- is.na(merged_CRISPR_df[["Entrez_ID"]])
  merged_CRISPR_df[["Entrez_ID_assignment"]] <- ifelse(!(are_ambiguous | are_NA), "Unambiguous", NA_character_)
  if (any(are_ambiguous)) {
    ambiguous_df <- AssignToGeneByNearbyTSS(merged_CRISPR_df[are_ambiguous, ], combined_TSS_df, prefix = "The gene symbol was ambiguous; ")
    merged_CRISPR_df[["Entrez_ID"]][are_ambiguous]            <- ifelse(is.na(ambiguous_df[["Entrez_ID"]]), merged_CRISPR_df[["Entrez_ID"]][are_ambiguous],   ambiguous_df[["Entrez_ID"]])
    merged_CRISPR_df[["Gene_symbol"]][are_ambiguous]          <- ifelse(is.na(ambiguous_df[["Entrez_ID"]]), merged_CRISPR_df[["Gene_symbol"]][are_ambiguous], ambiguous_df[["Gene_symbol"]])
    merged_CRISPR_df[["Combined_ID"]][are_ambiguous]          <- ifelse(is.na(ambiguous_df[["Entrez_ID"]]), merged_CRISPR_df[["Combined_ID"]][are_ambiguous], ambiguous_df[["Entrez_ID"]])
    merged_CRISPR_df[["Entrez_ID_assignment"]][are_ambiguous] <- ambiguous_df[["Assignment"]]
  }
  if (any(are_NA)) {
    NA_df <- AssignToGeneByNearbyTSS(merged_CRISPR_df[are_NA, ], combined_TSS_df, prefix = "No Entrez ID was found for the gene symbol; ")
    merged_CRISPR_df[["Entrez_ID"]][are_NA]                   <- ifelse(is.na(NA_df[["Entrez_ID"]]), merged_CRISPR_df[["Entrez_ID"]][are_NA],   NA_df[["Entrez_ID"]])
    merged_CRISPR_df[["Gene_symbol"]][are_NA]                 <- ifelse(is.na(NA_df[["Entrez_ID"]]), merged_CRISPR_df[["Gene_symbol"]][are_NA], NA_df[["Gene_symbol"]])
    merged_CRISPR_df[["Combined_ID"]][are_NA]                 <- ifelse(is.na(NA_df[["Entrez_ID"]]), merged_CRISPR_df[["Combined_ID"]][are_NA], NA_df[["Entrez_ID"]])
    merged_CRISPR_df[["Entrez_ID_assignment"]][are_NA]        <- NA_df[["Assignment"]]
  }

  # After assigning Entrez IDs to more sgRNAs (in the previous step), perform another merge with the data from GuideScan
  remerged_CRISPR_df <- MergeTSSandGuideScan(merged_CRISPR_df, guidescan_df, combined_TSS_df)
  for (my_column in names(remerged_CRISPR_df)) {
    merged_CRISPR_df[[my_column]] <- remerged_CRISPR_df[[my_column]]
  }
  merged_CRISPR_df <- MergeLocations(merged_CRISPR_df, combined_TSS_df)

  assign("merged_CRISPR_df_pre_duplicated", merged_CRISPR_df, envir = globalenv())


  if (!(minimal_version)) {
    ## Assign duplicated genes to a location
    merged_CRISPR_df <- FindBest0MMLocations(merged_CRISPR_df, parallel_mode = TRUE)
    merged_CRISPR_df <- ReassignTSS(merged_CRISPR_df, combined_TSS_df)
  }

  # Remove location data for sgRNAs that seem to have discrepant locations
  for (column in c("Chromosome", "Strand", "Start", "End", "PAM")) {
    merged_CRISPR_df[[column]][merged_CRISPR_df[["Discordant_locations"]]] <- NA
  }

  if (!(minimal_version)) {
    ### Calculate the "cut" location
    merged_CRISPR_df[["Cut_location"]] <- GetCutLocations(merged_CRISPR_df)

    ### Calculate the distance from the TSS
    merged_CRISPR_df[["Distance_from_TSS"]] <- GetDistanceFromTSS(merged_CRISPR_df)
  }

  # Eliminate duplicated sgRNAs
  merged_CRISPR_df <- ResolveDuplicates(merged_CRISPR_df,
                                        concatenate_columns = c("Sublibrary",
                                                                "hCRISPRa_v2_ID", "hCRISPRa_TSS_source",
                                                                "hCRISPRi_v2_ID", "hCRISPRi_TSS_source"
                                                                )
                                        )

  if (!(minimal_version)) {
    ####################################################################################################################
    ### The following section of code attempts to troubleshoot the issue of sgRNAs that were not found by GuideScan. ###
    ### Is it because the regions where these sgRNAs are located were not submitted to GuideScan?                    ###
    ### Or were they submitted, but the sgRNA was not present in GuideScan's database?                               ###
    ####################################################################################################################

    TSS_ranges_df <- data.frame(combined_TSS_df, TSSRangesForGuideScan(combined_TSS_df), stringsAsFactors = FALSE, row.names = NULL)
    TSS_ranges_df[["Region"]] <- TSSStringForGuideScan(combined_TSS_df)

    IDs_fac <- factor(merged_CRISPR_df[["Combined_ID"]], levels = unique(merged_CRISPR_df[["Combined_ID"]]))

    ### Create a new column that lists of all regions that were submitted to GuideScan for that gene
    regions_per_gene_list <- tapply(
      seq_len(nrow(merged_CRISPR_df)),
      IDs_fac,
      function(x) {
        are_this_ID <- TSS_ranges_df[["Combined_ID"]] == merged_CRISPR_df[["Combined_ID"]][[x[[1]]]]
        if (any(are_this_ID)) {
          my_result <- paste0(unique(TSS_ranges_df[["Region"]][are_this_ID]), collapse = "; ")
        } else {
          my_result <- NA_character_
        }
        rep.int(my_result, length(x))
      }
    )
    merged_CRISPR_df[["TSS_regions"]] <- unlist(regions_per_gene_list)

    ### Create a column that indicates whether a given genomic location was within the TSS regions searched by GuideScan
    searched_by_GuideScan_list <- tapply(
      seq_len(nrow(merged_CRISPR_df)),
      IDs_fac,
      function(x) {
        have_location <- !(is.na(merged_CRISPR_df[["Start"]][x]))
        are_this_ID <- TSS_ranges_df[["Combined_ID"]] == merged_CRISPR_df[["Combined_ID"]][[x[[1]]]]
        results_vec <- rep.int(NA_character_, length(x))
        if (any(have_location) && any(are_this_ID) && (any(merged_CRISPR_df[["Chromosome"]][x[have_location]] %in% TSS_ranges_df[["Chromosome"]][are_this_ID]))) {
          GRanges_object_sgRNAs <- RangesDfToGRangesObject(merged_CRISPR_df[x[have_location], ])
          GRanges_object_GuideScan <- RangesDfToGRangesObject(TSS_ranges_df[are_this_ID, ])
          overlap_Hits_object <- findOverlaps(GRanges_object_sgRNAs, GRanges_object_GuideScan, ignore.strand = TRUE)
          were_searched <- seq_len(sum(have_location)) %in% queryHits(overlap_Hits_object)
          results_vec[have_location] <- ifelse(were_searched, "Yes", "No")
        }
        return(results_vec)

      },
      simplify = FALSE
    )
    merged_CRISPR_df[["TSS_searched_by_GuideScan"]] <- unlist(searched_by_GuideScan_list)

    ### Further refine the data in the "TSS_searched_by_GuideScan" to reflect whether the sgRNA was located within the TSS search window for specifically this gene
    have_location <- !(is.na(merged_CRISPR_df[["Start"]]))
    GRanges_object_sgRNAs <- RangesDfToGRangesObject(merged_CRISPR_df[have_location, ])
    GRanges_object_GuideScan <- RangesDfToGRangesObject(TSS_ranges_df)
    overlap_Hits_object <- findOverlaps(GRanges_object_sgRNAs, GRanges_object_GuideScan, ignore.strand = TRUE)
    were_searched <- seq_len(sum(have_location)) %in% queryHits(overlap_Hits_object)

    submitted_to_GuideScan_vec <- rep.int(NA_character_, nrow(merged_CRISPR_df))
    submitted_to_GuideScan_vec[have_location] <- ifelse(were_searched, "Yes", "No")

    are_conflicting <- (merged_CRISPR_df[["TSS_searched_by_GuideScan"]] %in% "No") & (submitted_to_GuideScan_vec %in% "Yes")
    merged_CRISPR_df[["TSS_searched_by_GuideScan"]][are_conflicting] <- "Not this gene"

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################

    merged_CRISPR_df[["GuideScan_offtarget_category"]] <- GetOffTargetCategory(merged_CRISPR_df)

  }



  # Remove unnecessary columns
  remove_columns <- c("Hits_chromosome", "GuideScan_chromosome", "Hits_strand", "GuideScan_strand", "Hits_start",
                      "GuideScan_start", "Hits_end", "GuideScan_end",
                      "GuideScan_entrez_ID", "GuideScan_symbol", "gRNA", "gRNA_label", "Annotation"
                      )
  results_df <- merged_CRISPR_df[, !(colnames(merged_CRISPR_df) %in% remove_columns)]

  return(results_df)
}





# Miscellaneous functions -------------------------------------------------

CheckForInconsistentChromosomes <- function(CRISPR_df) {
  have_more_than_one_chr_vec <- sapply(
    unique(CRISPR_df[["Combined_ID"]]),
    function(x) {
      chromosomes_vec <- CRISPR_df[CRISPR_df[["Combined_ID"]] == x, "Chromosome"]
      chromosomes_vec <- chromosomes_vec[!(is.na(chromosomes_vec))]
      more_than_one <- length(unique(chromosomes_vec)) > 1
      return(more_than_one)
    })
  results_vec <- names(which(have_more_than_one_chr_vec))
  return(results_vec)
}




FindBest0MMLocations <- function(CRISPR_df, parallel_mode = TRUE, num_cores = NULL) {

  assign("original_CRISPR_df", CRISPR_df, envir = globalenv())

  use_TSS <- "Best_TSS" %in% names(CRISPR_df)
  are_not_controls <- CRISPR_df[["Is_control"]] %in% "No"

  combined_IDs_vec <- CRISPR_df[["Combined_ID"]][are_not_controls]
  combined_IDs_fac <- factor(combined_IDs_vec, levels = unique(combined_IDs_vec))

  # The following check is optional:
  stopifnot(identical(length(unique(combined_IDs_fac)),
                      length(rle(as.integer(combined_IDs_fac))[["lengths"]])
                      )
            )

  have_any_0MM <- tapply(CRISPR_df[["Num_0MM"]][are_not_controls] > 1,
                         combined_IDs_fac,
                         any
                         )

  split_df_list <- split(CRISPR_df[are_not_controls, ], combined_IDs_fac)

  strict_df <- CRISPR_df[, location_columns]
  colnames(strict_df) <- paste0(colnames(strict_df), "_strict")

  if (parallel_mode) {
    use_varlist <- c("split_df_list", "have_any_0MM", "use_TSS", "location_columns",
                     "Choose0MMLocation", "LocationStringToDf", "NoMatchForChromosome"
                     )
    if (use_TSS) {
      use_varlist <- c(use_varlist,
                       "GetBestTSSPositions", "FindSingleLocationFromTSS", "GetCutLocations"
                       )
    } else {
      use_varlist <- c(use_varlist,
                       "FindSingleLocationFromGene", "RangesDfToGRangesObject",
                       "human_genes_GRanges", "mcols", "seqnames",
                       "distanceToNearest", "GRanges", "IRanges"
                       )
    }

    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores() - 2
    }
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterExport(cl, varlist = use_varlist, envir = environment())
    lax_df_list <- parallel::parLapply(cl,
                                       split_df_list[have_any_0MM],
                                       function(x) Choose0MMLocation(x, use_TSS)
                                       )
    parallel::stopCluster(cl)
  } else {
    lax_df_list <- lapply(split_df_list[have_any_0MM], function(x) Choose0MMLocation(x, use_TSS))
  }

  lax_df_list <- lapply(split_df_list[have_any_0MM], function(x) Choose0MMLocation(x, use_TSS))
  split_df_list[have_any_0MM] <- lax_df_list

  assign("delete_split_df_list",    split_df_list,    envir = globalenv())
  assign("delete_are_not_controls", are_not_controls, envir = globalenv())
  assign("delete_CRISPR_df",        CRISPR_df,        envir = globalenv())

  results_df <- do.call(rbind.data.frame,
                        c(split_df_list,
                          list(CRISPR_df[!(are_not_controls), ],
                               stringsAsFactors = FALSE, make.row.names = FALSE
                        )))

  assign("delete_results_df", results_df, envir = globalenv())
  assign("delete_strict_df", strict_df, envir = globalenv())


  results_df <- cbind.data.frame(results_df, strict_df)
  results_df[["PAM"]] <- ifelse(is.na(results_df[["PAM"]]),
                                results_df[["PAM_strict"]],
                                results_df[["PAM"]]
                                )

  stopifnot(identical(nrow(results_df), nrow(CRISPR_df)))
  return(results_df)
}























