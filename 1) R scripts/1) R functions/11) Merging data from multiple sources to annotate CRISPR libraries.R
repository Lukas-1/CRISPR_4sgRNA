### 9th September 2019 ####




# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Functions for integrating with a genome search and GuideScan ------------

ExtendWithGenomeSearch <- function(CRISPR_df, search_df, allow_5pG = FALSE) {
  genome_matches <- match(toupper(CRISPR_df[, "sgRNA_sequence"]), toupper(search_df[, "Sequence"]))
  results_df <- data.frame(
    CRISPR_df,
    search_df[genome_matches, !(colnames(search_df) %in% c("Sequence", "Found"))],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  stopifnot(!(anyNA(results_df[, "Num_0MM"])))
  have_hit <- results_df[, "Num_0MM"] > 0
  if (allow_5pG) {
    have_hit <- ifelse(grepl("hCRISPRa-v2", CRISPR_df[, "Source"], fixed = TRUE) & (CRISPR_df[, "Is_control"] == "No"),
                       have_hit | ((results_df[, "Num_0MM"] == 0) & (results_df[, "Num_5G_MM"] > 0)),
                       have_hit
                       )
  }
  for (column_name in c("Chromosome", "Strand", "Start", "End")) {
    results_df[, column_name] <- ifelse(have_hit, results_df[, column_name], NA)
    colnames(results_df)[colnames(results_df) == column_name] <- paste0("Hits_", tolower(column_name))
  }
  return(results_df)
}



EntrezIDsToChromosomes <- function(entrez_vec) {
  unique_entrez_ID_strings <- unique(entrez_vec)
  expanded_entrez_IDs_df <- ExpandList(strsplit(unique_entrez_ID_strings, ", ", fixed = TRUE))
  expanded_entrez_IDs_df[, "Chromosome"] <- entrez_to_symbol_df[match(expanded_entrez_IDs_df[, "Value"], entrez_to_symbol_df[, "Entrez_ID"]), "Chromosome"]
  chromosomes_list <- split(expanded_entrez_IDs_df[, "Chromosome"], expanded_entrez_IDs_df[, "List_index"])
  chromosomes_list <- lapply(chromosomes_list, unique)
  chromosomes_vec <- vapply(chromosomes_list, function(x) if (length(x) > 1) NA_character_ else x, "")
  results_vec <- chromosomes_vec[match(entrez_vec, unique_entrez_ID_strings)]
  return(results_vec)
}



MergeTSSandGuideScan <- function(CRISPR_df, guidescan_df) {

  CRISPR_df[, "Entrez_chromosome"] <- EntrezIDsToChromosomes(CRISPR_df[, "Entrez_ID"])

  combined_TSS_CRISPRa_df[, "GuideScan_input"] <- TSSStringForGuideScan(combined_TSS_CRISPRa_df)
  tidy_guidescan_df <- TidyGuideScanColumns(guidescan_df)

  TSS_columns <- c("Strand", "Best_TSS", "First_TSS", "Last_TSS")

  IDs_per_region_list <- sapply(unique(tidy_guidescan_df[, "Region"]), function(x) {
    unique(combined_TSS_CRISPRa_df[combined_TSS_CRISPRa_df[, "GuideScan_input"] == x, c("Combined_ID", "Entrez_ID", "Gene_symbol", TSS_columns)])
  }, simplify = FALSE)
  guidescan_per_region_list <- split(tidy_guidescan_df, factor(tidy_guidescan_df[, "Region"], levels = unique(tidy_guidescan_df[, "Region"])))
  guidescan_per_region_extended_list <- sapply(names(guidescan_per_region_list), function(x) {
    num_reps_df <- IDs_per_region_list[[x]]
    gs_df <- guidescan_per_region_list[[x]]
    sub_list <- lapply(seq_len(nrow(num_reps_df)), function(x) data.frame(num_reps_df[rep.int(x, nrow(gs_df)), ], gs_df))
    do.call(rbind.data.frame, c(sub_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  }, simplify = FALSE)
  extended_guidescan_df <- do.call(rbind.data.frame, c(guidescan_per_region_extended_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  vec_for_matching_CRISPR               <- paste0(CRISPR_df[, "Combined_ID"], "__", toupper(CRISPR_df[, "sgRNA_sequence"]))
  vec_for_matching_GuideScan            <- paste0(extended_guidescan_df[, "Combined_ID"], "__", toupper(extended_guidescan_df[, "gRNA"]))
  vec_for_matching_CRISPR_chromosome    <- paste0(CRISPR_df[, "Entrez_chromosome"], "__", vec_for_matching_CRISPR)
  vec_for_matching_GuideScan_chromosome <- paste0(extended_guidescan_df[, "GuideScan_chromosome"], "__", vec_for_matching_GuideScan)

  matches_vec <- ifelse(is.na(CRISPR_df[, "Entrez_chromosome"]),
                        match(vec_for_matching_CRISPR, vec_for_matching_GuideScan),
                        match(vec_for_matching_CRISPR_chromosome, vec_for_matching_GuideScan_chromosome)
                        )

  colnames(extended_guidescan_df)[colnames(extended_guidescan_df) == "Entrez_ID"]   <- "GuideScan_entrez_ID"
  colnames(extended_guidescan_df)[colnames(extended_guidescan_df) == "Gene_symbol"] <- "GuideScan_symbol"
  colnames(extended_guidescan_df)[colnames(extended_guidescan_df) == "Strand"]      <- "Strand_of_TSS"
  colnames(extended_guidescan_df)[colnames(extended_guidescan_df) == "Region"]      <- "GuideScan_region"

  combined_df <- CRISPR_df
  for (column_name in colnames(extended_guidescan_df)[colnames(extended_guidescan_df) != "Combined_ID"]) {
    combined_df[, column_name] <- extended_guidescan_df[matches_vec, column_name]
  }
  return(combined_df)
}





# Functions for tidying and summarizing the data from GuideScan -----------

TidyGuideScanColumns <- function(guidescan_df) {
  results_df <- data.frame(
    guidescan_df[, c("Region", "gRNA")],
    "GuideScan_chromosome"  = guidescan_df[, "chromosome"],
    "GuideScan_strand"      = guidescan_df[, "strand"],
    "GuideScan_start"       = as.integer(guidescan_df[, "target_site_start_coordinate"]),
    "GuideScan_end"         = as.integer(guidescan_df[, "target_site_end_coordinate"]),
    "GuideScan_efficiency"  = as.numeric(ifelse(guidescan_df[, "cutting_efficiency_score"] == "*", NA_character_, guidescan_df[, "cutting_efficiency_score"])),
    "GuideScan_specificity" = as.numeric(guidescan_df[, "cutting_specificity_score"]),
    SplitOffTargetsSummary(guidescan_df[, "offtargets_summary"]),
    "GuideScan_Num_2or3MM"  = as.integer(guidescan_df[, "offtargets_sum"]),
    "Annotation"            = guidescan_df[, "annotation"],
    guidescan_df[, "gRNA_label", drop = FALSE],
    stringsAsFactors        = FALSE,
    row.names               = NULL
  )

  # Make the GuideScan locations consistent with the locations returned by Biostrings::matchPattern
  results_df[, "GuideScan_start"] <- ifelse(results_df[, "GuideScan_strand"] %in% "-",
                                            results_df[, "GuideScan_start"] + 3L,
                                            results_df[, "GuideScan_start"]
                                            )

  results_df[, "GuideScan_end"] <- ifelse(results_df[, "GuideScan_strand"] %in% "-",
                                          results_df[, "GuideScan_end"] + 1L,
                                          results_df[, "GuideScan_end"] - 2L
                                          )

  return(results_df)
}



SplitOffTargetsSummary <- function(off_targets_summary_vec) {
  offtargets_summary_splits <- strsplit(off_targets_summary_vec, "|", fixed = TRUE)
  results_df <- data.frame(
    "GuideScan_Num_2MM" = vapply(offtargets_summary_splits, function(x) if (all(is.na(x))) NA_integer_ else as.integer(sub("2:", "", x[[1]], fixed = TRUE)), integer(1)),
    "GuideScan_Num_3MM" = vapply(offtargets_summary_splits, function(x) if (all(is.na(x))) NA_integer_ else as.integer(sub("3:", "", x[[2]], fixed = TRUE)), integer(1)),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  return(results_df)
}



GetOffTargetCategory <- function(merged_CRISPR_df) {
  if ("Num_5G_MM" %in% colnames(merged_CRISPR_df)) {
    num_perfect_hits_vec <- rowSums(merged_CRISPR_df[, c("Num_0MM", "Num_5G_MM")])
  } else {
    num_perfect_hits_vec <- merged_CRISPR_df[, "Num_0MM"]
  }
  is_unspecific_vec <- (num_perfect_hits_vec > 1) | (merged_CRISPR_df[, "Num_1MM"] > 0)
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
    num_2or3MM <- merged_CRISPR_df[x, "GuideScan_Num_2or3MM"]
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

  if ("TSS_searched_by_GuideScan" %in% colnames(merged_CRISPR_df)) {
    were_not_scanned <- (results_fac %in% "Unknown") & (merged_CRISPR_df[, "TSS_searched_by_GuideScan"] %in% c("No", "Not this gene"))
    results_fac[were_not_scanned] <- "Not scanned"
  }
  return(results_fac)
}





# Functions for resolving ambiguous gene IDs or sgRNA locations -----------

AssignToGeneByNearbyTSS <- function(CRISPR_df, prefix = "") {
  # Depends on the data frame 'combined_TSS_CRISPRa_df' in the global environment

  GRanges_object_TSSs <- GRanges(
    seqnames = sub("chr", "", combined_TSS_CRISPRa_df[, "Chromosome"], fixed = TRUE),
    ranges   = IRanges(start = combined_TSS_CRISPRa_df[, "First_TSS"] - 1500L, end = combined_TSS_CRISPRa_df[, "Last_TSS"] + 1500L),
    strand   = combined_TSS_CRISPRa_df[, "Strand"]
  )
  combined_IDs <- unique(CRISPR_df[CRISPR_df[, "Is_control"] == "No", "Combined_ID"])

  new_entrezs_vec <- rep.int(NA_character_, nrow(CRISPR_df))
  assignment_vec <- rep.int(NA_character_, nrow(CRISPR_df))

  for (combined_ID in combined_IDs) {
    are_this_ID <- CRISPR_df[, "Combined_ID"] == combined_ID
    sub_df <- CRISPR_df[are_this_ID, ]
    mapped_entrezs <- unique(sub_df[, "GuideScan_entrez_ID"])
    mapped_entrezs <- mapped_entrezs[!(is.na(mapped_entrezs))]
    if (length(mapped_entrezs) == 1) {
      assignment <- "mapped using GuideScan"
    } else if (length(mapped_entrezs) == 0) {
      overlapping_entrezs_list <- lapply(which(!(is.na(sub_df[, "Start"]))), function(x) {
        this_GRange_object <- GRanges(
          seqnames = sub("chr", "", sub_df[x, "Chromosome"], fixed = TRUE),
          ranges   = IRanges(start = sub_df[x, "Start"], end = sub_df[x, "End"]),
          strand   = sub_df[x, "Strand"]
        )
        overlap_matches_df <- as.data.frame(findOverlaps(this_GRange_object, GRanges_object_TSSs, ignore.strand = TRUE))
        matched_entrezs <- unique(combined_TSS_CRISPRa_df[overlap_matches_df[, 2], "Entrez_ID"])
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
  ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] - 2L, ranges_df[, "Start"] + 3L)
}


LocationStringToDf <- function(location_char_vec) {

  chromosome_splits <- strsplit(location_char_vec, "(", fixed = TRUE)
  strand_splits <- strsplit(sapply(chromosome_splits, "[", 2), ")", fixed = TRUE)
  location_vec <- sapply(strand_splits, "[", 2)
  location_vec <- substr(location_vec, 2, nchar(location_vec))
  location_splits <- strsplit(location_vec, "-", fixed = TRUE)

  results_df <- data.frame(
    "Chromosome" = sapply(chromosome_splits, "[", 1),
    "Strand"     = sapply(strand_splits, "[", 1),
    "Start"      = as.integer(sapply(location_splits, "[", 1)),
    "End"        = as.integer(sapply(location_splits, "[", 2)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}


FindSingleLocationFromTSS <- function(location_vec, chromosome, TSS_location = NA) {
  results_df <- LocationStringToDf(location_vec)
  if (!(is.na(chromosome))) {
    are_this_chromosome <- results_df[, "Chromosome"] == chromosome
    stopifnot(any(are_this_chromosome))
    results_df <- results_df[are_this_chromosome, ]
  }
  if (is.na(TSS_location)) {
    results_df <- results_df[1, ]
  } else {
    cut_location_vec <- GetCutLocations(results_df)
    TSS_distances <- abs(TSS_location - cut_location_vec)
    results_df <- results_df[which.min(TSS_distances), ]
  }
  return(results_df)
}


FindSingleLocationFromGene <- function(location_vec, chromosome, gene_GRanges_object) {
  results_df <- LocationStringToDf(location_vec)
  if (!(is.na(chromosome))) {
    are_this_chromosome <- results_df[, "Chromosome"] == chromosome
    stopifnot(any(are_this_chromosome))
    results_df <- results_df[are_this_chromosome, ]
  }
  if (is.null(gene_GRanges_object)) {
    results_df <- results_df[1, ]
  } else {
    sgRNAs_GRanges_object <- GRanges(
      seqnames = results_df[, "Chromosome"],
      ranges   = IRanges(start = results_df[, "Start"], end = results_df[, "End"]),
      strand   = results_df[, "Strand"]
    )
    distance_vec <- mcols(distanceToNearest(sgRNAs_GRanges_object, gene_GRanges_object, ignore.strand = TRUE, select = "all"))[, 1]
    if (length(distance_vec) == 0) {
      results_df <- results_df[1, ] # i.e. gene_GRanges_object is on a different chromosome
    } else {
      stopifnot(length(distance_vec) == nrow(results_df))
      results_df <- results_df[which.min(distance_vec), ]
    }
  }
  return(results_df)
}



FindDuplicatedGenes <- function(CRISPR_df, fraction_cutoff = 0.8, absolute_cutoff = 20, min_to_map = 8) {
  # Depends on the gene models object 'human_genes_GRanges' in the global environment

  use_TSS <- "Best_TSS" %in% colnames(CRISPR_df)

  are_not_controls <- CRISPR_df[, "Is_control"] %in% "No"
  have_two_0MM <- CRISPR_df[, "Num_0MM"] == 2
  combined_IDs <- unique(CRISPR_df[are_not_controls & have_two_0MM, "Combined_ID"])

  results_df <- data.frame(
    "Is_duplicated_gene" = rep.int(FALSE, nrow(CRISPR_df)),
    "Chromosome"         = rep.int(NA_character_, nrow(CRISPR_df)),
    "Strand"             = rep.int(NA_character_, nrow(CRISPR_df)),
    "Start"              = rep.int(NA_integer_,   nrow(CRISPR_df)),
    "End"                = rep.int(NA_integer_,   nrow(CRISPR_df)),
    stringsAsFactors     = FALSE
  )

  for (combined_ID in combined_IDs) {

    are_this_gene <- CRISPR_df[, "Combined_ID"] == combined_ID
    this_gene_duplicates <- are_this_gene & have_two_0MM
    splits_list <- strsplit(CRISPR_df[this_gene_duplicates, "Locations_0MM"], "; ", fixed = TRUE)
    assign("delete_are_this_gene", are_this_gene, envir = globalenv())
    assign("original_splits_list", splits_list, envir = globalenv())
    assign("original_this_gene_duplicates", this_gene_duplicates, envir = globalenv())
    split_1_vec <- sapply(splits_list, "[", 1)
    split_2_vec <- sapply(splits_list, "[", 2)
    chromosome_1_vec <- sapply(strsplit(split_1_vec, "(", fixed = TRUE), "[", 1)
    chromosome_2_vec <- sapply(strsplit(split_2_vec, "(", fixed = TRUE), "[", 1)
    chromosome_combo_vec <- paste0(chromosome_1_vec, "_", chromosome_2_vec)
    most_common_combo <- names(sort(table(chromosome_combo_vec), decreasing = TRUE))[[1]]
    this_gene_duplicates[this_gene_duplicates] <- chromosome_combo_vec == most_common_combo
    is_duplicated_gene <- ((sum(this_gene_duplicates) / sum(are_this_gene)) > fraction_cutoff)

    if (!(is_duplicated_gene) && (sum(this_gene_duplicates) > 20)) {
      already_mapped <- !(is.na(CRISPR_df[are_this_gene, "Start"]))
      if (sum(already_mapped) < min_to_map) {
        gene_symbol <- unique(CRISPR_df[are_this_gene, "Gene_symbol"])
        if (length(gene_symbol) == 1) {
          gene_string <- paste0("symbol '", gene_symbol, "'")
        } else {
          gene_string <- paste0("combined ID '", combined_ID, "'")
        }
        message(paste0("The gene with the ", gene_string, " was counted as a duplicated gene, ",
                       "even though it did not reach the fraction cutoff of ", fraction_cutoff, ", ",
                       "because >", absolute_cutoff, " of its guides seem to be duplicated."
                       ))
        is_duplicated_gene <- TRUE
      }
    }

    if (is_duplicated_gene) {

      results_df[this_gene_duplicates, "Is_duplicated_gene"] <- TRUE

      splits_list <- strsplit(CRISPR_df[this_gene_duplicates, "Locations_0MM"], "; ", fixed = TRUE) # The split has to be done again, since this_gene_duplicates may have changed

      assign("delete_this_gene_duplicates", this_gene_duplicates, envir = globalenv())

      entrez_chromosome <- unique(CRISPR_df[this_gene_duplicates, "Entrez_chromosome"])
      if (entrez_chromosome %in% "chrX, chrY") {
        entrez_chromosome <- "chrX"
      }

      if (use_TSS) {

        best_TSS <- unique(GetBestTSSPositions(CRISPR_df[this_gene_duplicates, ]))

        assign("delete_combined_ID", combined_ID, envir = globalenv())
        assign("delete_this_gene_duplicates", this_gene_duplicates, envir = globalenv())
        assign("delete_best_TSS", best_TSS, envir = globalenv())
        assign("delete_entrez_chromosome", entrez_chromosome, envir = globalenv())
        assign("delete_CRISPR_sub_df", CRISPR_df[this_gene_duplicates, ], envir = globalenv())

        assign("delete_splits_list", splits_list, envir = globalenv())
        new_locations_df <- do.call(rbind.data.frame,
                                    c(lapply(splits_list, function(x) FindSingleLocationFromTSS(x, entrez_chromosome, best_TSS)),
                                      list(stringsAsFactors = FALSE, make.row.names = FALSE)
                                      )
                                    )

      } else {

        entrez_ID <- unique(CRISPR_df[this_gene_duplicates, "Entrez_ID"])

        if (is.na(entrez_ID)) {
          gene_GRanges_object <- NULL
        } else {
          are_this_entrez <- (mcols(human_genes_GRanges)[, "gene_id"] == entrez_ID) &
                             (seqnames(human_genes_GRanges) == entrez_chromosome)
          if (any(are_this_entrez)) {
            gene_GRanges_object <- human_genes_GRanges[mcols(human_genes_GRanges)[, "gene_id"] == entrez_ID]
          } else {
            gene_GRanges_object <- NULL
          }
        }

        new_locations_df <- do.call(rbind.data.frame,
                                    c(lapply(splits_list, function(x) FindSingleLocationFromGene(x, entrez_chromosome, gene_GRanges_object)),
                                      list(stringsAsFactors = FALSE, make.row.names = FALSE)
                                      )
                                    )

        assign("delete_splits_list", splits_list, envir = globalenv())
        assign("delete_new_locations_df", new_locations_df, envir = globalenv())
        assign("delete_entrez_ID", entrez_ID, envir = globalenv())
        assign("delete_gene_GRanges_object", gene_GRanges_object, envir = globalenv())

      }

      for (column_name in colnames(new_locations_df)) {
        results_df[this_gene_duplicates, column_name] <- new_locations_df[, column_name]
      }

    }
  }

  return(results_df)
}


ReplaceDuplicatedGenes <- function(CRISPR_df) {

  duplicated_genes_df <- FindDuplicatedGenes(CRISPR_df)
  are_to_replace <- duplicated_genes_df[, "Is_duplicated_gene"] & !(is.na(duplicated_genes_df[, "Start"]))

  location_columns <- c("Chromosome", "Strand", "Start", "End")

  for (column in location_columns) {
    CRISPR_df[are_to_replace, column] <- duplicated_genes_df[are_to_replace, column]
  }
  message(paste0(sum(are_to_replace), " sgRNAs from ", length(unique(CRISPR_df[are_to_replace, "Combined_ID"])),
                 ' duplicated genes were re-assigned to a single location.'
                 ))

  return(CRISPR_df)
}


# Functions for refining the genomic locations of sgRNAs ------------------


ReassignTSS <- function(merged_CRISPR_df) {
  # Assigns a TSS to all sgRNAs based on the gene ID they are associated with (also to those sgRNAs which were not found on GuideScan!)
  # Depends on the data frame 'combined_TSS_CRISPRa_df' in the global environment

  vec_for_matching_CRISPR <- paste0(merged_CRISPR_df[, "Combined_ID"], "__", merged_CRISPR_df[, "Chromosome"])
  vec_for_matching_TSS <- paste0(combined_TSS_CRISPRa_df[, "Combined_ID"], "__", combined_TSS_CRISPRa_df[, "Chromosome"])
  matches_vec <- match(vec_for_matching_CRISPR, vec_for_matching_TSS)
  for (column_name in c("Strand_of_TSS", "Best_TSS", "First_TSS", "Last_TSS")) {
    TSS_column_names <- c("Strand_of_TSS" = "Strand")
    if (column_name %in% names(TSS_column_names)) {
      TSS_column_name <- TSS_column_names[[column_name]]
    } else {
      TSS_column_name <- column_name
    }
    merged_CRISPR_df[, column_name] <- ifelse(is.na(merged_CRISPR_df[, column_name]),
                                              ifelse(is.na(merged_CRISPR_df[, "Chromosome"]), NA, combined_TSS_CRISPRa_df[matches_vec, TSS_column_name]),
                                              merged_CRISPR_df[, column_name]
                                              )
  }
  return(merged_CRISPR_df)
}




GetBestTSSPositions <- function(CRISPR_df) {
  ifelse(is.na(CRISPR_df[, "Best_TSS"]),
         ifelse(CRISPR_df[, "Strand_of_TSS"] == "+", CRISPR_df[, "First_TSS"], CRISPR_df[, "Last_TSS"]),
         CRISPR_df[, "Best_TSS"]
         )
}



MergeLocations <- function(merged_CRISPR_df) {

  ### Check for discrepancies between various methods of finding the location of an sgRNA
  discordant_mat <- cbind(
    "Chromosome_entrez_hits"      = !(is.na(merged_CRISPR_df[, "Entrez_chromosome"])) & !(is.na(merged_CRISPR_df[, "Hits_chromosome"]))      & (merged_CRISPR_df[, "Entrez_chromosome"] != merged_CRISPR_df[, "Hits_chromosome"]),
    "Chromosome_entrez_GuideScan" = !(is.na(merged_CRISPR_df[, "Entrez_chromosome"])) & !(is.na(merged_CRISPR_df[, "GuideScan_chromosome"])) & (merged_CRISPR_df[, "Entrez_chromosome"] != merged_CRISPR_df[, "GuideScan_chromosome"]),
    "Chromosome_hits_GuideScan"   = !(is.na(merged_CRISPR_df[, "Hits_chromosome"]))   & !(is.na(merged_CRISPR_df[, "GuideScan_chromosome"])) & (merged_CRISPR_df[, "Hits_chromosome"]   != merged_CRISPR_df[, "GuideScan_chromosome"]),
    "Strand"                      = !(is.na(merged_CRISPR_df[, "Hits_strand"]))       & !(is.na(merged_CRISPR_df[, "GuideScan_strand"]))     & (merged_CRISPR_df[, "Hits_strand"]       != merged_CRISPR_df[, "GuideScan_strand"]),
    "Start"                       = !(is.na(merged_CRISPR_df[, "Hits_start"]))        & !(is.na(merged_CRISPR_df[, "GuideScan_start"]))      & (merged_CRISPR_df[, "Hits_start"] != merged_CRISPR_df[, "GuideScan_start"]),
    "End"                         = !(is.na(merged_CRISPR_df[, "Hits_end"]))          & !(is.na(merged_CRISPR_df[, "GuideScan_end"]))        & (merged_CRISPR_df[, "Hits_end"]   != merged_CRISPR_df[, "GuideScan_end"])
  )
  merged_CRISPR_df[, "Is_discordant"] <- rowSums(discordant_mat) >= 1
  is_discordant_location <- rowSums(discordant_mat[, c("Chromosome_hits_GuideScan", "Strand", "Strand", "Start", "End")]) >= 1
  stopifnot(all(!(is_discordant_location), na.rm = TRUE))


  ### Assign each sgRNA to a location
  were_hits <- !(is.na(merged_CRISPR_df[, "Hits_start"]))
  stopifnot(identical(is.na(merged_CRISPR_df[, "Hits_start"]), is.na(merged_CRISPR_df[, "Hits_end"]))) ### DELETE THIS ###
  stopifnot(identical(is.na(merged_CRISPR_df[, "Hits_start"]), is.na(merged_CRISPR_df[, "Hits_strand"]))) ### DELETE THIS ###
  stopifnot(identical(is.na(merged_CRISPR_df[, "Hits_start"]), is.na(merged_CRISPR_df[, "Hits_chromosome"]))) ### DELETE THIS ###

  merged_CRISPR_df[, "Chromosome"] <- ifelse(were_hits,
                                             merged_CRISPR_df[, "Hits_chromosome"],
                                             ifelse(!(is.na(merged_CRISPR_df[, "Entrez_chromosome"])),
                                                    merged_CRISPR_df[, "Entrez_chromosome"],
                                                    merged_CRISPR_df[, "Entrez_chromosome"]
                                                    )

                                             )
  merged_CRISPR_df[, "Strand"] <- ifelse(were_hits, merged_CRISPR_df[, "Hits_strand"], merged_CRISPR_df[, "GuideScan_strand"])
  merged_CRISPR_df[, "Start"]  <- ifelse(were_hits, merged_CRISPR_df[, "Hits_start"],  merged_CRISPR_df[, "GuideScan_start"])
  merged_CRISPR_df[, "End"]    <- ifelse(were_hits, merged_CRISPR_df[, "Hits_end"],    merged_CRISPR_df[, "GuideScan_end"])

  merged_CRISPR_df <- ReassignTSS(merged_CRISPR_df)

  return(merged_CRISPR_df)
}





AdjustPositionColumns <- function(merged_CRISPR_df, guidescan_df, reorder_by_rank = TRUE, allow_5pG_MM = TRUE) {
  # Depends on the data frame 'combined_TSS_CRISPRa_df' in the global environment

  # Prepare for re-ordering the columns in a later step
  remove_columns <- c("Entrez_chromosome", "Hits_chromosome", "GuideScan_chromosome", "Hits_strand", "GuideScan_strand", "Hits_start",
                      "GuideScan_start", "Hits_end", "GuideScan_end"
                      )

  # Assign sgRNAs to their genomic locations
  merged_CRISPR_df <- MergeLocations(merged_CRISPR_df)

  # Resolve the issue of sgRNAs that could not be assigned to a single Entrez ID based on gene annotation alone (i.e. the gene symbol or Entrez ID provided)
  # (These will be assigned to a gene whose TSS is located within 1500 bp of the sgRNA)
  are_ambiguous <- grepl(",", merged_CRISPR_df[, "Entrez_ID"], fixed = TRUE)
  are_NA <- is.na(merged_CRISPR_df[, "Entrez_ID"])
  merged_CRISPR_df[, "Entrez_ID_assignment"] <- ifelse(!(are_ambiguous | are_NA), "Unambiguous", NA_character_)
  if (any(are_ambiguous)) {
    ambiguous_df <- AssignToGeneByNearbyTSS(merged_CRISPR_df[are_ambiguous, ], prefix = "The gene symbol was ambiguous; ")
    assign("delete_ambiguous_df", ambiguous_df, envir = globalenv())
    merged_CRISPR_df[are_ambiguous, "Entrez_ID"]            <- ifelse(is.na(ambiguous_df[, "Entrez_ID"]), merged_CRISPR_df[are_ambiguous, "Entrez_ID"],   ambiguous_df[, "Entrez_ID"])
    merged_CRISPR_df[are_ambiguous, "Gene_symbol"]          <- ifelse(is.na(ambiguous_df[, "Entrez_ID"]), merged_CRISPR_df[are_ambiguous, "Gene_symbol"], ambiguous_df[, "Gene_symbol"])
    merged_CRISPR_df[are_ambiguous, "Combined_ID"]          <- ifelse(is.na(ambiguous_df[, "Entrez_ID"]), merged_CRISPR_df[are_ambiguous, "Combined_ID"], ambiguous_df[, "Entrez_ID"])
    merged_CRISPR_df[are_ambiguous, "Entrez_ID_assignment"] <- ambiguous_df[, "Assignment"]
  }
  if (any(are_NA)) {
    NA_df <- AssignToGeneByNearbyTSS(merged_CRISPR_df[are_NA, ], prefix = "No Entrez ID was found for the gene symbol; ")
    assign("delete_NA_df", NA_df, envir = globalenv())
    merged_CRISPR_df[are_NA, "Entrez_ID"]                   <- ifelse(is.na(NA_df[, "Entrez_ID"]), merged_CRISPR_df[are_NA, "Entrez_ID"],   NA_df[, "Entrez_ID"])
    merged_CRISPR_df[are_NA, "Gene_symbol"]                 <- ifelse(is.na(NA_df[, "Entrez_ID"]), merged_CRISPR_df[are_NA, "Gene_symbol"], NA_df[, "Gene_symbol"])
    merged_CRISPR_df[are_NA, "Combined_ID"]                 <- ifelse(is.na(NA_df[, "Entrez_ID"]), merged_CRISPR_df[are_NA, "Combined_ID"], NA_df[, "Entrez_ID"])
    merged_CRISPR_df[are_NA, "Entrez_ID_assignment"]        <- NA_df[, "Assignment"]
  }

  # After assigning Entrez IDs to more sgRNAs (in the previous step), perform another merge with the data from GuideScan
  remerged_CRISPR_df <- MergeTSSandGuideScan(merged_CRISPR_df, guidescan_df)
  for (my_column in colnames(remerged_CRISPR_df)) {
    merged_CRISPR_df[, my_column] <- remerged_CRISPR_df[, my_column]
  }
  merged_CRISPR_df <- MergeLocations(merged_CRISPR_df)



  # For genes that are clearly duplicated (e.g. those that are located in the pseudo-autosomal regions of chromosomes X and Y),
  # arbitrarily assign the sgRNAs to the "first" of the two locations in the genome
  merged_CRISPR_df <- ReplaceDuplicatedGenes(merged_CRISPR_df)

  merged_CRISPR_df <- ReassignTSS(merged_CRISPR_df)



  # Remove location data for sgRNAs that seem to have discrepant locations
  assign("delete_pre_discrepant", merged_CRISPR_df, envir = globalenv())
  for (column in c("Chromosome", "Strand", "Start", "End", "PAM")) {
    merged_CRISPR_df[merged_CRISPR_df[, "Is_discordant"], column] <- NA
  }
  assign("delete_post_discrepant", merged_CRISPR_df, envir = globalenv())


  ### Calculate the "cut" location
  merged_CRISPR_df[, "Cut_location"] <- GetCutLocations(merged_CRISPR_df)



  ### Calculate the distance from the TSS
  TSS_position_vec <- GetBestTSSPositions(merged_CRISPR_df)

  have_neg_strand_TSS <- merged_CRISPR_df[, "Strand_of_TSS"] == "-"
  merged_CRISPR_df[, "Distance_from_TSS"] <- (merged_CRISPR_df[, "Cut_location"] - TSS_position_vec) * ifelse(have_neg_strand_TSS, -1L, 1L) +
                                              ifelse(have_neg_strand_TSS, +2L, -1L)



  # Eliminate duplicated sgRNAs
  if (legacy_mode) {
    merged_CRISPR_df <- ResolveDuplicates(merged_CRISPR_df, concatenate_columns = c("Sublibrary", "hCRISPRa_v2_ID"))
  } else {
    merged_CRISPR_df <- ResolveDuplicates(merged_CRISPR_df, concatenate_columns = c("Sublibrary", "hCRISPRa_v2_ID", "hCRISPRa_TSS_source"))
  }

  ####################################################################################################################
  ### The following section of code attempts to troubleshoot the issue of sgRNAs that were not found by GuideScan. ###
  ### Is it because the regions where these sgRNAs are located were not submitted to GuideScan?                    ###
  ### Or were they submitted, but the sgRNA was not present in GuideScan's database?                               ###
  ####################################################################################################################

  TSS_ranges_df <- data.frame(combined_TSS_CRISPRa_df, TSSRangesForGuideScan(combined_TSS_CRISPRa_df), stringsAsFactors = FALSE, row.names = NULL)
  TSS_ranges_df[, "Region"] <- TSSStringForGuideScan(combined_TSS_CRISPRa_df)

  IDs_fac <- factor(merged_CRISPR_df[, "Combined_ID"], levels = unique(merged_CRISPR_df[, "Combined_ID"]))

  ### Create a new column that lists of all regions that were submitted to GuideScan for that gene
  regions_per_gene_list <- tapply(
    seq_len(nrow(merged_CRISPR_df)),
    IDs_fac,
    function(x) {
      are_this_ID <- TSS_ranges_df[, "Combined_ID"] == merged_CRISPR_df[x[[1]], "Combined_ID"]
      if (any(are_this_ID)) {
        my_result <- paste0(unique(TSS_ranges_df[are_this_ID, "Region"]), collapse = "; ")
      } else {
        my_result <- NA_character_
      }
      rep.int(my_result, length(x))
    }
  )
  merged_CRISPR_df[, "TSS_regions"] <- unlist(regions_per_gene_list)

  ### Create a column that indicates whether a given genomic location was within the TSS regions searched by GuideScan
  searched_by_GuideScan_list <- tapply(
    seq_len(nrow(merged_CRISPR_df)),
    IDs_fac,
    function(x) {
      have_location <- !(is.na(merged_CRISPR_df[x, "Start"]))
      are_this_ID <- TSS_ranges_df[, "Combined_ID"] == merged_CRISPR_df[x[[1]], "Combined_ID"]
      results_vec <- rep.int(NA_character_, length(x))
      if (any(have_location) && any(are_this_ID) && (any(merged_CRISPR_df[x[have_location], "Chromosome"] %in% TSS_ranges_df[are_this_ID, "Chromosome"]))) {
        GRanges_object_sgRNAs <- GRanges(
          seqnames = sub("chr", "", merged_CRISPR_df[x[have_location], "Chromosome"], fixed = TRUE),
          ranges   = IRanges(start = merged_CRISPR_df[x[have_location], "Start"], end = merged_CRISPR_df[x[have_location], "End"]),
          strand   = merged_CRISPR_df[x[have_location], "Strand"]
        )
        GRanges_object_GuideScan <- GRanges(
          seqnames = sub("chr", "", TSS_ranges_df[are_this_ID, "Chromosome"], fixed = TRUE),
          ranges   = IRanges(start = TSS_ranges_df[are_this_ID, "Start"], end = TSS_ranges_df[are_this_ID, "End"]),
          strand   = TSS_ranges_df[are_this_ID, "Strand"]
        )
        overlap_matches_df <- as.data.frame(findOverlaps(GRanges_object_sgRNAs, GRanges_object_GuideScan, ignore.strand = TRUE))
        were_searched <- seq_len(sum(have_location)) %in% overlap_matches_df[, "queryHits"]
        results_vec[have_location] <- ifelse(were_searched, "Yes", "No")
      }
      return(results_vec)

    },
    simplify = FALSE
  )
  merged_CRISPR_df[, "TSS_searched_by_GuideScan"] <- unlist(searched_by_GuideScan_list)

  ### Further refine the data in the "TSS_searched_by_GuideScan" to reflect whether the sgRNA was located within the TSS search window for specifically this gene
  have_location <- !(is.na(merged_CRISPR_df[, "Start"]))
  GRanges_object_sgRNAs <- GRanges(
    seqnames = sub("chr", "", merged_CRISPR_df[have_location, "Chromosome"], fixed = TRUE),
    ranges   = IRanges(start = merged_CRISPR_df[have_location, "Start"], end = merged_CRISPR_df[have_location, "End"]),
    strand   = merged_CRISPR_df[have_location, "Strand"]
  )
  GRanges_object_GuideScan <- GRanges(
    seqnames = sub("chr", "", TSS_ranges_df[, "Chromosome"], fixed = TRUE),
    ranges   = IRanges(start = TSS_ranges_df[, "Start"], end = TSS_ranges_df[, "End"]),
    strand   = TSS_ranges_df[, "Strand"]
  )
  overlap_matches_df <- as.data.frame(findOverlaps(GRanges_object_sgRNAs, GRanges_object_GuideScan, ignore.strand = TRUE))
  were_searched <- seq_len(sum(have_location)) %in% overlap_matches_df[, "queryHits"]

  submitted_to_GuideScan_vec <- rep.int(NA_character_, nrow(merged_CRISPR_df))
  submitted_to_GuideScan_vec[have_location] <- ifelse(were_searched, "Yes", "No")

  are_conflicting <- (merged_CRISPR_df[, "TSS_searched_by_GuideScan"] %in% "No") & (submitted_to_GuideScan_vec %in% "Yes")
  merged_CRISPR_df[are_conflicting, "TSS_searched_by_GuideScan"] <- "Not this gene"

  ####################################################################################################################
  ####################################################################################################################
  ####################################################################################################################


  # Final steps
  merged_CRISPR_df[, "GuideScan_offtarget_category"] <- GetOffTargetCategory(merged_CRISPR_df)


  # Re-order the columns
  results_df <- merged_CRISPR_df[, !(colnames(merged_CRISPR_df) %in% remove_columns)]

  return(results_df)
}





# Miscellaneous functions -------------------------------------------------

CheckForInconsistentChromosomes <- function(CRISPR_df) {
  have_more_than_one_chr_vec <- sapply(
    unique(CRISPR_df[, "Combined_ID"]),
    function(x) {
      chromosomes_vec <- CRISPR_df[CRISPR_df[, "Combined_ID"] == x, "Chromosome"]
      chromosomes_vec <- chromosomes_vec[!(is.na(chromosomes_vec))]
      more_than_one <- length(unique(chromosomes_vec)) > 1
      return(more_than_one)
    })
  results_vec <- names(which(have_more_than_one_chr_vec))
  return(results_vec)
}






