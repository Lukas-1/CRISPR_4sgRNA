### 15th August 2021 ###



# Functions for preparing mutations_df ------------------------------------

FilterMutations <- function(extract_df, use_zmws, min_length = 19, max_length = 42) {

  are_sg <- extract_df[, "Feature"] %in% paste0("sg", 1:4)
  are_mutated <- extract_df[, "Category"] %in% "Mutation"
  are_eligible <- extract_df[, "ZMW"] %in% use_zmws

  mutations_df <- extract_df[are_sg & are_mutated & are_eligible, ]

  mutations_df[["Read_without_gaps"]] <- gsub("-", "", mutations_df[["Aligned_read"]], fixed = TRUE)

  names(mutations_df)[[which(names(mutations_df) == "Feature")]] <- "sg_number"
  mutations_df[["sg_number"]] <- as.integer(substr(mutations_df[["sg_number"]], 3, 3))

  sequence_lengths <- nchar(mutations_df[["Read_without_gaps"]])
  are_too_short <- sequence_lengths < min_length
  are_too_long <- sequence_lengths > max_length
  mutations_df <- mutations_df[!(are_too_short | are_too_long), ]
  row.names(mutations_df) <- NULL
  return(mutations_df)
}



AddPrefixToLast3 <- function(input_df, use_prefix) {
  column_names <- names(input_df)[(ncol(input_df) - 2):ncol(input_df)]
  column_names <- paste0(use_prefix, "_",
                         ifelse(substr(column_names, 2, 2) == toupper(substr(column_names, 2, 2)),
                                column_names,
                                paste0(tolower(substr(column_names, 1, 1)),
                                       substr(column_names, 2, nchar(column_names))
                                       )

                                )
                         )
  names(input_df)[(ncol(input_df) - 2):ncol(input_df)] <- column_names
  return(input_df)
}



FindHitsMutatedAndTemplate <- function(mut_df) {
  mut_df <- Add0MMHits(mut_df, "Read_without_gaps")
  mut_df <- AddPrefixToLast3(mut_df, "Mutated")

  mut_df <- Add0MMHits(mut_df, "Template")
  mut_df <- AddPrefixToLast3(mut_df, "Template")

  redundant_columns <- c("Category", "Is_correct", "Mostly_deleted",
                         "Is_contamination"
                         )
  mut_df <- mut_df[, !(names(mut_df) %in% redundant_columns)]
  return(mut_df)
}




# Functions for annotating mutations_df -----------------------------------

GetEntrezsList <- function(full_df) {
  entrezs_splits <- strsplit(full_df[["Affected_Entrez_IDs"]], ", ", fixed = TRUE)
  entrezs_split_splits <- split(entrezs_splits, full_df[["Index"]])
  entrezs_list <- lapply(entrezs_split_splits, function(x) as.integer(unique(unlist(x, use.names = FALSE))))
  return(entrezs_list)
}


GetNewTargets <- function(template_full_df, mutated_full_df) {

  template_entrezs <- GetEntrezsList(template_full_df)
  mutated_entrezs <- GetEntrezsList(mutated_full_df)

  only_mutated_entrezs <- mapply(setdiff, mutated_entrezs, template_entrezs)
  only_mutated_entrezs[is.na(only_mutated_entrezs)] <- list(integer(0))

  mutated_entrezs_splits <- strsplit(mutated_full_df[["Affected_Entrez_IDs"]], ", ", fixed = TRUE)
  mutated_entrezs_split_splits <- split(mutated_entrezs_splits, mutated_full_df[["Index"]])

  have_only_mutated <- mapply(function(x, y) x %in% y, mutated_entrezs_split_splits, only_mutated_entrezs)
  are_new_genes <- unlist(have_only_mutated, use.names = FALSE)

  are_intended <- (mutated_full_df[["Intended_Entrez_ID"]] == mutated_full_df[["Affected_Entrez_IDs"]]) %in% TRUE
  are_selected <- are_new_genes | are_intended

  are_selected_split <- split(are_selected, mutated_full_df[["Index"]])
  are_only_entry <- unlist(lapply(are_selected_split, function(x) {
     if (!(any(x))) {
       c(TRUE, rep(FALSE, length(x) - 1))
     } else {
       rep(FALSE, length(x))
     }
  }), use.names = FALSE)

  results_df <- mutated_full_df
  preserve_columns <- c(
    "Locus", "Index", "Intended_Entrez_ID", "Entrez_ID_available",
    "Intended_gene_symbol", "Num_loci"
  )
  other_columns <- setdiff(names(results_df), preserve_columns)
  for (other_column in other_columns) {
    results_df[[other_column]][are_only_entry] <- NA
  }
  are_to_retain <- are_selected | are_only_entry
  results_df <- results_df[are_to_retain, ]
  row.names(results_df) <- NULL
  return(results_df)
}



GetNumUnintended <- function(input_df) {
  targets_intended <- input_df[, "Loci_targeting_intended_gene"] >= 1
  num_unintended <- input_df[, "Num_affected_genes"] - as.integer(targets_intended)
  return(num_unintended)
}


RemoveIntendedGene <- function(char_vec) {
  assign("delete_char_vec", char_vec, envir = globalenv())
  are_present <- (char_vec != "") & !(is.na(char_vec))
  target_intended <- are_present & !(grepl(" is not ", char_vec, fixed = TRUE))
  only_unintended <- are_present & !(target_intended) & !(grepl("!$", char_vec))
  char_splits <- strsplit(char_vec[only_unintended], ": ", fixed = TRUE)
  results_vec <- character(length = length(char_vec))
  results_vec[target_intended] <- char_vec[are_present & target_intended]
  stopifnot(all(lengths(char_splits) == 2))
  results_vec[only_unintended] <- sapply(char_splits, "[[", 2)
  return(results_vec)
}



AnnotateMutations <- function(all_mut_df,
                              modalities_vec,
                              GeneLociFunction,
                              gene_loci_df
                              ) {

  ## Prepare the mutations data frame

  coord_columns <- c("Chromosome", "Strand", "Start", "End")
  for (column_name in coord_columns) {
    all_mut_df[[column_name]] <- NA
  }

  are_correct_modality <- all_mut_df[["Modality"]] %in% modalities_vec
  have_0MM <- all_mut_df[["Mutated_num_0MM"]] >= 1

  message(paste0("Out of a total of ", sum(are_correct_modality),
                 " mutated gRNAs (for ",
                 paste0(modalities_vec, collapse = " and "),
                 "), \n", sum(have_0MM & are_correct_modality),
                 " contained one or more perfect-match binding sites in the genome.\n",
                 "(Note: Incorrect gRNAs with >50% deleted base pairs or that",
                 " matched the guide RNA sequences\nof other wells perfectly",
                 " (i.e. contaminations) were excluded.)\n"
                 )
          )

  use_mut_df <- all_mut_df[are_correct_modality & have_0MM, ]
  row.names(use_mut_df) <- NULL



  ## Identify genes targeted by the mutated and template gRNAs

  orig_df      <- use_mut_df
  mut_df       <- use_mut_df
  names(orig_df)[names(orig_df) == "Template_loci_0MM"] <- "Locations_0MM"
  names(mut_df)[names(mut_df) == "Mutated_loci_0MM"]  <- "Locations_0MM"
  names(orig_df)[names(orig_df) == "Template_num_0MM"] <- "Num_0MM"
  names(mut_df)[names(mut_df) == "Mutated_num_0MM"]  <- "Num_0MM"

  template_list <- AlignSummaryDf(GeneLociFunction, orig_df, gene_loci_df)
  mutated_list  <- AlignSummaryDf(GeneLociFunction, mut_df, gene_loci_df)



  ## Identify genes that are UNIQUELY targeted by the mutated gRNAs

  new_full_df <- GetNewTargets(template_list[["full_df"]], mutated_list[["full_df"]])

  new_summary_df <- SummarizeFullDf(new_full_df)
  new_summary_df <- SummarizeSummaryDf(new_summary_df)

  stopifnot(nrow(new_summary_df) == max(mutated_list[["full_df"]][["Index"]]))

  are_valid <- !(is.na(use_mut_df[["Entrez_ID"]])) &
               !(grepl(",", use_mut_df[["Entrez_ID"]], fixed = TRUE))
  indices_vec <- rep(NA, nrow(use_mut_df))
  indices_vec[are_valid] <- seq_len(sum(are_valid))
  new_summary_df <- new_summary_df[indices_vec, ]
  row.names(new_summary_df) <- NULL

  stopifnot(nrow(new_summary_df) == nrow(use_mut_df))



  ## Merge the annotation data into a final data frame

  results_df <- use_mut_df[, !(names(use_mut_df) %in% coord_columns)]
  results_df[["Targets_intended_gene"]] <- mutated_list[["summary_df"]][, "Loci_targeting_intended_gene"] >= 1
  results_df[["Num_all_unintended_genes"]] <- GetNumUnintended(mutated_list[["summary_df"]])
  results_df[["Num_new_unintended_genes"]] <- GetNumUnintended(new_summary_df)

  results_df[["Entrez_ID_available"]] <- template_list[["summary_df"]][, "Entrez_ID_available"]

  results_df[["All_affected_gene_symbols"]]                <- mutated_list[["summary_df"]][, "Affected_gene_symbols"]
  results_df[["New_unintended_gene_symbols"]]              <- new_summary_df[, "Unintended_gene_symbols"]
  results_df[["Template_affected_gene_symbols"]]           <- template_list[["summary_df"]][, "Affected_gene_symbols"]
  results_df[["All_affected_gene_symbols_truncated"]]      <- mutated_list[["summary_df"]][, "Affected_gene_symbols_truncated"]
  results_df[["New_unintended_gene_symbols_truncated"]]    <- new_summary_df[, "Unintended_gene_symbols_truncated"]
  results_df[["Template_affected_gene_symbols_truncated"]] <- template_list[["summary_df"]][, "Affected_gene_symbols"]

  has_entrezs <- "Affected_Entrez_IDs" %in% names(mutated_list[["summary_df"]])
  if (has_entrezs) {
    results_df[["All_affected_Entrez_IDs"]]                <- mutated_list[["summary_df"]][, "Affected_Entrez_IDs"]
    results_df[["New_unintended_Entrez_IDs"]]              <- new_summary_df[, "Unintended_Entrez_IDs"]
    results_df[["Template_affected_Entrez_IDs"]]           <- template_list[["summary_df"]][, "Affected_Entrez_IDs"]
    results_df[["All_affected_Entrez_IDs_truncated"]]      <- mutated_list[["summary_df"]][, "Affected_Entrez_IDs_truncated"]
    results_df[["New_unintended_Entrez_IDs_truncated"]]    <- new_summary_df[, "Unintended_Entrez_IDs_truncated"]
    results_df[["Template_affected_Entrez_IDs_truncated"]] <- template_list[["summary_df"]][, "Affected_Entrez_IDs"]
  } else {
    results_df[["All_affected_gene_IDs"]]                <- mutated_list[["summary_df"]][, "Affected_gene_IDs"]
    results_df[["New_unintended_gene_IDs"]]              <- new_summary_df[, "Unintended_gene_IDs"]
    results_df[["Template_affected_gene_IDs"]]           <- template_list[["summary_df"]][, "Affected_gene_IDs"]
    results_df[["All_affected_gene_IDs_truncated"]]      <- mutated_list[["summary_df"]][, "Affected_gene_IDs_truncated"]
    results_df[["New_unintended_gene_IDs_truncated"]]    <- new_summary_df[, "Unintended_gene_IDs_truncated"]
    results_df[["Template_affected_gene_IDs_truncated"]] <- template_list[["summary_df"]][, "Affected_gene_IDs"]
  }


  ## Select and rename columns

  if (has_entrezs) {
    gene_ID_columns <- c("All_affected_Entrez_IDs",
                         "New_unintended_Entrez_IDs",
                         "Template_affected_Entrez_IDs"
                         )
    gene_IDs_truncated <- c("All_affected_Entrez_IDs_truncated",
                            "New_unintended_Entrez_IDs_truncated",
                            "Template_affected_Entrez_IDs_truncated"
                            )
  } else {
    gene_ID_columns <- c("All_affected_gene_IDs",
                         "New_unintended_gene_IDs",
                         "Template_affected_gene_IDs"
                         )
    gene_IDs_truncated <- c("All_affected_gene_IDs_truncated",
                            "New_unintended_gene_IDs_truncated",
                            "Template_affected_gene_IDs_truncated"
                            )
  }

  unintended_columns <- c("New_unintended_gene_symbols",
                          "New_unintended_gene_symbols_truncated",
                          gene_ID_columns[[2]],
                          gene_IDs_truncated[[2]]
                          )
  for (column_name in unintended_columns) {
    results_df[, column_name] <- RemoveIntendedGene(results_df[, column_name])
  }

  columns_in_order <- c(
    "Combined_ID",  "ZMW", "sg_number", "Modality",
    "Gene_symbol", "Entrez_ID", "TSS_number",

    "Mutated_num_0MM", "Mutated_loci_0MM",
    "Targets_intended_gene", "Num_all_unintended_genes", "Num_new_unintended_genes",

    "All_affected_gene_symbols", "New_unintended_gene_symbols",
    "Template_affected_gene_symbols",
    gene_ID_columns,
    "Entrez_ID_available",

    "Template", "Aligned_template", "Aligned_read", "Read_without_gaps",
    "Num_incorrect", "Num_missing",

    "All_affected_gene_symbols_truncated", "New_unintended_gene_symbols_truncated",
    "Template_affected_gene_symbols_truncated",
    gene_IDs_truncated
  )

  results_df <- results_df[, columns_in_order]

  rename_vec <- c(
    "Entrez_ID_available" = "Location_available_for_Entrez_ID"
  )
  for (column_name in names(rename_vec)) {
    names(results_df)[names(results_df) == column_name] <- rename_vec[[column_name]]
  }


  ## Report on the results

  message(paste0("\nOut these, ",
                 sum((results_df[["Num_all_unintended_genes"]] >= 1) %in% TRUE),
                 " mutated single guide RNAs targeted at least one gene.",
                 "\nHowever, only ",
                 sum((results_df[["Num_new_unintended_genes"]] >= 1) %in% TRUE),
                 " of these targeted any genes in addition\nto those already",
                 " targeted by the template (unmutated) guide RNA."
                 )
          )
  return(results_df)
}



ExportMutatedDf <- function(input_df, file_name) {
  truncated_columns <- grep("_truncated$", names(input_df), value = TRUE)
  untruncated_columns <- sub("_truncated$", "", truncated_columns)
  export_df <- input_df[, !(names(input_df) %in% untruncated_columns)]
  export_df[["Mutated_loci_0MM"]] <- TruncateLongEntries(export_df[["Mutated_loci_0MM"]])
  export_df[["Targets_intended_gene"]] <- as.integer(export_df[["Targets_intended_gene"]])
  export_df[["Location_available_for_Entrez_ID"]] <- ifelse(export_df[["Location_available_for_Entrez_ID"]],
                                                            "Yes", "No"
                                                            )
  # export_df[["Aligned_template"]] <- paste0("'", export_df[["Aligned_template"]])
  # export_df[["Aligned_read"]] <- paste0("'", export_df[["Aligned_read"]])
  new_order <- order((export_df[["Num_new_unintended_genes"]] >= 1) %in% TRUE,
                     decreasing = TRUE
                     )
  export_df <- export_df[new_order, ]
  for (i in seq_along(export_df)) {
    export_df[[i]] <- ifelse(is.na(export_df[[i]]), "", export_df[[i]])
  }
  for (i in which(sapply(export_df, is.character))) {
    export_df[[i]] <- ifelse(export_df[[i]] == "", " ", export_df[[i]])
  }

  write.table(export_df,
              file      = file.path(output_directory, paste0(file_name, ".tsv")),
              sep       = "\t",
              row.names = FALSE#,
              # quote     = FALSE
              )
  return(invisible(NULL))
}






