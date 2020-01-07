 ### 4th November 2019 ###




# Import packages and source code -----------------------------------------

library("parallel")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"

source(file.path(general_functions_directory, "09) Constants and settings.R"))
source(file.path(general_functions_directory, "12) Re-ordering sgRNAs based on their genomic location.R"))







# Define functions --------------------------------------------------------

NumberTSSs <- function(locations_vec, min_space = 1001L) {
  results_vec <- rep.int(NA_integer_, length(locations_vec))
  current_number <- 0L
  current_position <- -Inf
  for (i in seq_along(locations_vec)) {
    if (!(is.na(locations_vec[[i]]))) {
      if (locations_vec[[i]] >= (current_position + min_space)) {
        current_number <- current_number + 1L
      }
      results_vec[[i]] <- current_number
      current_position <- locations_vec[[i]]
    }
  }
  return(results_vec)
}



AllocateTSSs <- function(TSS_number_vec, original_TSS_vec, positions_vec, new_TSS_prefix = "N") {

  stopifnot(length(TSS_number_vec) == length(original_TSS_vec))
  allocated_vec <- rep.int(NA_character_, length(TSS_number_vec))
  TSS_numbers <- unique(TSS_number_vec)
  TSS_numbers <- TSS_numbers[!(is.na(TSS_numbers))]
  new_TSS_number <- 0L
  for (TSS_number in TSS_numbers) {
    are_this_number <- TSS_number_vec %in% TSS_number
    original_TSS <- unique(original_TSS_vec[are_this_number])
    original_TSS <- original_TSS[!(is.na(original_TSS))]
    if (length(original_TSS) == 0) {
      new_TSS_number <- new_TSS_number + 1L
      new_TSS_string <- paste0(new_TSS_prefix, new_TSS_number)
      allocated_vec[are_this_number] <- new_TSS_string
    } else if (length(original_TSS) > 1) {
      message("Multiple transcript IDs were found with the same TSS number!")
      are_eligible <- are_this_number & !(is.na(original_TSS_vec))
      for (i in which(are_this_number)) {
        nearest_TSS <- original_TSS_vec[are_eligible][[which.min(abs(positions_vec[are_eligible] - positions_vec[[i]]))]]
        allocated_vec[[i]] <- nearest_TSS
      }
    } else {
      allocated_vec[are_this_number] <- original_TSS
    }
  }
  results_vec <- rep.int(NA_character_, length(allocated_vec))
  unique_TSSs <- unique(allocated_vec)
  unique_TSSs <- unique_TSSs[!(is.na(unique_TSSs))]
  for (allocated_TSS in unique_TSSs) {
    are_this_TSS <- allocated_vec %in% allocated_TSS
    if (length(unique(TSS_number_vec[are_this_TSS])) > 1) {
      pasted_vec <- paste0(allocated_vec[are_this_TSS], TSS_number_vec[are_this_TSS])
      my_matches <- match(pasted_vec, unique(pasted_vec))
      results_vec[are_this_TSS] <- paste0(allocated_vec[are_this_TSS], paste0("_loc", my_matches))
    } else {
      results_vec[are_this_TSS] <- allocated_TSS
    }
  }
  return(results_vec)
}



SeparateByTSS <- function(CRISPR_sub_df) {

  MessageID(CRISPR_sub_df)

  reordered_list <- ReorderSubDfByLocation(CRISPR_sub_df)
  reordered_df <- reordered_list[["reordered_df"]]
  TSS_number_vec <- NumberTSSs(reordered_df[, "Cut_location"])
  allocated_TSS_vec <- AllocateTSSs(TSS_number_vec, reordered_df[, "hCRISPRa_v2_transcript"], reordered_df[, "Cut_location"])

  results_df <- data.frame(
    reordered_df,
    "TSS_number"     = TSS_number_vec,
    "Allocated_TSS"  = allocated_TSS_vec,
    stringsAsFactors = FALSE
  )
  new_order <- order(match(results_df[, "Allocated_TSS"], results_df[, "Allocated_TSS"]))

  result_df <- results_df[new_order, ]
  row.names(result_df) <- NULL
  return(results_df)
}




CountNumTSSs <- function(CRISPR_df, TSS_column = "hCRISPRa_v2_transcript") {
  IDs_fac <- factor(CRISPR_df[, "Combined_ID"], levels = unique(CRISPR_df[, "Combined_ID"]))
  num_transcripts_list <- tapply(
    seq_len(nrow(CRISPR_df)),
    IDs_fac,
    function(x) {
      transcripts <- CRISPR_df[x, TSS_column]
      transcripts <- unique(transcripts[!(is.na(transcripts))])
      num_transcripts <- length(transcripts)
      if (num_transcripts == 0) {
        num_transcripts <- 1L
      }
      results_vec <- rep.int(num_transcripts, length(x))
      return(results_vec)
    },
    simplify = FALSE
  )
  results_vec <- unlist(num_transcripts_list, use.names = FALSE)
  return(results_vec)
}



TSSCombinedIDs <- function(CRISPR_df, TSS_prefix = "T") {
  IDs_fac <- factor(CRISPR_df[, "Combined_ID"], levels = unique(CRISPR_df[, "Combined_ID"]))
  abbr_TSS_list <- tapply(
    seq_len(nrow(CRISPR_df)),
    IDs_fac,
    function(x) {
      results_vec <- rep.int(NA_character_, length(x))
      num_TSSs <- CRISPR_df[x[[1]], "Num_TSSs"]
      if (num_TSSs > 1) {
        TSS_IDs <- CRISPR_df[x, "Allocated_TSS"]
        are_standard <- TSS_IDs %in% c("P1", "P2", "P1P2", paste0("N", 1:30))
        results_vec[are_standard] <- TSS_IDs[are_standard]
        if (!(all(are_standard))) {
          unique_not_standard <- unique(TSS_IDs[!(are_standard)])
          not_standard_rename_vec <- paste0(TSS_prefix, seq_along(unique_not_standard))
          names(not_standard_rename_vec) <- unique_not_standard
          results_vec[!(are_standard)] <- not_standard_rename_vec[TSS_IDs[!(are_standard)]]
        }
      }
      return(results_vec)
    },
    simplify = FALSE
  )
  results_vec <- unlist(abbr_TSS_list, use.names = FALSE)
  return(results_vec)
}



AllocateTSSforAllGenes <- function(CRISPR_df, omit_optional_columns = FALSE, parallel_mode = TRUE, num_cores = NULL) {

  are_controls <- CRISPR_df[, "Is_control"] == "Yes"
  combined_IDs <- unique(CRISPR_df[!(are_controls), "Combined_ID"])

  if (parallel_mode) {
    if (is.null(num_cores)) {
      num_cores <- parallel::detectCores() - 2
    }
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterExport(cl,
                            varlist = c("SeparateByTSS", "MessageID", "ReorderSubDfByLocation",
                                        "NumberTSSs", "AllocateTSSs", "CRISPR_df"
                                        ),
                            envir = environment()
                            )
    reordered_df_list <- parallel::parLapply(cl,
                                             combined_IDs,
                                             function(x) SeparateByTSS(CRISPR_df[CRISPR_df[, "Combined_ID"] == x, , drop = FALSE])
                                             )
    parallel::stopCluster(cl)

  } else {
    reordered_df_list <- lapply(combined_IDs, function(x) SeparateByTSS(CRISPR_df[CRISPR_df[, "Combined_ID"] == x, , drop = FALSE]))
  }

  if (any(are_controls)) {
    controls_df <- CRISPR_df[are_controls, ]
    controls_df[, "TSS_number"] <- NA_integer_
    controls_df[, "Allocated_TSS"] <- NA_character_
  } else {
    controls_df <- NULL
  }
  results_df <- do.call(rbind.data.frame, c(reordered_df_list,
                                            list(controls_df),
                                            list(stringsAsFactors = FALSE, make.row.names = FALSE)
                                            )
                        )
  results_df[, "Num_TSSs"] <- CountNumTSSs(results_df, TSS_column = "Allocated_TSS")
  results_df[, "TSS_ID"] <- TSSCombinedIDs(results_df)
  results_df[, "AltTSS_ID"] <- ifelse(is.na(results_df[, "TSS_ID"]),
                                      results_df[, "Combined_ID"],
                                      paste0(results_df[, "Combined_ID"], "_", results_df[, "TSS_ID"])
                                      )
  if (omit_optional_columns) {
    results_df <- results_df[, !(names(results_df) %in% c("TSS_number", "Allocated_TSS", "TSS_ID"))]
  }
  return(results_df)
}

















