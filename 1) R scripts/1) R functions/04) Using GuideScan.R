### 26 July 2019 ###



# Define functions --------------------------------------------------------

TSSRangesForGuideScan <- function(TSS_df, TSS_range = 1202L, total_range = (1000L * 1000L - 2L)) {

  if (TSS_range > total_range) {
    stop("The TSSDfToVec function assumes that the search space around each TSS is smaller than the total permitted range!")
  }

  TSS_center_vec <- rowMeans(cbind(TSS_df[, "Last_TSS"], TSS_df[, "First_TSS"]))

  max_start_vec <- as.integer(ceiling(TSS_center_vec - (total_range / 2)))
  max_start_vec <- ifelse(max_start_vec < 0L, 0L, max_start_vec)
  max_end_vec <- as.integer(floor(TSS_center_vec + (total_range / 2)))

  preferred_start_vec <- TSS_df[, "First_TSS"] - (TSS_range / 2L)
  preferred_end_vec   <- TSS_df[, "Last_TSS"] + (TSS_range / 2L)

  exceed_range_vec <- (preferred_start_vec < max_start_vec) | (preferred_end_vec > max_end_vec)

  single_TSS_vec <- ifelse(is.na(TSS_df[, "Best_TSS"]), TSS_df[, "First_TSS"], TSS_df[, "Best_TSS"])

  single_start_vec <- single_TSS_vec - (total_range / 2L)
  single_end_vec   <- single_TSS_vec + (total_range / 2L)

  exceed_range_start_bounded <- exceed_range_vec & (preferred_start_vec > single_start_vec)
  exceed_range_end_bounded   <- exceed_range_vec & (preferred_end_vec < single_end_vec)
  exceed_range_centered      <- exceed_range_vec & !(exceed_range_start_bounded | exceed_range_end_bounded)

  results_start_vec <- ifelse(exceed_range_vec, NA_integer_, preferred_start_vec)
  results_end_vec   <- ifelse(exceed_range_vec, NA_integer_, preferred_end_vec)

  results_start_vec[exceed_range_end_bounded]   <- (preferred_end_vec - total_range)[exceed_range_end_bounded]
  results_end_vec[exceed_range_end_bounded]     <- preferred_end_vec[exceed_range_end_bounded]

  results_start_vec[exceed_range_start_bounded] <- preferred_start_vec[exceed_range_start_bounded]
  results_end_vec[exceed_range_start_bounded]   <- (preferred_start_vec + total_range)[exceed_range_start_bounded]

  results_start_vec[exceed_range_centered]      <- single_start_vec[exceed_range_centered]
  results_end_vec[exceed_range_centered]        <- single_end_vec[exceed_range_centered]

  results_df <- data.frame(
    "TSS_used"       = ifelse(exceed_range_vec, ifelse(is.na(TSS_df[, "Best_TSS"]), "First", "Highest FANTOM5 score"), "All"),
    "Is_centered"    = ifelse(!(exceed_range_vec) | exceed_range_centered, "Yes", "No"),
    "Start"          = results_start_vec,
    "End"            = results_end_vec,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  return(results_df)
}


TSSStringForGuideScan <- function(TSS_df, ...) {
  TSS_GuideScan_df <- TSSRangesForGuideScan(TSS_df, ...)
  results_vec <- paste0(TSS_df[, "Chromosome"], ":", TSS_GuideScan_df[, "Start"], "-", TSS_GuideScan_df[, "End"])
  return(results_vec)
}


sgRNAStringForGuideScan <- function(CRISPR_df) {
  start_vec <- ifelse(CRISPR_df[, "Strand"] == "+", CRISPR_df[, "Start"], CRISPR_df[, "Start"] - 3L)
  end_vec   <- ifelse(CRISPR_df[, "Strand"] == "+", CRISPR_df[, "End"] + 3L, CRISPR_df[, "End"])
  start_vec <- start_vec - 1L
  end_vec <- end_vec - 1L
  results_vec <- paste0(CRISPR_df[, "Chromosome"], ":", start_vec, "-", end_vec)
  return(results_vec)
}


ReadGuideScanOutput <- function(file_name, use_fread = FALSE) {
  # This function requires 'CRISPOR_files_directory' in the global workspace.
  # Unfortunately, use of data.table::fread does not lead to a significant
  # speedup for this function. The main delay is not caused by reading the
  # the file, but by the 'strsplit(string_vec, ",", fixed = TRUE)' operation,
  # which is very time-consuming.
  file_path <- file.path(GuideScan_files_directory, file_name)
  results_df <- read.csv(file = file_path, header = FALSE, row.names = NULL,
                         quote = "\"", stringsAsFactors = FALSE, comment.char = ""
                         )
  return(results_df)
}



GuideScanOutputToDf <- function(GuideScan_output_df) {

  are_title <- apply(GuideScan_output_df, 1, function(x) all(as.character(x[2:11]) == ""))
  file_number_vec <- rep.int(NA_integer_, nrow(GuideScan_output_df))
  file_number <- 0L
  for (i in seq_along(file_number_vec)) {
    if (are_title[[i]]) {
      file_number <- file_number + 1L
    }
    file_number_vec[[i]] <- file_number
  }

  guidescan_df_list <- split(GuideScan_output_df, file_number_vec)

  guidescan_df_list <- lapply(guidescan_df_list, function (x) {
    results_df <- x[3:nrow(x), , drop = FALSE]
    results_df <- data.frame("Region" = x[1, 1], results_df, stringsAsFactors = FALSE, row.names = NULL)
    colnames(results_df)[2:ncol(results_df)] <- as.character(x[2, ])
    return(results_df)
  })

  results_df <- do.call(rbind.data.frame, c(guidescan_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  colnames(results_df) <- gsub(" ", "_", colnames(results_df), fixed = TRUE)
  return(results_df)
}




MakeGuideScanFileNumbers <- function(raw_df) {
  are_title <- raw_df[, 2] == ""
  file_number_vec <- rep.int(NA_integer_, nrow(raw_df))
  file_number <- 0L
  for (i in seq_along(file_number_vec)) {
    if (are_title[[i]]) {
      file_number <- file_number + 1L
    }
    file_number_vec[[i]] <- file_number
  }
  return(file_number_vec)
}




BuildGuideScanDf <- function(raw_df, TSS_df, CRISPR_df) {

  ### Collect all sgRNAs in CRISPRa_df for a given region ###

  TSS_df[, "GuideScan_input"] <- TSSStringForGuideScan(TSS_df)
  combined_ID_for_region <- sapply(unique(TSS_df[, "GuideScan_input"]),
                                   function(x) unique(TSS_df[TSS_df[, "GuideScan_input"] == x, "Combined_ID"]),
                                   simplify = FALSE
                                   )
  sgRNAs_for_region <- lapply(combined_ID_for_region, function(x) unique(toupper(CRISPR_df[CRISPR_df[, "Combined_ID"] %in% x, "sgRNA_sequence"])))



  ### Split the rows of raw_df into the individual regions ###

  file_number_vec <- MakeGuideScanFileNumbers(raw_df)



  ### Generate a list of filtered GuidesScan data frames ###

  regions_vec <- tapply(seq_len(nrow(raw_df)), file_number_vec, function(x) raw_df[x[[1]], 1])
  regions_vec <- unname(regions_vec)

  guidescan_indices_list <- tapply(seq_len(nrow(raw_df)), file_number_vec, function(x) x[3:length(x)])

  sgRNAs_for_region_matches <- match(regions_vec, names(sgRNAs_for_region))

  guidescan_df_list <- lapply(seq_along(guidescan_indices_list), function(x) {
    if ((((x %% 100) == 0) && (x < 1000)) || (((x %% 1000)) == 0)) {
      message(paste0(x, " out of ", length(guidescan_indices_list), " locations have been processed."))
    }
    sub_df <- raw_df[guidescan_indices_list[[x]], ]
    are_present_sgRNAs <- sub_df[, 4] %in% sgRNAs_for_region[[sgRNAs_for_region_matches[[x]]]]
    if (any(are_present_sgRNAs)) {
      return(data.frame("Region" = regions_vec[[x]], sub_df[are_present_sgRNAs, ], stringsAsFactors = FALSE, row.names = NULL))
    } else {
      return(NULL)
    }
  })


  ### Assemble the filtered GuideScan output ###

  guidescan_all_genes_df <- do.call(rbind.data.frame, c(guidescan_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  colnames(guidescan_all_genes_df)[2:ncol(guidescan_all_genes_df)] <- gsub(" ", "_", raw_df[2, ], fixed = TRUE)

  return(guidescan_all_genes_df)
}



