### 9th November 2020 ###





# Define maps -------------------------------------------------------------

column_rename_vec <- c(
  "np" = "Num_full_passes",
  "rq" = "Read_quality",
  "qs" = "Clip_start",
  "qe" = "Clip_end",
  "bq" = "Barcode_combined_score",
  "bl" = "Barcode_1_sequence",
  "bt" = "Barcode_2_sequence",
  "ql" = "Barcode_1_quality",
  "qt" = "Barcode_2_quality",
  "RG" = "Read_group_ID"
)




# Define functions --------------------------------------------------------

ReadLimaReport <- function(file_path) {
  read.table(file_path,
             header = TRUE, sep = "\t",
             quote = "", stringsAsFactors = FALSE
             )
}



SplitCharVec <- function(char_vec, use_mode = "integer") {
  char_splits <- strsplit(char_vec, ",", fixed = TRUE)
  results_mat <- do.call(rbind, char_splits)
  mode(results_mat) <- use_mode
  return(results_mat)
}




TidyReportDf <- function(report_df) {

  clips_mat <- SplitCharVec(report_df[["ClipsCombined"]])
  colnames(clips_mat) <- paste0("Clip_", c("start", "end"))
  clips_mat[, 1] <- clips_mat[, 1] + 1L # Convert from 0-based to 1-based

  barcodes_mat <- SplitCharVec(report_df[["IdxsCombined"]])

  combo_IDs_vec <- paste0(report_df[["IdxLowestNamed"]], "--",
                          report_df[["IdxHighestNamed"]]
                          )
  results_df <- data.frame(
    "ZMW"                    = as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]]))),
    "Well_number"            = barcodes_to_wells_map[combo_IDs_vec],
    "Barcode_combined_score" = report_df[["ScoreCombined"]],
    "Barcode_score_lead"     = report_df[["ScoreLead"]],
    "Passed_filters"         = as.logical(report_df[["PassedFilters"]]),
    "Clipped_read_length"    = clips_mat[, 2] - clips_mat[, 1] + 1L,
    clips_mat,
    stringsAsFactors = FALSE
  )
  orientation_vec <- ifelse(is.na(results_df[["Well_number"]]),
                            NA,
                            barcodes_mat[, 1] > barcodes_mat[, 2]
                            )
  results_df[["Barcodes_orientation_fwd"]] <- orientation_vec
  return(results_df)
}



IntegrateReportDf <- function(sam_df, report_df, ccs5_report_df = NULL) {

  tidied_report_df <- TidyReportDf(report_df)

  ZMW_matches <- match(sam_df[["ZMW"]], tidied_report_df[["ZMW"]])

  matched_report_df <- tidied_report_df[ZMW_matches, ]

  if (!("Barcode_1_length" %in% names(sam_df))) {
    for (column_name in c("Sequence", "Quality")) {
      vec_1 <- substr(sam_df[[column_name]], 1, matched_report_df[["Clip_start"]] - 1)
      vec_2 <- substr(sam_df[[column_name]], matched_report_df[["Clip_end"]] + 1, nchar(sam_df[[column_name]]))
      new_column_names <- paste0("Barcode_", 1:2, "_", tolower(column_name))
      sam_df[[new_column_names[[1]]]] <- vec_1
      sam_df[[new_column_names[[2]]]] <- vec_2
    }
    sam_df[["Barcode_1_length"]] <- nchar(vec_1)
    sam_df[["Barcode_2_length"]] <- nchar(vec_2)
  }

  duplicated_columns <- intersect(names(sam_df), names(matched_report_df))
  for (column_name in duplicated_columns) {
    stopifnot(identical(sam_df[[column_name]], matched_report_df[[column_name]]))
  }
  matched_report_df <- matched_report_df[, !(names(matched_report_df) %in% duplicated_columns)]

  results_df <- data.frame(
    sam_df,
    matched_report_df,
    stringsAsFactors = FALSE
  )

  if (!(is.null(ccs5_report_df))) {
    ccs5_zmws <- as.integer(substr(ccs5_report_df[["ZMW"]], 22, nchar(ccs5_report_df[["ZMW"]])))
    results_df[["Pass_CCS5"]] <- results_df[["ZMW"]] %in% ccs5_zmws
  }

  all_columns <- c(
    "ZMW", "Well_number",
    "Clipped_read_length",
    "Passed_filters", "Pass_CCS5",
    "Num_full_passes", "Read_quality", "Barcode_combined_score", "Barcode_score_lead",
    "Clip_start", "Clip_end", "Barcodes_orientation_fwd",
    "Barcode_1_length", "Barcode_2_length",
    "Barcode_1_sequence", "Barcode_2_sequence", "Barcode_1_quality", "Barcode_2_quality",
    "Sequence", "Quality"
  )
  results_df <- results_df[, all_columns]
  return(results_df)
}





ProcessSAM <- function(SAM_file_path) {

  sam_vec <- as.character(readLines(SAM_file_path))
  tab_search_list <- gregexpr("\t", sam_vec, fixed = TRUE)

  num_tabs <- max(lengths(tab_search_list))
  are_header <- lengths(tab_search_list) < num_tabs
  last_header <- max(which(are_header))
  stopifnot(all(are_header[seq_len(last_header)]))
  message(paste0("The first ", last_header, " lines seemed to belong to the ",
                 "header and were skipped!"
                 )
          )
  use_indices <- seq(from = last_header + 1, to = length(sam_vec))
  tab_search_list <- tab_search_list[use_indices]
  sam_vec <- sam_vec[use_indices]

  sequence_vec <- vapply(seq_along(sam_vec), function(x) {
    substr(sam_vec[[x]], tab_search_list[[x]][[9]] + 1, tab_search_list[[x]][[10]] - 1)
  }, "")
  sequence_lengths <- nchar(sequence_vec)
  quality_vec <- vapply(seq_along(sam_vec), function(x) {
    substr(sam_vec[[x]], tab_search_list[[x]][[10]] + 1, tab_search_list[[x]][[10]] + sequence_lengths[[x]])
  }, "")
  last_part_vec <- vapply(seq_along(sam_vec), function(x) {
    substr(sam_vec[[x]], tab_search_list[[x]][[10]] + sequence_lengths[[x]] + 2, nchar(sam_vec[[x]]))
  }, "")

  last_part_list <- strsplit(last_part_vec, "\t", fixed = TRUE)
  last_part_vec_list <- lapply(seq_len(unique(lengths(last_part_list))), function(x) {
    sapply(last_part_list, "[[", x)
  })

  last_part_vec_list_list <- lapply(last_part_vec_list, function(x) {
    split_list <- strsplit(x, ":", fixed = TRUE)
    var_name <- unique(sapply(split_list, "[[", 1))
    results_vec <- sapply(split_list, function(x) x[[length(x)]])
    var_type <- unique(sapply(split_list, "[[", 2))
    if (var_type == "i") {
      results_vec <- as.integer(results_vec)
    } else if (var_type == "f") {
      results_vec <- as.numeric(results_vec)
    }
    results_list <- list(
      "name" = var_name,
      "vector" = results_vec
    )
    return(results_list)
  })

  last_part_df <- do.call(data.frame,
                          c(lapply(last_part_vec_list_list, function(x) x[["vector"]]),
                            list(stringsAsFactors = FALSE)
                            )
                          )
  names(last_part_df) <- vapply(last_part_vec_list_list, function(x) x[["name"]], "")

  last_part_df[["sn"]] <- sub("^f,", "", last_part_df[["sn"]])
  sn_mat <- SplitCharVec(last_part_df[["sn"]], "numeric")
  colnames(sn_mat) <- paste0("sn_", c("a", "c", "g", "t"))

  if ("bx" %in% names(last_part_df)) {
    last_part_df[["bx"]] <- sub("^i,", "", last_part_df[["bx"]])
    last_part_df[["bc"]] <- sub("^S,", "", last_part_df[["bc"]])
    barcode_mat <- cbind(SplitCharVec(last_part_df[["bc"]]),
                         SplitCharVec(last_part_df[["bx"]])
                         )
    barcode_mat[, 1:2] <- barcode_mat[, 1:2] + 1L # convert from 0-based to 1-based
    colnames(barcode_mat) <- paste0("Barcode_",
                                    c(c("lowest", "highest"), paste0(1:2, "_length"))
                                    )
    last_part_df[["qs"]] <- last_part_df[["qs"]] + 1L # convert from 0-based to 1-based
    last_part_df <- data.frame(last_part_df, barcode_mat, stringsAsFactors = FALSE)
  }

  omit_columns <- c("sn", "bc", "bx")
  check_uniqueness_columns <- c("RG", "cx")
  for (column_name in check_uniqueness_columns) {
    if (length(unique(last_part_df[[column_name]])) == 1) {
      omit_columns <- c(omit_columns, column_name)
    }
  }

  last_part_df <- data.frame(
    last_part_df[, !(names(last_part_df) %in% omit_columns)],
    sn_mat,
    stringsAsFactors = FALSE
  )
  are_to_rename <- names(last_part_df) %in% names(column_rename_vec)
  names(last_part_df)[are_to_rename] <- column_rename_vec[names(last_part_df)[are_to_rename]]

  results_df <- data.frame(
    "ZMW"      = last_part_df[["zm"]],
    "Sequence" = sequence_vec,
    "Quality"  = quality_vec,
    last_part_df[, names(last_part_df) != "zm"],
    stringsAsFactors = FALSE
  )
  return(results_df)
}


