### 26th July 2021 ###



# Define functions --------------------------------------------------------

TidySequencesDf <- function(CRISPR_df) {
  if ("AltTSS_ID" %in% names(CRISPR_df)) {
    ID_column <- "AltTSS_ID"
  } else {
    ID_column <- "Combined_ID"
    CRISPR_df[["TSS_number"]] <- NA
  }
  CRISPR_df[["Target_ID"]] <- CRISPR_df[[ID_column]]
  use_columns <- c("Plate_string", "Well_number", "Rank", "sgRNA_sequence",
                   "Target_ID", "Entrez_ID", "Gene_symbol", "TSS_number",
                   "Modality"
                   )
  results_df <- CRISPR_df[, use_columns]
  results_df[["Plate_string"]] <- sub("_tf", "_", results_df[["Plate_string"]], fixed = TRUE)
  results_df[["Plate_string"]] <- sub("_sg[1-4]$", "", results_df[["Plate_string"]])
  results_df[["Plate_string"]] <- toupper(results_df[["Plate_string"]])
  results_df[["Plate_string"]] <- sub("+", "plus", results_df[["Plate_string"]], fixed = TRUE)
  names(results_df) <- c("Plate_name", "Well_number", "Sg_number", "Sequence",
                         use_columns[5:length(use_columns)]
                         )
  results_df[["Sequence"]] <- toupper(results_df[["Sequence"]])
  return(results_df)
}



CreateControlsDf <- function(controls_Excel_df, tidy_CRISPRko_df) {
  controls_mat <- as.matrix(controls_Excel_df[15:(15 + 16 - 1), 3:25])
  controls_mat <- cbind(controls_mat, rep(NA, 16))
  colnames(controls_mat) <- NULL
  wells_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
  are_present <- !(is.na(controls_mat))
  control_wells <- wells_mat[are_present]

  are_controls <- (tidy_CRISPRko_df[["Plate_name"]] %in% "HO_38") &
                  (tidy_CRISPRko_df[["Well_number"]] %in% control_wells)
  control_plate_df <- tidy_CRISPRko_df[are_controls, ]
  control_plate_df[["Plate_name"]] <- "Intctrl"
  row.names(control_plate_df) <- NULL
  return(control_plate_df)
}



MakeLibraryDf <- function(input_plates_df, tidy_CRISPR_df) {
  sequences_list <- lapply(1:4, function(x) tidy_CRISPR_df[tidy_CRISPR_df[, "Sg_number"] %in% x, "Sequence"])
  sequences_mat <- do.call(cbind, sequences_list)
  colnames(sequences_mat) <- paste0("Sequence_sg", 1:4)

  use_columns <- setdiff(names(tidy_CRISPR_df), c("Sequence", "Sg_number"))
  library_df <- data.frame("Combined_ID" = NA,
                           "Plate_number" = NA,
                           tidy_CRISPR_df[tidy_CRISPR_df[, "Sg_number"] %in% 1, use_columns],
                           sequences_mat,
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )

  plate_matches <- match(library_df[, "Plate_name"], input_plates_df[["Plate_name"]])
  plate_numbers <- plates_df[plate_matches, "Plate_number"]
  plate_number_width <- max(nchar(as.character(plate_numbers)))
  well_names <- paste0("Plate", formatC(plate_numbers, width = plate_number_width, flag = "0"),
                       "_Well", formatC(library_df[, "Well_number"], width = 3, flag = "0")
                       )
  library_df[["Combined_ID"]] <- well_names
  library_df[["Plate_number"]] <- plate_numbers
  return(library_df)
}





