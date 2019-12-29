### 30th October 2019 ###




# Import packages and source code -----------------------------------------

library("Biostrings")





# Define functions --------------------------------------------------------

OnesAndZeros <- function(my_vec) {
  current_value <- my_vec[[1]]
  current_result <- 0L
  result_vec <- rep.int(NA_integer_, length(my_vec))
  for (i in seq_along(my_vec)) {
    if (current_value != my_vec[[i]]) {
      current_value <- my_vec[[i]]
      current_result <- abs(current_result - 1L)
    }
    result_vec[[i]] <- current_result
  }
  return(result_vec)
}



MoveAfterColumn <- function(my_df, after_this_column, column_to_move) {
  stopifnot(all(c(after_this_column, column_to_move) %in% colnames(my_df)))
  column_index <- match(after_this_column, colnames(my_df))
  before_columns <- colnames(my_df)[seq_len(column_index - 1)]
  if (column_index != ncol(my_df)) {
    after_columns <- colnames(my_df)[(column_index + 1):ncol(my_df)]
  } else {
    after_columns <- c()
  }
  before_columns <- before_columns[before_columns != column_to_move]
  after_columns <- after_columns[after_columns != column_to_move]
  all_columns <- c(before_columns, after_this_column, column_to_move, after_columns)
  return(my_df[, all_columns])
}


FormatForExcel <- function(my_df,
                           remove_columns            = NULL,
                           probability_to_percentage = FALSE,
                           convert_excluded_to_3     = TRUE,
                           add_primers               = FALSE,
                           allow_curated             = FALSE
                           ) {

  # Depends on 'source_abbreviations_vec' in the global environment

  is_CRISPRa <- "Calabrese_rank" %in% colnames(my_df)

  ones_and_zeros_vec <- OnesAndZeros(my_df[, "Combined_ID"])
  if (convert_excluded_to_3) {
    are_to_be_excluded <- !(MeetCriteria(my_df, allow_curated = allow_curated))
    ones_and_zeros_vec[are_to_be_excluded] <- 2L
  }

  are_controls <- my_df[, "Is_control"] == "Yes"
  if (any(are_controls)) {
    if (is_CRISPRa) {
      my_df[are_controls, "hCRISPRa_v2_transcript"] <- NA_character_
    }
    my_df[are_controls, "Gene_symbol"]            <- my_df[are_controls, "Combined_ID"]
    my_df[are_controls, "Original_symbol"]        <- NA_character_
    my_df[are_controls, "Rank"]                   <- NA_integer_
  }

  my_df[, "Num_overlaps"] <- ifelse(is.na(my_df[, "Num_overlaps"]),
                                    ifelse(my_df[, "Rank"] %in% 1:4, ifelse(my_df[, "Spacing"] %in% 10, "<10bp", "n.d."), NA_character_),
                                    paste0(my_df[, "Num_overlaps"], "|", my_df[, "Spacing"], "bp")
                                    )
  if (add_primers) {
    assign("delete_my_df_2", my_df, envir = globalenv())
    my_df[, "Sequence_with_primers"] <- AddPrimers(my_df)
    my_df <- MoveAfterColumn(my_df, "PAM", "Sequence_with_primers")
  }

  if (all(c("Exon_number_Brunello", "Exon_number_TKOv3", "Exon_number_GPP") %in% colnames(my_df))) {
    Brunello_exon_vec <- my_df[, "Exon_number_Brunello"]
    TKOv3_exon_vec <- my_df[, "Exon_number_TKOv3"]
    GPP_exon_vec <- my_df[, "Exon_number_GPP"]
    my_df[, "Exon_number_Brunello"] <- vapply(seq_len(nrow(my_df)), function(x) {
      exon_vec <- unique(c(Brunello_exon_vec[[x]], TKOv3_exon_vec[[x]], GPP_exon_vec[[x]]))
      exon_vec <- exon_vec[!(is.na(exon_vec))]
      return(paste0(exon_vec, collapse = " | "))
    }, "")
    colnames(my_df)[colnames(my_df) == "Exon_number_Brunello"] <- "Exon_number"
    my_df <- my_df[, colnames(my_df) != "Exon_number_TKOv3"]
    my_df <- my_df[, colnames(my_df) != "Exon_number_GPP"]
  }

  if ("CRISPOR_Graf_status" %in% colnames(my_df)) {
    my_df[, "CRISPOR_Graf_status"] <- ifelse(my_df[, "CRISPOR_Graf_status"] == "GrafOK", "OK", my_df[, "CRISPOR_Graf_status"])
  }

  if (!(is.null(remove_columns))) {
    my_df <- my_df[, !(colnames(my_df) %in% remove_columns)]
  }

  SNP_ID_column <- grep("_SNP_IDs_", colnames(my_df), fixed = TRUE)
  SNP_AF_column <- grep("_SNP_AF_(max|sum)_", colnames(my_df))
  if ((length(SNP_ID_column) == 1) && (length(SNP_AF_column) == 1)) {
    my_df[is.na(my_df[, SNP_AF_column]), SNP_ID_column] <- NA_character_
    my_df[is.na(my_df[, SNP_ID_column]), SNP_AF_column] <- NA_real_
  }

  if (is_CRISPRa) {
    my_df[, "Calabrese_rank"] <- gsub("/", " or ", my_df[, "Calabrese_rank"], fixed = TRUE)
    my_df[, "Calabrese_rank"] <- sub(" or ", ", ", my_df[, "Calabrese_rank"], fixed = TRUE)
  }
  if (probability_to_percentage) {
    for (column_index in grep("_AF_(sum|max)_", colnames(my_df), fixed = TRUE)) {
      my_df[, column_index] <- ifelse(is.na(my_df[, column_index]), NA_character_, paste0(my_df[, column_index] * 100, "%"))
    }
  }
  for (i in seq_len((ncol(my_df) - 7))) {
    my_df[, i] <- ifelse(is.na(my_df[, i]), "", as.character(my_df[, i]))
  }
  my_df <- AbbreviateColumns(my_df)
  for (i in (ncol(my_df) - 6):ncol(my_df)) {
    my_df[, i] <- ifelse(is.na(my_df[, i]), " ", as.character(my_df[, i]))
  }

  have_multiple_sources <- grepl(", ", my_df[, "Source"], fixed = TRUE)

  for (source in names(source_abbreviations_vec)) {
    my_df[have_multiple_sources, "Source"] <- sub(source, source_abbreviations_vec[[source]], my_df[have_multiple_sources, "Source"], fixed = TRUE)
  }

  my_df[, "Color"] <- ones_and_zeros_vec + 1L

  assign("delete_my_df", my_df, envir = globalenv())
  return(my_df)
}





DfToTSV <- function(CRISPR_df,file_name, remove_columns = full_omit_columns, probability_to_percentage = FALSE,
                    add_primers = FALSE, allow_curated = FALSE
                    ) {
  write.table(FormatForExcel(CRISPR_df, remove_columns = remove_columns, probability_to_percentage = probability_to_percentage,
                             add_primers = add_primers, allow_curated = allow_curated
                             ),
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )
  return(invisible(NULL))
}



AbbreviateColumns <- function(CRISPR_df) {
  for (column_name in c("Locations_0MM", "Locations_1MM", "Sequences_1MM")) {
    my_lengths <- lengths(strsplit(CRISPR_df[, column_name], "; ", fixed = TRUE))
    max_entries <- 10
    are_too_long <- my_lengths > max_entries
    if (any(are_too_long)) {
      CRISPR_df[are_too_long, column_name] <- "too long"
      message(paste0(sum(are_too_long), " rows in the ", column_name, " column were ommitted, because they were too long (>", max_entries, " entries)."))
    }
  }
  return(CRISPR_df)
}


primer_sequences <- list(
  "sg1" = c("ttgtggaaaggacgaaacaccG", "GTTTAAGAGCTAAGCTG"),
  "sg2" = c("cttggagaaaagccttgtttG", "GTTTGAGAGCTAAGCAGA"),
  "sg3" = c("gtatgagaccactctttcccG", "GTTTCAGAGCTAAGCACA"),
  "sg4" = c("gccgcttgggtacctcG", "GTTTCAGAGCTACAGCAGAAAT")
)


AddPrimers <- function(CRISPR_df) {
  results_vec <- rep.int(NA_character_, nrow(CRISPR_df))
  for (i in 1:4) {
    are_this_rank <- CRISPR_df[, "Rank"] %in% i
    new_sequences <- paste0(primer_sequences[[i]][[1]], CRISPR_df[are_this_rank, "sgRNA_sequence"], primer_sequences[[i]][[2]])
    if (i == 4) {
      new_sequences <- Biostrings::reverseComplement(Biostrings::DNAStringSet(new_sequences))
    }
    results_vec[are_this_rank] <- new_sequences
  }
  return(results_vec)
}


AreCompleteTranscripts <- function(CRISPR_df) {
  are_incomplete <- rep.int(NA, nrow(CRISPR_df))
  unique_TSS_IDs <- unique(CRISPR_df[CRISPR_df[, "Is_control"] == "No", "AltTSS_ID"])
  for (unique_TSS_ID in unique_TSS_IDs) {
    are_this_TSS <- CRISPR_df[, "AltTSS_ID"] %in% unique_TSS_ID
    is_complete <- all(1:4 %in% CRISPR_df[are_this_TSS, "Rank"])
    are_incomplete[are_this_TSS] <- is_complete
  }
  return(are_incomplete)
}



FormatFixedWidthInteger <- function(integer_vec) {
  integer_width <- max(nchar(as.character(as.integer(integer_vec))))
  result <- formatC(integer_vec, width = integer_width, flag = "0")
  return(result)
}














