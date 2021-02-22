### 30th October 2019 ###




# Import packages and source code -----------------------------------------

library("Biostrings")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R"))




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
  stopifnot(all(c(after_this_column, column_to_move) %in% names(my_df)))
  column_index <- match(after_this_column, names(my_df))
  before_columns <- names(my_df)[seq_len(column_index - 1)]
  if (column_index != ncol(my_df)) {
    after_columns <- names(my_df)[(column_index + 1):ncol(my_df)]
  } else {
    after_columns <- c()
  }
  before_columns <- before_columns[before_columns != column_to_move]
  after_columns <- after_columns[after_columns != column_to_move]
  all_columns <- c(before_columns, after_this_column, column_to_move, after_columns)
  return(my_df[, all_columns])
}



RoundNumericColumns <- function(my_df, num_digits = 7) {
  are_floats <- vapply(seq_along(my_df), function(x) is.double(my_df[[x]]) , logical(1))
  for (i in which(are_floats)) {
    my_df[[i]] <- round(my_df[[i]], digits = num_digits)
  }
  return(my_df)
}




MakeIDs <- function(CRISPR_df) {
  paste0(CRISPR_df[["Entrez_ID"]], "__",
         toupper(CRISPR_df[["sgRNA_sequence"]])
         )

}


AddOtherTargets <- function(CRISPR_df, targets_sgRNA_df, targets_4sg_df = NULL) {
  stopifnot(nrow(CRISPR_df) == nrow(targets_sgRNA_df))
  if ("Unintended_Entrez_IDs_truncated" %in% names(targets_sgRNA_df)) {
    gene_ID_column <- "Unintended_Entrez_IDs_truncated"
  } else {
    gene_ID_column <- "Unintended_gene_IDs_truncated"
  }
  CRISPR_df[["Other_target_Entrez_IDs"]] <- targets_sgRNA_df[[gene_ID_column]]
  CRISPR_df[["Other_target_symbols"]] <- targets_sgRNA_df[["Unintended_gene_symbols_truncated"]]
  if (!(is.null(targets_4sg_df))) {
    CRISPR_df[["Other_Entrez_IDs_4sg"]] <- NA
    CRISPR_df[["Other_symbols_4sg"]] <- NA
    are_top4 <- CRISPR_df[["Rank"]] %in% 1:4
    matches_vec <- match(CRISPR_df[["Entrez_ID"]][are_top4],
                         targets_4sg_df[["Intended_Entrez_ID"]]
                         )
    CRISPR_df[["Other_Entrez_IDs_4sg"]][are_top4] <- targets_4sg_df[["Unintended_gene_IDs_truncated"]][matches_vec]
    CRISPR_df[["Other_symbols_4sg"]][are_top4] <- targets_4sg_df[["Unintended_gene_symbols_truncated"]][matches_vec]
  }
  return(CRISPR_df)
}



AddMainTSS <- function(CRISPR_df) {
  main_TSS_vec <- GetMainTSS(CRISPR_df)
  are_main_TSS <- (main_TSS_vec == CRISPR_df[["AltTSS_ID"]]) %in% TRUE
  CRISPR_df[["Is_main_TSS"]] <- ifelse(CRISPR_df[["Num_TSSs"]] == 1,
                                       "Single",
                                       ifelse(are_main_TSS,
                                              "Main",
                                              ifelse(is.na(CRISPR_df[["TSS_ID"]]), NA, "Other")
                                              )
                                       )
  return(CRISPR_df)
}



FormatForExcel <- function(my_df,
                           remove_columns             = NULL,
                           probability_to_percentage  = FALSE,
                           convert_excluded_to_3      = TRUE,
                           convert_controls_to_4      = TRUE,
                           add_primers                = FALSE,
                           allow_curated              = FALSE,
                           add_padding_between_plates = FALSE
                           ) {
  # Requires the object 'source_abbreviations_vec' in the global environment

  is_CRISPRko <- "Exon_number_GPP" %in% names(my_df)

  if (is_CRISPRko) {
    source_abbreviations_vec <- c(
      "Brunello" = "Bru",
      "TKOv3"    = "TKOv3",
      "GPP"      = "GPP"
    )
  } else {
    if ("Entrez_source_Dolcetto" %in% names(my_df)) {
      source_abbreviations_vec <- c(
        "Curated"       = "Cur",
        "Dolcetto"      = "Dol",
        "hCRISPRi-v2.1" = "hC-v2",
        "hCRISPRi-v2.0" = "hC-v2",
        "GPP"           = "GPP"
      )
    } else {
      source_abbreviations_vec <- c(
        "Curated"     = "Cur",
        "Calabrese"   = "Cal",
        "hCRISPRa-v2" = "hC-v2",
        "GPP"         = "GPP"
      )
    }
  }


  ones_and_zeros_vec <- OnesAndZeros(my_df[["Combined_ID"]])

  assign("delete_my_df", my_df, envir = globalenv())
  assign("delete_first_ones_and_zeros_vec", ones_and_zeros_vec, envir = globalenv())
  if (convert_excluded_to_3) {
    are_to_be_excluded <- !(MeetCriteria(my_df, allow_curated = allow_curated))
    ones_and_zeros_vec[are_to_be_excluded] <- 2L
  }

  my_df <- RoundNumericColumns(my_df, num_digits = 5)

  are_controls <- my_df[["Is_control"]] == "Yes"
  if (any(are_controls)) {
    if (convert_controls_to_4) {
      are_valid_controls <- are_controls & (my_df[["Num_0MM"]] == 0) & (my_df[["Num_1MM"]] == 0) &
                            !(grepl("TTTT", my_df[["sgRNA_sequence"]], ignore.case = TRUE))
      ones_and_zeros_vec[are_valid_controls] <- 3L
    }
    if (!(is_CRISPRko)) {
      transcript_column <- grep("_v2_transcript", colnames(my_df), fixed = TRUE)
      my_df[[transcript_column]][are_controls] <- NA_character_
    }
    my_df[["Gene_symbol"]][are_controls] <- my_df[["Combined_ID"]][are_controls]
    my_df[["Original_symbol"]][are_controls] <- ifelse(grepl("TKOv3", my_df[["Source"]][are_controls], fixed = TRUE), # The TKOv3 library targets EGFP, luciferase, and LacZ
                                                       my_df[["Original_symbol"]][are_controls],
                                                       NA_character_
                                                       )
    if (!(add_primers)) {
      my_df[["Rank"]][are_controls] <- NA_integer_
    }
    if ("GuideScan_offtarget_category" %in% names(my_df)) {
      my_df[["GuideScan_offtarget_category"]][are_controls] <- NA_character_
    }
  }

  my_df[["Num_overlaps"]] <- ifelse(is.na(my_df[["Num_overlaps"]]),
                                    ifelse((my_df[["Rank"]] %in% 1:4) & (my_df[["Is_control"]] == "No"),
                                           ifelse(my_df[["Spacing"]] %in% 12, "<8bp", "n.d."),
                                           NA_character_
                                           ),
                                    paste0(my_df[["Num_overlaps"]], "|", my_df[["Spacing"]], "bp")
                                    )
  if (add_primers) {
    my_df[["Sequence_with_primers"]] <- AddPrimers(my_df)
    my_df <- MoveAfterColumn(my_df, "PAM", "Sequence_with_primers")
  }

  if (all(c("Exon_number_Brunello", "Exon_number_TKOv3", "Exon_number_GPP") %in% names(my_df))) {
    GPP_exon_vec <- my_df[["Exon_number_GPP"]]
    Brunello_exon_vec <- my_df[["Exon_number_Brunello"]]
    TKOv3_exon_vec <- my_df[["Exon_number_TKOv3"]]
    my_df[["Exon_number_Brunello"]] <- vapply(seq_len(nrow(my_df)), function(x) {
      exon_vec <- unique(c(GPP_exon_vec[[x]], Brunello_exon_vec[[x]], TKOv3_exon_vec[[x]])) # This should match the ordering in the "Source" column
      exon_vec <- exon_vec[!(is.na(exon_vec))]
      return(paste0(exon_vec, collapse = " | "))
    }, "")
    names(my_df)[names(my_df) == "Exon_number_Brunello"] <- "Exon_number"
    my_df <- my_df[, !(names(my_df) %in% c("Exon_number_TKOv3", "Exon_number_GPP"))]
  }

  if ("CRISPOR_Graf_status" %in% names(my_df)) {
    my_df[["CRISPOR_Graf_status"]] <- ifelse(my_df[["CRISPOR_Graf_status"]] == "GrafOK", "OK", my_df[["CRISPOR_Graf_status"]])
  }

  SNP_ID_column <- grep("_SNP_IDs_", names(my_df), fixed = TRUE)
  SNP_AF_column <- grep("_SNP_AF_(max|sum)_", names(my_df))
  if ((length(SNP_ID_column) == 1) && (length(SNP_AF_column) == 1)) {
    my_df[[SNP_ID_column]][is.na(my_df[[SNP_AF_column]])] <- NA_character_
    my_df[[SNP_AF_column]][is.na(my_df[[SNP_ID_column]])] <- NA_real_
  }

  if (!(is_CRISPRko)) {
    rank_column <- grep("^(Calabrese|Dolcetto|Caprano|Dolomiti)_", colnames(my_df))
    # This is to prevent the cell value "1/2/3" from being converted into a date by Excel
    my_df[[rank_column]] <- gsub("/", " or ", my_df[[rank_column]], fixed = TRUE)
    my_df[[rank_column]] <- sub(" or ", ", ", my_df[[rank_column]], fixed = TRUE)
  }
  if (probability_to_percentage) {
    for (column_index in grep("_AF_(sum|max)_", names(my_df), fixed = TRUE)) {
      my_df[[column_index]] <- ifelse(is.na(my_df[[column_index]]), NA_character_, paste0(my_df[[column_index]] * 100, "%"))
    }
  }
  my_df[["Locations_0MM"]] <- TruncateLongEntries(my_df[["Locations_0MM"]])

  have_multiple_sources <- grepl(", ", my_df[["Source"]], fixed = TRUE)
  for (source in names(source_abbreviations_vec)) {
    my_df[["Source"]][have_multiple_sources] <- sub(source,
                                                    source_abbreviations_vec[[source]],
                                                    my_df[["Source"]][have_multiple_sources],
                                                    fixed = TRUE
                                                    )
  }
  if (is_CRISPRko) {
    my_df[["Source"]] <- ifelse(my_df[["Source"]] == "GPP, Bru, TKOv3",
                                "GPP, Bru, tk3",
                                my_df[["Source"]]
                                )
  }

  my_df[["Color"]] <- ones_and_zeros_vec + 1L

  if (add_padding_between_plates) {
    my_df <- AddPaddingToPlates(my_df)
    my_df[["Color"]][is.na(my_df[["Color"]])] <- 5L
  }

  if (!(is.null(remove_columns))) {
    my_df <- my_df[, !(names(my_df) %in% remove_columns)]
  }
  for (i in seq_len((ncol(my_df) - 6))) {
    my_df[[i]] <- ifelse(is.na(my_df[[i]]), "", as.character(my_df[[i]]))
  }
  fill_indices <- (ncol(my_df) - 5):ncol(my_df)
  if (is_CRISPRko) {
    fill_columns <- c("Other_Entrez_IDs_4sg", "Other_target_symbols",
                      "Other_symbols_4sg", "Original_symbol"
                      )
    fill_indices <- union(fill_indices, which(names(my_df) %in% fill_columns))
  }
  for (i in fill_indices) {
    my_df[[i]] <- ifelse(is.na(my_df[[i]]) | (my_df[[i]] %in% ""),
                         " ",
                         as.character(my_df[[i]])
                         )
  }

  assign("delete_my_df", my_df, envir = globalenv())
  assign("delete_ones_and_zeros_vec_99", ones_and_zeros_vec, envir = globalenv())

  return(my_df)
}





DfToTSV <- function(CRISPR_df,
                    file_name,
                    remove_columns             = full_omit_columns,
                    probability_to_percentage  = FALSE,
                    add_primers                = FALSE,
                    allow_curated              = FALSE,
                    add_padding_between_plates = FALSE
                    ) {
  formatted_df <- FormatForExcel(CRISPR_df,
                                 remove_columns             = remove_columns,
                                 probability_to_percentage  = probability_to_percentage,
                                 add_primers                = add_primers,
                                 allow_curated              = allow_curated,
                                 add_padding_between_plates = add_padding_between_plates
                                 )
  # Requires the objects 'full_omit_columns' and 'file_output_directory' in the global environment
  write.table(formatted_df,
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )
  return(invisible(formatted_df))
}



# AbbreviateColumns <- function(CRISPR_df, max_entries = 10) {
#   for (column_name in c("Locations_0MM", "Locations_1MM", "Sequences_1MM")) {
#     my_lengths <- lengths(strsplit(CRISPR_df[[column_name]], "; ", fixed = TRUE))
#     are_too_long <- my_lengths > max_entries
#     if (any(are_too_long)) {
#       CRISPR_df[[column_name]][are_too_long] <- "too long"
#       message(paste0(sum(are_too_long), " rows in the ",
#                      column_name, " column were omitted, because they were too long (>",
#                      max_entries, " entries)."
#                      )
#               )
#     }
#   }
#   return(CRISPR_df)
# }


primer_sequences <- list(
  "sg1" = c("ttgtggaaaggacgaaacaccG", "GTTTAAGAGCTAAGCTG"),
  "sg2" = c("cttggagaaaagccttgtttG", "GTTTGAGAGCTAAGCAGA"),
  "sg3" = c("gtatgagaccactctttcccG", "GTTTCAGAGCTAAGCACA"),
  "sg4" = c("gccgcttgggtacctcG", "GTTTCAGAGCTACAGCAGAAAT")
)


AddPrimers <- function(CRISPR_df) {
  results_vec <- rep.int(NA_character_, nrow(CRISPR_df))
  for (i in 1:4) {
    are_this_rank <- CRISPR_df[["Rank"]] %in% i
    new_sequences <- paste0(primer_sequences[[i]][[1]], CRISPR_df[["sgRNA_sequence"]][are_this_rank], primer_sequences[[i]][[2]])
    if (i == 4) {
      new_sequences <- Biostrings::reverseComplement(Biostrings::DNAStringSet(new_sequences))
    }
    results_vec[are_this_rank] <- new_sequences
  }
  return(results_vec)
}



FormatFixedWidthInteger <- function(integer_vec) {
  integer_width <- max(nchar(as.character(as.integer(integer_vec))))
  result <- formatC(integer_vec, width = integer_width, flag = "0")
  return(result)
}




AddPaddingToPlates <- function(CRISPR_df) {
  stopifnot("Plate_number" %in% colnames(CRISPR_df))
  CheckThatFactorIsInOrder(CRISPR_df[["Plate_number"]])
  padding_mat <- matrix(nrow = 10,
                        ncol = ncol(CRISPR_df),
                        dimnames = list(NULL, colnames(CRISPR_df))
                        )
  CRISPR_plates_df_list <- split(CRISPR_df, CRISPR_df[["Plate_number"]])
  num_plates <- length(CRISPR_plates_df_list)
  use_seq <- seq_len(num_plates - 1)
  CRISPR_plates_df_list <- c(unlist(lapply(use_seq, function(x) list(CRISPR_plates_df_list[[x]], padding_mat)),
                                    use.names = FALSE, recursive = FALSE
                                    ),
                             CRISPR_plates_df_list[num_plates]
                             )
  results_df <- do.call(rbind.data.frame,
                        c(CRISPR_plates_df_list,
                          list(make.row.names = FALSE, stringsAsFactors = FALSE)
                          )
                        )
  return(results_df)
}




ExportSharedDf <- function(shared_df, file_name) {
  if ("AltTSS_ID" %in% names(shared_df)) {
    shared_df[["AltTSS_ID"]] <- ifelse(shared_df[["AltTSS_ID"]] == shared_df[["Entrez_ID"]],
                                       " ",
                                       shared_df[["AltTSS_ID"]]
                                       )
  }
  write.table(shared_df,
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )

  return(invisible(NULL))
}








