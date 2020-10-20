### 20th October 2020 ###





# Import packages and source code -----------------------------------------

library("Biostrings")




# Define feature positions ------------------------------------------------

features_list <- list(
  "column_barcode" = c(1, 10),
  "column_primer"  = c(11, 30),

  "promoter1_hU6"  = c(177, 426),
  "sg1"            = c(427, 446),
  "tracrRNA1"      = c(447, 532),

  "promoter2_mU6"  = c(850, 1165),
  "sg2"            = c(1166, 1185),
  "tracrRNA2"      = c(1186, 1273),

  "promoter3_hH1"  = c(1281, 1504),
  "sg3"            = c(1505, 1524),
  "tracrRNA3"      = c(1525, 1612),

  "promoter4_h7SK" = c(1620, 1863),
  "sg4"            = c(1864, 1883),
  "tracrRNA4"      = c(1884, 1969),

  "row_primer"     = c(2216, 2235),
  "row_barcode"    = c(2236, 2245)
)





# Define functions --------------------------------------------------------

AdjustForsgRNALength <- function(features_df, sg_vec) {
  stopifnot(length(sg_vec) == 4)
  sg_lengths <- nchar(sg_vec)
  results_mat <- as.matrix(features_df[, c("Start", "End")])
  if (!(all(sg_lengths == 20))) {
    for (i in 1:4) {
      if (sg_lengths[[i]] != 20) {
        sg_row <- match(paste0("sg", i), features_df[, "Feature"])
        addend <- -20L + sg_lengths[[i]]
        results_mat[sg_row, "End"] <- results_mat[sg_row, "End"] + addend
        results_mat[(sg_row + 1):nrow(results_mat), ] <- results_mat[(sg_row + 1):nrow(results_mat), ] + addend
      }
    }
  }
  return(results_mat)
}



ExtractAlignedSequences <- function(use_sl7 = TRUE) {

  stopifnot(all(c("features_df", "features_mat_list", "barcoded_plasmids") %in% ls(envir = globalenv())))

  if (use_sl7) {
    alignments_df <- sl7_alignments_df
    ccs_list <- sl7_ccs3_ccs
  } else {
    alignments_df <- sl9_alignments_df
    ccs_list <- sl9_ccs3_ccs
  }

  qualities_vec <- as.character(ccs_list[["qual"]])
  ccs_zmws  <- as.integer(substr(ccs_list[["qname"]], 22, nchar(ccs_list[["qname"]]) - 4))

  num_features <- nrow(features_df)
  features_vec <- features_df[["Feature"]]

  well_df_list <- lapply(seq_len(384), function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- alignments_df[["Well_number"]] == well_number

    aligned_plasmid_vec <- alignments_df[["Aligned_plasmid"]][are_this_well]
    aligned_read_vec <- alignments_df[["Aligned_read"]][are_this_well]

    aligned_plasmid_char_list <- strsplit(aligned_plasmid_vec, "")
    aligned_read_char_list <- strsplit(aligned_read_vec, "")

    this_well_zmws <- alignments_df[["ZMW"]][are_this_well]
    zmw_matches <- match(this_well_zmws, ccs_zmws)
    unaligned_qual_vec <- qualities_vec[zmw_matches]
    unaligned_qual_char_list <- strsplit(unaligned_qual_vec, "")

    num_reads <- length(this_well_zmws)

    extracted_mat_list <- lapply(seq_len(num_reads), function(x) {

      plasmid_char_numbers <- cumsum(aligned_plasmid_char_list[[x]] != "-")

      are_gaps <- aligned_read_char_list[[x]] == "-"
      read_char_numbers <- cumsum(!(are_gaps))
      aligned_qual_char_vec <- unaligned_qual_char_list[[x]][read_char_numbers]
      aligned_qual_char_vec[are_gaps] <- " "

      original_read_length <- length(unaligned_qual_char_list[[x]])
      stopifnot(max(read_char_numbers) == original_read_length)
      stopifnot(sum(!(are_gaps)) == original_read_length)
      # stopifnot(original_numbers[[length(original_numbers)]] == 2245)

      feature_indices <- lapply(features_indices_list[[well_number]],
                                function(x) which(plasmid_char_numbers %in% x)
                                )

      aligned_templates <- vapply(feature_indices,
                                  function(y) paste0(aligned_plasmid_char_list[[x]][y], collapse = ""),
                                  ""
                                  )
      extracted_sequences <- vapply(feature_indices,
                                    function(y) paste0(aligned_read_char_list[[x]][y], collapse = ""),
                                    ""
                                    )
      extracted_qualities <- vapply(feature_indices,
                                    function(x) paste0(aligned_qual_char_vec[x], collapse = ""),
                                    ""
                                    )

      extracted_mat <- cbind(
        "Aligned_template" = aligned_templates,
        "Aligned_read"     = extracted_sequences,
        "Quality"          = extracted_qualities
      )
      return(extracted_mat)
    })
    well_mat <- do.call(rbind, extracted_mat_list)

    well_df <- data.frame(
      "Well_number" = well_number,
      "ZMW"         = this_well_zmws,
      "Feature"     = rep(features_vec, times = num_reads),
      "Template"    = rep(features_templates_list[[well_number]], times = num_reads),
      well_mat,
      stringsAsFactors = FALSE
    )
    return(well_df)
  })
  message("Collating the final data frame...")
  results_df <- do.call(rbind.data.frame, c(well_df_list, stringsAsFactors = FALSE, make.row.names = FALSE))
  results_df <- ProcessExtractedDf(results_df)
  return(results_df)
}



ProcessExtractedDf <- function(extracted_df) {

  quality_int_list <- as(PhredQuality(extracted_df[["Quality"]]), "IntegerList")
  quality_int_list <- lapply(quality_int_list, function(x) x[x != -1L])

  mean_quality_vec <- vapply(quality_int_list,
                             function(x) if (length(x) == 0) NA_real_ else mean(x),
                             numeric(1)
                             )

  sequence_splits <- strsplit(extracted_df[["Aligned_read"]], "", fixed = TRUE)

  num_bases_lost <- vapply(sequence_splits, function(x) sum(x == "-"), integer(1))

  sequence_lengths <- nchar(extracted_df[["Aligned_read"]])

  results_df <- data.frame(
    extracted_df,
    "Is_correct"     = extracted_df[["Template"]] == extracted_df[["Aligned_read"]],
    "Mean_quality"   = mean_quality_vec,
    "Num_missing"    = num_bases_lost,
    "Mostly_deleted" = num_bases_lost >= (sequence_lengths / 2)
  )
  return(results_df)
}


