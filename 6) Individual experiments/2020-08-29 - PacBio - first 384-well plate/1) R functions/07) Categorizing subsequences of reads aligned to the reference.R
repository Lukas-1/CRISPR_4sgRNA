### 20th October 2020 ###





# Import packages and source code -----------------------------------------

library("Biostrings")




# Define feature positions ------------------------------------------------

features_list <- list(
  "column_barcode" = c(1, 10),
  "column_primer"  = c(11, 30),

  "promoter1_hU6"  = c(177, 426),
  "sg1"            = c(427, 446),
  "sg1_cr1"        = c(427, 532),
  "tracrRNA1"      = c(447, 532),

  "TpR_DHFR"       = c(606, 842),

  "promoter2_mU6"  = c(850, 1165),
  "sg2"            = c(1166, 1185),
  "sg2_cr2"        = c(1166, 1273),
  "tracrRNA2"      = c(1186, 1273),

  "promoter3_hH1"  = c(1281, 1504),
  "sg3"            = c(1505, 1524),
  "sg3_cr3"        = c(1505, 1612),
  "tracrRNA3"      = c(1525, 1612),

  "promoter4_h7SK" = c(1620, 1863),
  "sg4"            = c(1864, 1883),
  "sg4_cr4"        = c(1864, 1969),
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



ExtractAlignedSequences <- function(ccs_df, alignments_df, wells_vec = seq_len(384)) {

  stopifnot(all(c("features_df", "features_mat_list", "barcoded_plasmids") %in% ls(envir = globalenv())))

  num_features <- nrow(features_df)
  features_vec <- features_df[["Feature"]]

  well_df_list <- lapply(wells_vec, function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- alignments_df[["Well_number"]] == well_number

    aligned_plasmid_vec <- alignments_df[["Aligned_plasmid"]][are_this_well]
    aligned_read_vec <- alignments_df[["Aligned_read"]][are_this_well]

    aligned_plasmid_char_list <- strsplit(aligned_plasmid_vec, "")
    aligned_read_char_list <- strsplit(aligned_read_vec, "")

    this_well_zmws <- alignments_df[["ZMW"]][are_this_well]
    zmw_matches <- match(this_well_zmws, ccs_df[["ZMW"]])
    unaligned_qual_vec <- ccs_df[["Quality"]][zmw_matches]
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
      "ZMW"         = rep(this_well_zmws, each = num_features),
      "Feature"     = rep(features_vec, times = num_reads),
      "Template"    = rep(features_templates_list[[well_number]], times = num_reads),
      well_mat,
      stringsAsFactors = FALSE
    )
    return(well_df)
  })
  message("Collating the final data frame...")
  results_df <- do.call(rbind.data.frame, c(well_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  message("Processing the final data frame...")
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

  category_vec <- ifelse(results_df[["Is_correct"]],
                         "Correct",
                         ifelse(results_df[["Mostly_deleted"]],
                                "Deletion",
                                "Mutation"
                                )
                         )


  ## Check for contaminations

  contamination_possible <- !(results_df[["Is_correct"]] | results_df[["Mostly_deleted"]])

  sg_list <- lapply(1:4, function(x) sg_sequences_df[[paste0("Sequence_sg", x)]])

  all_guides_vec <- toupper(unlist(sg_list))

  tracrRNA_matches <- match(paste0("tracrRNA", 1:4), features_df[["Feature"]])
  tracrRNA_sequences <- features_templates_list[[1]][tracrRNA_matches]
  all_tracRNAs_vec <- unlist(lapply(tracrRNA_sequences, function(x) paste0(all_guides_vec, x)))

  are_contamination_sg <- rep(FALSE, nrow(results_df))
  for (i in 1:4) {
    are_this_sg <- extracted_df[["Feature"]] == paste0("sg", i)
    are_eligible <- are_this_sg & contamination_possible
    are_contamination_sg[are_eligible] <- extracted_df[["Aligned_read"]][are_eligible] %in% all_guides_vec
  }
  are_contamination_sg_cr <- rep(FALSE, nrow(results_df))
  for (i in 1:4) {
    are_this_sg_cr <- extracted_df[["Feature"]] == paste0("sg", i, "_cr", i)
    are_eligible <- are_this_sg_cr & contamination_possible
    are_contamination_sg_cr[are_eligible] <- extracted_df[["Aligned_read"]][are_eligible] %in% all_tracRNAs_vec
  }
  stopifnot(!(any(are_contamination_sg & are_contamination_sg_cr)))
  are_contamination <- are_contamination_sg | are_contamination_sg_cr
  stopifnot(all(category_vec[are_contamination] == "Mutation"))
  category_vec[are_contamination] <- "Contamination"


  ## Check for small insertions that may or may not be considered to fall within
  ## the sgRNA + tracrRNA region
  wells_vec <- unique(extracted_df[["Well_number"]])
  insertion_possible <- !(results_df[["Is_correct"]] | results_df[["Mostly_deleted"]] | are_contamination)
  have_flanking_insertion <- rep(FALSE, nrow(results_df))
  for (include_tracrRNA in c(FALSE, TRUE)) {
    for (i in 1:4) {
      if (include_tracrRNA) {
        this_feature <- paste0("sg", i, "_cr", i)
      } else {
        this_feature <- paste0("sg", i)
      }
      message(paste0("Checking for flanking insertions for ", this_feature, "..."))
      are_this_feature <- extracted_df[["Feature"]] == this_feature
      for (well_number in wells_vec) {
        are_this_well <-(extracted_df[["Well_number"]] == well_number) & are_this_feature
        well_template <- unique(extracted_df[["Template"]][are_this_well])
        are_selected <- are_this_well & insertion_possible
        have_flanking_insertion[are_selected] <- grepl(well_template,
                                                       extracted_df[["Aligned_read"]][are_selected],
                                                       fixed = TRUE
                                                       )
      }
    }
  }

  category_vec[have_flanking_insertion] <- "Flanking insertion"

  results_df[["Category"]] <- category_vec
  return(results_df)
}


