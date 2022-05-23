## 2022-04-08



# Load packages and source code -------------------------------------------

library("Biostrings")
library("ShortRead")



# Define functions --------------------------------------------------------

GetMeanQuality <- function(qualities, rescale = TRUE) {
  if (!(identical("PhredQuality", as.character(class(qualities))))) {
    qualities <- Biostrings::PhredQuality(qualities)
  }
  mean_qualities <- ShortRead::alphabetScore(qualities) / Biostrings::width(qualities)
  if (rescale) {
    mean_qualities <- mean_qualities / 93 * 100
  }
  return(mean_qualities)
}


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


FeaturesListToDf <- function(features_list) {
  features_mat <- do.call(rbind, features_list)
  mode(features_mat) <- "integer"
  rownames(features_mat) <- NULL
  colnames(features_mat) <- c("Start", "End")
  features_df <- data.frame("Feature" = names(features_list),
                            features_mat,
                            stringsAsFactors = FALSE
                            )
  return(features_df)
}


ExtractAlignedSequences <- function(align_df) {

  stopifnot(all(c("features_df", "features_indices_list") %in% ls(envir = globalenv())))

  features_vec <- features_df[["Feature"]]

  aligned_plasmid_char_list <- strsplit(align_df[["Aligned_ref"]], "")
  aligned_read_char_list <- strsplit(align_df[["Aligned_read"]], "")

  unaligned_qual_vec <- align_df[["Read_quality"]]
  are_fwd <- align_df[["Orientation_fwd"]]

  unaligned_qual_vec[!(are_fwd)] <- reverse(unaligned_qual_vec[!(are_fwd)])

  unaligned_qual_char_list <- strsplit(unaligned_qual_vec, "")

  num_reads <- nrow(align_df)

  extracted_mat_list <- lapply(seq_len(num_reads), function(x) {

    plasmid_char_numbers <- cumsum(aligned_plasmid_char_list[[x]] != "-")

    are_gaps <- aligned_read_char_list[[x]] == "-"
    read_char_numbers <- cumsum(!(are_gaps))
    read_char_numbers[read_char_numbers == 0] <- NA
    aligned_qual_char_vec <- unaligned_qual_char_list[[x]][read_char_numbers]
    aligned_qual_char_vec[are_gaps] <- " "
    stopifnot(!(anyNA(aligned_qual_char_vec)))

    original_read_length <- length(unaligned_qual_char_list[[x]])
    stopifnot(max(read_char_numbers, na.rm = TRUE) == original_read_length)
    stopifnot(sum(!(are_gaps)) == original_read_length)

    feature_indices <- lapply(features_indices_list,
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

  results_df <- data.frame(
    "Feature" = rep(features_vec, times = num_reads),
    well_mat,
    stringsAsFactors = FALSE
  )
  return(results_df)
}



ExtractAlignedSgRNAs <- function(alignments_df) {

  ## Prepare data
  alignments_df[, "Mean_quality"] <- GetMeanQuality(alignments_df[, "Read_sequence"])
  features_df <- FeaturesListToDf(features_list)
  for (column in c("Start", "End")) {
    features_df[[column]] <- features_df[[column]] - 10L
  }
  features_df <- features_df[features_df[, "Feature"] %in% paste0("sg", 1:4), ]
  features_indices_list <- lapply(seq_len(nrow(features_df)),
                                  function(y) seq(from = features_df[y, "Start"], to = features_df[y, "End"])
                                  )

  assign("features_df", features_df, envir = globalenv())
  assign("features_indices_list", features_indices_list, envir = globalenv())

  ## Extract aligned sequences
  num_reads <- nrow(alignments_df)
  reads_per_chunk <- 10000
  num_chunks <- ceiling(num_reads / reads_per_chunk)
  chunks_vec <- rep(seq_len(num_chunks), each = reads_per_chunk)[seq_len(num_reads)]
  chunks_list <- vector(mode = "list", length = num_chunks)
  first_vec <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[1]]))
  last_vec  <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[length(x)]]))
  chunk_numbers <- format(seq_len(num_chunks))

  for (i in seq_len(num_chunks)) {
    are_this_chunk <- chunks_vec == i
    message("Processing chunk #", chunk_numbers[[i]], " of ",
            chunk_numbers[[length(chunk_numbers)]],  " (extracting reads ",
            first_vec[[i]], " to ", last_vec[[i]], ")..."
            )
    sub_df <- ExtractAlignedSequences(alignments_df[are_this_chunk, ])
    chunks_list[[i]] <- sub_df
  }

  extracted_df <- do.call(rbind.data.frame,
                          c(chunks_list,
                            stringsAsFactors = FALSE,
                            make.row.names = FALSE
                            )
                          )
  extracted_df[, "Sg_number"] <- as.integer(sub("sg", "", extracted_df[, "Feature"], fixed = TRUE))

  extracted_df_list <- split(extracted_df[, !(names(extracted_df) %in% c("Feature", "Sg_number"))],
                             extracted_df[, "Sg_number"]
                             )

  for (i in 1:4) {
    names(extracted_df_list[[i]]) <- paste0(names(extracted_df_list[[i]]), "_sg", i)
  }

  extracted_df <- data.frame(
    alignments_df[, c("Orientation_fwd", "Score_fwd", "Score_rev", "Mean_quality")],
    extracted_df_list[[1]], extracted_df_list[[2]],
    extracted_df_list[[3]], extracted_df_list[[4]]
  )

  return(extracted_df)
}
