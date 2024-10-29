### 2024-08-11


# Load packages and source code -------------------------------------------

library("Rsamtools")
library("GenomicAlignments")
library("ShortRead")




# Define functions --------------------------------------------------------

ExtractCoords <- function(GAlignments_object, use_coords, columns_postfix = NULL) {

  stopifnot((length(use_coords) == 2) && is.numeric(use_coords))

  sequence_name <- seqlevels(GAlignments_object)
  stopifnot(length(sequence_name) == 1)

  if (is.null(names(GAlignments_object))) {
    names(GAlignments_object) <- mcols(GAlignments_object)[, "qname"]
  }
  stopifnot(!(anyDuplicated(names(GAlignments_object))))

  use_roi <- GRanges(sequence_name, IRanges(start = use_coords[[1]], end = use_coords[[2]]))
  mapped_ranges <- mapToAlignments(use_roi, GAlignments_object)

  matches_vec <- match(as.character(seqnames(mapped_ranges)), names(GAlignments_object))

  sequences_vec <- as.character(subseq(mcols(GAlignments_object)[, "seq"][matches_vec],
                                       start(mapped_ranges), end(mapped_ranges)
                                       ))

  phred_qualities <- subseq(mcols(GAlignments_object)[, "qual"][matches_vec],
                            start(mapped_ranges), end(mapped_ranges)
                            )
  mean_qualities <- GetMeanQuality(phred_qualities)

  results_df <- data.frame(
    "Name"         = mcols(GAlignments_object)[, "qname"][matches_vec],
    "Sequence"     = sequences_vec,
    "Quality"      = as.character(phred_qualities),
    "Mean_quality" = mean_qualities,
    stringsAsFactors = FALSE
  )
  if (!(is.null(columns_postfix))) {
    names(results_df) <- paste0(names(results_df), columns_postfix)
  }
  return(results_df)
}



ReadBamData <- function(bam_paths) {
  param1 <- ScanBamParam(flag = scanBamFlag(isSupplementaryAlignment = FALSE,
                                            isSecondaryAlignment = FALSE,
                                            isMinusStrand = FALSE
                                            ),
                         what = c("qname", "flag", "seq", "qual")
                         )

  bam_views <- BamViews(bamPaths = bam_paths)
  bam_data <- readGAlignments(bam_views, param = param1, use.names = FALSE)

  merge_bam_string <- paste0(paste0("bam_data[[", seq_along(bam_paths), "]]"), collapse = ", ")
  merge_bam_string <- paste0("bam_data <- c(", merge_bam_string, ")")

  eval(parse(text = merge_bam_string))
  return(bam_data)
}



BamExtractGuides <- function(bam_data, meta_df) {

  ### Extract sgRNAs
  sg_df_list <- lapply(1:4, function(x) {
    use_sg <- paste0("sg", x)
    message("Extracting ", use_sg, "...")
    ExtractCoords(bam_data, features_list[[use_sg]] - 10, columns_postfix = paste0("_", use_sg))
  })

  all_reads <- unique(unlist(lapply(sg_df_list, function(x) x[, 1])))
  file_numbers <- as.integer(substr(all_reads, 2, 2))
  read_numbers <- as.integer(substr(all_reads, 4, nchar(all_reads)))
  use_order <- order(file_numbers, read_numbers)
  all_reads <- all_reads[use_order]
  # all_read_numbers <- as.integer(substr(all_reads, 4, nchar(all_reads)))
  # all_read_numbers <- sort(all_read_numbers)

  sg_df_list <- lapply(sg_df_list, function(x) {
    matches_vec <- match(all_reads, x[, 1])
    x[matches_vec, 2:ncol(x)]
  })

  ### Combine data
  num_per_file <- table(meta_df[, "Part_number"])
  demuxed_read_IDs <- paste0("p", meta_df[, "Part_number"], "_",
                             unlist(lapply(num_per_file, seq_len), use.names = FALSE)
                             )
  matches_vec <- match(all_reads, demuxed_read_IDs)

  extracted_df <- data.frame("File_number" = meta_df[, "File_number"][matches_vec],
                             "Read_number" = read_numbers[use_order],
                             meta_df[matches_vec, names(meta_df) != "File_number"],
                             lapply(sg_df_list, "[", 3),
                             lapply(sg_df_list, "[", 1),
                             lapply(sg_df_list, "[", 2),
                             stringsAsFactors = FALSE,
                             row.names = NULL
                             )

  return(extracted_df)
}


