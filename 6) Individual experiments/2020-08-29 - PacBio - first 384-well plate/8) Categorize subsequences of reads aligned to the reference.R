### 15th October 2020 ###



# Import packages and source code -----------------------------------------

library("Rsamtools")

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "2) Functions for analyzing reads.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "3) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed.RData"))





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



ExtractAlignedSequences <- function(use_ccs3 = TRUE, use_sl7 = TRUE) {

  stopifnot(all(c("features_df", "features_mat_list", "barcoded_plasmids") %in% ls(envir = globalenv())))


  if (use_ccs3) {
    if (use_sl7) {
      ccs_list <- sl7_ccs3_ccs
      report_df <- sl7_ccs3_report_df
    } else {
      ccs_list <- sl9_ccs3_ccs
      report_df <- sl9_ccs3_report_df
    }
  } else {
    if (use_sl7) {
      ccs_list <- sl7_ccs5_ccs
      report_df <- sl7_ccs5_report_df
    } else {
      ccs_list <- sl9_ccs5_ccs
      report_df <- sl9_ccs5_report_df
    }
  }

  ccs_zmws <- as.integer(substr(ccs_list[["qname"]],  22, nchar(ccs_list[["qname"]])  - 4))
  ccs_well_numbers <- GetWellNumbers(report_df)

  stopifnot(identical(ccs_zmws, as.integer(substr(report_df[[1]], 22, nchar(report_df[[1]])))))

  num_features <- nrow(features_df)

  well_barcodes_df_list <- lapply(seq_len(384), function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- ccs_well_numbers %in% well_number
    ccs_seq <- ccs_list[["seq"]][are_this_well]
    ccs_qual <- ccs_list[["qual"]][are_this_well]

    plasmid <- DNAStringSet(barcoded_plasmids[[well_number]])
    row_template <- row_bc_vec[[well_number]]
    column_template <- column_bc_vec[[well_number]]

    fwd_alignments <- pairwiseAlignment(ccs_seq, plasmid, type = "global")
    rev_alignments <- pairwiseAlignment(reverseComplement(ccs_seq), plasmid, type = "global")
    are_forward_vec <- score(fwd_alignments) > score(rev_alignments)

    extracted_list <- lapply(seq_len(sum(are_this_well)), function(x) {

      if (are_forward_vec[[x]]) {
        use_alignment <- fwd_alignments[x]
      } else {
        use_alignment <- rev_alignments[x]
      }

      aligned_plasmid_seq <- alignedSubject(use_alignment)
      aligned_plasmid_chars <- strsplit(as.character(aligned_plasmid_seq), "")[[1]]

      aligned_read_seq <- alignedPattern(use_alignment)
      aligned_read_chars <- strsplit(as.character(aligned_read_seq), "")[[1]]

      original_numbers <- cumsum(aligned_plasmid_chars != "-")
      # stopifnot(max(original_numbers) == 2245)
      # stopifnot(original_numbers[[length(original_numbers)]] == 2245)

      extracted_sequences <- vapply(features_indices_list[[well_number]],
                                    function(x) paste0(aligned_read_chars[original_numbers %in% x], collapse = ""),
                                    ""
                                    )
      are_the_same <- vapply(seq_len(num_features),
                             function(x) extracted_sequences[[x]] == features_templates_list[[well_number]][[x]],
                             logical(1)
                             )

      are_mostly_lost <- vapply(extracted_sequences,
                                function(x) sum(strsplit(x, "") == "-") > (nchar(x) / 2),
                                logical(1)
                                )

      results_list <- list(
        "extracted_sequences" = extracted_sequences,
        "are_the_same"        = are_the_same,
        "are_mostly_lost"     = are_mostly_lost
      )

      return(results_list)
    })
    extracted_sequences_mat <- do.call(rbind, lapply(extracted_list, function(x) x[["extracted_sequences"]]))
    are_the_same_mat <- do.call(rbind, lapply(extracted_list, function(x) x[["are_the_same"]]))
    are_mostly_lost_mat <- do.call(rbind, lapply(extracted_list, function(x) x[["are_mostly_lost"]]))
    results_list <- list("extracted_sequences_mat" = extracted_sequences_mat,
                         "are_the_same_mat"        = are_the_same_mat,
                         "are_mostly_lost_mat"     = are_mostly_lost_mat
                         )
    return(results_list)
  })
  extracted_sequences_mat <- do.call(rbind, lapply(results_list, function(x) x[["extracted_sequences_mat"]]))
  are_the_same_mat <- do.call(rbind, lapply(results_list, function(x) x[["are_the_same_mat"]]))
  extracted_sequences_mat <- do.call(rbind, lapply(results_list, function(x) x[["are_mostly_lost_mat"]]))
  colnames(extracted_sequences_mat) <- features_df[["Features"]]
  colnames(are_the_same_mat)        <- features_df[["Features"]]
  colnames(extracted_sequences_mat) <- features_df[["Features"]]
  results_list <- list("extracted_sequences_mat" = extracted_sequences_mat,
                       "are_the_same_mat"        = are_the_same_mat,
                       "are_mostly_lost_mat"     = are_mostly_lost_mat
                       )

  return(results_list)
}






# Define all important features, and generate per-well coordinates --------

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

features_mat <- do.call(rbind, features_list)
mode(features_mat) <- "integer"
rownames(features_mat) <- NULL
colnames(features_mat) <- c("Start", "End")
features_df <- data.frame("Feature" = names(features_list),
                          features_mat,
                          stringsAsFactors = FALSE
                          )

features_mat_list <- lapply(seq_len(384), function(x) {
  sg_vec <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  AdjustForsgRNALength(features_df, sg_vec)
})

features_indices_list <- lapply(features_mat_list, function(x) {
  lapply(seq_len(nrow(x)), function(y) seq(from = x[y, "Start"], to = x[y, "End"]))
})


barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

features_templates_list <- lapply(seq_len(384), function(x) {
  use_mat <- features_mat_list[[x]]
  vapply(seq_len(nrow(use_mat)),
         function(y) substr(barcoded_plasmids[[x]], use_mat[y, "Start"], use_mat[y, "End"]),
         ""
         )
})








# Extract barcodes --------------------------------------------------------

sl7_ccs3_alignments_df <- ExtractAlignedSequences(use_sl7 = TRUE, use_ccs3 = TRUE)
sl7_ccs5_alignments_df <- ExtractAlignedSequences(use_sl7 = TRUE, use_ccs3 = FALSE)
sl9_ccs3_alignments_df <- ExtractAlignedSequences(use_sl7 = FALSE, use_ccs3 = TRUE)
sl9_ccs5_alignments_df <- ExtractAlignedSequences(use_sl7 = FALSE, use_ccs3 = FALSE)






# Save data ---------------------------------------------------------------

save(list = paste0(c("sl7_ccs3", "sl7_ccs5", "sl9_ccs3", "sl9_ccs5"), "_barcodes_df"),
     file = file.path(R_objects_directory, "7) Extract barcode sequences and quality scores.RData")
     )











