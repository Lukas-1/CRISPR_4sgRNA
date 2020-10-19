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
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed - ccs3.RData"))
load(file.path(R_objects_directory, "7) Perform pairwise alignments with the reference sequence.RData"))





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
    alignments_list <- sl7_alignments_list
  } else {
    alignments_list <- sl9_alignments_list
  }

  num_features <- nrow(features_df)

  well_list <- lapply(seq_len(384), function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    well_alignments <- alignments_list[[well_number]][["alignments"]]

    aligned_plasmid_vec <- as.character(alignedSubject(well_alignments))
    aligned_read_vec <- as.character(alignedPattern(well_alignments))

    aligned_plasmid_char_list <- strsplit(aligned_plasmid_vec, "")
    aligned_read_char_list <- strsplit(aligned_read_vec, "")

    extracted_list <- lapply(seq_along(well_alignments), function(x) {

      original_numbers <- cumsum(aligned_plasmid_char_list[[x]] != "-")
      # stopifnot(max(original_numbers) == 2245)
      # stopifnot(original_numbers[[length(original_numbers)]] == 2245)

      extracted_sequences <- vapply(features_indices_list[[well_number]],
                                    function(y) paste0(aligned_read_char_list[[x]][original_numbers %in% y], collapse = ""),
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
  extracted_sequences_mat <- do.call(rbind, lapply(well_list, function(x) x[["extracted_sequences_mat"]]))
  are_the_same_mat <- do.call(rbind, lapply(well_list, function(x) x[["are_the_same_mat"]]))
  are_mostly_lost_mat <- do.call(rbind, lapply(well_list, function(x) x[["are_mostly_lost_mat"]]))
  colnames(extracted_sequences_mat) <- features_df[["Features"]]
  colnames(are_the_same_mat)        <- features_df[["Features"]]
  colnames(are_mostly_lost_mat)     <- features_df[["Features"]]
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






# Extract sequences -------------------------------------------------------

sl7_alignments_df <- ExtractAlignedSequences(use_sl7 = TRUE)
sl9_alignments_df <- ExtractAlignedSequences(use_sl7 = FALSE)






# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_alignments_df"),
     file = file.path(R_objects_directory, "9) Extract barcode sequences and quality scores.RData")
     )











