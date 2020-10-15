### 5th October 2020 ###




# Import packages and source code -----------------------------------------

library("Rsamtools")
library("DECIPHER")





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
file_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory      <- file.path(file_directory, "2) Input")
R_objects_directory       <- file.path(file_directory, "3) R objects")
file_output_directory     <- file.path(file_directory, "5) Output")
four_seq_output_directory <- file.path(file_output_directory, "4 sequences")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed.RData"))
load(file.path(R_objects_directory, "9) Process demultiplexed PacBio reads.RData"))





# Define functions --------------------------------------------------------


GetReadForZMW <- function(reads_list, zmw) {
  all_zmws <- as.integer(substr(reads_list[["qname"]], 22, nchar(reads_list[["qname"]]) - 4))
  zmw_index <- match(zmw, all_zmws)
  read_sequence <- reads_list[["seq"]][zmw_index]
  return(read_sequence)
}





Get4Sequences <- function(use_ccs3 = TRUE, use_sl7 = TRUE) {

  if (use_ccs3) {
    if (use_sl7) {
      ccs_list <- sl7_ccs3_ccs
      lima_list <- sl7_ccs3_lima
      reads_df <- sl7_ccs3_df_list[["individual_reads_df"]]
      output_folder <- "SmrtLink7_CCS3"
    } else {
      ccs_list <- sl9_ccs3_ccs
      lima_list <- sl9_ccs3_lima
      reads_df <- sl9_ccs3_df_list[["individual_reads_df"]]
      output_folder <- "SmrtLink7_CCS5"
    }
  } else {
    if (use_sl7) {
      ccs_list <- sl7_ccs5_ccs
      lima_list <- sl7_ccs5_lima
      reads_df <- sl7_ccs5_df_list[["individual_reads_df"]]
      output_folder <- "SmrtLink9_CCS3"
    } else {
      ccs_list <- sl9_ccs5_ccs
      lima_list <- sl9_ccs5_lima
      reads_df <- sl9_ccs5_df_list[["individual_reads_df"]]
      output_folder <- "SmrtLink9_CCS3"
    }
  }

  lima_zmws <- as.integer(substr(lima_list[["qname"]], 22, nchar(lima_list[["qname"]]) - 4))
  ccs_zmws  <- as.integer(substr(ccs_list[["qname"]],  22, nchar(ccs_list[["qname"]])  - 4))

  well_sequences_list <- lapply(seq_len(384), function(well_number) {

    message(paste0("Processing well #", well_number, "..."))

    are_this_well <- reads_df[["Well_number"]] == well_number
    this_well_zmws <- reads_df[["ZMW"]][are_this_well]

    lima_matches <- match(this_well_zmws, lima_zmws)
    stopifnot(!(anyNA(lima_matches)))
    lima_seq <- lima_list[["seq"]][lima_matches]
    lima_qual <- lima_list[["qual"]][lima_matches]

    ccs_matches <- match(this_well_zmws, ccs_zmws)
    stopifnot(!(anyNA(ccs_matches)))
    ccs_seq <- ccs_list[["seq"]][ccs_matches]
    ccs_qual <- ccs_list[["qual"]][ccs_matches]

    well_dir <- file.path(four_seq_output_directory,
                          output_folder,
                          paste0("Well", formatC(well_number, flag = "0", width = 3))
                          )
    dir.create(well_dir)

    qual_list <- lapply(seq_along(this_well_zmws), function(x) {

      # fwd_patterns <- c(ccs_seq[x], lima_seq[x])
      # rev_patterns <- reverseComplement(fwd_patterns)
      # fwd_align <- pairwiseAlignment(fwd_patterns[1], plasmid, type = "global")
      # rev_align <- pairwiseAlignment(rev_patterns[2], plasmid, type = "global")
      # is_forward <- sum(score(fwd_align)) > sum(score(rev_align))

      # sub_seq <- c(DNAStringSet(c(barcoded_plasmids_NNN[[well_number]],
      #                             barcoded_plasmids[[well_number]]
      #                           )),
      #              ccs_seq[x], lima_seq[x] # if (is_forward) fwd_patterns else rev_patterns
      #              )
      # sub_qual <- c(PhredQuality(".")[c(1, 1)], ccs_qual[x], lima_qual[x])
      # result_object <- QualityScaledBStringSet(sub_seq, sub_qual)
      # names(result_object) <- paste0(this_well_zmws[[x]],
      #                                c("_plasmid_NNN",
      #                                  "_plasmid",
      #                                  "_ccs",
      #                                  "_demux"
      #                                  )
      #                                )

      sub_seq <- c(ccs_seq[x], lima_seq[x])
      sub_qual <- c(ccs_qual[x], lima_qual[x])
      result_object <- QualityScaledBStringSet(sub_seq, sub_qual)
      names(result_object) <- paste0(this_well_zmws[[x]],
                                     c("_ccs",
                                       "_demux"
                                       )
                                     )
      # widths <- width(sub_seq)
      assign("delete_sub_seq", sub_seq, envir = globalenv())
      assign("delete_sub_qual", sub_qual, envir = globalenv())
      stop()
      # stopifnot((widths[[1]] - 20) == widths[[2]])
      # ccs_subseq <- sub_seq[[1]][11:(widths[[1]] - 10)]
      # stopifnot(identical(as.character(ccs_subseq), as.character(sub_seq[[2]])))
      file_path <- file.path(well_dir, paste0(this_well_zmws[[x]], ".fastq"))
      writeQualityScaledXStringSet(result_object, filepath = file_path)
      return(result_object)
    })
    names(qual_list) <- this_well_zmws
    return(qual_list)
  })
  return(well_sequences_list)
}















AlignmentsForZMW <- function(zmw, use_ccs3 = TRUE, use_sl7 = TRUE, display_alignment = TRUE) {

  if (!(is.integer(zmw))) {
    zmw <- as.integer(zmw)
  }

  if (use_ccs3) {
    if (use_sl7) {
      ccs_list <- sl7_ccs3_ccs
      lima_list <- sl7_ccs3_lima
      reads_df <- sl7_ccs3_df_list[["individual_reads_df"]]
    } else {
      ccs_list <- sl9_ccs3_ccs
      lima_list <- sl9_ccs3_lima
      reads_df <- sl9_ccs3_df_list[["individual_reads_df"]]
    }
  } else {
    if (use_sl7) {
      ccs_list <- sl7_ccs5_ccs
      lima_list <- sl7_ccs5_lima
      reads_df <- sl7_ccs5_df_list[["individual_reads_df"]]
    } else {
      ccs_list <- sl9_ccs5_ccs
      lima_list <- sl9_ccs5_lima
      reads_df <- sl9_ccs5_df_list[["individual_reads_df"]]
    }
  }

  zmw_index <- match(zmw, reads_df[["ZMW"]])
  if (is.na(zmw_index)) {
    stop(paste0("The zmw'", zmw, "' was not found!"))
  }
  use_well <- reads_df[["Well_number"]][zmw_index]

  plasmid <- DNAStringSet(barcoded_plasmids[[use_well]])

  consensus_sequence <- GetReadForZMW(ccs_list, zmw)
  demultiplexed_sequence <- GetReadForZMW(lima_list, zmw)
  assign("delete_seq",
         c("plasmid"       = plasmid,
           "ccs"           = as.character(consensus_sequence),
           "demultiplexed" = as.character(demultiplexed_sequence)
         ),
         envir = globalenv()
         )
  fwd_patterns <- c(consensus_sequence, demultiplexed_sequence)
  rev_patterns <- reverseComplement(fwd_patterns)

  fwd_align <- pairwiseAlignment(fwd_patterns, plasmid, type = "global")
  rev_align <- pairwiseAlignment(rev_patterns, plasmid, type = "global")

  if (sum(score(fwd_align)) > sum(score(rev_align))) {
    is_forward <- TRUE
    message("Forward orientation")
    align_result <- fwd_align
  } else {
    is_forward <- FALSE
    message("Reverse orientation")
    align_result <- rev_align
  }

  assign("delete_align_result", align_result, envir = globalenv())

  if (display_alignment) {
    plasmid_NNN <- DNAStringSet(barcoded_plasmids_NNN[[use_well]])
    if (is_forward) {
      use_patterns <- fwd_patterns
    } else {
      use_patterns <- rev_patterns
    }
    use_patterns <- c(plasmid_NNN, use_patterns)
    assign("delete_seq", c(plasmid, use_patterns), envir = globalenv())
    names(use_patterns) <- c("Plasmid_NNN", "Consensus", "Demultiplexed")
    names(plasmid) <- "Plasmid"
    display_align <- pairwiseAlignment(use_patterns, plasmid, type = "global")
    # ExamineAlignment(display_align)
  }
  return(align_result)
}




ExamineAlignment <- function(alignment) {

  assign("delete_alignment", alignment, envir = globalenv())

  plasmid <- alignedSubject(alignment)
  if (length(plasmid) == 3) {
    plasmid <- plasmid[2:3]
  }
  plasmid <- unique(plasmid)
  stopifnot(length(plasmid) == 1)

  reads <- alignedPattern(alignment)
  actual_lengths <- vapply(strsplit(as.character(reads), ""), function(x) sum(x != "-"), integer(1))
  assign("delete_actual_lengths", actual_lengths, envir = globalenv())
  stopifnot(identical(actual_lengths[[2]] - 20L, actual_lengths[[3]]))

  all_sequences <- c(plasmid, reads)

  BrowseSeqs(all_sequences)
  return(all_sequences)
}






# Add barcodes to plasmid sequences ---------------------------------------

barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)
barcoded_plasmids_NNN <- paste0(column_bc_vec, plasmid_string, row_bc_vec)






# Try stuff ---------------------------------------------------------------

# sl7_ccs3_4seqs <- Get4Sequences(use_sl7 = TRUE, use_ccs3 = TRUE)
# sl7_ccs5_4seqs <- Get4Sequences(use_sl7 = TRUE, use_ccs3 = FALSE)
# sl9_ccs3_4seqs <- Get4Sequences(use_sl7 = FALSE, use_ccs3 = TRUE)
# sl9_ccs5_4seqs <- Get4Sequences(use_sl7 = FALSE, use_ccs3 = FALSE)









# Try stuff ---------------------------------------------------------------


# alg <- AlignmentsForZMW(5112773)
#
# newseq <- unique(c(alignedPattern(alg), alignedSubject(alg)))
#
#
# AlignmentsForZMW(4194619, try_reverse = FALSE)
#
#
# head(sl7_ccs3_df_list[["individual_reads_df"]])
#
#
# pairwiseAlignment(as.character(sl7_ccs3_ccs[[2]][[1]]), barcoded_plasmids[[1]])
#

alignment <- AlignmentsForZMW(4194619)

alignment <- AlignmentsForZMW(5112773)

AlignmentsForZMW(5112773)

alignment <- AlignmentsForZMW(32965024)

BrowseSeqs(alignment)



alignment <- AlignmentsForZMW(59244630)



show_seq <- c(alignedPattern(alignment[1]), alignedSubject(alignment[1]))



BrowseSeqs(show_seq)

