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
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed - ccs3.RData"))





# Define functions --------------------------------------------------------

ExtractAlignedSequences <- function(use_sl7 = TRUE) {

  barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

  if (use_sl7) {
    ccs_list <- sl7_ccs3_ccs
    report_df <- sl7_ccs3_report_df
  } else {
    ccs_list <- sl9_ccs3_ccs
    report_df <- sl9_ccs3_report_df
  }

  ccs_zmws <- as.integer(substr(ccs_list[["qname"]],  22, nchar(ccs_list[["qname"]])  - 4))
  ccs_well_numbers <- GetWellNumbers(report_df)

  stopifnot(identical(ccs_zmws, as.integer(substr(report_df[[1]], 22, nchar(report_df[[1]])))))

  alignments_list <- lapply(seq_len(384), function(well_number) {

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

    meta_df <- data.frame("ZMW"             = ccs_zmws[are_this_well],
                          "Orientation_fwd" = are_forward_vec,
                          "Score_fwd"       = score(fwd_alignments),
                          "Score_rev"       = score(rev_alignments),
                          stringsAsFactors = FALSE
                          )

    new_indices <- c(which(are_forward_vec), which(!(are_forward_vec)))

    combined_alignments <- c(fwd_alignments[are_forward_vec],
                             rev_alignments[!(are_forward_vec)]
                             )
    combined_alignments <- combined_alignments[order(new_indices)]

    results_list <- list(
      "meta_df" = meta_df,
      "alignments" = combined_alignments
    )
    return(results_list)
  })
  return(alignments_list)
}




# Extract barcodes --------------------------------------------------------

sl7_alignments_list <- ExtractAlignedSequences(use_sl7 = TRUE)
sl9_alignments_list <- ExtractAlignedSequences(use_sl7 = FALSE)




# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_alignments_list"),
     file = file.path(R_objects_directory, "7) Perform pairwise alignments with the reference sequence.RData")
     )








