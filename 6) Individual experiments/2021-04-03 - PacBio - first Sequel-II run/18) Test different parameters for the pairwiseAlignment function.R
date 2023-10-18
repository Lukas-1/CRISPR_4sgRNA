### 14th May 2021 ###




# Import packages and source code -----------------------------------------

library("Rsamtools")
library("DECIPHER")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))






# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "04) Create reference sequences for each well - constant sequences.RData"))
load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Define functions --------------------------------------------------------

AlignmentsForZMW <- function(zmw, gap_opening = 10, display_alignment = TRUE) {

  stopifnot(all(c("sg_sequences_df", "ccs_df") %in% ls(envir = globalenv())))

  if (!(is.integer(zmw))) {
    zmw <- as.integer(zmw)
  }

  ccs_index <- match(zmw, ccs_df[["ZMW"]])
  if (is.na(ccs_index)) {
    stop(paste0("The zmw'", zmw, "' was not found!"))
  }
  combined_ID <- ccs_df[["Combined_ID"]][[ccs_index]]
  library_index <- match(combined_ID, sg_sequences_df[["Combined_ID"]])

  reference_sequence <- DNAStringSet(sg_sequences_df[library_index, "Barcoded_plasmid"])
  ccs_sequence <- DNAStringSet(ccs_df[ccs_index, "Sequence"])

  rev_ccs_sequence <- reverseComplement(ccs_sequence)

  fwd_align <- pairwiseAlignment(ccs_sequence, reference_sequence,
                                 type = "global", gapOpening = gap_opening
                                 )
  rev_align <- pairwiseAlignment(rev_ccs_sequence, reference_sequence,
                                 type = "global", gapOpening = gap_opening
                                 )

  if (sum(score(fwd_align)) > sum(score(rev_align))) {
    is_forward <- TRUE
    message("Forward orientation")
    align_result <- fwd_align
  } else {
    is_forward <- FALSE
    message("Reverse orientation")
    align_result <- rev_align
  }

  if (display_alignment) {
    plasmid_NNN <- DNAStringSet(sg_sequences_df[library_index, "Barcoded_plasmid_NNN"])
    if (is_forward) {
      use_pattern <- ccs_sequence
    } else {
      use_pattern <- rev_ccs_sequence
    }
    use_patterns <- c(plasmid_NNN, use_pattern)
    names(use_patterns) <- c("Plasmid_NNN", "CCS")
    names(reference_sequence) <- "Plasmid"
    display_align <- pairwiseAlignment(use_patterns, reference_sequence,
                                       type = "global", gapOpening = gap_opening
                                       )
    ExamineAlignment(display_align)
  }
  return(align_result)
}



ExamineAlignment <- function(alignment) {

  plasmid <- alignedSubject(alignment)
  if (length(plasmid) == 3) {
    plasmid <- plasmid[2:3]
  }
  plasmid <- unique(plasmid)
  # stopifnot(length(plasmid) == 1)

  reads <- alignedPattern(alignment)
  actual_lengths <- vapply(strsplit(as.character(reads), ""), function(x) sum(x != "-"), integer(1))

  all_sequences <- c(plasmid, reads)

  BrowseSeqs(all_sequences)
  return(all_sequences)
}





# Prepare lists of ZMWs ---------------------------------------------------

ccs_df[["Passed_filters"]] <- ccs_df[["Plate_passed_filters"]] &
                              (ccs_df[["Well_passed_filters"]] %in% TRUE)

ccs3_zmws <- GetCCS3_ZMWs(ccs_df)
ccs5_zmws <- GetCCS5_ZMWs(ccs_df)
ccs7_zmws <- GetCCS7_ZMWs(ccs_df)

reads_df <- ccs3_df_list[["individual_reads_df"]]

pass_bc   <- reads_df[["Passes_barcode_filters"]] == 1
pass_read <- reads_df[["Passes_read_quality"]] == 1
pass_sg   <- reads_df[["Passes_sg_quality"]] == 1

passing_bc_zmws   <- reads_df[["ZMW"]][pass_bc]
passing_read_zmws <- reads_df[["ZMW"]][pass_bc & pass_read]




# Find problematic ZMWs ---------------------------------------------------

use_df <- ccs7_df_list[["individual_reads_df"]]

four_contam_search <- (use_df[, "Num_contaminating_guides"] >= 4)
four_contam_align  <- (rowSums(use_df[, paste0("sg", 1:4, "_cr", 1:4, "_category")] == "Contamination") == 4)

search_greater_align <- four_contam_search & !(four_contam_align)
align_greater_search <- four_contam_align & !(four_contam_search)

head(use_df[search_greater_align, c("ZMW", "Combined_ID")])
head(use_df[align_greater_search, c("ZMW", "Combined_ID")])





# Try stuff ---------------------------------------------------------------

head(intersect(ccs7_zmws, passing_read_zmws))

show_seq <- AlignmentsForZMW(97978568,  gap_opening = 30)
show_seq <- AlignmentsForZMW(106171035, gap_opening = 30)
show_seq <- AlignmentsForZMW(131269636, gap_opening = 30)
show_seq <- AlignmentsForZMW(83952293,  gap_opening = 30)

show_seq <- AlignmentsForZMW(131072176, gap_opening = 30)
show_seq <- AlignmentsForZMW(66062756, gap_opening = 30)
show_seq <- AlignmentsForZMW(66062756, gap_opening = 100)




BrowseSeqs(show_seq)

