### 2023-10-18


# Load packages and source code -------------------------------------------

library("Rsamtools")
library("GenomicAlignments")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_pacbio_dir      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate", "1) R functions")
source(file.path(first_pacbio_dir, "02) Analyzing reads.R")) # For GetMeanQuality
source(file.path(first_pacbio_dir, "07) Categorizing subsequences of reads aligned to the reference.R"))



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2023-10-05 - prepooled vs postpooled - Nanopore")
rdata_dir   <- file.path(project_dir, "03_R_objects")
import_dir  <- file.path(project_dir, "04_intermediate_files", "2_input_to_R")



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

  sequences_vec <- as.character(subseq(mcols(bam_data)[, "seq"][matches_vec],
                                       start(mapped_ranges), end(mapped_ranges)
                                       ))

  phred_qualities <- subseq(mcols(bam_data)[, "qual"][matches_vec],
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



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_demuxed_meta_df.RData"))



# Read in data ------------------------------------------------------------

bam_file <- file.path(import_dir, "all_reads_output.bam")

param1 <- ScanBamParam(flag = scanBamFlag(isSupplementaryAlignment = FALSE,
                                          isSecondaryAlignment = FALSE,
                                          isMinusStrand = FALSE
                                          ),
                       what = c("qname", "flag", "seq", "qual")
                       )
bam_data <- readGAlignments(bam_file, param = param1, use.names = FALSE)




# Explore metadata --------------------------------------------------------

table(strand(bam_data))
table(mcols(bam_data)[, "flag"])
seqlevels(bam_data)



# Extract sgRNAs ----------------------------------------------------------

sg_df_list <- lapply(1:4, function(x) {
  use_sg <- paste0("sg", x)
  message("Extracting ", use_sg, "...")
  ExtractCoords(bam_data, features_list[[use_sg]] - 10, columns_postfix = paste0("_", use_sg))
})

all_reads <- unique(unlist(lapply(sg_df_list, function(x) x[, 1])))
all_read_numbers <- as.integer(substr(all_reads, 6, nchar(all_reads)))
all_read_numbers <- sort(all_read_numbers)

sg_df_list <- lapply(sg_df_list, function(x) {
  reads_vec <- as.integer(substr(x[, 1], 6, nchar(x[, 1])))
  matches_vec <- match(all_read_numbers, reads_vec)
  x[matches_vec, 2:ncol(x)]
})



# Combine data ------------------------------------------------------------

num_demuxed_reads <- nrow(demuxed_meta_df)
matches_vec <- match(all_read_numbers, seq_len(num_demuxed_reads))

extracted_df <- data.frame("Read_number" = all_read_numbers,
                           demuxed_meta_df[matches_vec, ],
                           lapply(sg_df_list, "[", 3),
                           lapply(sg_df_list, "[", 1),
                           lapply(sg_df_list, "[", 2),
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )



# Save data ---------------------------------------------------------------

save(list = c("extracted_df", "num_demuxed_reads"),
     file = file.path(rdata_dir, "03_extract_aligned_sgRNAs_from_SAM_file.RData")
     )

