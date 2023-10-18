## 2023-09-28


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))



# Define functions --------------------------------------------------------

ProcessSamples <- function(all_reads_list) {
  path_list <- strsplit(names(all_reads_list), "/", fixed = TRUE)
  file_names <- vapply(path_list, function(x) x[[length(x)]], "")
  file_splits <- strsplit(file_names, "_", fixed = TRUE)
  file_names <- vapply(file_splits, function(x) paste0(x[2:length(x)], collapse = "_"), "")
  file_names <- sub("_001.fastq.gz", "", file_names, fixed = TRUE)
  sample_names <- sub("_R[12]$", "", file_names)
  sample_name_splits <- strsplit(sample_names, "(?<=[0-9])-", perl = TRUE)
  results_df <- data.frame(
    "Sample_name"    = sapply(sample_name_splits, "[[", 2),
    "Sample_number"  = as.integer(sapply(sample_name_splits, "[[", 1)),
    "Num_reads"      = vapply(all_reads_list, length, integer(1)),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  return(results_df)
}



MakeFastqDf <- function(all_reads_list) {

  ## Prepare sample names
  names_df <- ProcessSamples(all_reads_list)
  indices_list <- split(seq_along(all_reads_list),
                        factor(names_df[, "Num_reads"],
                               levels = unique(names_df[, "Num_reads"])
                               )
                        )

  ## Assemble the data frame
  df_list <- lapply(indices_list, function(x) {
    data.frame(
      "Sample_number"  = names_df[, "Sample_number"][[x[[1]]]],
      "Sample_name"    = names_df[, "Sample_name"][[x[[1]]]],
      "Sequence_sg2"   = substr(as.character(ShortRead::sread(all_reads_list[[x[[1]]]])), 1, 20),
      "Sequence_sg3"   = substr(as.character(ShortRead::sread(all_reads_list[[x[[2]]]])), 1, 20),
      "Quality_sg2"    = substr(as.character(as(Biostrings::quality(all_reads_list[[x[[1]]]]), "PhredQuality")), 1, 20),
      "Quality_sg3"    = substr(as.character(as(Biostrings::quality(all_reads_list[[x[[2]]]]), "PhredQuality")), 1, 20),
      stringsAsFactors = FALSE
    )
  })
  fastq_df <- do.call(rbind.data.frame, c(df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  fastq_df[, "Mean_quality_sg2"] <- GetMeanQuality(fastq_df[, "Quality_sg2"])
  fastq_df[, "Mean_quality_sg3"] <- GetMeanQuality(fastq_df[, "Quality_sg3"])
  return(fastq_df)
}





