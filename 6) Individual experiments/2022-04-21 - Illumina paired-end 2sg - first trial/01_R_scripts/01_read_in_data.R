## 2022-04-21


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir    <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
project_dir           <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir       <- file.path(project_dir, "01_R_scripts", "R_functions")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))
source(file.path(R_functions_dir, "06_read_level_QC.R"))



# Define folder paths -----------------------------------------------------

input_dir   <- file.path(project_dir, "02_input_data")
rdata_dir   <- file.path(project_dir, "03_R_objects")
reads_dir   <- file.path(input_dir, "Raw reads")



# Read in FASTQ files -----------------------------------------------------

all_reads_list_run1 <- sapply(list.files(path.expand(file.path(reads_dir, "20220408")), full.names = TRUE),
                              function(x) ShortRead::readFastq(x),
                              simplify = FALSE
                              )
all_reads_list_run2 <- sapply(list.files(path.expand(file.path(reads_dir, "20220412")), full.names = TRUE),
                              function(x) ShortRead::readFastq(x),
                              simplify = FALSE
                              )


# Define functions --------------------------------------------------------

ProcessSamples <- function(all_reads_list) {
  path_list <- strsplit(names(all_reads_list), "/", fixed = TRUE)
  file_names <- vapply(path_list, function(x) x[[length(x)]], "")
  file_splits <- strsplit(file_names, "_", fixed = TRUE)
  file_names <- vapply(file_splits, function(x) paste0(x[3:length(x)], collapse = "_"), "")
  file_names <- sub(".fastq.gz", "", file_names, fixed = TRUE)
  sample_names <- sub("_R[12]$", "", file_names)
  sample_name_splits <- strsplit(sample_names, "-", fixed = TRUE)
  sample_numbers <- as.integer(sapply(sample_name_splits, "[[", 1))
  sample_names <- sapply(sample_name_splits, "[[", 2)
  results_df <- data.frame(
    "Sample_name"    = sample_names,
    "Sample_number"  = sample_numbers,
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
      "Sequence_sg1"   = as.character(ShortRead::sread(all_reads_list[[x[[1]]]])),
      "Sequence_sg2"   = as.character(ShortRead::sread(reverseComplement(all_reads_list[[x[[2]]]]))),
      "Quality_sg1"    = as.character(as(Biostrings::quality(all_reads_list[[x[[1]]]]), "PhredQuality")),
      "Quality_sg2"    = as.character(as(Biostrings::quality(reverseComplement(all_reads_list[[x[[2]]]])), "PhredQuality")),
      stringsAsFactors = FALSE
    )
  })
  fastq_df <- do.call(rbind.data.frame, c(df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  fastq_df[, "Mean_quality_sg1"] <- GetMeanQuality(fastq_df[, "Quality_sg1"])
  fastq_df[, "Mean_quality_sg2"] <- GetMeanQuality(fastq_df[, "Quality_sg2"])

  fastq_df[, "Sequence_sg1"] <- substr(fastq_df[, "Sequence_sg1"], 1, 19)
  fastq_df[, "Sequence_sg2"] <- substr(fastq_df[, "Sequence_sg2"],
                                       nchar(fastq_df[, "Sequence_sg2"]) - 18,
                                       nchar(fastq_df[, "Sequence_sg2"])
                                       )
  fastq_df[, "Quality_sg1"] <- substr(fastq_df[, "Quality_sg1"], 1, 19)
  fastq_df[, "Quality_sg2"] <- substr(fastq_df[, "Quality_sg2"],
                                      nchar(fastq_df[, "Quality_sg2"]) - 18,
                                      nchar(fastq_df[, "Quality_sg2"])
                                      )
  return(fastq_df)
}



# Count the number of reads per sample ------------------------------------

run1_samples_df <- unique(ProcessSamples(all_reads_list_run1), MARGIN = 1)
run2_samples_df <- unique(ProcessSamples(all_reads_list_run2), MARGIN = 1)

samples_df <- data.frame(
  run1_samples_df[, 1:2],
  "Num_reads_run1" = run1_samples_df[, "Num_reads"],
  "Num_reads_run2" = run2_samples_df[, "Num_reads"],
  stringsAsFactors = FALSE
)
samples_df[, "Num_reads"] <- rowSums(samples_df[, c("Num_reads_run1", "Num_reads_run2")])




# Assemble fastq data frames ----------------------------------------------

run1_fastq_df <- MakeFastqDf(all_reads_list_run1)
run2_fastq_df <- MakeFastqDf(all_reads_list_run2)




# Split fastq data frames -------------------------------------------------

run2_chunk1_fastq_df <- run2_fastq_df[1:45000000, ]
run2_chunk2_fastq_df <- run2_fastq_df[45000001:nrow(run2_fastq_df), ]
row.names(run2_chunk2_fastq_df) <- NULL




# Save data ---------------------------------------------------------------

save(samples_df, file = file.path(rdata_dir, "01_read_in_data_num_reads.RData"))
save(run1_fastq_df, file = file.path(rdata_dir, "01_read_in_data_run1.RData"))
save(run2_chunk1_fastq_df, file = file.path(rdata_dir, "01_read_in_data_run2_chunk1.RData"))
save(run2_chunk2_fastq_df, file = file.path(rdata_dir, "01_read_in_data_run2_chunk2.RData"))






