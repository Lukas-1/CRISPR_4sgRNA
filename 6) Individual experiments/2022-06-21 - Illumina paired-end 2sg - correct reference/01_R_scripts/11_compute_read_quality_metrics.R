## 2022-10-25


# Load packages and source code -------------------------------------------

library("ShortRead")
CRISPR_root_directory        <- "~/CRISPR_4sgRNA"
experiments_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir     <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_illumina_functions_dir <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")
source(file.path(first_illumina_functions_dir, "06_read_level_QC.R"))




# Define folder paths -----------------------------------------------------

first_rdata_dir <- file.path(first_illumina_trial_dir, "03_R_objects")
project_dir     <- file.path(experiments_directory, "2022-06-21 - Illumina paired-end 2sg - correct reference")
rdata_dir       <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "01_read_in_data_num_reads.RData"))
load(file.path(first_rdata_dir, "01_read_in_data_run1.RData"))
load(file.path(first_rdata_dir, "01_read_in_data_run2_chunk1.RData"))
load(file.path(first_rdata_dir, "01_read_in_data_run2_chunk2.RData"))



# Define functions --------------------------------------------------------

ProcessSamples <- function(all_reads_list) {
  results_df <- data.frame(
    "Sample_name"    = names(all_reads_list),
    "Sample_number"  = match(names(all_reads_list), unique(names(all_reads_list))),
    "Num_reads"      = vapply(all_reads_list, length, integer(1)),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  return(results_df)
}


SimulateShortReadQ <- function(sg_number) {
  use_columns <- c("Sample_number", paste0(c("Sequence_sg", "Quality_sg"), sg_number))
  combined_df <- rbind.data.frame(run1_fastq_df[, use_columns],
                                  run2_chunk1_fastq_df[, use_columns],
                                  run2_chunk2_fastq_df[, use_columns],
                                  stringsAsFactors = FALSE,
                                  make.row.names = FALSE
                                  )
  reads_list <- split(combined_df[, names(combined_df) != "Sample_number"],
                      combined_df[, "Sample_number"]
                      )
  reads_list <- lapply(reads_list, function(x) {
    ShortReadQ(sread   = DNAStringSet(x[, paste0("Sequence_sg", sg_number)]),
               quality = BStringSet(x[, paste0("Quality_sg", sg_number)])
               )
  })
  sample_names <- sub("_2sg", "", samples_df[, "Sample_name"], fixed = TRUE)
  sample_names <- sub("_off", "", sample_names, fixed = TRUE)
  names(reads_list) <- sample_names
  return(reads_list)
}



# Simulate FASTQ files ----------------------------------------------------

all_reads <- c(SimulateShortReadQ(1), SimulateShortReadQ(2))
all_reads <- all_reads[order(match(names(all_reads), names(all_reads)))]



# Count the number of reads per sample ------------------------------------

samples_df <- ExtendSamplesDf(ProcessSamples(all_reads))
samples_df <- unique(samples_df[order(samples_df[, "Rank"]), 1:6], MARGIN = 1)



# Compute per-base mean quality scores ------------------------------------

base_qual_mat <- MeanBaseQuality(all_reads)



# Compute density curves --------------------------------------------------

reads_options <- list("Read 1" = 1L, "Read 2" = 2L, "Both" = 1:2)
sequence_qual_densities <- lapply(reads_options, function(x) QualityDensities(all_reads, x))
GC_content_densities <- lapply(reads_options, function(x) GC_Densities(all_reads, x))



# Save data ---------------------------------------------------------------

save(list = c("samples_df", "base_qual_mat",
              "sequence_qual_densities", "GC_content_densities"
              ),
     file = file.path(rdata_dir, "11_compute_read_quality_metrics.RData")
     )


