## 2022-10-21


# Load packages and source code -------------------------------------------

CRISPR_root_directory        <- "~/CRISPR_4sgRNA"
experiments_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir     <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_illumina_functions_dir <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")
project_dir                  <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
R_functions_dir              <- file.path(project_dir, "01_R_scripts", "R_functions")
source(file.path(first_illumina_functions_dir, "06_read_level_QC.R"))
source(file.path(R_functions_dir, "01_reading_in_data.R"))




# Define folder paths -----------------------------------------------------

input_dir <- file.path(project_dir, "02_input_data")
rdata_dir <- file.path(project_dir, "03_R_objects")
reads_dir <- file.path(input_dir, "Raw reads")




# Read in FASTQ files -----------------------------------------------------

all_files <- list.files(path.expand(reads_dir), full.names = TRUE)
all_reads <- sapply(all_files, function(x) ShortRead::readFastq(x), simplify = FALSE)



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





