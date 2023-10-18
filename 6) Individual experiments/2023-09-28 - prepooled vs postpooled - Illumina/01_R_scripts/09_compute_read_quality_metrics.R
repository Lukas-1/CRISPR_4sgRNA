## 2023-10-07


# Load packages and source code -------------------------------------------

CRISPR_root_directory        <- "~/CRISPR_4sgRNA"
experiments_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir     <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_illumina_functions_dir <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")
project_dir                  <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")
R_functions_dir              <- file.path(project_dir, "01_R_scripts", "R_functions")
source(file.path(first_illumina_functions_dir, "06_read_level_QC.R"))
source(file.path(R_functions_dir, "01_reading_in_data.R"))



# Define folder paths -----------------------------------------------------

input_dir   <- file.path(project_dir, "02_input")
rdata_dir   <- file.path(project_dir, "03_R_objects")



# Read in FASTQ files -----------------------------------------------------

all_files <- list.files(path.expand(input_dir), full.names = TRUE)
all_reads <- sapply(all_files, function(x) ShortRead::readFastq(x), simplify = FALSE)
all_reads <- lapply(all_reads, function(x) {
  reads_object <- subseq(ShortRead::sread(x), 1, 20)
  quality_object <- subseq(as(Biostrings::quality(x), "PhredQuality"), 1, 20)
  ShortReadQ(sread = reads_object , quality = quality_object)
})



# Construct samples_df ----------------------------------------------------

samples_df <- ProcessSamples(all_reads)
samples_df[, "Sample_name"] <- sub("-pool", "pool", samples_df[, "Sample_name"], fixed = TRUE)

sample_splits <- strsplit(samples_df[, "Sample_name"], "_", fixed = TRUE)
replicate_numbers <- as.integer(sub("^R", "", sapply(sample_splits, "[[", 3)))
timepoint_numbers <- match(sapply(sample_splits, "[[", 2),
                           c("T0", "T12")
                           )
read_numbers <- rep(1:2, times = nrow(samples_df) / 2)
conditions_vec <- c("Prepool" = "Pre-pooled", "Postpool" = "Post-pooled")[sapply(sample_splits, "[[", 1)]
conditions_numbers <- match(conditions_vec, c("Pre-pooled", "Post-pooled"))
samples_df <- data.frame(
  samples_df,
  "Condition" = c("Prepool" = "Pre-pooled", "Postpool" = "Post-pooled")[sapply(sample_splits, "[[", 1)],
  "Timepoint" = sapply(sample_splits, "[[", 2),
  "Replicate" = replicate_numbers,
  "Read"      = read_numbers,
  "Rank"      = order(order(conditions_numbers, timepoint_numbers, replicate_numbers, read_numbers)),
  row.names   = NULL,
  stringsAsFactors = FALSE
)
samples_df[, "Long_name"] <- paste0(sub("-", "", samples_df[, "Condition"], fixed = TRUE),
                                    "_", samples_df[, "Timepoint"],
                                    "_rep", samples_df[, "Replicate"],
                                    "_read", samples_df[, "Read"]
                                    )



# Compute per-base mean quality scores ------------------------------------

base_qual_mat <- MeanBaseQuality(all_reads, samples_df = samples_df)



# Compute density curves --------------------------------------------------

reads_options <- list("Read 1" = 1L, "Read 2" = 2L, "Both" = 1:2)
sequence_qual_densities <- lapply(reads_options, function(x) QualityDensities(all_reads, x, samples_df))
GC_content_densities <- lapply(reads_options, function(x) GC_Densities(all_reads, x, samples_df))



# Finish samples_df -------------------------------------------------------

samples_df <- unique(samples_df[order(samples_df[, "Rank"]), 1:6], MARGIN = 1)
row.names(samples_df) <- NULL



# Save data ---------------------------------------------------------------

save(list = c("samples_df", "base_qual_mat",
              "sequence_qual_densities", "GC_content_densities"
              ),
     file = file.path(rdata_dir, "09_compute_read_quality_metrics.RData")
     )





