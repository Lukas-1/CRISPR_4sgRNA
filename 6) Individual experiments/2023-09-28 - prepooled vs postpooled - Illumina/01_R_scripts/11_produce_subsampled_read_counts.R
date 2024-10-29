## 2023-12-13


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
source(file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions", "04_obtaining_counts_for_2sg_plasmids.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_parts_combined.RData"))



# Define functions --------------------------------------------------------

AllSamplesCountsWithSubsampling <- function(input_df,
                                            use_fractions   = c(1, 0.75, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025, 0.01, 0.005, 0.001),
                                            num_repetitions = 3
                                            ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  set.seed(1)

  results_list <- lapply(use_fractions, function(use_fraction) {

    message(paste0("\nUsing a subsampling ratio of ", use_fraction, "..."))
    num_reads <- round(nrow(input_df) * use_fraction)
    if (use_fraction == 1) {
      reps_vec <- 1
    } else {
      reps_vec <- seq_len(num_repetitions)
    }
    rep_list <- lapply(reps_vec, function(i) {

      sample_indices <- sample(seq_len(nrow(input_df)), size = num_reads)
      subsampled_reads_df <- input_df[sample_indices, ]
      row.names(subsampled_reads_df) <- NULL

      message(paste0("\nProcessing repetition #", i, "..."))

      counts_mat <- AllSamplesCounts(subsampled_reads_df,
                                     only_0MM = FALSE,
                                     no_template_switch = TRUE
                                     )

      return(counts_mat)
    })
    names(rep_list) <- paste0("rep", reps_vec)
    return(rep_list)
  })
  names(results_list) <- paste0(use_fractions * 100, "% sampled")
  return(results_list)
}




# Modify lumi_df ----------------------------------------------------------

names(lumi_df) <- sub("_sg2", "_sg1", names(lumi_df), fixed = TRUE)
names(lumi_df) <- sub("_sg3", "_sg2", names(lumi_df), fixed = TRUE)

lumi_df <- lumi_df[, c("Sample_name", "Sample_number",
                       "Plasmid_sg1", "Plasmid_sg2",
                       "Has_template_switch", "Num_MM_sg1", "Num_MM_sg2",
                       "Num_matched_sgRNAs"
                       )]
gc()



# Obtain read counts from subsamples of reads -----------------------------

subsampled_counts_mat_list <- AllSamplesCountsWithSubsampling(lumi_df)



# Adjust column names and fix the column order ----------------------------

samples_in_order <- c(
  "Prepool_T0_R1", "Prepool_T0_R2", "Prepool_T12_R1", "Prepool_T12_R2",
  "Postpool_T0_R1", "Postpool_T0_R2", "Postpool_T12_R1", "Postpool_T12_R2"
)

subsampled_counts_mat_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    colnames(y) <- sub("-", "", colnames(y), fixed = TRUE)
    colnames(y) <- sub("_S[0-8]$", "", colnames(y))
    y <- y[, samples_in_order]
    return(y)
  })
})



# Save data ---------------------------------------------------------------

save(list = "subsampled_counts_mat_list",
     file = file.path(rdata_dir, "11_produce_subsampled_read_counts.RData")
     )

