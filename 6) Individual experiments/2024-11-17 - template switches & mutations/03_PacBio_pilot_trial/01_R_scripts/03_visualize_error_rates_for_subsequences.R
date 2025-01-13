## 2024-11-18


# Load packages and source code -------------------------------------------

root_dir      <- "~/CRISPR_4sgRNA"
exper_dir     <- file.path(root_dir, "6) Individual experiments")
illumina_dir  <- file.path(exper_dir, "2022-04-21 - Illumina paired-end 2sg - first trial")
lumi_func_dir <- file.path(illumina_dir, "01_R_scripts", "R_functions")
project_dir   <- file.path(exper_dir, "2024-11-17 - template switches & mutations")

source(file.path(lumi_func_dir, "01_violin_swarm_plots.R")) # For RepositionByGroups
source(file.path(project_dir, "01_R_functions", "03_visualizing_error_rates_for_subsequences.R"))


# Define paths ------------------------------------------------------------

plasmids_dir <- file.path(project_dir, "02_PacBio_plasmids", "02_R_objects")
sub_dir      <- file.path(project_dir, "03_PacBio_pilot_trial")
rdata_dir    <- file.path(sub_dir, "02_R_objects")
output_dir   <- file.path(sub_dir, "03_output")


# Load data ---------------------------------------------------------------

load(file.path(plasmids_dir, "03_compute_error_rates__subsequences.RData"))
load(file.path(rdata_dir, "01_extract_and_categorize_subsequences__features_df.RData"))
load(file.path(rdata_dir, "02_compute_error_rates_for_subsequences.RData"))



# Add comparison data from the original plasmids --------------------------

sg_pairs <- names(subsequence_error_mat_list)[2:7]

full_reads_deletions_df_list <- lapply(full_reads_deletions_df_list, function(x) {
  x[, "Fraction_reference"] <- subsequence_error_mat_list[["Full"]][, "Fraction_deleted"]
  x
})

full_reads_errors_df_list <- lapply(full_reads_errors_df_list, function(x) {
  x[, "Fraction_reference"] <- subsequence_error_mat_list[["Full"]][, "Fraction_incorrect"]
  x
})

all_reads_deletions_df_list <- lapply(1:6, function(x) {
  use_df <- all_reads_deletions_df_list[[x]]
  use_df[, "Fraction_reference"] <- subsequence_error_mat_list[[sg_pairs[[x]]]][, "Fraction_deleted"]
  use_df
})
names(all_reads_deletions_df_list) <- sg_pairs

all_reads_errors_df_list <- lapply(1:6, function(x) {
  use_df <- all_reads_errors_df_list[[x]]
  use_df[, "Fraction_reference"] <- subsequence_error_mat_list[[sg_pairs[[x]]]][, "Fraction_incorrect"]
  use_df
})
names(all_reads_errors_df_list) <- sg_pairs




# Display error rates -----------------------------------------------------

PDF_width <- 7
PDF_height <- 8

for (create_PDF in c(FALSE, TRUE)) {

  if (create_PDF) {
    pdf(file.path(output_dir, "Dumbell plots - a) error rates - only fully mapped reads.pdf"),
        width = PDF_width, height = PDF_height
        )
  }
  ErrorDumbBells(full_reads_errors_df_list[["sg1_sg2"]], 1, 2)
  ErrorDumbBells(full_reads_errors_df_list[["sg2_sg3"]], 2, 3)
  ErrorDumbBells(full_reads_errors_df_list[["sg3_sg4"]], 3, 4)
  ErrorDumbBells(full_reads_errors_df_list[["sg1_sg3"]], 1, 3)
  ErrorDumbBells(full_reads_errors_df_list[["sg2_sg4"]], 2, 4)
  ErrorDumbBells(full_reads_errors_df_list[["sg1_sg4"]], 1, 4)
  dev.off()


  if (create_PDF) {
    pdf(file.path(output_dir, "Dumbell plots - b) deletion rates - only fully mapped reads.pdf"),
        width = PDF_width, height = PDF_height
        )
  }
  ErrorDumbBells(full_reads_deletions_df_list[["sg1_sg2"]], 1, 2, x_axis_label = "Deletion rate", x_upper_limit = 5)
  ErrorDumbBells(full_reads_deletions_df_list[["sg2_sg3"]], 2, 3, x_axis_label = "Deletion rate", x_upper_limit = 5)
  ErrorDumbBells(full_reads_deletions_df_list[["sg3_sg4"]], 3, 4, x_axis_label = "Deletion rate", x_upper_limit = 5)
  ErrorDumbBells(full_reads_deletions_df_list[["sg1_sg3"]], 1, 3, x_axis_label = "Deletion rate", x_upper_limit = 5)
  ErrorDumbBells(full_reads_deletions_df_list[["sg2_sg4"]], 2, 4, x_axis_label = "Deletion rate", x_upper_limit = 5)
  ErrorDumbBells(full_reads_deletions_df_list[["sg1_sg4"]], 1, 4, x_axis_label = "Deletion rate", x_upper_limit = 5)
  dev.off()


  if (create_PDF) {
    pdf(file.path(output_dir, "Dumbell plots - c) error rates - all reads.pdf"),
        width = PDF_width, height = PDF_height
        )
  }
  ErrorDumbBells(all_reads_errors_df_list[["sg1_sg2"]])
  ErrorDumbBells(all_reads_errors_df_list[["sg2_sg3"]])
  ErrorDumbBells(all_reads_errors_df_list[["sg3_sg4"]])
  ErrorDumbBells(all_reads_errors_df_list[["sg1_sg3"]])
  ErrorDumbBells(all_reads_errors_df_list[["sg2_sg4"]])
  ErrorDumbBells(all_reads_errors_df_list[["sg1_sg4"]])
  if (create_PDF) {
    dev.off()
  }


  if (create_PDF) {
    pdf(file.path(output_dir, "Dumbell plots - d) deletion rates - all reads.pdf"),
        width = PDF_width, height = PDF_height
        )
  }
  ErrorDumbBells(all_reads_deletions_df_list[["sg1_sg2"]], x_axis_label = "Deletion rate")
  ErrorDumbBells(all_reads_deletions_df_list[["sg2_sg3"]], x_axis_label = "Deletion rate")
  ErrorDumbBells(all_reads_deletions_df_list[["sg3_sg4"]], x_axis_label = "Deletion rate")
  ErrorDumbBells(all_reads_deletions_df_list[["sg1_sg3"]], x_axis_label = "Deletion rate")
  ErrorDumbBells(all_reads_deletions_df_list[["sg2_sg4"]], x_axis_label = "Deletion rate")
  ErrorDumbBells(all_reads_deletions_df_list[["sg1_sg4"]], x_axis_label = "Deletion rate")
  if (create_PDF) {
    dev.off()
  }

}



