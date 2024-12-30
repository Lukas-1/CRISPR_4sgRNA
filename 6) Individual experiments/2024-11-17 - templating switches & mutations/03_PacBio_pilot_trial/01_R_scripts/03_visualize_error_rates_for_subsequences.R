## 2024-11-18


# Load packages and source code -------------------------------------------

root_dir      <- "~/CRISPR_4sgRNA"
exper_dir     <- file.path(root_dir, "6) Individual experiments")
illumina_dir  <- file.path(exper_dir, "2022-04-21 - Illumina paired-end 2sg - first trial")
lumi_func_dir <- file.path(illumina_dir, "01_R_scripts", "R_functions")
project_dir   <- file.path(exper_dir, "2024-11-17 - templating switches & mutations")

source(file.path(lumi_func_dir, "01_violin_swarm_plots.R")) # For RepositionByGroups
source(file.path(project_dir, "01_R_functions", "03_visualizing_error_rates_for_subsequences.R"))


# Define paths ------------------------------------------------------------

sub_dir    <- file.path(project_dir, "03_PacBio_pilot_trial")
rdata_dir  <- file.path(sub_dir, "02_R_objects")
output_dir <- file.path(sub_dir, "03_output")


# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "01_extract_and_categorize_sequences__features_df.RData"))
load(file.path(rdata_dir, "02_compute_error_rates.RData"))



# Display error rates -----------------------------------------------------

PDF_width <- 6.5
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
    pdf(file.path(output_dir, "Dumbell plots - b) error rates - all reads.pdf"),
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
    pdf(file.path(output_dir, "Dumbell plots - c) deletion rates - all reads.pdf"),
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



