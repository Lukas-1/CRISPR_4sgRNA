## 2022-02-03


# Load packages and source code -------------------------------------------

library("RColorBrewer")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

rdata_dir   <- file.path(project_dir, "03_R_objects")
figures_dir <- file.path(project_dir, "04_output_data", "Figures")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_combine_aligned_reads.RData"))



# Prepare data ------------------------------------------------------------

alignments_df[, "Read_length"] <- nchar(alignments_df[, "Read_sequence"])



# Draw a histogram of read lengths ----------------------------------------

ReadLengthsHistogram(alignments_df[, "Read_length"])


pdf(file = file.path(figures_dir, "PDFs", "Read length histogram.pdf"),
    width = 5, height = 4
    )
ReadLengthsHistogram(alignments_df[, "Read_length"])
dev.off()


png(filename = file.path(figures_dir, "PNGs", "Read length histogram.png"),
    width = 5, height = 4, units = "in", res = 600
    )
ReadLengthsHistogram(alignments_df[, "Read_length"])
dev.off()





