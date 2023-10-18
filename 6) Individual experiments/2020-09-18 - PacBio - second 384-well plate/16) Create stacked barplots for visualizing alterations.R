### 21st October 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))
source(file.path(R_functions_directory, "11) Creating stacked barplots for visualizing alterations.R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")
file_output_directory  <- file.path(plate2_directory, "3) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")
manuscript_directory   <- file.path(plots_output_directory, "For the manuscript")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
load(file.path(p2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))
load(file.path(p2_R_objects_directory, "10) Process demultiplexed reads - with subsampling.RData"))





# Exclude problematic wells, for the time being ---------------------------

sg_sequences_df[["Empty_well"]] <- ifelse(sg_sequences_df[["Well_number"]] == 2,
                                          TRUE, FALSE
                                          )





# Export individual graphics ----------------------------------------------

ExportAlterationsForManuscript(sl7_ccs7_df_list[["filtered_summary_df"]],
                               "Colony-picked controls"
                               )


pdf(file = file.path(manuscript_directory, "Sand chart - 54 colony-picked controls - second 384-well plate - CCS7 (filtered).pdf"),
    width = 3.8, height = 6.7
    )
DrawReorderedSandPlots(sl7_ccs7_df_list[["filtered_summary_df"]],
                       main_title = "Genes re-ordered by the percentage of correct sgRNAs",
                       exclude_blocks = 2
                       )
dev.off()

pdf(file = file.path(manuscript_directory, "Sand chart - 51 genes - second 384-well plate - CCS7 (filtered).pdf"),
    width = 3.8, height = 6.7
    )
DrawReorderedSandPlots(sl7_ccs7_df_list[["filtered_summary_df"]],
                       main_title = "Genes re-ordered by the percentage of correct sgRNAs",
                       exclude_blocks = 1
                       )
dev.off()




# Set up loop -------------------------------------------------------------

AllSmrtLinkVersionsForOnePlate(DrawAlterationBarplot,
                               "Stacked barplots",
                               "Alteration barplots"
                               )

DrawAllSubsampledPlotsForOnePlate(DrawAlterationBarplot,
                                  "Stacked barplots",
                                  "Alteration barplots"
                                  )




