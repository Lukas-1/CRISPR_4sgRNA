### 5th August 2021 ###



# Import packages and source code -----------------------------------------

library("devEMF")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))
source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions
source(file.path(R_functions_directory, "20) Summarizing data across wells.R"))
source(file.path(R_functions_directory, "23) Systematically comparing different read quality thresholds.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
file_output_directory    <- file.path(s2r2_directory, "5) Output", "Figures", "Explore read quality cut-offs")
PNGs_output_directory    <- file.path(s2r2_directory, "5) Output", "PNGs", "Explore read quality cut-offs")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - plates_analysis_list.RData"))
load(file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - contam_df.RData"))
load(file.path(s2r2_R_objects_directory, "10) Identify and characterize deletions.RData"))
load(file.path(s2r2_R_objects_directory, "23) Summarize data across wells - plate selections.RData"))





# Create the 384-well-plate "distance list" -------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)





# Summarize wells for different read quality cut-offs ---------------------

summary_list_list <- SummarizeWellsForAllCutoffs(plates_analysis_list,
                                                 filter_cross_plate = TRUE
                                                 )




# Modify plate selections -------------------------------------------------

plate_selection_titles_list <- c(
  list("Colony-picked (second run)" = "Single-colony picked controls"),
  plate_selection_titles_list
)
plate_selection_prefixes <- c("00", plate_selection_prefixes)




# Produce example plots ---------------------------------------------------

ViolinBoxAllCutoffs(summary_list_list, "Count_all_4", "All plates (second run)")
ViolinBoxAllCutoffs(summary_list_list, "Count_no_contam_all_4", "All plates (second run)")

ViolinBoxAllCutoffs(summary_list_list, "Count_all_4", "All plates (second run)",
                    filter_mean_quality = TRUE
                    )

ViolinBoxAllCutoffs(summary_list_list, "Num_contaminated_reads", "All plates (second run)")

ViolinBoxAllCutoffs(summary_list_list, "Num_contaminated_reads",
                    "Colony-picked (second run)", use_y_limits = c(0, 0.1)
                    )

ViolinBoxAllCutoffs(summary_list_list, "Num_reads_with_deletions_exceeding_20bp", "All plates (second run)")
ViolinBoxAllCutoffs(summary_list_list, "Num_reads_with_deletions_spanning_tracrRNAs", "All plates (second run)")
ViolinBoxAllCutoffs(summary_list_list, "Count_no_contam_sg4_cr4", "All plates (second run)")



# Export .emf files for selected plots ------------------------------------

use_lwd <- 0.8
use_cex <- 0.8
thesis_dir <- file.path(s2r2_directory, "5) Output", "Figures", "Thesis")


devEMF::emf(file = file.path(thesis_dir, paste0("Quality cutoffs - 1) Single colony-derived controls.emf")),
            width = 7.1, height = 3.8, emfPlus = FALSE, coordDPI = 1500
            )
par(mai = c(0.8, 0.8, 0.6, 0.3), lwd = use_lwd, cex = use_cex)
ViolinBoxAllCutoffs(summary_list_list,
                    "Count_all_4",
                    "Colony-picked (second run)",
                    custom_title     = "Single colony-derived controls",
                    y_axis_label     = "All 4 sgRNAs and tracrRNAs correct",
                    title_line       = 1.4,
                    large_gap_lines  = 2.5,
                    bold_read_counts = TRUE,
                    set_mar          = FALSE,
                    embed_PNG        = TRUE,
                    transparent_PNG  = FALSE,
                    brewer_palette   = "Purples",
                    point_cex        = 0.6,
                    use_wex          = 0.7,
                    violin_wex       = 0.5,
                    sina_plot        = TRUE,
                    box_lwd          = 1.25,
                    med_lwd          = 3,
                    grid_lwd         = 0.7
                    )
dev.off()



devEMF::emf(file = file.path(thesis_dir, paste0("Quality cutoffs - 2) APPEAL plasmids.emf")),
            width = 7.1, height = 3.8, emfPlus = FALSE, coordDPI = 1500
            )
par(mai = c(0.8, 0.8, 0.6, 0.3), lwd = use_lwd, cex = use_cex)
ViolinBoxAllCutoffs(summary_list_list,
                    "Count_all_4",
                    "All plates (second run)",
                    custom_title     = "Plasmids cloned by APPEAL",
                    y_axis_label     = "All 4 sgRNAs and tracrRNAs correct",
                    title_line       = 1.4,
                    large_gap_lines  = 2.5,
                    bold_read_counts = TRUE,
                    set_mar          = FALSE,
                    embed_PNG        = TRUE,
                    transparent_PNG  = FALSE,
                    point_cex        = 0.4,
                    use_wex          = 0.7,
                    violin_wex       = 0.7,
                    sina_plot        = TRUE,
                    box_lwd          = 1.25,
                    med_lwd          = 4,
                    grid_lwd         = 0.7
                    )
dev.off()





# Export violin/box plots -------------------------------------------------

DrawAllCutoffComparisons(export_PNGs = FALSE)






