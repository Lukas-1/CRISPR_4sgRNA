### 3rd August 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory        <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions
source(file.path(R_functions_directory, "20) Summarizing data across wells.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
file_output_directory    <- file.path(s2r2_directory, "5) Output", "Figures", "Summaries across wells")
PNGs_output_directory    <- file.path(s2r2_directory, "5) Output", "PNGs", "Summaries across wells")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Define plate selections -------------------------------------------------

are_CRISPRa  <- grepl("HA", plates_df[["Plate_name"]], fixed = TRUE)
are_CRISPRko <- grepl("HO", plates_df[["Plate_name"]], fixed = TRUE)

are_library <- are_CRISPRa | are_CRISPRko

plate_selection_list <- list(
  "All plates (second run)"    = plates_df[["Plate_name"]][are_library],
  "Colony-picked (second run)" = plates_df[["Plate_name"]][plates_df[["Colony_picked"]]],
  "CRISPRa (second run)"       = plates_df[["Plate_name"]][are_CRISPRa],
  "CRISPRko (second run)"      = plates_df[["Plate_name"]][are_CRISPRko]
)

plate_selection_titles_list <- list(
  "All plates (second run)"    = "All plates (2nd PacBio run)",
  "CRISPRa (second run)"       = "CRISPRa plates (2nd PacBio run)",
  "CRISPRko (second run)"      = "CRISPRko plates (2nd PacBio run)"
)




# Create labels for individual plates -------------------------------------

plate_labels <- paste0("Plate #", plates_df[["Plate_number"]], " \u2013 ", plates_df[["Plate_name"]])
names(plate_labels) <- plates_df[["Plate_name"]]
plate_labels <- plate_labels[!(plates_df[["Colony_picked"]])]
are_single_plates <- c(rep(FALSE, length(plate_selection_titles_list)),
                       rep(TRUE, length(plate_labels))
                       )
plate_selection_titles_list <- c(plate_selection_titles_list, as.list(plate_labels))





# Prepare file name strings for plate selections --------------------------

plate_matches <- match(names(plate_selection_titles_list)[are_single_plates],
                       plates_df[["Plate_name"]]
                       )
plate_numbers <- plates_df[["Plate_number"]][plate_matches]
plate_selection_prefixes <- c(paste0("0", letters[seq_len(sum(!(are_single_plates)))]),
                              formatC(plate_numbers, width = 2, flag = "0")
                              )



# Produce example plots ---------------------------------------------------

use_df <- ccs7_df_list[["filtered_summary_df"]]

SummaryBoxPlot(use_df, "All plates (second run)")

LollipopPlot(use_df,
             "All plates (second run)",
             c("Num_contaminated_reads",
               "Num_contaminated_reads_aligned",
               "Num_cross_plate_contaminated"
               ),
             use_y_limits = c(0, 0.05)
             )

LollipopPlot(use_df,
             "All plates (second run)",
             c("Num_contaminated_reads"),
             use_y_limits = c(0, 0.05)
             )

LollipopPlot(use_df,
             "All plates (second run)",
             c("Num_reads_with_sgRNA_deletion",
               "Num_reads_with_deletions_exceeding_20bp",
               "Num_reads_with_deletions_spanning_tracrRNAs",
               "Num_reads_with_deletions_spanning_promoters"
               )
             )




# Export lollipop plots and violin/box plots ------------------------------

DrawAllLollipopsAndViolins(export_PNGs = FALSE)





# Save data ---------------------------------------------------------------

save(list = c("plate_selection_list",
              "plate_selection_titles_list",
              "plate_selection_prefixes"
              ),
     file = file.path(s2r2_R_objects_directory, "23) Summarize data across wells - plate selections.RData")
     )



