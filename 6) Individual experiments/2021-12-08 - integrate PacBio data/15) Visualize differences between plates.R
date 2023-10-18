### 31st December 2021 ###



# Import packages and source code -----------------------------------------

library("svglite")
CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "01) Define titles and labels.R"))
source(file.path(p1_R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust
source(file.path(s2r1_R_functions_directory, "04) Visualizing differences between plates.R"))
source(file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial",
                 "01_R_scripts", "R_functions", "01_violin_swarm_plots.R" # For StartEmbedPNG and StopEmbedPNG
                 ))



# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")
file_output_directory    <- file.path(s2rI_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Compare plates")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Compare plates")
thesis_plots_directory   <- file.path(file_output_directory, "Figures", "Manuscript", "Thesis")



# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Make preparations -------------------------------------------------------

titles_list <- c(list("Count_total" = "Number of reads per well"),
                 titles_list
                 )



# Draw example plots ------------------------------------------------------

ComparePlates(ccs7_df_list[["filtered_summary_df"]], "Count_total",
              use_cex = 0.075, side_space = -2
              )




# Draw ridgeline plots ----------------------------------------------------

use_summary_df <- ccs3_df_list[["original_summary_df"]]

CRISPRa_plate_numbers  <- plates_df[grepl("^HA", plates_df[, "Plate_name"]), "Plate_number"]
CRISPRko_plate_numbers <- plates_df[grepl("^HO", plates_df[, "Plate_name"]), "Plate_number"]

CRISPRa_summary_df  <- use_summary_df[use_summary_df[, "Plate_number"] %in% CRISPRa_plate_numbers, ]
CRISPRko_summary_df <- use_summary_df[use_summary_df[, "Plate_number"] %in% CRISPRko_plate_numbers, ]



PlateRidgelines(CRISPRa_summary_df,
                fill_colors   = c("#a83e76", "#f49aca"),
                border_colors = c("#621e42", "#892a5c"),
                show_title    = "T.gonfio library",
                use_alpha     = 0.75,
                zebra_pattern = TRUE,
                zebra_palify  = 0.3
                )

PlateRidgelines(CRISPRko_summary_df,
                fill_colors   = c("#225faa", "#9ac2f4"),
                border_colors = c("#113055", "#1d5395"),
                num_slots     = length(CRISPRa_plate_numbers),
                show_title    = "T.spiezzo library",
                use_alpha     = 0.75,
                zebra_pattern = TRUE
                )


PlateRidgelines(CRISPRa_summary_df,
                fill_colors    = c("#dd3c92", "#f49aca"),
                border_colors  = c("#f7d4e7", "#f8edf2"),
                show_title     = "T.gonfio library",
                use_alpha      = 0.75,
                zebra_pattern  = TRUE,
                zebra_palify   = 0.3,
                reorder_plates = TRUE
                )

PlateRidgelines(CRISPRko_summary_df,
                fill_colors   = c("#225faa", "#8cbaf2"),
                border_colors = c("#d4e4f7", "#d1e3fa"),
                num_slots     = length(CRISPRa_plate_numbers),
                show_title    = "T.spiezzo library",
                use_alpha     = 0.75,
                zebra_pattern = TRUE,
                reorder_plates = TRUE
                )



svg_width <- 2.7
svg_height <- 5
use_mai <- c(0.8, 0.3, 0.3, 0.4)

svg_path <- file.path(thesis_plots_directory, "Ridgeline plot - T.gonfio.svg")
svglite::svglite(svg_path, width = svg_width, height = svg_height, bg = "transparent")
par(mai = use_mai, cex = 0.6, lwd = 0.6)
PlateRidgelines(CRISPRa_summary_df,
                fill_colors      = c("#d74291", "#f49aca"),
                border_colors    = c("#f7d4e7", "#f8edf2"),
                show_title       = "T.gonfio library",
                use_alpha        = 0.75,
                zebra_pattern    = TRUE,
                zebra_palify     = 0.3,
                embed_PNG        = FALSE,
                PNG_lwd          = 0.8,
                ridge_height     = 2.5,
                plate_label_cex  = 0.5,
                grid_lwd         = 0.7,
                for_powerpoint   = TRUE
                )
dev.off()


svglite::svglite(file.path(thesis_plots_directory, "Ridgeline plot - T.spiezzo.svg"),
    width = svg_width, height = svg_height, bg = "transparent"
    )
par(mai = use_mai, cex = 0.6, lwd = 0.6)
PlateRidgelines(CRISPRko_summary_df,
                fill_colors     = c("#225faa", "#75acf0"),
                border_colors   = c("#d4e4f7", "#d1e3fa"),
                num_slots       = length(CRISPRa_plate_numbers),
                show_title      = "T.spiezzo library",
                use_alpha       = 0.75,
                zebra_pattern   = TRUE,
                embed_PNG       = FALSE,
                PNG_lwd         = 0.8,
                label_on_right  = TRUE,
                ridge_height    = 2.5,
                plate_label_cex = 0.5,
                grid_lwd        = 0.7
                )
dev.off()




# Export beeswarm plots ---------------------------------------------------

DrawAllPlateComparisons(use_cex          = 0.175,
                        beeswarm_spacing = 0.3,
                        beeswarm_corral  = "omit",
                        side_space       = -3.5,
                        use_width        = 22,
                        use_height       = 6.5,
                        export_PNGs      = TRUE,
                        exclude_controls = TRUE
                        )



