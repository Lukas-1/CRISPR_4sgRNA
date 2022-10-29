### 29th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
file_output_directory    <- file.path(s2r2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Schematics of a 384-well plate")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Schematics")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# General preparation -----------------------------------------------------

sg_sequences_df[["Empty_well"]] <- FALSE



# Export all plates -------------------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]#[order(plates_df[["Plate_rank"]])]

DrawSchematicsForAllPlates()




# Export .emf graphics from a selected plate ------------------------------

library("devEMF")
thesis_dir <- file.path(file_output_directory, "Figures", "Thesis")

summary_df <- ccs7_df_list[["filtered_summary_df"]]
summary_df <- summary_df[summary_df[, "Plate_number"] %in% 33, ]
row.names(summary_df) <- NULL


par("cex" = 2)
emf(file.path(thesis_dir, "Plate HA_21 - schematic - Num_contaminated_reads.emf"),
    width = 7.5, height = 6, emfPlus = FALSE
    )
WellLayoutBarplot(summary_df,
                  "Num_contaminated_reads",
                  sg_sequences_df,
                  label_rows_with_letters = TRUE
                  )
CenterText("Percentage of reads with contaminations",
           y = par("usr")[[4]] + diff(grconvertY(c(0, 1.3), from = "lines", to = "user")),
           use_cex = par("cex") * 1.3
           )
dev.off()



