### 4th April 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
plate1_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "15) Examining the Hamming distances of barcodes.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))




# Calculate distances -----------------------------------------------------

distances_df <- GetDistances(plates_df[["Barcode_sequence"]])



# Display results ---------------------------------------------------------

title_label <- " pairs of plate barcodes "
PlotDistances(distances_df, title_label)




# Export results ----------------------------------------------------------

pdf(file = file.path(plots_output_directory,
                     "Check that barcode distances are random.pdf"
                     ),
    width = 7.5, height = 6.5
    )
par(oma = c(0.5, 1, 0.25, 1))
PlotDistances(distances_df, title_unit = title_label)
dev.off()





