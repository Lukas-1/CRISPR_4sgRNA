### 8th September 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "15) Examining the Hamming distances of barcodes.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
R_objects_directory    <- file.path(file_directory, "3) R objects")

file_output_directory  <- file.path(file_directory, "5) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))




# Calculate distances -----------------------------------------------------

row_distances_df <- GetDistances(row_barcodes)
column_distances_df <- GetDistances(column_barcodes)




# Prepare plots -----------------------------------------------------------

row_title <- " pairs of row barcodes "
column_title <- " pairs of column barcodes "



# Display results ---------------------------------------------------------

PlotDistances(row_distances_df, title_unit = row_title)
PlotDistances(column_distances_df, title_unit = column_title)




# Export results ----------------------------------------------------------

pdf(file = file.path(plots_output_directory,
                     "Check that barcode distances are random.pdf"
                     ),
    width = 7.5, height = 6.5
    )
par(oma = c(0.5, 1, 0.25, 1))
PlotDistances(row_distances_df, title_unit = row_title)
PlotDistances(column_distances_df, title_unit = column_title)
dev.off()





