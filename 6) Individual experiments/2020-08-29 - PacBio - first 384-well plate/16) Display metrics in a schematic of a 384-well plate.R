### 19th September 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

file_input_directory   <- file.path(file_directory, "2) Input")
file_output_directory  <- file.path(file_directory, "5) Output")
R_objects_directory    <- file.path(file_directory, "3) R objects")
plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "11) Process demultiplexed PacBio reads.RData"))





# Display metrics in the layout of a 384-well plate -----------------------

DrawAllSchematicsForOnePlate()

BarPlotPanel(sl7_ccs3_df_list[["original_summary_df"]],
             "Longest_subsequence",
             sg_sequences_df,
             indicate_homologies = TRUE
             )







# Export selected graphics as .emf files ----------------------------------

library("devEMF")
manuscript_directory <- file.path(plots_output_directory, "For the manuscript")

summary_df <- sl7_ccs3_df_list[["filtered_summary_df"]]


emf(file.path(manuscript_directory, "Thesis", "Count_at_least_2.emf"),
    width = 7.5, height = 6, emfPlus = FALSE, coordDPI = 1500
    )
WellLayoutBarplot(summary_df,
                  "Count_at_least_2",
                  sg_sequences_df,
                  label_rows_with_letters = TRUE, label_cex = 0.8,
                  label_space_adjust = 0.1, label_font = 1
                  )
CenterText("Percentage of reads with at least 2 sgRNAs that are entirely correct",
           y = par("usr")[[4]] + diff(grconvertY(c(0, 1.3), from = "lines", to = "user")),
           use_cex = par("cex") * 1.15
           )
dev.off()



emf(file.path(manuscript_directory, "Thesis", "Num_contaminated_reads.emf"),
    width = 7.5, height = 6, emfPlus = FALSE, coordDPI = 1500
    )
WellLayoutBarplot(summary_df,
                  "Num_contaminated_reads",
                  sg_sequences_df,
                  label_rows_with_letters = TRUE, label_cex = 0.8,
                  label_space_adjust = 0.1, label_font = 1
                  )
CenterText("Percentage of reads with contaminations (sgRNAs from other wells)",
           y = par("usr")[[4]] + diff(grconvertY(c(0, 1.3), from = "lines", to = "user")),
           use_cex = par("cex") * 1.15
           )
dev.off()






