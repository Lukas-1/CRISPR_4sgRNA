### 31st May 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")
source(file.path(R_functions_directory, "14) Drawing sigmoid curves (acceptable wells vs. percentage sampled).R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Sub-sampled", "Sigmoid curves")





# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed reads - with subsampling.RData"))





# Draw plots --------------------------------------------------------------

filtered_results_mat <- DrawFourCurves(subsampled_list, "filtered_summary_df")

summary_df_names <- c("original_summary_df", "filtered_summary_df", "filtered_gRNAs_df")

for (df_i in seq_along(summary_df_names)) {
  file_name <- paste0("Sigmoid curves - ",
                      c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[df_i]]
                      )

  pdf(file = file.path(plots_output_directory, paste0(file_name, ".pdf")),
      height = 7, width  = 7
      )
  old_mar <- par(mar = c(5, 5.5, 5, 2.5))

  data_mat <- DrawFourCurves(subsampled_list, summary_df_names[[df_i]])

  write.table(data_mat,
              file      = file.path(plots_output_directory, "Tables", paste0(file_name, ".tsv")),
              sep       = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE
              )

  par(old_mar)
  dev.off()
}



