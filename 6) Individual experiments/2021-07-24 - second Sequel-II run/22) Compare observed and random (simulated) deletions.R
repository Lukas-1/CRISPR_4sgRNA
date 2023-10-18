### 3rd August 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "22) Comparing observed and simulated deletions.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
file_output_directory    <- file.path(s2r2_directory, "5) Output", "Figures", "Compare simulated and observed deletions")
PNGs_directory           <- file.path(s2r2_directory, "5) Output", "PNGs", "Compare simulated and observed deletions")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "21) Simulate random deletions (to compare with observed deletions).RData"))





# Explore the metrics -----------------------------------------------------

colMeans(one_read_simul_mat) * 100
colSums(one_read_del_df[, colnames(one_read_simul_mat)])  / nrow(one_read_del_df) * 100

colMeans(all_reads_simul_mat) * 100
colSums(all_reads_del_df[, colnames(one_read_simul_mat)]) / nrow(all_reads_del_df) * 100





# Draw plots --------------------------------------------------------------

SimVsObsDeletions(one_read_simul_mat, one_read_del_df)
SimVsObsDeletions(all_reads_simul_mat, all_reads_del_df)
SimVsObsDeletions(one_read_simul_mat[, c(1, 2, 4)], one_read_del_df)


DrawAllSimVsObsPlots(export_PNGs = FALSE)





