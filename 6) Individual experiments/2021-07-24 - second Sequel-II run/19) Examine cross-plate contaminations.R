### 31st July 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(s2r1_R_functions_directory, "05) Examining cross-plate contaminations.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
output_directory         <- file.path(s2r2_directory, "5) Output", "Tables", "Aligned contaminations")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Integrate data on thresholds passed -------------------------------------

contam_df <- AddFilterPassed(contam_df)




# Explore the sources of contaminations -----------------------------------

table(contam_df[["Are_cross_plate"]])

table(contam_df[["Are_cross_plate"]][contam_df[["Filter_passed"]] != "None"])
table(contam_df[["Are_cross_plate"]][contam_df[["Filter_passed"]] == "CCS7"])

head(sort(table(contam_df[["Reference_ID"]][contam_df[["Are_cross_plate"]]]), decreasing = TRUE))
head(sort(table(contam_df[["Contaminating_ID"]][contam_df[["Are_cross_plate"]]]), decreasing = TRUE))

table(contam_df[["Reference_plate_number"]][contam_df[["Are_cross_plate"]]])

length(unique(contam_df[["ZMW"]][contam_df[["Are_cross_plate"]]]))
length(unique(contam_df[["ZMW"]]))

length(unique(contam_df[["ZMW"]][contam_df[["Filter_passed"]] != "None"]))
length(unique(contam_df[["ZMW"]][contam_df[["Are_cross_plate"]] & (contam_df[["Filter_passed"]] != "None")]))




# Prepare for export ------------------------------------------------------

export_contam_df <- MakeExportContamDf(contam_df)




# Export data -------------------------------------------------------------

write.table(export_contam_df,
            file      = file.path(output_directory, "Aligned contaminations.tsv"),
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE
            )



# Save data ---------------------------------------------------------------

save(list = "contam_df",
     file = file.path(s2r2_R_objects_directory, "19) Check for cross-plate contaminations.RData")
     )



