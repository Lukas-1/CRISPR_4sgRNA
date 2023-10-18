### 26th May 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
sql2_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_functions_directory <- file.path(sql2_directory, "1) R functions")
source(file.path(sql2_R_functions_directory, "05) Examining cross-plate contaminations.R"))




# Define folder paths -----------------------------------------------------

sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
output_directory         <- file.path(sql2_directory, "5) Output", "Tables", "Aligned contaminations")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




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
     file = file.path(sql2_R_objects_directory, "19) Examine cross-plate contaminations.RData")
     )



