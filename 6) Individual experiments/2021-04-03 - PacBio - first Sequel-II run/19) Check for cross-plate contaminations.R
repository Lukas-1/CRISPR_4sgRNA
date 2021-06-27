### 26th May 2021 ###




# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
output_directory         <- file.path(sql2_directory, "5) Output", "Tables", "Aligned contaminations")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Integrate data on thresholds passed -------------------------------------

threshold_vec <- ifelse(contam_df[["ZMW"]] %in% ccs7_df_list[["individual_reads_df"]][["ZMW"]],
                        "CCS7",
                        ifelse(contam_df[["ZMW"]] %in% ccs5_df_list[["individual_reads_df"]][["ZMW"]],
                               "CCS5",
                               ifelse(contam_df[["ZMW"]] %in% ccs3_df_list[["individual_reads_df"]][["ZMW"]],
                                      "CCS3",
                                      "None"
                                      )
                               )
                        )
contam_df[["Filter_passed"]] <- threshold_vec




# Re-order the columns ----------------------------------------------------

move_after_column <- "Are_cross_plate"
insert_column <- "Filter_passed"
move_after_index <- which(names(contam_df) == move_after_column)
are_before <- seq_along(contam_df) < move_after_index
are_after <- seq_along(contam_df) > move_after_index
are_current <- names(contam_df) == insert_column
columns_reordered <- c(names(contam_df)[are_before],
                       move_after_column, insert_column,
                       names(contam_df)[are_after & !(are_current)]
                       )
contam_df <- contam_df[, columns_reordered]



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

export_contam_df <- contam_df
export_contam_df <- export_contam_df[order(export_contam_df[["Are_cross_plate"]], decreasing = TRUE), ]
export_contam_df <- export_contam_df[export_contam_df[["Filter_passed"]] != "None", ]

row.names(export_contam_df) <- NULL

export_contam_df[["Are_cross_plate"]] <- ifelse(export_contam_df[["Are_cross_plate"]], "Yes", "No")
names(export_contam_df)[names(export_contam_df) == "Are_cross_plate"] <- "Cross_plate"




# Export data -------------------------------------------------------------

write.table(export_contam_df,
            file      = file.path(output_directory, "Aligned contaminations.tsv"),
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE
            )



# Save data ---------------------------------------------------------------

save(list = "contam_df",
     file = file.path(sql2_R_objects_directory, "19) Check for cross-plate contaminations.RData")
     )



