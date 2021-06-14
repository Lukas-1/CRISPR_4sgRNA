### 26th May 2021 ###




# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
output_directory         <- file.path(sql2_directory, "5) Output", "Tables", "Aligned contaminations")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Look for cross-plate contaminations -------------------------------------

are_cross_plate <- rep(NA, nrow(extracted_df))

guides_mat <- as.matrix(sg_sequences_df[, paste0("Sequence_sg", 1:4)])
all_guides <- as.character(guides_mat)

are_sg <- extracted_df[["Feature"]] %in% paste0("sg", 1:4)
are_sg_contam <- extracted_df[["Is_contamination"]] & are_sg

for (plate_number in unique(sg_sequences_df[["Plate_number"]])) {
  are_eligible <- (extracted_df[["Plate_number"]] == plate_number) & are_sg_contam
  this_plate_guides <- as.character(guides_mat[sg_sequences_df[["Plate_number"]] == plate_number, ])
  other_plate_guides <- setdiff(all_guides, this_plate_guides)
  sub_cross_plate <- extracted_df[["Aligned_read"]][are_eligible] %in% other_plate_guides
  are_cross_plate[are_eligible] <- sub_cross_plate
}




# Explore contaminated wells ----------------------------------------------

table(are_cross_plate[are_sg_contam])

table(extracted_df[["Combined_ID"]][are_cross_plate %in% TRUE])
table(extracted_df[["Plate_number"]][are_cross_plate %in% TRUE])

length(unique(extracted_df[["ZMW"]][are_cross_plate %in% TRUE]))
length(unique(extracted_df[["ZMW"]][are_sg_contam]))





# Create a helper data frame for finding contaminations -------------------

long_df_list <- lapply(1:4,
                       function(x) data.frame(sg_sequences_df[, c("Combined_ID", "Plate_number", "Well_number")],
                                              "Sg_number" = x,
                                              "Sequence"  = guides_mat[, x],
                                              "Num_occurrences" = sg_sequences_df[, paste0("Num_occurrences_sg", x)],
                                              stringsAsFactors = FALSE,
                                              row.names = NULL
                                              )
                       )
long_df <- do.call(rbind.data.frame, c(long_df_list, stringsAsFactors = FALSE))

new_order <- order(long_df[, "Plate_number"], long_df[, "Well_number"])
long_df <- long_df[new_order, ]
row.names(long_df) <- NULL




# Integrate data on the contaminating well --------------------------------

contam_matches <- match(extracted_df[["Aligned_read"]][are_sg_contam], long_df[["Sequence"]])

contam_df <- long_df[contam_matches, ]
row.names(contam_df) <- NULL
names(contam_df) <- paste0("Contaminating_",
                           tolower(substr(names(contam_df), 1, 1)),
                           substr(names(contam_df), 2, nchar(names(contam_df)))
                           )

redundant_columns <- c("Num_missing", "Mostly_deleted", "Is_correct",
                       "Category", "Is_contamination", "Aligned_template"
                       )
contam_df <- data.frame(extracted_df[are_sg_contam, !(names(extracted_df) %in% redundant_columns)],
                        "Are_cross_plate" = are_cross_plate[are_sg_contam],
                        contam_df[, names(contam_df) != "Contaminating_sequence"],
                        stringsAsFactors = FALSE,
                        row.names = NULL
                        )




# Tidy the data frame on contaminations -----------------------------------

contam_df[["Reference_sg_number"]] <- as.integer(substr(contam_df[["Feature"]], 3, 3))
contam_df <- contam_df[, names(contam_df) != "Feature"]

contam_df[["Filter_passed"]] <- NA

rename_vec <- c(
                                    "ZMW",
                                    "Are_cross_plate",
                                    "Filter_passed",
  "Combined_ID"                   = "Reference_ID",
  "Contaminating_combined_ID"     = "Contaminating_ID",
  "Plate_number"                  = "Reference_plate_number",
                                    "Contaminating_plate_number",
  "Well_number"                   = "Reference_well_number",
                                    "Contaminating_well_number",
                                    "Reference_sg_number",
                                    "Contaminating_sg_number",
  "Template"                      = "Reference_sequence",
  "Aligned_read"                  = "Contaminating_sequence",
                                    "Quality",
                                    "Mean_quality",
  "Contaminating_num_occurrences" = "Contaminating_sequence_occurrences_in_library"
)

for (column_name in setdiff(names(rename_vec), "")) {
  names(contam_df)[names(contam_df) == column_name] <- rename_vec[[column_name]]
}

contam_df <- contam_df[, order(match(names(contam_df), rename_vec))]

contam_df[["Mean_quality"]] <- contam_df[["Mean_quality"]] / 93





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





# Explore the sources of contaminations -----------------------------------

table(contam_df[["Contaminating_ID"]][contam_df[["Are_cross_plate"]]])

table(contam_df[["Are_cross_plate"]][contam_df[["Filter_passed"]] != "None"])
table(contam_df[["Are_cross_plate"]][contam_df[["Filter_passed"]] == "CCS7"])
table(contam_df[["Are_cross_plate"]])


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
     file = file.path(sql2_R_objects_directory, "17) Check for cross-plate contaminations.RData")
     )




# are_identical <- vapply(names(contam_df), function(x) {
#   identical(contam_df[[x]], new_contam_df[[x]])
# }, logical(1))
#

















