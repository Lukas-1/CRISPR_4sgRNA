### 14th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory     <- file.path(RData_directory, "7) Mouse - CRISPRa")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "Mouse - CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "06) Read in gene lists.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))






# Assign membrane sgRNAs to plates ----------------------------------------

are_membrane <- (merged_replaced_CRISPRa_df[["Combined_ID"]] %in% membrane_het_df[["Entrez_ID"]]) |
                (merged_replaced_CRISPRa_df[["Is_control"]] %in% "Yes")
membrane_CRISPRa_df <- merged_replaced_CRISPRa_df[are_membrane, ]





# Add the control guides to the membrane plates ---------------------------

membrane_4sg_df <- AllocateAllGuidesToPlates(membrane_CRISPRa_df,
                                             list("Membrane" = membrane_het_df[["Entrez_ID"]]),
                                             num_control_wells = 6,
                                             reorder_df = FALSE
                                             )




# Re-order the data frame according to the plate layout -------------------

membrane_4sg_df[["Old_order"]] <- seq_len(nrow(membrane_4sg_df))
membrane_4sg_reordered_df <- ReorderPlates(membrane_4sg_df)





# Move the control guides to the last plate -------------------------------

are_controls <- membrane_4sg_reordered_df[["Is_control"]] == "Yes"
last_plate_number <- max(membrane_4sg_reordered_df[["Plate_number"]])
max_well_last_plate <- max(membrane_4sg_reordered_df[["Well_number"]][(membrane_4sg_reordered_df[["Plate_number"]] == last_plate_number) & !(are_controls)])
control_wells_seq <- rep(seq(from = max_well_last_plate + 1L,
                             to = max_well_last_plate + sum(are_controls) / 4,
                             ),
                         each = 4
                         )
membrane_4sg_reordered_df[["Well_number"]][are_controls] <- control_wells_seq
membrane_4sg_reordered_df[["Plate_number"]][are_controls] <- last_plate_number
membrane_4sg_reordered_df <- AssignPlateStrings(membrane_4sg_reordered_df, use_prefix = "mmb_")



# Restore the original plate order ----------------------------------------

membrane_4sg_df <- membrane_4sg_reordered_df[order(membrane_4sg_reordered_df[["Old_order"]]), ]
row.names(membrane_4sg_df) <- NULL
membrane_4sg_reordered_df <- membrane_4sg_reordered_df[, names(membrane_4sg_reordered_df) != "Old_order"]
membrane_4sg_df <- membrane_4sg_df[, names(membrane_4sg_df) != "Old_order"]



# Export the membrane heterogeneity plate layouts -------------------------

membrane_4sg_df <- membrane_4sg_df[, names(membrane_4sg_df) != "Sublibrary_4sg"]
membrane_4sg_reordered_df <- membrane_4sg_reordered_df[, names(membrane_4sg_reordered_df) != "Sublibrary_4sg"]

ExportPlates(membrane_4sg_df, "Membrane_4sg_original_order", sub_folder = "Plate layout - membrane")
ExportPlates(membrane_4sg_reordered_df, "Membrane_4sg_reordered", sub_folder = "Plate layout - membrane")

for (i in 1:4) {
  use_df <- membrane_4sg_reordered_df[membrane_4sg_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("Membrane_4sg_reordered_sg", i),
               sub_folder = "Plate layout - membrane",
               add_padding_between_plates = TRUE
               )
}




# Save data ---------------------------------------------------------------

save(list = c("membrane_4sg_df", "membrane_4sg_reordered_df"),
     file = file.path(CRISPRa_RData_directory, "27) Distribute sgRNAs for the whole genome onto plates.RData")
     )


