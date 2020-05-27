### 14th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRi")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRi_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))

load(file.path(CRISPRi_RData_directory, "24) Export libraries to .tsv files - vacuolation_entrezs.RData"))





# Assign vacuolation sgRNAs to plates -------------------------------------

are_vacuolation <- (merged_replaced_CRISPRi_df[["Combined_ID"]] %in% vacuolation_entrezs) |
                   (merged_replaced_CRISPRi_df[["Is_control"]] %in% "Yes")
vacuolation_CRISPRi_df <- merged_replaced_CRISPRi_df[are_vacuolation, ]





# Add the control guides to vacuolation plate 2 ---------------------------

vac_4sg_df <- AllocateAllGuidesToPlates(vacuolation_CRISPRi_df,
                                        list("Vacuolation" = vacuolation_entrezs),
                                        num_control_wells = 10,
                                        reorder_df = FALSE
                                        )

vac_4sg_df <- AssignPlateStrings(vac_4sg_df, use_prefix = "vac")




# Re-order the data frame according to the plate layout -------------------

vac_4sg_reordered_df <- ReorderPlates(vac_4sg_df)





# Move the control guides to vacuolation plate 2 --------------------------

are_controls <- vac_4sg_reordered_df[["Is_control"]] == "Yes"
max_well_plate_2 <- max(vac_4sg_reordered_df[["Well_number"]][(vac_4sg_reordered_df[["Plate_number"]] == 2) & !(are_controls)])
control_wells_seq <- seq(from = max_well_plate_2 + 1L,
                         to = max_well_plate_2 + sum(are_controls),
                         )
vac_4sg_reordered_df[["Well_number"]][are_controls] <- control_wells_seq
vac_4sg_reordered_df[["Plate_number"]][are_controls] <- 2L




# Export the plate layouts ------------------------------------------------

vac_4sg_df <- vac_4sg_df[, colnames(vac_4sg_df) != "Sublibrary_4sg"]
vac_4sg_reordered_df <- vac_4sg_reordered_df[, colnames(vac_4sg_reordered_df) != "Sublibrary_4sg"]

ExportPlates(vac_4sg_df, "Vacuolation_4sg_original_order", sub_folder = "Plate layouts")
ExportPlates(vac_4sg_reordered_df, "Vacuolation_4sg_reordered", sub_folder = "Plate layouts")

for (i in 1:4) {
  use_df <- vac_4sg_reordered_df[vac_4sg_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("Vacuolation_4sg_reordered_sg", i),
               sub_folder = "Plate layouts",
               add_padding_between_plates = TRUE
               )
}







# Save data ---------------------------------------------------------------

save(list = c("vac_4sg_df", "vac_4sg_reordered_df"),
     file = file.path(CRISPRi_RData_directory, "27) Distribute sgRNAs for the whole genome onto plates.RData")
     )





