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

load(file.path(CRISPRi_RData_directory, "20) Create a gene-based summary of the human genome - vacuolation_entrezs.RData"))






# Assign all guides to plates ---------------------------------------------

sg4_df <- AllocateAllGuidesToPlates(merged_replaced_CRISPRi_df,
                                    sublibraries_entrezs_list = sublibraries_all_entrezs_list,
                                    num_control_wells = 96,
                                    reorder_df = FALSE
                                    )
sg4_df[["Old_order"]] <- seq_len(nrow(sg4_df))

sg4_reordered_df <- ReorderPlates(sg4_df)
sg4_reordered_df <- RenumberPlatesContinuously(sg4_reordered_df)
sg4_reordered_df <- AssignPlateStrings(sg4_reordered_df)

sg4_df <- sg4_reordered_df[order(sg4_reordered_df[["Old_order"]]), ]

sg4_df <- sg4_df[, colnames(sg4_df) != "Old_order"]
sg4_reordered_df <- sg4_reordered_df[, colnames(sg4_reordered_df) != "Old_order"]




# Examine the sub-library allocation --------------------------------------

are_selected <- (sg4_df[["Is_control"]] == "No") &
                 !(duplicated(sg4_df[["Combined_ID"]]))

table(sg4_df[["Sublibrary_4sg"]][are_selected])

are_selected_TSSs <- (sg4_df[["Is_control"]] == "No") &
                     !(duplicated(sg4_df[["AltTSS_ID"]]))

sum(are_selected_TSSs) # Number of TSSs targeted by the library





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


vac_4sg_df[(vac_4sg_df[["Rank"]] %in% 1) & (vac_4sg_df[["Is_control"]] %in% "Yes"), ]





# Re-order the data frame according to the plate layout -------------------

vac_4sg_df[["Old_order"]] <- seq_len(nrow(vac_4sg_df))
vac_4sg_reordered_df <- ReorderPlates(vac_4sg_df)





# Move the control guides to vacuolation plate 2 --------------------------

are_controls <- vac_4sg_reordered_df[["Is_control"]] == "Yes"
max_well_plate_2 <- max(vac_4sg_reordered_df[["Well_number"]][(vac_4sg_reordered_df[["Plate_number"]] == 2) & !(are_controls)])
control_wells_seq <- rep(seq(from = max_well_plate_2 + 1L,
                             to = max_well_plate_2 + sum(are_controls) / 4,
                             ),
                         each = 4
                         )
vac_4sg_reordered_df[["Well_number"]][are_controls] <- control_wells_seq
vac_4sg_reordered_df[["Plate_number"]][are_controls] <- 2L
vac_4sg_reordered_df <- AssignPlateStrings(vac_4sg_reordered_df, use_prefix = "vac")




# Restore the original plate order ----------------------------------------

vac_4sg_df <- vac_4sg_reordered_df[order(vac_4sg_reordered_df[["Old_order"]]), ]
row.names(vac_4sg_df) <- NULL
vac_4sg_reordered_df <- vac_4sg_reordered_df[, colnames(vac_4sg_reordered_df) != "Old_order"]
vac_4sg_df <- vac_4sg_df[, colnames(vac_4sg_df) != "Old_order"]




# Export the vacuolation plate layouts ------------------------------------

vac_4sg_df <- vac_4sg_df[, colnames(vac_4sg_df) != "Sublibrary_4sg"]
vac_4sg_reordered_df <- vac_4sg_reordered_df[, colnames(vac_4sg_reordered_df) != "Sublibrary_4sg"]

ExportPlates(vac_4sg_df, "Vacuolation_4sg_original_order", sub_folder = "Plate layout - vacuolation")
ExportPlates(vac_4sg_reordered_df, "Vacuolation_4sg_reordered", sub_folder = "Plate layout - vacuolation")

for (i in 1:4) {
  use_df <- vac_4sg_reordered_df[vac_4sg_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("Vacuolation_4sg_reordered_sg", i),
               sub_folder = "Plate layout - vacuolation",
               add_padding_between_plates = TRUE
               )
}



# Export the whole-genome plate layouts -----------------------------------

use_sub_folder <- "Plate layout - all genes"

ExportPlates(sg4_df, "All_sublibraries_original_order", sub_folder = use_sub_folder)
ExportPlates(sg4_reordered_df, "All_sublibraries_reordered", sub_folder = use_sub_folder)

for (i in 1:4) {
  use_df <- sg4_reordered_df[sg4_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("4sg_reordered_sg", i),
               sub_folder = use_sub_folder,
               add_padding_between_plates = TRUE
               )
}





# Save data ---------------------------------------------------------------

save(list = c("sg4_reordered_df", "sg4_df",
              "vac_4sg_df", "vac_4sg_reordered_df"
              ),
     file = file.path(CRISPRi_RData_directory, "27) Distribute sgRNAs for the whole genome onto plates.RData")
     )





