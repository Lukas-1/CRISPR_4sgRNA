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
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "Mouse - CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "06) Read in gene lists.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))






# Create an empty SNP column ----------------------------------------------

merged_replaced_CRISPRa_df[[preferred_AF_max_column]] <- NA_real_






# Assign all guides to plates ---------------------------------------------

sg4_by_gene_df <- AllocateAllGuidesToPlates(merged_replaced_CRISPRa_df,
                                            sublibraries_entrezs_list = list("All" = collected_entrez_IDs),
                                            num_control_wells = 96,
                                            reorder_df = FALSE
                                            )
sg4_by_gene_df[["Old_order"]] <- seq_len(nrow(sg4_by_gene_df))

sg4_by_well_df <- ReorderPlates(sg4_by_gene_df)
sg4_by_well_df <- RenumberPlatesContinuously(sg4_by_well_df)
sg4_by_well_df <- AssignPlateStrings(sg4_by_well_df)

sg4_by_gene_df <- sg4_by_well_df[order(sg4_by_well_df[["Old_order"]]), ]

sg4_by_gene_df <- sg4_by_gene_df[, names(sg4_by_gene_df) != "Old_order"]
sg4_by_well_df <- sg4_by_well_df[, names(sg4_by_well_df) != "Old_order"]






# Assign membrane sgRNAs to plates ----------------------------------------

are_membrane <- (merged_replaced_CRISPRa_df[["Combined_ID"]] %in% membrane_het_df[["Entrez_ID"]]) |
                (merged_replaced_CRISPRa_df[["Is_control"]] %in% "Yes")
membrane_CRISPRa_df <- merged_replaced_CRISPRa_df[are_membrane, ]

membrane_4sg_by_gene_df <- AllocateAllGuidesToPlates(membrane_CRISPRa_df,
                                                     list("Membrane" = membrane_het_df[["Entrez_ID"]]),
                                                     num_control_wells = 6,
                                                     reorder_df = FALSE
                                                     )




# Re-order the membrane sgRNAs according to the plate layout --------------

membrane_4sg_by_gene_df[["Old_order"]] <- seq_len(nrow(membrane_4sg_by_gene_df))
membrane_4sg_by_well_df <- ReorderPlates(membrane_4sg_by_gene_df)





# Move the membrane control guides to the last plate ----------------------

are_controls <- membrane_4sg_by_well_df[["Is_control"]] == "Yes"
last_plate_number <- max(membrane_4sg_by_well_df[["Plate_number"]])
are_non_control_last <- (membrane_4sg_by_well_df[["Plate_number"]] == last_plate_number) & !(are_controls)
max_well_last_plate <-
max(membrane_4sg_by_well_df[["Well_number"]][are_non_control_last])
control_wells_seq <- rep(seq(from = max_well_last_plate + 1L,
                             to = max_well_last_plate + sum(are_controls) / 4,
                             ),
                         each = 4
                         )
membrane_4sg_by_well_df[["Well_number"]][are_controls] <- control_wells_seq
membrane_4sg_by_well_df[["Plate_number"]][are_controls] <- last_plate_number
membrane_4sg_by_well_df <- AssignPlateStrings(membrane_4sg_by_well_df, use_prefix = "mmb_")




# Restore the original plate order ----------------------------------------

membrane_4sg_by_gene_df <- membrane_4sg_by_gene_df[order(membrane_4sg_by_gene_df[["Old_order"]]), ]
row.names(membrane_4sg_by_gene_df) <- NULL
membrane_4sg_by_well_df <- membrane_4sg_by_well_df[, names(membrane_4sg_by_well_df) != "Old_order"]
membrane_4sg_by_gene_df <- membrane_4sg_by_gene_df[, names(membrane_4sg_by_gene_df) != "Old_order"]





# Export all shared / duplicated sgRNAs -----------------------------------

shared_sgRNAs_df <- SharedsgRNAsDf(sg4_by_gene_df)

use_sub_folder <- "Plate layout - all genes"

ExportSharedDf(shared_sgRNAs_df,
               file.path(use_sub_folder,
                         "Duplicated mouse CRISPRa sgRNAs (shared between genes)"
                         )
               )





# Export the membrane heterogeneity plate layouts -------------------------

membrane_4sg_by_gene_df <- membrane_4sg_by_gene_df[, names(membrane_4sg_by_gene_df) != "Sublibrary_4sg"]
membrane_4sg_by_well_df <- membrane_4sg_by_well_df[, names(membrane_4sg_by_well_df) != "Sublibrary_4sg"]

export_columns <- setdiff(export_columns, preferred_AF_max_column)

ExportPlates(membrane_4sg_by_gene_df, "Membrane_4sg_by_gene", sub_folder = "Plate layout - membrane")
ExportPlates(membrane_4sg_by_well_df, "Membrane_4sg_by_well", sub_folder = "Plate layout - membrane")

for (i in 1:4) {
  use_df <- membrane_4sg_by_gene_df[membrane_4sg_by_gene_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("PD_sg", i),
               sub_folder = "Plate layout - membrane",
               add_padding_between_plates = TRUE
               )
}




# Export the whole-genome plate layouts -----------------------------------

ExportPlates(sg4_by_gene_df,
             "All_sublibraries_ordered_by_gene",
             use_sub_folder,
             add_primers = FALSE
             )
ExportPlates(sg4_by_well_df,
             "All_sublibraries_ordered_by_well",
             use_sub_folder,
             add_primers = FALSE
             )

for (i in 1:4) {
  use_df <- sg4_by_well_df[sg4_by_well_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("4sg_ordered_by_well_sg", i),
               sub_folder = use_sub_folder,
               add_padding_between_plates = TRUE
               )
}




# Remove the empty SNP column ---------------------------------------------

membrane_4sg_by_gene_df <- membrane_4sg_by_gene_df[, names(membrane_4sg_by_gene_df) != preferred_AF_max_column]
membrane_4sg_by_well_df <- membrane_4sg_by_well_df[, names(membrane_4sg_by_well_df) != preferred_AF_max_column]




# Save data ---------------------------------------------------------------

save(list = c("sg4_by_gene_df", "sg4_by_well_df",
              "membrane_4sg_by_well_df", "membrane_4sg_by_gene_df"
              ),
     file = file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData")
     )


