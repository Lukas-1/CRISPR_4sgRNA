### 14th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
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




# Identify the main TSS by its inclusion in the Caprano library -----------

TSS_IDs_fac <- factor(merged_replaced_CRISPRa_df[["AltTSS_ID"]],
                      levels = unique(merged_replaced_CRISPRa_df[["AltTSS_ID"]])
                      )
CheckThatFactorIsInOrder(TSS_IDs_fac)
include_Caprano <- tapply(merged_replaced_CRISPRa_df[["Source"]],
                          TSS_IDs_fac,
                          function(x) any(grepl("Caprano", x, fixed = TRUE))
                          )
merged_replaced_CRISPRa_df[["Is_Caprano_TSS"]] <- rep(include_Caprano,
                                                      as.integer(table(TSS_IDs_fac))
                                                      )


gene_IDs_fac <- factor(merged_replaced_CRISPRa_df[["Combined_ID"]],
                       levels = unique(merged_replaced_CRISPRa_df[["Combined_ID"]])
                       )
CheckThatFactorIsInOrder(gene_IDs_fac)

are_main <- tapply(merged_replaced_CRISPRa_df[["Is_Caprano_TSS"]],
                   gene_IDs_fac,
                   function(x) if (!(any(x))) rep(TRUE, length(x)) else x,
                   simplify = FALSE
                   )

merged_replaced_CRISPRa_df[["Use_as_main_TSS"]] <- unlist(are_main)




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





# Assign synthetic lethality sgRNAs to plates -----------------------------

# No controls are needed

are_synthetic_lethal <- (merged_replaced_CRISPRa_df[["Combined_ID"]] %in% prions_lethality_df[["Entrez_ID"]])
lethal_CRISPRa_df <- merged_replaced_CRISPRa_df[are_synthetic_lethal, ]

are_synthetic_lethal <- (sg4_by_gene_df[["Combined_ID"]] %in% prions_lethality_df[["Entrez_ID"]])
lethal_4sg_df <- sg4_by_gene_df[are_synthetic_lethal, ]
stopifnot(identical(unique(as.integer(table(lethal_4sg_df[["AltTSS_ID"]]))), 4L))
TSS_splits <- lapply(split(lethal_4sg_df[["AltTSS_ID"]], lethal_4sg_df[["Entrez_ID"]]), unique)
matches_vec <- match(prions_lethality_df[["Entrez_ID"]], names(TSS_splits))
prions_lethality_df[["Num_TSSs"]] <- lengths(TSS_splits)[matches_vec]

are_first_plate <- cumsum(prions_lethality_df[["Num_TSSs"]]) <= 96

prions_lethality_df[["Plate_number"]] <- ifelse(are_first_plate, 1L, 2L)

lethal_4sg_by_gene_df <- AllocateAllGuidesToPlates(lethal_CRISPRa_df,
                                                   list("Lethality" = prions_lethality_df[["Entrez_ID"]]),
                                                   num_control_wells = 0,
                                                   reorder_df = FALSE
                                                   )

matches_vec <- match(lethal_4sg_by_gene_df[["Entrez_ID"]], prions_lethality_df[["Entrez_ID"]])
new_order <- order(prions_lethality_df[["Plate_number"]][matches_vec],
                   lethal_4sg_by_gene_df[["Well_number"]]
                   )

lethal_4sg_by_gene_df <- lethal_4sg_by_gene_df[new_order, ]
row.names(lethal_4sg_by_gene_df) <- NULL
lethal_4sg_by_gene_df[["Well_number"]] <- rep(seq_len(nrow(lethal_4sg_by_gene_df) / 4), each = 4)

matches_vec <- match(lethal_4sg_by_gene_df[["Entrez_ID"]], prions_lethality_df[["Entrez_ID"]])
categories_fac <- prions_lethality_df[["Category"]][matches_vec]
categories_vec <- sub(" overlapped top hits", "", as.character(categories_fac), fixed = TRUE)
categories_vec <- gsub(" ", "_", categories_vec, fixed = TRUE)
lethal_4sg_by_gene_df[["Plate_string"]] <- paste0(categories_vec, "__well",
                                                  lethal_4sg_by_gene_df[["Well_number"]],
                                                  paste0("__sg", lethal_4sg_by_gene_df[["Rank"]])
                                                  )

lethal_4sg_by_well_df <- lethal_4sg_by_gene_df
new_order <- order(as.integer(lethal_4sg_by_well_df[["Entrez_ID"]]))
lethal_4sg_by_gene_df <- lethal_4sg_by_well_df[new_order, ]
row.names(lethal_4sg_by_gene_df) <- NULL





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



# Export the prion synthetic lethality plate layouts ----------------------

lethal_4sg_by_gene_df <- lethal_4sg_by_gene_df[, names(lethal_4sg_by_gene_df) != "Sublibrary_4sg"]
lethal_4sg_by_well_df <- lethal_4sg_by_well_df[, names(lethal_4sg_by_well_df) != "Sublibrary_4sg"]

ExportPlates(lethal_4sg_by_gene_df, "Prions_synthetic_lethal_4sg_by_gene", sub_folder = "Plate layout - prion synthetic lethality")
ExportPlates(lethal_4sg_by_well_df, "Prions_synthetic_lethal_4sg_by_well", sub_folder = "Plate layout - prion synthetic lethality")

for (i in 1:4) {
  use_df <- lethal_4sg_by_well_df[lethal_4sg_by_well_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("Prions_lethal_", i),
               sub_folder = "Plate layout - prion synthetic lethality",
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
lethal_4sg_by_gene_df <- lethal_4sg_by_gene_df[, names(lethal_4sg_by_gene_df) != preferred_AF_max_column]
lethal_4sg_by_well_df <- lethal_4sg_by_well_df[, names(lethal_4sg_by_well_df) != preferred_AF_max_column]




# Save data ---------------------------------------------------------------

save(list = c("sg4_by_gene_df", "sg4_by_well_df",
              "membrane_4sg_by_well_df", "membrane_4sg_by_gene_df",
              "lethal_4sg_by_gene_df", "lethal_4sg_by_well_df"
              ),
     file = file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData")
     )


