### 14th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R")) # For GetMainTSS




# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR_4sgRNA"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory     <- file.path(RData_directory, "2) CRISPRa")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")
previous_versions_directory <- file.path(RData_directory, "5) Previous versions of the library")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "17) Read in additional gene lists.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRa_RData_directory, "24) Find all TSSs targeted by each sgRNA - summary data.RData"))
load(file.path(CRISPRa_RData_directory, "25) Export libraries to .tsv files - TF_v1_CRISPRa_df.RData"))





# Ensure that the guide RNAs are ordered by rank --------------------------

new_order <- order(match(merged_replaced_CRISPRa_df[, "AltTSS_ID"], merged_replaced_CRISPRa_df[, "AltTSS_ID"]),
                   merged_replaced_CRISPRa_df[, "Rank"]
                   )
merged_replaced_CRISPRa_df <- merged_replaced_CRISPRa_df[new_order, ]
row.names(merged_replaced_CRISPRa_df) <- NULL





# Add data on other (unintended) targeted TSSs ----------------------------

merged_replaced_CRISPRa_df <- AddOtherTargets(merged_replaced_CRISPRa_df, TSS_targets_df)




# Add data on the main TSS ------------------------------------------------

merged_replaced_CRISPRa_df <- AddMainTSS(merged_replaced_CRISPRa_df)




# Assign the sgRNAs to plates ---------------------------------------------

legacy_PD_4sg_entrezs <- setdiff(PD_4sg_entrezs, c("51142", "11315", "9842")) # This is how the full library was ordered

sg4_by_well_df <- AllocateAllGuides_v2(merged_replaced_CRISPRa_df,
                                       sublibraries_entrezs_list  = sublibraries_all_entrezs_list,
                                       previous_version_CRISPR_df = TF_v1_CRISPRa_df,
                                       candidate_entrezs          = legacy_PD_4sg_entrezs
                                       )
sg4_by_gene_df <- RestoreOriginalOrder(sg4_by_well_df)






#  Combine the TF library with the rest of the 4sg library ----------------

are_TF_controls <- TF_v1_CRISPRa_df[["Is_control"]] %in% "Yes"
TF_v1_CRISPRa_df[["AltTSS_ID"]][are_TF_controls] <- TF_v1_CRISPRa_df[["Combined_ID"]][are_TF_controls]
TF_v1_CRISPRa_df <- RenumberPlatesContinuously(TF_v1_CRISPRa_df)

full_sg4_list <- MergeTFWithRest(sg4_by_well_df, TF_v1_CRISPRa_df)

full_4sg_by_well_df <- full_sg4_list[["full_4sg_by_well_df"]]
full_4sg_by_gene_df <- full_sg4_list[["full_4sg_by_gene_df"]]
rm(full_sg4_list)



# Examine the sub-library allocation --------------------------------------

sg4_allTFs_df <- AllocateAllGuidesToPlates(merged_replaced_CRISPRa_df,
                                           sublibraries_entrezs_list = sublibraries_all_entrezs_list,
                                           num_control_wells = 96
                                           )

are_unique_gene <- !(duplicated(sg4_allTFs_df[["Combined_ID"]]))
table(sg4_allTFs_df[["Sublibrary_4sg"]][are_unique_gene], useNA = "ifany")

are_controls <- sg4_allTFs_df[["Is_control"]] == "Yes"
are_selected_TSSs <- !(are_controls) & !(duplicated(sg4_allTFs_df[["AltTSS_ID"]]))

sum(are_selected_TSSs) # Number of TSSs targeted by the library

table(sg4_allTFs_df[["Num_TSSs"]][!(are_controls) & are_unique_gene] >= 2) # Number of genes with multiple TSSs

table(full_4sg_by_well_df[["Is_control"]]) / 4 # Number of controls





# Assign just the PD genes to plates --------------------------------------

are_PD <- (merged_replaced_CRISPRa_df[["Combined_ID"]] %in% PD_all_entrezs)
PD_CRISPRa_df <- merged_replaced_CRISPRa_df[are_PD, ]

PD_4sg_by_gene_df <- AllocateAllGuidesToPlates(PD_CRISPRa_df,
                                               list("PD" = PD_all_entrezs),
                                               num_control_wells = 0,
                                               reorder_df = FALSE
                                               )

PD_4sg_by_gene_df <- AssignPlateStrings(PD_4sg_by_gene_df, use_prefix = "PD_")
PD_4sg_by_gene_df[["Plate_string"]] <- sub("_1_", "_", PD_4sg_by_gene_df[["Plate_string"]], fixed = TRUE)

PD_4sg_by_well_df <- ReorderPlates(PD_4sg_by_gene_df)

PD_4sg_by_gene_df[["Sublibrary_4sg"]] <- NULL
PD_4sg_by_well_df[["Sublibrary_4sg"]] <- NULL




# Export all shared / duplicated sgRNAs -----------------------------------

shared_sgRNAs_df <- SharedsgRNAsDf(full_4sg_by_gene_df)

ExportSharedDf(shared_sgRNAs_df,
               file.path("4sg plate layout (complete)",
                         "Duplicated CRISPRa sgRNAs (shared between genes)"
                         )
               )




# Export the plate layouts for the complete 4sg library -------------------

ExportPlates(full_4sg_by_gene_df,
             "All_sublibraries_ordered_by_gene",
             "4sg plate layout (complete)",
             add_primers = FALSE
             )
ExportPlates(full_4sg_by_well_df,
             "All_sublibraries_ordered_by_well",
             "4sg plate layout (complete)",
             add_primers = FALSE
             )




# Export the plate layouts for the 4sg library (excluding TFs) ------------

ExportPlates(sg4_by_gene_df, "4sg_non_TF_ordered_by_gene", "4sg plate layout (non-TF)")
ExportPlates(sg4_by_well_df, "4sg_non_TF_ordered_by_well", "4sg plate layout (non-TF)")

for (i in 1:4) {
  use_df <- sg4_by_well_df[sg4_by_well_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("4sg_non_TF_sg", i),
               "4sg plate layout (non-TF)",
               add_padding_between_plates = TRUE
               )
}




# Export the PD plate layout ----------------------------------------------

ExportPlates(PD_4sg_by_gene_df, "PD_ordered_by_gene", sub_folder = "PD plate layout")
ExportPlates(PD_4sg_by_well_df, "PD_ordered_by_well", sub_folder = "PD plate layout")

for (i in 1:4) {
  use_df <- PD_4sg_by_well_df[PD_4sg_by_well_df[["Rank"]] %in% i, ]
  ExportPlates(use_df, paste0("PD_sg", i), sub_folder = "PD plate layout")
}




# Export sparse data / only the most important columns --------------------

export_columns <- c("Sublibrary_4sg", "Plate_number", "Well_number",
                    "Entrez_ID", "Gene_symbol", "TSS_ID",
                    "sgRNA_sequence"
                    )

ExportPlates(full_4sg_by_well_df,
             "CRISPRa_4sg_by_well",
             sub_folder  = "4sg plate layout (complete)/Sparse",
             add_colors  = FALSE,
             add_primers = FALSE
             )




# Save data ---------------------------------------------------------------

save(list = c("full_4sg_by_gene_df", "full_4sg_by_well_df",
              "sg4_by_gene_df", "sg4_by_well_df",
              "shared_sgRNAs_df"
              ),
     file = file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData")
     )

save(list = c("PD_4sg_by_gene_df", "PD_4sg_by_well_df"),
     file = file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates - PD genes.RData")
     )




