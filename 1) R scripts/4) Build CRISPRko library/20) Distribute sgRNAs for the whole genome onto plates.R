### 14th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R")) # For GetMinEntrez
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory    <- file.path(RData_directory, "3) CRISPRko")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko")
previous_versions_directory <- file.path(RData_directory, "5) Previous versions of the library")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "17) Read in additional gene lists.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))
load(file.path(CRISPRko_RData_directory, "16) Find all genes targeted by each sgRNA - summary data.RData"))
load(file.path(CRISPRko_RData_directory, "17) Export libraries to .tsv files - TF_v1_CRISPRko_df.RData"))






# Add data on other (unintended) targeted genes ---------------------------

merged_CRISPRko_df <- AddOtherTargets(merged_CRISPRko_df,
                                      guides_CDS_df,
                                      deletions_CDS_df
                                      )




# Assign the sgRNAs to plates ---------------------------------------------

legacy_PD_4sg_entrezs <- setdiff(PD_4sg_entrezs, c("51142", "11315", "9842")) # This is how the full library was ordered

sg4_by_well_df <- AllocateAllGuides_v2(merged_CRISPRko_df,
                                         sublibraries_entrezs_list  = sublibraries_all_entrezs_list,
                                         previous_version_CRISPR_df = TF_v1_CRISPRko_df,
                                         candidate_entrezs          = legacy_PD_4sg_entrezs
                                         )
sg4_by_gene_df <- RestoreOriginalOrder(sg4_by_well_df)





#  Combine the TF library with the rest of the 4sg library ----------------

full_sg4_list <- MergeTFWithRest(sg4_by_well_df, TF_v1_CRISPRko_df)

full_4sg_by_well_df <- full_sg4_list[["full_4sg_by_well_df"]]
full_4sg_by_gene_df <- full_sg4_list[["full_4sg_by_gene_df"]]
rm(full_sg4_list)





# Examine the sub-library allocation --------------------------------------

sg4_allTFs_df <- AllocateAllGuidesToPlates(merged_CRISPRko_df,
                                           sublibraries_entrezs_list = sublibraries_all_entrezs_list,
                                           num_control_wells = 96
                                           )

are_selected <- (sg4_allTFs_df[["Is_control"]] == "No") &
                (sg4_allTFs_df[["Rank"]] == 1)

table(sg4_allTFs_df[["Sublibrary_4sg"]][are_selected])





# Assign just the PD genes to plates --------------------------------------

are_PD <- (merged_CRISPRko_df[["Combined_ID"]] %in% PD_all_entrezs) #|
          #(merged_replaced_CRISPRi_df[["Is_control"]] %in% "Yes")
PD_CRISPRko_df <- merged_CRISPRko_df[are_PD, ]

PD_4sg_df <- AllocateAllGuidesToPlates(PD_CRISPRko_df,
                                       list("PD" = PD_all_entrezs),
                                       num_control_wells = 0,
                                       reorder_df = FALSE
                                       )
PD_4sg_df <- AssignPlateStrings(PD_4sg_df, use_prefix = "PD_")
PD_4sg_df[["Plate_string"]] <- sub("_1_", "_", PD_4sg_df[["Plate_string"]], fixed = TRUE)

PD_4sg_reordered_df <- ReorderPlates(PD_4sg_df)

PD_4sg_df[["Sublibrary_4sg"]] <- NULL
PD_4sg_reordered_df[["Sublibrary_4sg"]] <- NULL





# Export the plate layouts for the complete 4sg library -------------------

ExportPlates(full_4sg_by_gene_df, "All_sublibraries_ordered_by_gene", "4sg plate layout (complete)")
ExportPlates(full_4sg_by_well_df, "All_sublibraries_ordered_by_well", "4sg plate layout (complete)")





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

ExportPlates(PD_4sg_df, "PD_ordered_by_gene", sub_folder = "PD plate layout")
ExportPlates(PD_4sg_reordered_df, "PD_ordered_by_well", sub_folder = "PD plate layout")

for (i in 1:4) {
  use_df <- PD_4sg_reordered_df[PD_4sg_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df, paste0("PD_sg", i), sub_folder = "PD plate layout")
}





# Save data ---------------------------------------------------------------

save(list = c("full_4sg_by_gene_df", "full_4sg_by_well_df",
              "sg4_by_gene_df", "sg4_by_well_df"
              ),
     file = file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData")
     )

save(list = c("PD_4sg_df", "PD_4sg_reordered_df"),
     file = file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates - PD genes.RData")
     )








