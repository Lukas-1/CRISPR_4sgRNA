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
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory     <- file.path(RData_directory, "2) CRISPRa")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")
previous_versions_directory <- file.path(RData_directory, "5) Previous versions of the library")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "18) Read in additional gene lists.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(previous_versions_directory, "01) CRISPRa transcription factor sub-library (1st version) - TF_v1_CRISPRa_df.RData"))







# Assign the sgRNAs to plates ---------------------------------------------

sg4_reordered_df <- AllocateAllGuides_v2(merged_replaced_CRISPRa_df,
                                         sublibraries_entrezs_list  = sublibraries_all_entrezs_list,
                                         previous_version_CRISPR_df = TF_v1_CRISPRa_df,
                                         candidate_entrezs          = PD_4sg_entrezs
                                         )
sg4_df <- RestoreOriginalOrder(sg4_reordered_df)




# Examine the sub-library allocation --------------------------------------

sg4_allTFs_df <- AllocateAllGuidesToPlates(merged_replaced_CRISPRa_df,
                                           sublibraries_entrezs_list = sublibraries_all_entrezs_list,
                                           num_control_wells = 96
                                           )

are_selected <- (sg4_allTFs_df[["Is_control"]] == "No") &
                 !(duplicated(sg4_allTFs_df[["Combined_ID"]]))

table(sg4_allTFs_df[["Sublibrary_4sg"]][are_selected])

are_selected_TSSs <- (sg4_allTFs_df[["Is_control"]] == "No") &
                     !(duplicated(sg4_allTFs_df[["AltTSS_ID"]]))

sum(are_selected_TSSs) # Number of TSSs targeted by the library



# Export the plate layouts ------------------------------------------------

ExportPlates(sg4_df, "All_sublibraries_original_order")
ExportPlates(sg4_reordered_df, "All_sublibraries_reordered")

for (i in 1:4) {
  use_df <- sg4_reordered_df[sg4_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df, paste0("4sg_reordered_sg", i), add_padding_between_plates = TRUE)
}




# Save data ---------------------------------------------------------------

save(list = c("sg4_df", "sg4_reordered_df"),
     file = file.path(CRISPRa_RData_directory, "27) Distribute sgRNAs for the whole genome onto plates.RData")
     )









