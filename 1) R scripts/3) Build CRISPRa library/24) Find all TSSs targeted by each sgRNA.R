### 18th October 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))




# Find all TSSs targeted by each sgRNA ------------------------------------

# The gene symbol for MEMO1 translates to the wrong Entrez ID.
# Correcting this avoids an error resulting from a discrepancy between the
# number of genes targeted analyzing gene symbols vs. Entrez IDs.

are_to_replace <- (merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% "7795") &
                  (merged_replaced_CRISPRa_df[["Gene_symbol"]] %in% "MEMO1")
merged_replaced_CRISPRa_df[["Entrez_ID"]][are_to_replace] <- "51072"

nearby_list <- AlignSummaryDf(FindNearbyTSSs,
                              merged_replaced_CRISPRa_df,
                              all_TSS_df
                              )

nearby_protein_list <- AlignSummaryDf(FindNearbyTSSs,
                                      merged_replaced_CRISPRa_df,
                                      all_TSS_df,
                                      only_protein_coding = TRUE
                                      )

nearby_main_TSS_list <- AlignSummaryDf(FindNearbyTSSs,
                                       merged_replaced_CRISPRa_df,
                                       all_TSS_df,
                                       only_best_TSS = TRUE
                                       )

nearby_main_protein_list <- AlignSummaryDf(FindNearbyTSSs,
                                           merged_replaced_CRISPRa_df,
                                           all_TSS_df,
                                           only_protein_coding = TRUE,
                                           only_best_TSS = TRUE
                                           )




# Prepare for saving data frames ------------------------------------------

TSS_targets_df              <- nearby_list[["summary_df"]]
TSS_protein_targets_df      <- nearby_protein_list[["summary_df"]]
main_TSS_targets_df         <- nearby_main_TSS_list[["summary_df"]]
main_TSS_protein_targets_df <- nearby_main_protein_list[["summary_df"]]


TSS_targets_full_df              <- nearby_list[["full_df"]]
TSS_protein_targets_full_df      <- nearby_protein_list[["full_df"]]
main_TSS_targets_full_df         <- nearby_main_TSS_list[["full_df"]]
main_TSS_protein_targets_full_df <- nearby_main_protein_list[["full_df"]]





# Save data ---------------------------------------------------------------

save(list = c("TSS_targets_df", "TSS_protein_targets_df", "main_TSS_targets_df", "main_TSS_protein_targets_df"),
     file = file.path(CRISPRa_RData_directory, "24) Find all TSSs targeted by each sgRNA - summary data.RData")
     )
save(list = c("TSS_targets_full_df", "TSS_protein_targets_full_df", "main_TSS_targets_full_df", "main_TSS_protein_targets_full_df"),
     file = file.path(CRISPRa_RData_directory, "24) Find all TSSs targeted by each sgRNA - full data.RData")
     )



