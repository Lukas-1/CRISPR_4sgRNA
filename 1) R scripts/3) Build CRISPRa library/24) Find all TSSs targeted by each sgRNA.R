### 18th October 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))




# Find all TSSs targeted by each sgRNA ------------------------------------

have_single_entrez <- !(is.na(merged_replaced_CRISPRa_df[["Entrez_ID"]])) &
                      !(grepl(",", merged_replaced_CRISPRa_df[["Entrez_ID"]], fixed = TRUE)) &
                      (!(merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% "7795")) # The gene symbol for MEMO1 translates to the wrong Entrez ID

indices_vec <- rep(NA, nrow(merged_replaced_CRISPRa_df))
indices_vec[have_single_entrez] <- seq_len(sum(have_single_entrez))

nearby_list <- FindNearbyTSSs(merged_replaced_CRISPRa_df[have_single_entrez, ],
                              all_TSS_df
                              )

TSS_targets_df <- nearby_list[["summary_df"]]
TSS_targets_df <- TSS_targets_df[indices_vec, ]
row.names(TSS_targets_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = "TSS_targets_df",
     file = file.path(CRISPRa_RData_directory, "24) Find all TSSs targeted by each sgRNA.RData")
     )





