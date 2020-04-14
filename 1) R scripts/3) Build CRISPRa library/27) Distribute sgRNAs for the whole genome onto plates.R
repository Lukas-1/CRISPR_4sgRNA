### 11th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # for FormatFixedWidthInteger
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))





# Do everything -----------------------------------------------------------

all_df <- AssignAllGuides(merged_replaced_CRISPRa_df, sublibraries_all_entrezs_list, reorder_df = TRUE)


unique_well_IDs <- paste0(all_df[["Plate_ID"]], "__", all_df[["Well_number"]])

split_seqs <- split(toupper(all_df[["sgRNA_sequence"]]), unique_well_IDs)

num_sequences_goo <- vapply(split_seqs, LongestSharedSubsequence, integer(1))

table(all_df[["Plate_ID"]])

goo
num_occurrences_vec <- table(unique_well_IDs)[unique_well_IDs]


# Pick control guides -----------------------------------------------------

# merged_replaced_CRISPRa_df <- AddRandomized4sgControls(merged_replaced_CRISPRa_df)
#
# table(merged_replaced_CRISPRa_df[["Source"]][!(is.na(merged_replaced_CRISPRa_df[["Control_group_4sg"]]))])
#
# merged_replaced_CRISPRa_df <- AssignControlsToPlates(merged_replaced_CRISPRa_df)
#
# merged_replaced_CRISPRa_df[["Plate_ID"]][!(is.na(merged_replaced_CRISPRa_df[["Control_group_4sg"]]))]
#
# goo
#
#
#
# # Exclude most genes that are not protein-coding --------------------------
#
# all_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% unlist(sublibraries_all_entrezs_list, use.names = FALSE), ]
#
# all_CRISPRa_df <- AssignControlsToPlates(all_CRISPRa_df)
#
#
#
#
# # Include only valid 4sg combinations -------------------------------------
#
# are_top4_mat <- CRISPRaAreTop4Mat(all_CRISPRa_df)
#
# ShowProblematicGuides(all_CRISPRa_df, are_top4_mat)
#
# are_valid_chosen <- are_top4_mat[, "Are_chosen_4sg"] & are_top4_mat[, "Have_valid_guides"]
#
# all_CRISPRa_df <- all_CRISPRa_df[are_valid_chosen, ]
# rownames(all_CRISPRa_df) <- NULL
#
#
#
#
#
# # Allocate guides to plates -----------------------------------------------
#
# all_CRISPRa_df <- AddSublibrary(all_CRISPRa_df, sublibraries_all_entrezs_list)
# all_CRISPRa_df <- AssignToPlates(all_CRISPRa_df)
#
#
#
#
#
#











