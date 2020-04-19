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
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))





# Assign the sgRNAs to plates ---------------------------------------------

sg4_df <- AssignAllGuides(merged_replaced_CRISPRa_df, sublibraries_all_entrezs_list, reorder_df = FALSE)
sg4_reordered_df <- ReorderPlates(sg4_df)






# Export the plate layouts ------------------------------------------------

ExportPlates(sg4_df, "All_sublibraries_original_order")
ExportPlates(sg4_reordered_df, "All_sublibraries_reordered")



# sub_df <- merged_replaced_CRISPRa_df[(merged_replaced_CRISPRa_df[["Gene_symbol"]] %in% "NAT1"), ]
#
# sub_df <- sub_df[sub_df[["hCRISPRa_v2_rank"]] %in% as.character(1:4), ]
#
# sub_df <- sub_df[sub_df[["hCRISPRa_v2_transcript"]] %in% "P1", ]
#
# LongestSharedSubsequence(sub_df[["sgRNA_sequence"]])


length(unique(Calabrese_df[["Combined_ID"]][sg4_df[["Is_control"]] == "No"]))

