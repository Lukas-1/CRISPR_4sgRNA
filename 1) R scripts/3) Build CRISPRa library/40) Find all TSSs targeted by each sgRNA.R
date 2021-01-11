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





# Define functions --------------------------------------------------------




# Do stuff ----------------------------------------------------------------

# example_entrezs <- unique(merged_replaced_CRISPRa_df[["Entrez_ID"]])[1:1000]
example_entrezs <- unlist(sublibraries_all_entrezs_list, use.names = FALSE)
are_example_genes <- merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% example_entrezs


summary_df <- FindNearbyTSSs(merged_replaced_CRISPRa_df[are_example_genes, ],
                             all_TSS_df
                             )




