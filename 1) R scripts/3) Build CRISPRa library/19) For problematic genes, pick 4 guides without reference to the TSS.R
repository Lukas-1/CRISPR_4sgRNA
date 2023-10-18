### 1st March 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))
source(file.path(general_functions_directory, "25) Choosing guides without reference to the TSS for problematic genes.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRa_RData_directory, "18) Pick 4 guides per TSS.RData"))





# Pick a new set of 4 guides for genes with invalid 4sgs ------------------

merged_replaced_CRISPRa_df <- ReplaceUnspacedGuides(merged_replaced_CRISPRa_df)





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "19) For problematic genes, pick 4 guides without reference to the TSS.RData"
                      )
     )


