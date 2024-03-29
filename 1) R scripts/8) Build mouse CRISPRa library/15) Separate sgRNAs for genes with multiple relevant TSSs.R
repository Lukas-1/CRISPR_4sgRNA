### 9th October 2020 ###





# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "14) Find the nearest genes to each sgRNA.RData"))





# Find combinations of non-overlapping sgRNAs -----------------------------

merged_replaced_CRISPRa_df <- AllocateTSSforAllGenes(merged_replaced_CRISPRa_df,
                                                     parallel_mode = TRUE,
                                                     tolerate_divergent_chromosomes = TRUE
                                                     )





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData")
     )







