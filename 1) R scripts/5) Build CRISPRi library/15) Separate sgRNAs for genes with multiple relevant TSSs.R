### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "14) Find overlaps between sgRNA sequences and genetic polymorphisms.RData"))





# Find combinations of non-overlapping sgRNAs -----------------------------

merged_replaced_CRISPRi_df <- AllocateTSSforAllGenes(merged_replaced_CRISPRi_df,
                                                     parallel_mode = TRUE,
                                                     tolerate_divergent_chromosomes = TRUE
                                                     )





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRi_df",
     file = file.path(CRISPRi_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData")
     )







