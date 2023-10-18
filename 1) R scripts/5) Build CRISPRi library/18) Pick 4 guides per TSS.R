### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "17) Integrate the output from CRISPOR.RData"))





# Rank sgRNAs -------------------------------------------------------------

merged_replaced_CRISPRi_df <- RankCRISPRDfByTSS(merged_replaced_CRISPRi_df)





# Find combinations of non-overlapping sgRNAs -----------------------------

merged_replaced_CRISPRi_df <- PrioritizeNonOverlapping(merged_replaced_CRISPRi_df,
                                                       ID_column = "AltTSS_ID",
                                                       parallel_mode = TRUE,
                                                       tolerate_divergent_chromosomes = TRUE
                                                       )



# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRi_df",
     file = file.path(CRISPRi_RData_directory,
                      "18) Pick 4 guides per TSS.RData"
                      )
     )




