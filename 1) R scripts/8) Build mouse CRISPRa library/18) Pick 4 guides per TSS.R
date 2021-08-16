### 10th October 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData"))





# Create an empty SNP column ----------------------------------------------

merged_replaced_CRISPRa_df[[preferred_AF_max_column]] <- NA_real_



# Rank sgRNAs -------------------------------------------------------------

merged_replaced_CRISPRa_df <- RankCRISPRDfByTSS(merged_replaced_CRISPRa_df)





# Find combinations of non-overlapping sgRNAs -----------------------------

merged_replaced_CRISPRa_df <- PrioritizeNonOverlapping(merged_replaced_CRISPRa_df,
                                                       ID_column = "AltTSS_ID",
                                                       parallel_mode = FALSE,
                                                       tolerate_divergent_chromosomes = TRUE
                                                       )



# Remove the empty SNP column ---------------------------------------------

merged_replaced_CRISPRa_df <- merged_replaced_CRISPRa_df[, names(merged_replaced_CRISPRa_df) != preferred_AF_max_column]




# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "18) Pick 4 guides per TSS.RData"
                      )
     )




