### 5th November 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "10) Find overlaps between CRISPRko sgRNA sequences and genetic polymorphisms.RData"))






# Rank sgRNAs -------------------------------------------------------------

merged_CRISPRko_df <- RankCRISPRDf(merged_CRISPRko_df, ID_column = "Combined_ID")





# Find combinations of non-overlapping sgRNAs -----------------------------

merged_CRISPRko_df <- PrioritizeNonOverlapping(merged_CRISPRko_df, ID_column = "Combined_ID", parallel_mode = TRUE)






# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory,
                      "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"
                      )
     )



