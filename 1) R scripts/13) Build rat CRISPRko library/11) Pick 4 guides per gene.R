### 21st November 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "12) Rat - CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "09) Integrate the output from CRISPOR.RData"))






# Create an empty SNP column ----------------------------------------------

merged_CRISPRko_df[[preferred_AF_max_column]] <- NA_real_




# Rank sgRNAs -------------------------------------------------------------

merged_CRISPRko_df <- RankCRISPRDf(merged_CRISPRko_df, ID_column = "Combined_ID")





# Find combinations of non-overlapping sgRNAs -----------------------------

merged_CRISPRko_df <- PrioritizeNonOverlapping(merged_CRISPRko_df,
                                               ID_column = "Combined_ID",
                                               parallel_mode = FALSE,
                                               tolerate_divergent_chromosomes = FALSE
                                               )




# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory,
                      "11) Pick 4 guides per gene.RData"
                      )
     )


