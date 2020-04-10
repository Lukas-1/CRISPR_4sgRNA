### 8th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(CRISPRi_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRi-v2 - CRISPRi_df.RData"))





# Find all unique identifiers in CRISPRi_df -------------------------------

unique_CRISPRi_df <- CRISPRi_df[!(duplicated(CRISPRi_df[, c("Entrez_ID", "Original_symbol")])), ]
row.names(unique_CRISPRi_df) <- NULL




# Map genes to combined TSS summaries -------------------------------------

combined_TSS_CRISPRi_df <- FilterAndCombineEntrezSymbol(combined_TSS_summary_df,
                                                        unique_CRISPRi_df[["Entrez_ID"]],
                                                        unique_CRISPRi_df[["Gene_symbol"]],
                                                        unique_CRISPRi_df[["Combined_ID"]]
                                                        )



# Save data ---------------------------------------------------------------

save("combined_TSS_CRISPRi_df",
     file = file.path(CRISPRi_RData_directory, "03) Map CRISPR libraries to TSS data.RData")
     )




