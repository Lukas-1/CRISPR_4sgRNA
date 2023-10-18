### 25th July 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "05) Compile TSS (transcription start site) data.RData"))
load(file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from mCRISPRa-v2 - CRISPRa_df.RData"))





# Find all unique identifiers in CRISPRa_df -------------------------------

unique_CRISPRa_df <- CRISPRa_df[!(duplicated(CRISPRa_df[, c("Entrez_ID", "Original_symbol")])), ]
row.names(unique_CRISPRa_df) <- NULL





# Map genes to combined TSS summaries -------------------------------------

combined_TSS_CRISPRa_df <- FilterAndCombineEntrezSymbol(combined_TSS_summary_df,
                                                        unique_CRISPRa_df[["Entrez_ID"]],
                                                        unique_CRISPRa_df[["Gene_symbol"]],
                                                        unique_CRISPRa_df[["Combined_ID"]]
                                                        )





# Save data ---------------------------------------------------------------

save(list = "combined_TSS_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData")
     )




