### 25th July 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRa-v2 - CRISPRa_df.RData"))





# Find all unique identifiers in CRISPRa_df -------------------------------

unique_CRISPRa_df <- CRISPRa_df[!(duplicated(CRISPRa_df[, c("Entrez_ID", "Original_symbol")])), ]
row.names(unique_CRISPRa_df) <- NULL





# # Map genes to FANTOM5 TSS summaries --------------------------------------
#
# FANTOM5_CRISPRa_df <- FilterAndCombineEntrezSymbol(FANTOM5_summary_df,
#                                                    unique_CRISPRa_df[["Entrez_ID"]],
#                                                    unique_CRISPRa_df[["Gene_symbol"]],
#                                                    unique_CRISPRa_df[["Combined_ID"]]
#                                                    )
#
#
#
#
# # Map genes to BioMart TSS summaries --------------------------------------
#
# BioMart_CRISPRa_df <- FilterAndCombineEntrezSymbol(BioMart_summary_df,
#                                                    unique_CRISPRa_df[["Entrez_ID"]],
#                                                    unique_CRISPRa_df[["Gene_symbol"]],
#                                                    unique_CRISPRa_df[["Combined_ID"]]
#                                                    )
#



# Map genes to combined TSS summaries -------------------------------------

combined_TSS_CRISPRa_df <- FilterAndCombineEntrezSymbol(combined_TSS_summary_df,
                                                        unique_CRISPRa_df[["Entrez_ID"]],
                                                        unique_CRISPRa_df[["Gene_symbol"]],
                                                        unique_CRISPRa_df[["Combined_ID"]]
                                                        )





# Save data ---------------------------------------------------------------

save(list = c("combined_TSS_CRISPRa_df"), #, "FANTOM5_CRISPRa_df", "BioMart_CRISPRa_df"),
     file = file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData")
     )




