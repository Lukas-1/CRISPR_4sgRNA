### 31st March 2020 ###




# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(general_RData_directory, "17) Compile the information on gene type.RData"))






# Do stuff ----------------------------------------------------------------

head(combined_TSS_df)
head(FANTOM5_df)
head(BioMart_tidied_df)

grep(",", FANTOM5_df[["Entrez_ID"]], value = TRUE)




















