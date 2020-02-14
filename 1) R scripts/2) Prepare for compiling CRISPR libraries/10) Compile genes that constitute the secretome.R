### 13 February 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")
secretome_directory     <- file.path(CRISPR_input_directory, "Sublibraries", "Uhlen - Science Sig 2019")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))




# Read in data ------------------------------------------------------------

secretome_df <- as.data.frame(read_excel(file.path(secretome_directory, "2019 - The human secretome - Data S2.xlsx")),
                              stringsAsFactors = FALSE, check.names = FALSE
                              )



# Modify the data frame ---------------------------------------------------

colnames(secretome_df) <- gsub(" ", "_", colnames(secretome_df), fixed = TRUE)
colnames(secretome_df)[colnames(secretome_df) == "Ensembl_gene_id"] <- "Ensembl_gene_ID"
colnames(secretome_df)[[2]] <- "Gene_symbol"

secretome_df[["UniProt_accession"]] <- ifelse(secretome_df[["UniProt_accession"]] == "NA",
                                              NA_character_,
                                              secretome_df[["UniProt_accession"]]
                                              )




# Map secretome genes to Entrez IDs ---------------------------------------

ensembl_mappings_df <- MapEnsemblIDs(secretome_df)






# Examine secretome genes with problematic mappings -----------------------

View(ensembl_mappings_df[ensembl_mappings_df[["Are_problematic"]], ])






# Select columns for the secretome data frame -----------------------------

secretome_df <- data.frame(
  ensembl_mappings_df["Combined_ID"],
  "Entrez_ID"      = ensembl_mappings_df[["Consensus_entrez"]],
  "Gene_symbol"    = ensembl_mappings_df[["Consensus_symbol"]],
  ensembl_mappings_df["Original_symbol"],
  secretome_df[, c("Ensembl_gene_ID", "UniProt_accession", "Annotated_category")],
  "Supercategory" = ifelse(secretome_df[, "Annotated_category"] %in% c("Intracellular or membrane-bound", "Secreted to blood"),
                           secretome_df[, "Annotated_category"],
                           "Locally secreted"
                           ),
  stringsAsFactors = FALSE,
  row.names        = NULL
)




# Save data ---------------------------------------------------------------

save(list = "secretome_df",
     file = file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData")
     )




















