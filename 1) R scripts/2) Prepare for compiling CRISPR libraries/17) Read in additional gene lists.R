### 19th April 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
general_RData_directory <- file.path(RData_directory, "1) General")
gene_lists_directory    <- file.path(CRISPR_root_directory, "2) Input data", "Gene lists")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))





# Read in data ------------------------------------------------------------

PD_path <- file.path(gene_lists_directory, "PD_collaborators", "PD-GWAS list.xlsx")
PD_candidates_df <- data.frame(read_excel(PD_path), stringsAsFactors = FALSE)




# Collect Entrez IDs from our collaborators' PD GWAS study ----------------

## Correct the input gene symbols
PD_symbols_vec <- PD_candidates_df[[1]]
PD_symbols_vec[[42]] <- "CHCHD2"
PD_symbols_vec[[which(PD_symbols_vec == "PARK7/DJ1")]] <- "PARK7"
PD_symbols_vec[[which(PD_symbols_vec == "PLEKHMI")]]   <- "PLEKHM1"


## Map to Entrez IDs
PD_mapped_df <- MapToEntrezs(symbols_vec = PD_symbols_vec)


## Manually assign ambiguous Entrez IDs
are_ambiguous <- grepl(",", PD_mapped_df[["Entrez_ID"]], fixed = TRUE)
PD_mapped_df[are_ambiguous, ]
PD_mapped_df[["Entrez_ID"]][[which(are_ambiguous & (PD_mapped_df[["Original_symbol"]] %in% "QARS"))]] <- "5859"
PD_entrezs_vec <- PD_mapped_df[["Entrez_ID"]]

## Identify genes that are unavailable in our library
are_missing <- is.na(PD_entrezs_vec)
are_protein_coding <- PD_entrezs_vec %in% collected_entrez_IDs
are_included <- PD_entrezs_vec %in% unlist(sublibraries_all_entrezs_list)


## Examine missing genes
identical(are_protein_coding, are_included)
table(are_missing)
table(are_included)


## Examine duplicated genes
num_occurrences <- table(PD_entrezs_vec)[PD_entrezs_vec]
are_duplicated <- (num_occurrences >= 2) %in% TRUE
PD_mapped_df[are_duplicated, ]


## Define gene lists
PD_all_entrezs <- PD_entrezs_vec[!(are_missing)]
PD_4sg_entrezs <- PD_entrezs_vec[are_included]




# Save data ---------------------------------------------------------------

save(list = c("PD_all_entrezs", "PD_4sg_entrezs", "PD_mapped_df"),
     file = file.path(general_RData_directory, "17) Read in additional gene lists.RData")
     )











