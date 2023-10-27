## 2022-06-21


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts", "1) R functions")
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))



# Define folder paths -----------------------------------------------------

experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir  <- file.path(experiments_directory, "2022-06-21 - Illumina paired-end 2sg - correct reference")
rdata_dir    <- file.path(project_dir, "03_R_objects")
input_dir    <- file.path(project_dir, "02_input_data")
library_path <- file.path(input_dir, "20200513_library_1_2_unbalanced.csv")




# Read in the CRISPR-off library ------------------------------------------

CRISPRoff_df <- read.csv(library_path, stringsAsFactors = FALSE)[, -1]




# Process CRISPRoff 2sg library genes -------------------------------------

symbols_vec <- sapply(strsplit(CRISPRoff_df[, "sgID_AB"], "_", fixed = TRUE), "[[", 1)
all_entrezs_df <- SymbolToEntrezDf(symbols_vec)

CRISPRoff_df[, "Entrez_IDs"] <- all_entrezs_df[, "Entrez_ID"]
CRISPRoff_df[, "Gene_symbols"] <- all_entrezs_df[, "Gene_symbol"]
CRISPRoff_df[, "Combined_ID"] <- ifelse(is.na(CRISPRoff_df[, "Entrez_IDs"]),
                                        CRISPRoff_df[, "gene"],
                                        CRISPRoff_df[, "Entrez_IDs"]
                                        )




# Process plasmid IDs -----------------------------------------------------

are_NT <- CRISPRoff_df[, "gene"] == "negative_control"
NT_IDs <- paste0("NT_", gsub("non-targeting_", "", CRISPRoff_df[, "sgID_AB"][are_NT], fixed = TRUE))

num_occurrences <- table(symbols_vec)[symbols_vec]

plasmid_IDs_vec <- ifelse(are_NT,
                          NT_IDs,
                          ifelse(num_occurrences == 1,
                                 symbols_vec,
                                 paste0(symbols_vec, "_", CRISPRoff_df[, "transcript"])
                                 )
                          )

num_plasmid_occurrences <- table(plasmid_IDs_vec)[plasmid_IDs_vec]
plasmid_IDs_vec <- ifelse(num_plasmid_occurrences == 2,
                          CRISPRoff_df[, "gene"],
                          plasmid_IDs_vec
                          )

stopifnot(!(anyDuplicated(plasmid_IDs_vec)))

CRISPRoff_df[, "Plasmid_ID"] <- gsub(",", "/", plasmid_IDs_vec, fixed = TRUE)




# Explore the annotated CRISPRoff library ---------------------------------

CRISPRoff_df[!(are_NT) & is.na(CRISPRoff_df[, "Entrez_IDs"]), ]
CRISPRoff_df[grepl("," , CRISPRoff_df[, "Entrez_IDs"], fixed = TRUE), ]




# Perform a manual fix in preparation for disambiguating Entrez IDs -------

# The gene symbol for MEMO1 translates to the wrong Entrez ID.
# Correcting this avoids an error resulting from a discrepancy between the
# number of genes targeted analyzing gene symbols vs. Entrez IDs.

are_to_replace <- (CRISPRoff_df[, "Entrez_IDs"] %in% "7795") &
                  (CRISPRoff_df[, "Gene_symbols"] %in% "MEMO1")
CRISPRoff_df[, "Entrez_IDs"][are_to_replace] <- "51072"




# Save data ---------------------------------------------------------------

save(list = "CRISPRoff_df",
     file = file.path(rdata_dir, "02_annotate_dJR072_CRISPRoff_library.RData")
     )

