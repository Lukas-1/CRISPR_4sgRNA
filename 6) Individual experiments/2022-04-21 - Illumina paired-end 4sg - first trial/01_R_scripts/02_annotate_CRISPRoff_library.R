## 2022-04-21


# Load packages and source code -------------------------------------------

library("readxl")

CRISPR_root_directory <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts", "1) R functions")
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))



# Define folder paths -----------------------------------------------------

experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir  <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 4sg - first trial")
rdata_dir    <- file.path(project_dir, "03_R_objects")
input_dir    <- file.path(project_dir, "02_input_data")
library_path <- file.path(input_dir, "2021 - Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing - Table S3.xlsx")





# Read in the CRISPR-off library ------------------------------------------

CRISPRoff_df <- data.frame(read_excel(library_path, skip = 3), stringsAsFactors = FALSE)




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

same_transcripts <- mapply(identical, CRISPRoff_df[, "transcript_A"], CRISPRoff_df[, "transcript_B"])

transcripts_vec <- mapply(function(x, y) paste0(unique(c(x, y)), collapse = "/"),
                          CRISPRoff_df[, "transcript_A"],
                          CRISPRoff_df[, "transcript_B"]
                          )

plasmid_IDs_vec <- ifelse(are_NT,
                          NT_IDs,
                          ifelse(num_occurrences == 1,
                                 symbols_vec,
                                 paste0(symbols_vec, "_", transcripts_vec)
                                 )
                          )

num_plasmid_occurrences <- table(plasmid_IDs_vec)[plasmid_IDs_vec]
plasmid_IDs_vec <- ifelse(num_plasmid_occurrences == 2,
                          CRISPRoff_df[, "gene"],
                          plasmid_IDs_vec
                          )

stopifnot(!(any(duplicated(plasmid_IDs_vec))))

CRISPRoff_df[, "Plasmid_ID"] <- gsub(",", "/", plasmid_IDs_vec, fixed = TRUE)




# Explore the annotated CRISPRoff library ---------------------------------

CRISPRoff_df[!(are_NT) & is.na(CRISPRoff_df[, "Entrez_IDs"]), ]
CRISPRoff_df[grepl("," , CRISPRoff_df[, "Entrez_IDs"], fixed = TRUE), ]




# Save data ---------------------------------------------------------------

save(list = "CRISPRoff_df",
     file = file.path(rdata_dir, "02_annotate_CRISPRoff_library.RData")
     )

