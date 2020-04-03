### 31st March 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))




# Read in data ------------------------------------------------------------

CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
FANTOM5_input_directory <- file.path(CRISPR_input_directory, "Human genome", "FANTOM5_liftover")
Ensembl_input_directory <- file.path(CRISPR_input_directory, "Human genome", "Ensembl")
HGNC_input_directory    <- file.path(CRISPR_input_directory, "Human genome", "HGNC")



FANTOM5_ann_df <- read.table(file.path(FANTOM5_input_directory, "hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, fill = TRUE, check.names = FALSE, comment.char = ""
                             )

FANTOM5_bed_df <- read.table(file.path(FANTOM5_input_directory, "hg38_liftover_CAGE_peaks_phase1and2.bed"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                             )

# The BioMart file was downloaded from https://www.ensembl.org/biomart/martview
BioMart_df     <- read.table(file.path(CRISPR_input_directory, "Human genome", "Ensembl", "BioMart_human_2020-03-25_mart_export.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                             )

HGNC_df <- read.table(file.path(HGNC_input_directory, "Custom_HGNC_download_including_NCBI_IDs__2020_03_11.tsv"),
                      sep = "\t", quote = "", fill = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE, header = TRUE, row.names = NULL
                      )



# Modify the data frames --------------------------------------------------

HGNC_df[["Entrez_ID"]] <- ifelse(is.na(HGNC_df[["NCBI Gene ID"]]),
                                 HGNC_df[["NCBI Gene ID(supplied by NCBI)"]],
                                 HGNC_df[["NCBI Gene ID"]]
                                 )

HGNC_matches <- match(FANTOM5_ann_df[["HGNC/MGI_ID"]], HGNC_df[["HGNC ID"]])

FANTOM5_ann_df[["HGNC_entrez"]] <- as.character(HGNC_df[["Entrez_ID"]][HGNC_matches])

FANTOM5_ann_df[["FANTOM5_entrez"]] <- ifelse(FANTOM5_ann_df[["GeneID"]] == "",
                                        NA_character_,
                                        FANTOM5_ann_df[["GeneID"]]
                                        )

FANTOM5_ann_df[["FANTOM5_entrez"]] <- gsub(" ", ", ", FANTOM5_ann_df[["FANTOM5_entrez"]], fixed = TRUE)

FANTOM5_ann_df[["Consensus_entrez"]] <- ifelse(is.na(FANTOM5_ann_df[["FANTOM5_entrez"]]),
                                               FANTOM5_ann_df[["HGNC_entrez"]],
                                               ifelse(grepl(",", FANTOM5_ann_df[["FANTOM5_entrez"]], fixed = TRUE),
                                                      ifelse(is.na(FANTOM5_ann_df[["HGNC_entrez"]]),
                                                             FANTOM5_ann_df[["FANTOM5_entrez"]],
                                                             FANTOM5_ann_df[["HGNC_entrez"]]
                                                             ),
                                                      FANTOM5_ann_df[["HGNC_entrez"]]
                                                      )
                                               )



# Check for diverging annotations -----------------------------------------

FANTOM5_ann_columns <- c(
  "Distance", "FANTOM5_entrez", "HGNC_entrez", "Consensus_entrez",
  "Gene_symbol", "Gene_synonyms", "Gene_name",
  "HGNC/MGI_ID", "Gene_source", "UniProt_ID"
)

are_identical <- mapply(identical, FANTOM5_ann_df[["FANTOM5_entrez"]], FANTOM5_ann_df[["HGNC_entrez"]])

View(FANTOM5_ann_df[!(are_identical), FANTOM5_ann_columns])


table(grepl(",", FANTOM5_ann_df[["Consensus_entrez"]], fixed = TRUE))
table(grepl(",", FANTOM5_ann_df[["FANTOM5_entrez"]], fixed = TRUE))

View(FANTOM5_ann_df[grepl(",", FANTOM5_ann_df[["Consensus_entrez"]], fixed = TRUE), FANTOM5_ann_columns])







# View the data frames ----------------------------------------------------



View(FANTOM5_ann_df[, FANTOM5_ann_columns])





# Try stuff ---------------------------------------------------------------

head(FANTOM5_ann_df)
head(FANTOM5_ann_df)
head(BioMart_df)




are_identical <- mapply(identical, FANTOM5_ann_df[["HGNC/MGI_ID"]], FANTOM5_ann_df[["Gene_source"]])




View(FANTOM5_ann_df[!(are_identical), FANTOM5_ann_columns])




View(FANTOM5_ann_df[, FANTOM5_ann_columns])



FANTOM5_ann_df[["HGNC_Entrez_ID"]] <- HGNC_df[["NCBI Gene ID"]]



are_identical <- mapply(identical, HGNC_df[["NCBI Gene ID"]], HGNC_df[["NCBI Gene ID(supplied by NCBI)"]])
View(HGNC_df[!(are_identical) & !(is.na(HGNC_df[["NCBI Gene ID"]])), HGNC_columns])

all(is.na(HGNC_df[["NCBI Gene ID"]][!(are_identical)]))














