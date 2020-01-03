### 30 October 2019 ###






# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))






# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")






# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# on 21 July 2019
NCBI_Hs_info_df <- read.table(file.path(CRISPR_input_directory, "Human genome", "NCBI", "Homo_sapiens.gene_info"),
                              sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                              fill = TRUE, check.names = FALSE, comment.char = ""
                              )



# Downloaded from ftp://ftp.ensembl.org/pub/release-85/gtf/homo_sapiens/Homo_sapiens.GRCh38.85.gtf.gz
# on 14 October 2019
Ensembl_Hs_entrez_df <- read.table(file.path(CRISPR_input_directory, "Human genome", "ENSEMBL", "Homo_sapiens.GRCh38.98.entrez.tsv"),
                                   sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                                   check.names = FALSE
                                   )





# Collect Entrez IDs from various sources ---------------------------------

Ensembl_Hs_have_entrez <- !(is.na(suppressWarnings(as.integer(Ensembl_Hs_entrez_df[, "xref"]))))
Ensembl_Hs_are_protein_coding <- Ensembl_Hs_entrez_df[, "protein_stable_id"] != "-"

collected_entrez_IDs <- unique(c(as.character(NCBI_Hs_info_df[NCBI_Hs_info_df[, "type_of_gene"] == "protein-coding", "GeneID"]),
                                 Ensembl_Hs_entrez_df[Ensembl_Hs_have_entrez & Ensembl_Hs_are_protein_coding, "xref"]
                                 )
                               )





# Save data ---------------------------------------------------------------

save(list = "collected_entrez_IDs",
     file = file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData")
     )



