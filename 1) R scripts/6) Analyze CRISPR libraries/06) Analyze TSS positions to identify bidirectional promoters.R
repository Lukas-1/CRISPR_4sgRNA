### 9th March 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the human genome.RData"))






# Identify bidirectional promoters ----------------------------------------

all_genes_TSS_df <- PrepareTSSDf(all_TSS_df)

CRISPRa_genes <- sgRNAs_overview_df[["Entrez_ID"]][sgRNAs_overview_df[["In_4sg_library"]] %in% "Yes"]
CRISPRa_genes_TSS_df <- all_genes_TSS_df[all_genes_TSS_df[["Entrez_ID"]] %in% CRISPRa_genes, ]
row.names(CRISPRa_genes_TSS_df) <- NULL











