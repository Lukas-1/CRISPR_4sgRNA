### 7th July 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "27) Creating 4sg plasmid sequences.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")
plasmid_output_directory <- file.path(file_output_directory, "Plasmids")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))





# Export vectors ----------------------------------------------------------

export_genes <- c(
  "TXNIP", "KIF1A", "C19orf48"
)


export_genes <- c(
  "ADGRG6", "ADGRD1", "ADGRG1", "Control_TF1"
)

for (gene_symbol in export_genes) {
  ExportVectorsForGene(gene_symbol, full_4sg_by_gene_df)
}



















