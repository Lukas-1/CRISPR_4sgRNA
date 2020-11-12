### 27th June 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "27) Creating 4sg plasmid sequences.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory  <- file.path(RData_directory, "4) CRISPRi")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRi")
plasmid_output_directory <- file.path(file_output_directory, "Plasmids")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "27) Distribute sgRNAs for the whole genome onto plates.RData"))





# Export vectors ----------------------------------------------------------

export_genes <- c(
  "RGS18", "LSM4", "PML", "MED31", "ZNF441",
  "SERPINB2", "FAAH2", "LSM2", "SMN1", "PIKFYVE"
)

export_genes <- c(
  "TXNIP", "KIF1A", "C19orf48"
)


export_genes <- c(
  c(paste0("LSM", c(1, 3, 5, 6, 7, 8)), "PRPF31")
)

export_genes %in% sg4_df[["Gene_symbol"]]

for (gene_symbol in export_genes) {
  ExportVectorsForGene(gene_symbol, sg4_df)
}

ExportVectorsForGene("Control_1", vac_4sg_reordered_df)
ExportVectorsForGene("Control_2", vac_4sg_reordered_df)
ExportVectorsForGene("Control_3", vac_4sg_reordered_df)

















