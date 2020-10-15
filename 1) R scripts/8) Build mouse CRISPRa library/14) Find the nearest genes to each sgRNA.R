### 12th October 2020 ###



# Import packages and source code -----------------------------------------

library("TxDb.Mmusculus.UCSC.mm10.knownGene")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRa_RData_directory, "13) Integrate the output from GuideScan for individual sgRNA locations.RData"))




# Load global variables ---------------------------------------------------

mouse_genes_mm10_GRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
rm(human_genes_GRanges)





# Identify sgRNAs that have a defined location ----------------------------

are_mapped <- !(is.na(merged_replaced_CRISPRa_df[["Start"]]))
are_mitochondrial <- merged_replaced_CRISPRa_df[["Chromosome"]] %in% "chrM" # The TxDb.Mmusculus.UCSC.mm10.knownGene object does not contain mitochondrial genes

are_eligible <- are_mapped & !(are_mitochondrial)

mapped_indices <- rep(NA_integer_, length(are_eligible))
mapped_indices[are_eligible] <- seq_len(sum(are_eligible))





# Find the nearest genes (0MM locations only) -----------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End")
nearest_columns <- c("Nearest_Entrez_IDs", "Nearest_symbols", "Distance")


nearest_genes_df <- FindNearestGenes(merged_replaced_CRISPRa_df[are_eligible, location_columns],
                                     mouse_genes_mm10_GRanges,
                                     is_mouse = TRUE
                                     )[, nearest_columns]
colnames(nearest_genes_df)[[3]] <- "Nearest_gene_distance"





# Merge the data frame ----------------------------------------------------

merged_replaced_CRISPRa_df <- data.frame(merged_replaced_CRISPRa_df,
                                         nearest_genes_df[mapped_indices, ],
                                         stringsAsFactors = FALSE,
                                         row.names = NULL
                                         )



# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "14) Find the nearest genes to each sgRNA.RData")
     )



