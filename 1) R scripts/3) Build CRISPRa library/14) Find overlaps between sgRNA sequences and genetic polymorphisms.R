### 26 July 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "08) Checking for overlaps with genetic polymorphisms.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "05) Compile data on genetic polymorphisms.RData"))
load(file.path(CRISPRa_RData_directory, "13) Integrate the output from GuideScan for individual sgRNA locations.RData"))





# Identify sgRNAs that have a defined location ----------------------------

are_mapped <- !(is.na(merged_replaced_CRISPRa_df[["Start"]]))

mapped_indices <- rep(NA_integer_, length(are_mapped))
mapped_indices[are_mapped] <- seq_len(sum(are_mapped))





# Search for overlaps between sgRNAs and genetic polymorphisms ------------

sgRNA_polymorphisms_df <- AllPolymorphisms(merged_replaced_CRISPRa_df[are_mapped, ])





# Find the nearest genes (0MM locations only) -----------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End")
nearest_columns <- c("Nearest_Entrez_IDs", "Nearest_symbols", "Distance")

nearest_genes_df <- FindNearestGenes(merged_replaced_CRISPRa_df[are_mapped, location_columns])[, nearest_columns]
names(nearest_genes_df)[[3]] <- "Nearest_gene_distance"




# Merge the data frame ----------------------------------------------------

merged_replaced_CRISPRa_df <- data.frame(merged_replaced_CRISPRa_df,
                                         nearest_genes_df[mapped_indices, ],
                                         sgRNA_polymorphisms_df[mapped_indices, ],
                                         stringsAsFactors = FALSE,
                                         row.names = NULL
                                         )



# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "14) Find overlaps between sgRNA sequences and genetic polymorphisms.RData")
     )








