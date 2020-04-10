### 13th November 2019 ###





# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "14) Find overlaps between sgRNA sequences and genetic polymorphisms.RData"))





# Find combinations of non-overlapping sgRNAs -----------------------------

# CRISPR_sub_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Gene_symbol"]] %in% "TRIM74", ]
# SeparateByTSS(CRISPR_sub_df)[, c("Source", "Gene_symbol", "Entrez_ID", "hCRISPRa_v2_transcript", "hCRISPRa_v2_ID", "TSS_number", "Allocated_TSS")]

merged_replaced_CRISPRa_df <- AllocateTSSforAllGenes(merged_replaced_CRISPRa_df,
                                                     parallel_mode = TRUE,
                                                     tolerate_divergent_chromosomes = TRUE
                                                     )





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData")
     )



# try_CRISPR_sub_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% "100288966", ]
#
#
# merged_replaced_CRISPRa_df[is.na(merged_replaced_CRISPRa_df[["TSS_number"]]), use_columns]
#
# table(is.na(merged_replaced_CRISPRa_df[["TSS_number"]]))
#
# table(is.na(merged_replaced_CRISPRa_df[["Allocated_TSS"]]))
#
#
# use_columns <- c("Source", "Gene_symbol", "Entrez_ID", "hCRISPRa_v2_transcript",
#                  "sgRNA_sequence", "PAM", "Chromosome", "Entrez_chromosome",
#                  "Best_TSS", "Start", #"Locations_0MM",
#                  "TSS_number", "Allocated_TSS", "Num_TSSs",
#                  "TSS_ID", "AltTSS_ID"
#                  )
#
# merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Gene_symbol"]] %in% "GOLGA6L9", use_columns]








