### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRi-v2 - CRISPRi_df.RData"))
load(file.path(CRISPRi_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))
load(file.path(CRISPRi_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData"))
load(file.path(CRISPRi_RData_directory, "09) Filter the output from GuideScan for TSS regions after fixing 5'G substitutions.RData"))
load(file.path(CRISPRi_RData_directory, "10) Find matches for sgRNA sequences after fixing 5'G substitutions - replaced_genome_search_df.RData"))






# Integrate the results of a genome search for sgRNA sequences ------------

replaced_CRISPRi_df <- replaced_merged_CRISPRi_df[, names(CRISPRi_df)]

extended_replaced_CRISPRi_df <- ExtendWithGenomeSearch(replaced_CRISPRi_df, replaced_genome_search_df)





# Flip lower-/upper-case to indicate 5' G replacements --------------------

first_letter_vec <- substr(extended_replaced_CRISPRi_df[["sgRNA_sequence"]], 1, 1)
first_letter_cases_vec  <- GetCases(first_letter_vec)
second_letter_cases_vec <- GetCases(substr(extended_replaced_CRISPRi_df[["sgRNA_sequence"]], 2, 2))
are_not_discordant_vec  <- first_letter_cases_vec == second_letter_cases_vec

flipped_vec <- ifelse(are_not_discordant_vec & (replaced_merged_CRISPRi_df[["Exchanged_5pG"]] %in% "Yes"),
                      paste0(ifelse(first_letter_cases_vec == "upper", tolower(first_letter_vec), toupper(first_letter_vec)),
                             substr(extended_replaced_CRISPRi_df[["sgRNA_sequence"]], 2, nchar(extended_replaced_CRISPRi_df[["sgRNA_sequence"]]))
                             ),
                      extended_replaced_CRISPRi_df[["sgRNA_sequence"]]
                      )
extended_replaced_CRISPRi_df[["sgRNA_sequence"]] <- flipped_vec




# Merge with GuideScan data -----------------------------------------------

full_merged_replaced_CRISPRi_df <- MergeTSSandGuideScan(extended_replaced_CRISPRi_df, replaced_guidescan_all_genes_df, combined_TSS_CRISPRi_df)
head(full_merged_replaced_CRISPRi_df[is.na(full_merged_replaced_CRISPRi_df[["Hits_start"]]) & !(is.na(full_merged_replaced_CRISPRi_df[["GuideScan_start"]])), ])

merged_replaced_CRISPRi_df <- AdjustPositionColumns(full_merged_replaced_CRISPRi_df, replaced_guidescan_all_genes_df, combined_TSS_CRISPRi_df)




# Check for inconsistent chromosome mappings ------------------------------

inconsistent_IDs <- CheckForInconsistentChromosomes(merged_replaced_CRISPRi_df)
message(paste0("The following combined IDs had sgRNAs that mapped to more than one chromosome: ",
               paste0(inconsistent_IDs, collapse = ", "), "!"
               )
        )

selected_gene_columns <- c("Source", "Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "Chromosome")
unique(merged_replaced_CRISPRi_df[merged_replaced_CRISPRi_df[["Combined_ID"]] %in% inconsistent_IDs, selected_gene_columns])





# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRi_df",
     file = file.path(CRISPRi_RData_directory, "11) Refine the genomic locations of sgRNA sequences after fixing 5'G substitutions.RData")
     )


