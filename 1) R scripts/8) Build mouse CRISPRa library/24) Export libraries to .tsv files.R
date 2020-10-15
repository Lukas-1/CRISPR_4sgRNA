### 12th October 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory <- file.path(RData_directory, "7) Mouse - CRISPRa")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "Mouse - CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Read in gene lists.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))






# Create an empty SNP column ----------------------------------------------

merged_replaced_CRISPRa_df[[preferred_AF_max_column]] <- NA_real_





# Re-arrange the columns --------------------------------------------------

selected_columns <- c("Combined_ID", "Entrez_ID", "Gene_symbol", "Original_entrez",
                      "Original_symbol",

                      "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS", "Num_TSSs",

                      "Rank",

                      "Original_rank",
                      "Best_combination_rank",

                      "Num_overlaps",  "Overlaps_tolerance", "Spacing",

                      "Source", "mCRISPRa_v2_transcript", "Is_control",
                      "Entrez_source_Caprano", "Entrez_source_mCRISPRa_v2",

                      "sgRNA_sequence", "PAM", "Original_PAM",

                      "Caprano_rank", "GPP_rank", "mCRISPRa_v2_rank",
                      "Predicted_score", "Empirical_score", "Off_target_stringency",
                      "Sublibrary", "mCRISPRa_v2_ID", "mCRISPRa_TSS_source",
                      "Best_TSS", "First_TSS", "Last_TSS", "Strand_of_TSS",

                      "Chromosome", "Strand",
                      "Start", "End", "Cut_location", "Distance_from_TSS",
                      "TSS_searched_by_GuideScan", "TSS_regions",

                      "GuideScan_efficiency", "CRISPOR_Doench_efficacy",
                      "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score", "CRISPOR_Graf_status",

                      "GuideScan_specificity", "CRISPOR_3MM_specificity",
                      "CRISPOR_4MM_specificity", "CRISPOR_CFD_specificity", "CRISPOR_MIT_specificity",

                      "Num_0MM", "Num_1MM", "Num_5G_MM",
                      "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM",
                      "GuideScan_offtarget_category",
                      "CRISPOR_Num_0MM", "CRISPOR_Num_1MM", "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
                      "CRISPOR_off_target_count", "CRISPOR_Num_2or3MM",

                      "Nearest_Entrez_IDs", "Nearest_symbols", "Nearest_gene_distance",

                      "PAM_0MM", "PAM_1MM",

                      "Locations_0MM",

                      "Locations_1MM", "Sequences_1MM",

                      preferred_AF_max_column

                      # "Entrez_nearest_0MM", "Symbol_nearest_0MM",
                      # "Entrez_nearest_1MM", "Symbol_nearest_1MM"
                      )

merged_replaced_CRISPRa_df <- merged_replaced_CRISPRa_df[, selected_columns]






# Make adjustments to the 5' G-substituted library ------------------------

merged_replaced_CRISPRa_df[["Num_1MM"]] <- as.integer(rowSums(merged_replaced_CRISPRa_df[, c("Num_5G_MM", "Num_1MM")]))





# Look at the mCRISPRa-v2 "off-target stringency" annotations -------------

show_columns <- c("Combined_ID", "Entrez_ID", "Gene_symbol",
                  "Caprano_rank", "mCRISPRa_v2_rank", "Predicted_score", "Empirical_score",
                  "Off_target_stringency", "Best_TSS", "Chromosome", "Strand",
                  "GuideScan_specificity", "Num_0MM", "Num_5G_MM", "Num_1MM",
                  "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM", "GuideScan_offtarget_category"
                  )

head(merged_replaced_CRISPRa_df[(merged_replaced_CRISPRa_df[["Off_target_stringency"]] > 0) %in% TRUE, show_columns])







# Count the genes with multiple transcripts in mCRISPRa-v2 ----------------

num_transcripts_vec <- tapply(merged_replaced_CRISPRa_df[["mCRISPRa_v2_transcript"]],
                              factor(merged_replaced_CRISPRa_df[["Combined_ID"]], levels = unique(merged_replaced_CRISPRa_df[["Combined_ID"]])),
                              function(x) length(unique(x[!(is.na(x))]))
                              )
table(num_transcripts_vec > 1)





# Subset data / define sublibraries ---------------------------------------

membrane_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% membrane_het_df[["Entrez_ID"]], ]





# Check for identical sub-sequences ---------------------------------------

top_4_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Rank"]] %in% 1:4, ]
have_homologies <- CheckForIdenticalSubsequences(top_4_df, 9)

unique(top_4_df[have_homologies, c("Combined_ID", "Gene_symbol", "AltTSS_ID")])




# Check for duplicated sgRNAs with the same TSS_ID ------------------------

ID_paste_columns <- c("AltTSS_ID", "sgRNA_sequence", "mCRISPRa_v2_transcript")

sgRNA_ID_vec        <- do.call(paste, c(as.list(merged_replaced_CRISPRa_df[, ID_paste_columns[1:2]]), list(sep = "__")))
sgRNA_strict_ID_vec <- do.call(paste, c(as.list(merged_replaced_CRISPRa_df[, ID_paste_columns]), list(sep = "__")))
num_occurrences_lax <- table(sgRNA_ID_vec)[sgRNA_ID_vec]
num_occurrences_strict <- table(sgRNA_strict_ID_vec)[sgRNA_strict_ID_vec]

show_columns_TSS <- c(
  "sgRNA_ID", "PAM", "Original_PAM", "Entrez_ID", "Nearest_Entrez_IDs",
  "Gene_symbol", "Original_symbol", "Nearest_symbols",
  "Source", "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS",
  "mCRISPRa_v2_transcript", "Num_TSSs",
  "Rank", "Original_rank", "GPP_rank", "mCRISPRa_v2_rank"
)

duplicated_sgRNA_IDs_df <- merged_replaced_CRISPRa_df
duplicated_sgRNA_IDs_df[["sgRNA_ID"]] <- sgRNA_ID_vec
are_duplicated <- num_occurrences_lax > 1
duplicated_sgRNA_IDs_df <- duplicated_sgRNA_IDs_df[are_duplicated, ]
new_order <- order(match(duplicated_sgRNA_IDs_df[["AltTSS_ID"]], duplicated_sgRNA_IDs_df[["AltTSS_ID"]]),
                   duplicated_sgRNA_IDs_df[["sgRNA_ID"]]
                   )
duplicated_sgRNA_IDs_df <- duplicated_sgRNA_IDs_df[new_order, ]
num_occurrences_strict_short <- num_occurrences_strict[are_duplicated][new_order]
head(duplicated_sgRNA_IDs_df[, show_columns_TSS])

duplicated_sgRNA_IDs_df[num_occurrences_strict_short > 1, show_columns_TSS]






# Select columns to export ------------------------------------------------

omit_columns <- c("Combined_ID", "Sublibrary", "mCRISPRa_v2_ID", "Original_PAM",

                  "TSS_searched_by_GuideScan", "TSS_regions",
                  "Best_combination_rank", "Original_rank",

                  "TSS_number", "Num_TSSs", "Allocated_TSS",
                  "mCRISPRa_v2_transcript", "AltTSS_ID",
                  "Spacing", "Overlaps_tolerance",

                  "Num_5G_MM",

                  "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score",
                  "CRISPOR_MIT_specificity",
                  "CRISPOR_Graf_status",
                  "CRISPOR_Num_2or3MM", "CRISPOR_off_target_count",

                  "Is_control", "Entrez_ID_assignment",
                  "Nearest_Entrez_IDs", "Nearest_symbols", "Nearest_gene_distance",
                  # "Entrez_nearest_0MM", "Symbol_nearest_0MM", "Entrez_nearest_1MM", "Symbol_nearest_1MM",
                  "PAM_0MM", "PAM_1MM",
                  "Best_TSS", "First_TSS", "Last_TSS",
                  "mCRISPRa_TSS_source", "Strand_of_TSS",
                  "Start", "End", "GuideScan_Num_2or3MM", "Original_entrez",
                  "Entrez_source_mCRISPRa_v2", "Entrez_source_Caprano", "Off_target_stringency",

                  preferred_AF_max_column
                  )





# Write CRISPRa sgRNA libraries to disk -----------------------------------

DfToTSV(membrane_CRISPRa_df, "CRISPRa_membrane_heterogeneity", remove_columns = omit_columns)

DfToTSV(merged_replaced_CRISPRa_df, "CRISPRa_all_genes", remove_columns = omit_columns)










