### 9th September 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory     <- file.path(RData_directory, "2) CRISPRa")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")
previous_versions_directory <- file.path(RData_directory, "5) Previous versions of the library")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData")) # for candidates_CRISPRa_df
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRa_RData_directory, "21) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRa_RData_directory, "23) Allocate transcription factor sgRNAs to plates.RData"))
load(file.path(previous_versions_directory, "01) CRISPRa transcription factor sub-library (1st version) - TF_v1_CRISPRa_df.RData"))







# Re-arrange the columns --------------------------------------------------

selected_columns <- c("Combined_ID", "Entrez_ID", "Gene_symbol", "Original_entrez",
                      "Original_symbol",

                      "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS", "Num_TSSs",

                      "Rank",

                      "Original_rank",
                      "Best_combination_rank",

                      "Num_overlaps",  "Overlaps_tolerance", "Spacing",

                      "Source", "hCRISPRa_v2_transcript", "Is_control",
                      "Entrez_source_Calabrese", "Entrez_source_hCRISPRa_v2",

                      "sgRNA_sequence", "PAM", "Original_PAM",

                      "Calabrese_rank", "GPP_rank", "hCRISPRa_v2_rank",
                      "Predicted_score", "Empirical_score", "Off_target_stringency",
                      "Sublibrary", "hCRISPRa_v2_ID", "hCRISPRa_TSS_source",
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

                      SNP_column_names,

                      "Nearest_Entrez_IDs", "Nearest_symbols", "Nearest_gene_distance",

                      "PAM_0MM", "PAM_1MM",

                      "Locations_0MM", #"Entrez_nearest_0MM", "Symbol_nearest_0MM",

                      "Locations_1MM", "Sequences_1MM"
                      # "Entrez_nearest_1MM", "Symbol_nearest_1MM"
                      )

merged_replaced_CRISPRa_df <- merged_replaced_CRISPRa_df[, selected_columns]

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[, c("Plate_number", "Well_number", selected_columns)]






# Make adjustments to the 5' G-substituted library ------------------------

merged_replaced_CRISPRa_df[["Num_1MM"]] <- as.integer(rowSums(merged_replaced_CRISPRa_df[, c("Num_5G_MM", "Num_1MM")]))





# Look at the hCRISPRa-v2 "off-target stringency" annotations -------------

show_columns <- c("Combined_ID", "Entrez_ID", "Gene_symbol",
                  "Calabrese_rank", "hCRISPRa_v2_rank", "Predicted_score", "Empirical_score",
                  "Off_target_stringency", "Best_TSS", "Chromosome", "Strand",
                  "GuideScan_specificity", "Num_0MM", "Num_5G_MM", "Num_1MM",
                  "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM", "GuideScan_offtarget_category"
                  )

head(merged_replaced_CRISPRa_df[(merged_replaced_CRISPRa_df[["Off_target_stringency"]] > 0) %in% TRUE, show_columns])







# Count the genes with multiple transcripts in hCRISPRa-v2 ----------------

num_transcripts_vec <- tapply(merged_replaced_CRISPRa_df[["hCRISPRa_v2_transcript"]],
                              factor(merged_replaced_CRISPRa_df[["Combined_ID"]], levels = unique(merged_replaced_CRISPRa_df[["Combined_ID"]])),
                              function(x) length(unique(x[!(is.na(x))]))
                              )
table(num_transcripts_vec > 1)





# Subset data / define sublibraries ---------------------------------------

is_curated <- grepl("^Curated", merged_replaced_CRISPRa_df[["Source"]])
replaced_curated_CRISPRa_df <- merged_replaced_CRISPRa_df[is_curated, ]

merged_replaced_candidates_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% candidates_CRISPRa_df[["Combined_ID"]], ]

replaced_TF_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% TF_overview_df[["Combined_ID"]], ]

secretome_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% secretome_df[["Combined_ID"]], ]




# Check for missing GuideScan data ----------------------------------------

top_4_TF_df <- replaced_TF_CRISPRa_df[replaced_TF_CRISPRa_df[["Rank"]] %in% 1:4, ]

unique(top_4_TF_df[is.na(top_4_TF_df[["GuideScan_specificity"]]), c("Combined_ID", "Gene_symbol", "AltTSS_ID")])
table(top_4_TF_df[["GuideScan_specificity"]] < 0.2, useNA = "ifany")






# Check for identical sub-sequences ---------------------------------------

top_4_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Rank"]] %in% 1:4, ]
have_homologies <- CheckForIdenticalSubsequences(top_4_df, 9)

unique(top_4_df[have_homologies, c("Combined_ID", "Gene_symbol", "AltTSS_ID")])

for (num_bases in 8:9) {
  unique_TSS_columns <- c("Combined_ID", "Gene_symbol", "AltTSS_ID")
  total_num <- nrow(unique(top_4_TF_df[, unique_TSS_columns]))
  homologous_TF_df <- top_4_TF_df[CheckForIdenticalSubsequences(top_4_TF_df, num_bases), ]
  unique_homologous_TF_df <- unique(homologous_TF_df[, unique_TSS_columns])
  message(paste0("For subsequences of length ", num_bases, ", ", nrow(unique_homologous_TF_df),
                 " problematic 4 sg combinations were found (out of a total of ",
                 total_num, ")!"
                 )
          )
  print(unique_homologous_TF_df)
  message("")
}





# Check for sgRNAs that would have been wrongly excluded ------------------
# ... due to single-nucleotide polymorphisms at the N(GG) position within the PAM!

overlap_with_SNP <- (top_4_TF_df[["all22_SNP_AF_max_Kaviar"]] > 0.01) %in% TRUE
overlap_with_indeterminate <- !(overlap_with_SNP) & ((top_4_TF_df[["all23_SNP_AF_max_Kaviar"]] > 0.01) %in% TRUE)
table(overlap_with_indeterminate)
table(overlap_with_SNP)
table((replaced_TF_CRISPRa_df[["all22_SNP_AF_max_Kaviar"]] > 0.01) %in% TRUE, useNA = "ifany")

SNP_indeterminate_columns <- c("Gene_symbol", "Rank", "Chromosome", "Cut_location", "sgRNA_sequence",
                               "all23_SNP_IDs_vcf", "all22_SNP_AF_max_Kaviar", "all23_SNP_AF_max_Kaviar"
                               )

top_4_TF_df[overlap_with_indeterminate, SNP_indeterminate_columns]
top_4_TF_df[overlap_with_SNP, SNP_indeterminate_columns]




# Check for duplicated sgRNAs with the same TSS_ID ------------------------

ID_paste_columns <- c("AltTSS_ID", "sgRNA_sequence", "hCRISPRa_v2_transcript")

sgRNA_ID_vec        <- do.call(paste, c(as.list(merged_replaced_CRISPRa_df[, ID_paste_columns[1:2]]), list(sep = "__")))
sgRNA_strict_ID_vec <- do.call(paste, c(as.list(merged_replaced_CRISPRa_df[, ID_paste_columns]), list(sep = "__")))
num_occurrences_lax <- table(sgRNA_ID_vec)[sgRNA_ID_vec]
num_occurrences_strict <- table(sgRNA_strict_ID_vec)[sgRNA_strict_ID_vec]

show_columns_TSS <- c(
  "sgRNA_ID", "PAM", "Original_PAM", "Entrez_ID", "Nearest_Entrez_IDs", "Gene_symbol", "Original_symbol", "Nearest_symbols",
  "Source", "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS", "hCRISPRa_v2_transcript", "Num_TSSs",
  "Rank", "Original_rank", "GPP_rank", "hCRISPRa_v2_rank"
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






# Check for duplicated sgRNAs among the transcription factors -------------

duplicated_df <- FindSharedsgRNAs(replaced_TF_CRISPRa_df)

duplicated_df_for_export <- duplicated_df

duplicated_df_for_export[["Num_all"]] <- ifelse(is.na(duplicated_df_for_export[["Num_all"]]),
                                                "",
                                                duplicated_df_for_export[["Num_all"]]
                                                )
duplicated_df_for_export[["AltTSS_ID"]] <- ifelse(duplicated_df_for_export[["AltTSS_ID"]] == duplicated_df_for_export[["Entrez_ID"]],
                                                  " ",
                                                  duplicated_df_for_export[["AltTSS_ID"]]
                                                  )

write.table(duplicated_df_for_export,
            file = file.path(file_output_directory, "CRISPRa_duplicated_sgRNAs_transcription_factors.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
            )





# Select columns to export ------------------------------------------------

omit_columns <- c("Combined_ID", "Sublibrary", "hCRISPRa_v2_ID", "Original_PAM",

                  "TSS_searched_by_GuideScan", "TSS_regions",
                  "Best_combination_rank", "Original_rank",

                  "TSS_number", "Num_TSSs", "Allocated_TSS", "hCRISPRa_v2_transcript", "AltTSS_ID",
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
                  "hCRISPRa_TSS_source", "Strand_of_TSS",
                  "Start", "End", "GuideScan_Num_2or3MM", "Original_entrez",
                  "Entrez_source_hCRISPRa_v2", "Entrez_source_Calabrese", "Off_target_stringency"
                  )
omit_SNP_columns <- grep("SNP", names(merged_replaced_CRISPRa_df), value = TRUE)
omit_SNP_columns <- setdiff(omit_SNP_columns, c(preferred_rsID_column, preferred_AF_max_column))
full_omit_columns <- c(omit_columns, omit_SNP_columns)





# Re-construct the TF sub-library in its original order -------------------

ID_paste_all_columns <- c(ID_paste_columns, c("Original_PAM", "hCRISPRa_v2_rank", "Original_symbol"))

TF_sgRNA_plates_df_copy <- TF_sgRNA_plates_df
TF_sgRNA_plates_df_copy[["Combined_ID"]][TF_sgRNA_plates_df_copy[["Is_control"]] == "Yes"] <- "Control"

randomized_guide_IDs <- do.call(paste, c(as.list(TF_sgRNA_plates_df_copy[, ID_paste_all_columns]),    list(sep = "__")))
original_guide_IDs   <- do.call(paste, c(as.list(merged_replaced_CRISPRa_df[, ID_paste_all_columns]), list(sep = "__")))

original_matches <- match(randomized_guide_IDs, original_guide_IDs)

if (any(is.na(original_matches))) {
  stop("Not all sgRNAs in TF_sgRNA_plates_df were found in the full library!")
}
TF_sgRNA_original_order_df <- TF_sgRNA_plates_df[order(original_matches), ]





# Write CRISPRa sgRNA libraries to disk -----------------------------------

DfToTSV(replaced_TF_CRISPRa_df, "CRISPRa_transcription_factors")
DfToTSV(secretome_CRISPRa_df, "CRISPRa_secretome")


TF_folder_name <- "TF library plate layout"
for (i in 1:4) {
  subset_df <- TF_sgRNA_plates_df[TF_sgRNA_plates_df[["Rank"]] %in% i, ]
  file_name <- paste0(file.path(TF_folder_name, "CRISPRa_TF_randomized_sg"), i)
  DfToTSV(subset_df, file_name, add_primers = TRUE)
}
DfToTSV(TF_sgRNA_plates_df,         file.path(TF_folder_name, "CRISPRa_TF_randomized_all_4_guides"),     add_primers = TRUE)
DfToTSV(TF_sgRNA_original_order_df, file.path(TF_folder_name, "CRISPRa_TF_original_order_all_4_guides"), add_primers = TRUE)

DfToTSV(replaced_curated_CRISPRa_df, file.path("Candidate genes", "CRISPRa_16_sgRNAs"), allow_curated = TRUE)
DfToTSV(merged_replaced_candidates_CRISPRa_df, file.path("Candidate genes", "CRISPRa_20_trial_genes"))
DfToTSV(merged_replaced_CRISPRa_df, "CRISPRa_all_genes")

DfToTSV(merged_replaced_CRISPRa_df, "CRISPRa_all_genes_all_SNP_databases", remove_columns = omit_columns)





# Write changed wells to disk ---------------------------------------------

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[TF_sgRNA_plates_df[["Is_control"]] == "No", ]
old_TF_sgRNA_plates_df <- TF_v1_CRISPRa_df[TF_v1_CRISPRa_df[["Is_control"]] == "No", ]

old_TF_sgRNA_plates_df[["Sequences_1MM"]] <- TruncateLongEntries(old_TF_sgRNA_plates_df[["Sequences_1MM"]])
old_TF_sgRNA_plates_df[["Locations_1MM"]] <- TruncateLongEntries(old_TF_sgRNA_plates_df[["Locations_1MM"]])

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[, c("Plate_number", "Well_number", selected_columns)]
old_TF_sgRNA_plates_df <- old_TF_sgRNA_plates_df[, c("Plate_number", "Well_number", intersect(selected_columns, colnames(old_TF_sgRNA_plates_df)))]

are_new      <- !(TF_sgRNA_plates_df[["sgRNA_sequence"]] %in% old_TF_sgRNA_plates_df[["sgRNA_sequence"]])
are_obsolete <- !(old_TF_sgRNA_plates_df[["sgRNA_sequence"]] %in% TF_sgRNA_plates_df[["sgRNA_sequence"]])
DfToTSV(TF_sgRNA_plates_df[are_new, ],
        file.path(TF_folder_name, "CRISPRa_TF_randomized_all_4_guides__new_guides"),
        add_primers = TRUE
        )
DfToTSV(old_TF_sgRNA_plates_df[are_obsolete, ],
        file.path(TF_folder_name, "CRISPRa_TF_randomized_all_4_guides__obsolete_guides"),
        add_primers = TRUE
        )





