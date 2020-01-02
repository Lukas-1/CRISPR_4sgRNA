### 30th September 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library.RData"))
load(file.path(CRISPRko_RData_directory, "14) Allocate sgRNAs to plates.RData"))





# Define lookup maps ------------------------------------------------------

source_abbreviations_vec <- c(
  "Brunello" = "Bru",
  "TKOv3"    = "TKOv3",
  "GPP"      = "GPP"
)





# Re-arrange the columns --------------------------------------------------

rearranged_column_names <- c(
  "Combined_ID", "Entrez_ID", "Gene_symbol", "Original_entrez", "Original_symbol",

  "Entrez_source_Brunello", "Entrez_source_TKOv3", "Is_control",
  "Exon_number_Brunello", "Exon_number_TKOv3", "Exon_number_GPP", "Transcript_ID", "Genomic_sequence_ID",

  "Original_order",
  "Rank", "Original_rank", "Best_combination_rank",
  "Num_overlaps",  "Overlaps_tolerance", "Spacing",

  "Source",

  "sgRNA_sequence", "PAM", "Original_PAM",
  "sgRNA_context_sequence", "Original_cut_position", "Original_orientation",
  "GPP_rank", "TKOv3_ID",

  "Chromosome", "Strand", "Start", "End", "Cut_location",

  "Rule_set_2_score",

  "GuideScan_efficiency", "CRISPOR_Doench_efficacy",
  "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score", "CRISPOR_Graf_status",

  "GuideScan_specificity", "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity",
  "CRISPOR_CFD_specificity", "CRISPOR_MIT_specificity",

  "Num_0MM", "Num_1MM",
  "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_Num_2or3MM",
  "GuideScan_offtarget_category",
  "CRISPOR_Num_0MM", "CRISPOR_Num_1MM", "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
  "CRISPOR_off_target_count", "CRISPOR_Num_2or3MM",

  SNP_column_names, #[SNP_column_names %in% colnames(merged_CRISPRko_df)],

  "PAM_0MM", "PAM_1MM",

  "Locations_0MM", "Entrez_overlapping_0MM", "Symbol_overlapping_0MM",
  "Locations_1MM", "Sequences_1MM", "Entrez_overlapping_1MM", "Symbol_overlapping_1MM"
)

merged_CRISPRko_df <- merged_CRISPRko_df[, rearranged_column_names]

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[, c("Plate_number", "Well_number", rearranged_column_names)]




# Select columns to export ------------------------------------------------

omit_columns <- c("Combined_ID",
                  "Is_control",
                  "Original_order",
                  "Original_PAM", "sgRNA_context_sequence", "Original_cut_position", "Original_orientation",
                  "Entrez_overlapping_0MM", "Symbol_overlapping_0MM", "Entrez_overlapping_1MM", "Symbol_overlapping_1MM", "PAM_0MM", "PAM_1MM",
                  "Start", "End", "GuideScan_Num_2or3MM", "Original_entrez", "Entrez_source_Brunello", "Entrez_source_TKOv3", "Off_target_stringency",

                  # "CRISPOR_Num_0MM", "CRISPOR_Num_1MM",  "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
                  "CRISPOR_Num_2or3MM", "CRISPOR_off_target_count",

                  "Rule_set_2_score",
                  "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score", "CRISPOR_MIT_specificity",

                  "TKOv3_ID",
                  "Original_rank", "Best_combination_rank",
                  "Spacing", "Overlaps_tolerance"
                  )

omit_SNP_columns <- grep("SNP", colnames(merged_CRISPRko_df), value = TRUE)
omit_SNP_columns <- setdiff(omit_SNP_columns, c(preferred_rsID_column, preferred_AF_max_column))
full_omit_columns <- c(omit_columns, omit_SNP_columns)





# Subset data / define sublibraries ---------------------------------------

merged_TF_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]






# Check for sgRNAs that violate the criteria of Graf et al. ---------------

TF_are_top4 <- merged_TF_CRISPRko_df[, "Rank"] %in% 1:4
TF_violate_Graf <- merged_TF_CRISPRko_df[, "CRISPOR_Graf_status"] != "GrafOK"

table(merged_TF_CRISPRko_df[TF_are_top4, "CRISPOR_Graf_status"])

merged_TF_CRISPRko_df[(TF_are_top4 & TF_violate_Graf) %in% TRUE, "CRISPOR_Doench_efficacy"]






# Write CRISPRko sgRNA libraries to disk ----------------------------------

TF_folder_name <- "TF library plate layout"
for (i in 1:4) {
  subset_df <- TF_sgRNA_plates_df[TF_sgRNA_plates_df[, "Rank"] %in% i, ]
  file_name <- paste0(file.path(TF_folder_name, "CRISPRko_TF_randomized_sg"), i)
  DfToTSV(subset_df, file_name, add_primers = TRUE)
}
DfToTSV(TF_sgRNA_plates_df, file.path(TF_folder_name, "CRISPRko_TF_randomized_all_4_guides"), add_primers = TRUE)

DfToTSV(merged_TF_CRISPRko_df, "CRISPRko_transcription_factors")
DfToTSV(merged_CRISPRko_df, "CRISPRko_all_genes")









