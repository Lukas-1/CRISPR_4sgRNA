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





# Define lookup maps ------------------------------------------------------

source_abbreviations_vec <- c(
  "Brunello" = "Bru",
  "TKOv3"    = "TKOv3",
  "GPP"      = "GPP"
)





# Re-arrange the columns --------------------------------------------------

SNP_column_names <- c("sgRNA_SNP_IDs_vcf",     "sgRNA_SNP_AFs_1kGenomes", "sgRNA_SNP_AF_max_1kGenomes", "sgRNA_SNP_AF_sum_1kGenomes",
                      "sgRNA_SNP_AFs_TOPMED",  "sgRNA_SNP_AF_max_TOPMED", "sgRNA_SNP_AF_sum_TOPMED",
                      "sgRNA_SNP_AFs_Kaviar",  "sgRNA_SNP_AF_max_Kaviar", "sgRNA_SNP_AF_sum_Kaviar",
                      "sgRNA_SNP_IDs_1kG_ph1", "sgRNA_SNP_AFs_1kG_ph1",   "sgRNA_SNP_AF_max_1kG_ph1", "sgRNA_SNP_AF_sum_1kG_ph1",
                      "sgRNA_SNP_IDs_1kG_ph3", "sgRNA_SNP_AFs_1kG_ph3",   "sgRNA_SNP_AF_max_1kG_ph3", "sgRNA_SNP_AF_sum_1kG_ph3",
                      "sgRNA_SNP_IDs_gnomAD",  "sgRNA_SNP_AFs_gnomAD",    "sgRNA_SNP_AF_max_gnomAD",  "sgRNA_SNP_AF_sum_gnomAD",

                      "PAM_SNP_IDs_vcf",     "PAM_SNP_AFs_1kGenomes", "PAM_SNP_AF_max_1kGenomes", "PAM_SNP_AF_sum_1kGenomes",
                      "PAM_SNP_AFs_TOPMED",  "PAM_SNP_AF_sum_TOPMED", "PAM_SNP_AF_max_TOPMED",
                      "PAM_SNP_AFs_Kaviar",  "PAM_SNP_AF_sum_Kaviar", "PAM_SNP_AF_max_Kaviar",
                      "PAM_SNP_IDs_1kG_ph1", "PAM_SNP_AFs_1kG_ph1",   "PAM_SNP_AF_max_1kG_ph1", "PAM_SNP_AF_sum_1kG_ph1",
                      "PAM_SNP_IDs_1kG_ph3", "PAM_SNP_AFs_1kG_ph3",   "PAM_SNP_AF_max_1kG_ph3", "PAM_SNP_AF_sum_1kG_ph3",
                      "PAM_SNP_IDs_gnomAD",  "PAM_SNP_AFs_gnomAD",    "PAM_SNP_AF_max_gnomAD",  "PAM_SNP_AF_sum_gnomAD",

                      "all23_SNP_IDs_vcf",     "all23_SNP_AFs_1kGenomes", "all23_SNP_AF_max_1kGenomes", "all23_SNP_AF_sum_1kGenomes",
                      "all23_SNP_AFs_TOPMED",  "all23_SNP_AF_sum_TOPMED", "all23_SNP_AF_max_TOPMED",
                      "all23_SNP_AFs_Kaviar",  "all23_SNP_AF_sum_Kaviar", "all23_SNP_AF_max_Kaviar",
                      "all23_SNP_IDs_1kG_ph1", "all23_SNP_AFs_1kG_ph1",   "all23_SNP_AF_max_1kG_ph1", "all23_SNP_AF_sum_1kG_ph1",
                      "all23_SNP_IDs_1kG_ph3", "all23_SNP_AFs_1kG_ph3",   "all23_SNP_AF_max_1kG_ph3", "all23_SNP_AF_sum_1kG_ph3",
                      "all23_SNP_IDs_gnomAD",  "all23_SNP_AFs_gnomAD",    "all23_SNP_AF_max_gnomAD",  "all23_SNP_AF_sum_gnomAD",

                      "all22_SNP_AF_max_1kGenomes", "all22_SNP_AF_max_TOPMED",  "all22_SNP_AF_max_Kaviar",
                      "all22_SNP_AF_max_1kG_ph1",   "all22_SNP_AF_max_1kG_ph3", "all22_SNP_AF_max_gnomAD"
                      )


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
omit_SNP_columns <- omit_SNP_columns[!(omit_SNP_columns %in% c(preferred_rsID_column, preferred_AF_max_column))]
full_omit_columns <- c(omit_columns, omit_SNP_columns)





# Subset data / define sublibraries ---------------------------------------

merged_TF_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]






# Check for sgRNAs that violate the criteria of Graf et al. ---------------

TF_are_top4 <- merged_TF_CRISPRko_df[, "Rank"] %in% 1:4
TF_violate_Graf <- merged_TF_CRISPRko_df[, "CRISPOR_Graf_status"] != "GrafOK"

table(merged_TF_CRISPRko_df[TF_are_top4, "CRISPOR_Graf_status"])

merged_TF_CRISPRko_df[(TF_are_top4 & TF_violate_Graf) %in% TRUE, "CRISPOR_Doench_efficacy"]





# Write data to disk ------------------------------------------------------

DfToTSV(merged_CRISPRko_df, "CRISPRko_all_genes")
DfToTSV(merged_TF_CRISPRko_df, "CRISPRko_transcription_factors")









