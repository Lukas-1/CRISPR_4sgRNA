### 23rd November 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "10) Rat - General")
CRISPRko_RData_directory    <- file.path(RData_directory, "12) Rat - CRISPRko")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "Rat - CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Read in gene lists.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Create an empty SNP column ----------------------------------------------

merged_CRISPRko_df[[preferred_AF_max_column]] <- NA_real_





# Re-arrange the columns --------------------------------------------------

rearranged_column_names <- c(
  "Combined_ID", "Entrez_ID", "Entrez_overlapping_0MM", "Gene_symbol", "Symbol_overlapping_0MM", "Original_entrez", "Original_symbol",

  "Is_control",
  "Exon_number_GPP", "Transcript_ID",

  "Rank", "Original_rank", "Best_combination_rank",
  "Num_overlaps",  "Overlaps_tolerance", "Spacing",

  "Source",

  "sgRNA_sequence", "PAM", "Original_PAM",
  "sgRNA_context_sequence", "Original_cut_position", "Original_orientation",
  "GPP_rank",

  "Chromosome", "Strand", "Start", "End", "Cut_location",

  "Rule_set_2_score",

  "GuideScan_efficiency", "CRISPOR_Doench_efficacy",
  "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score", "CRISPOR_Graf_status",

  "GuideScan_specificity", "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity",
  "CRISPOR_CFD_specificity", "CRISPOR_MIT_specificity",

  "Num_0MM", "Num_1MM",
  "GuideScan_Num_2MM", "GuideScan_Num_3MM",
  "CRISPOR_Num_0MM", "CRISPOR_Num_1MM", "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
  "CRISPOR_off_target_count", "CRISPOR_Num_2or3MM",

  "PAM_0MM", "PAM_1MM",

  "Locations_0MM", #"Entrez_overlapping_0MM", "Symbol_overlapping_0MM",
  "Locations_1MM", "Sequences_1MM", "Entrez_overlapping_1MM", "Symbol_overlapping_1MM",

  preferred_AF_max_column
)

merged_CRISPRko_df <- merged_CRISPRko_df[, rearranged_column_names]





# Select columns to export ------------------------------------------------

omit_columns <- c("Combined_ID",
                  "Is_control",
                  "Original_order",
                  "Original_PAM", "sgRNA_context_sequence", "Original_cut_position", "Original_orientation",
                  # "Entrez_overlapping_0MM", "Symbol_overlapping_0MM",
                  "Entrez_overlapping_1MM", "Symbol_overlapping_1MM", "PAM_0MM", "PAM_1MM",
                  "Start", "End", "GuideScan_Num_2or3MM", "Original_entrez",
                  "Entrez_source_Brunello", "Entrez_source_TKOv3", "Off_target_stringency",

                  "CRISPOR_Num_2or3MM", "CRISPOR_off_target_count",

                  "Rule_set_2_score",
                  "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score", "CRISPOR_MIT_specificity",

                  "TKOv3_ID",
                  "Original_rank", "Best_combination_rank",
                  "Spacing", "Overlaps_tolerance",

                  "GuideScan_specificity", "GuideScan_efficiency",
                  "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_offtarget_category",

                  preferred_AF_max_column
                  )




# Subset data / define sublibraries ---------------------------------------

pharma_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% pharma_df[["Entrez_ID"]], ]






# Check for sgRNAs that violate the criteria of Graf et al. ---------------

are_top4 <- merged_CRISPRko_df[["Rank"]] %in% 1:4
violate_Graf <- merged_CRISPRko_df[["CRISPOR_Graf_status"]] != "GrafOK"

table(merged_CRISPRko_df[["CRISPOR_Graf_status"]][are_top4])

merged_CRISPRko_df[["CRISPOR_Doench_efficacy"]][(are_top4 & violate_Graf) %in% TRUE]






# Write CRISPRko sgRNA libraries to disk ----------------------------------

DfToTSV(merged_CRISPRko_df, "CRISPRko_all_genes", remove_columns = omit_columns)









