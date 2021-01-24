### 30th September 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory    <- file.path(RData_directory, "3) CRISPRko")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko")
previous_versions_directory <- file.path(RData_directory, "5) Previous versions of the library")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))
load(file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRko_RData_directory, "15) Allocate transcription factor sgRNAs to plates.RData"))
load(file.path(CRISPRko_RData_directory, "16) Find all genes targeted by each sgRNA - summary data.RData"))
load(file.path(previous_versions_directory, "02) CRISPRko transcription factor sub-library (1st version) - TF_v1_CRISPRko_df.RData"))





# Add data on other (unintended) targeted genes ---------------------------

merged_CRISPRko_df <- AddOtherTargets(merged_CRISPRko_df,
                                      guides_CDS_df,
                                      deletions_CDS_df
                                      )

matches_vec <- match(MakeIDs(TF_sgRNA_plates_df),
                     MakeIDs(merged_CRISPRko_df)
                     )
TF_sgRNA_plates_df <- AddOtherTargets(TF_sgRNA_plates_df,
                                      guides_CDS_df[matches_vec, ],
                                      deletions_CDS_df
                                      )





# Re-arrange the columns --------------------------------------------------

rearranged_column_names <- c(
  "Combined_ID",
  "Entrez_ID", "Other_target_Entrez_IDs", "Other_Entrez_IDs_4sg", "Original_entrez",
  "Gene_symbol", "Other_target_symbols", "Other_symbols_4sg", "Original_symbol",

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

  SNP_column_names,

  "PAM_0MM", "PAM_1MM",

  "Locations_0MM", #"Entrez_overlapping_0MM", "Symbol_overlapping_0MM",
  "Locations_1MM", "Sequences_1MM", "Entrez_overlapping_1MM", "Symbol_overlapping_1MM"
)

merged_CRISPRko_df <- merged_CRISPRko_df[, rearranged_column_names]

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[, c("Plate_number", "Well_number", rearranged_column_names)]




# Select columns to export ------------------------------------------------

omit_columns <- c("Combined_ID",
                  "Is_control",
                  "Original_order",
                  "Original_PAM", "sgRNA_context_sequence", "Original_cut_position", "Original_orientation",
                  # "Entrez_overlapping_0MM", "Symbol_overlapping_0MM",
                  "Entrez_overlapping_1MM", "Symbol_overlapping_1MM", "PAM_0MM", "PAM_1MM",
                  "Start", "End", "GuideScan_Num_2or3MM", "Original_entrez", "Entrez_source_Brunello", "Entrez_source_TKOv3", "Off_target_stringency",

                  # "CRISPOR_Num_0MM", "CRISPOR_Num_1MM",  "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
                  "CRISPOR_Num_2or3MM", "CRISPOR_off_target_count",

                  "Rule_set_2_score",
                  "CRISPOR_Moreno_Mateos", "CRISPOR_out_of_frame", "CRISPOR_lindel_score", "CRISPOR_MIT_specificity",

                  "TKOv3_ID",
                  "Original_rank", "Best_combination_rank",
                  "Spacing", "Overlaps_tolerance"
                  )

omit_SNP_columns <- grep("SNP", names(merged_CRISPRko_df), value = TRUE)
omit_SNP_columns <- setdiff(omit_SNP_columns, c(preferred_rsID_column, preferred_AF_max_column))
full_omit_columns <- c(omit_columns, omit_SNP_columns)





# Subset data / define sublibraries ---------------------------------------

merged_TF_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% TF_overview_df[["Combined_ID"]], ]

secretome_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% secretome_df[["Combined_ID"]], ]





# Check for sgRNAs that violate the criteria of Graf et al. ---------------

TF_are_top4 <- merged_TF_CRISPRko_df[["Rank"]] %in% 1:4
TF_violate_Graf <- merged_TF_CRISPRko_df[["CRISPOR_Graf_status"]] != "GrafOK"

table(merged_TF_CRISPRko_df[["CRISPOR_Graf_status"]][TF_are_top4])

merged_TF_CRISPRko_df[["CRISPOR_Doench_efficacy"]][(TF_are_top4 & TF_violate_Graf) %in% TRUE]





# Re-construct the TF sub-library in its original order -------------------

randomized_guide_IDs <- paste0(ifelse(TF_sgRNA_plates_df[["Is_control"]] == "Yes", "Control", TF_sgRNA_plates_df[["Combined_ID"]]),
                               "__", TF_sgRNA_plates_df[["sgRNA_sequence"]]
                               )
original_guide_IDs <- paste0(ifelse(merged_CRISPRko_df[["Is_control"]] == "Yes", "Control", merged_CRISPRko_df[["Combined_ID"]]),
                             "__", merged_CRISPRko_df[["sgRNA_sequence"]]
                             )
original_matches <- match(randomized_guide_IDs, original_guide_IDs)
if (any(is.na(original_matches))) {
  stop("Not all sgRNAs in TF_sgRNA_plates_df were found in the full library!")
}
TF_sgRNA_original_order_df <- TF_sgRNA_plates_df[order(original_matches), ]





# Write CRISPRko sgRNA libraries to disk ----------------------------------

TF_folder_name <- "TF sub-library (updated algorithm)"
for (i in 1:4) {
  subset_df <- TF_sgRNA_plates_df[TF_sgRNA_plates_df[["Rank"]] %in% i, ]
  file_name <- paste0(file.path(TF_folder_name, "CRISPRko_TF_sg"), i)
  DfToTSV(subset_df, file_name, add_primers = TRUE)
}
DfToTSV(TF_sgRNA_original_order_df, file.path(TF_folder_name, "CRISPRko_TF_ordered_by_gene"), add_primers = TRUE)
DfToTSV(TF_sgRNA_plates_df, file.path(TF_folder_name, "CRISPRko_TF_ordered_by_well"), add_primers = TRUE)

DfToTSV(merged_TF_CRISPRko_df, "CRISPRko_transcription_factors")
DfToTSV(secretome_CRISPRko_df, "CRISPRko_secretome")

DfToTSV(merged_CRISPRko_df, "CRISPRko_all_genes")




# Update the original transcription factor data frame ---------------------

TF_deletions_CDS_df <- FindOverlapsWithDeletions(TF_v1_CRISPRko_df[TF_v1_CRISPRko_df[["Is_control"]] == "No", ],
                                                 CDS_or_exon_locations_df
                                                 )[["summary_df"]]

TF_guides_CDS_df <- AlignSummaryDf(FindOverlappingGenes,
                                   TF_v1_CRISPRko_df,
                                   CDS_or_exon_locations_df
                                   )[["summary_df"]]

TF_v1_CRISPRko_df <- AddOtherTargets(TF_v1_CRISPRko_df,
                                     TF_guides_CDS_df,
                                     TF_deletions_CDS_df
                                     )



# Write changed wells to disk ---------------------------------------------

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[TF_sgRNA_plates_df[["Is_control"]] == "No", ]
old_TF_sgRNA_plates_df <- TF_v1_CRISPRko_df[TF_v1_CRISPRko_df[["Is_control"]] == "No", ]

TF_sgRNA_plates_df <- TF_sgRNA_plates_df[, c("Plate_number", "Well_number", rearranged_column_names)]
old_TF_sgRNA_plates_df <- old_TF_sgRNA_plates_df[, c("Plate_number", "Well_number", rearranged_column_names)]

are_new      <- !(TF_sgRNA_plates_df[["sgRNA_sequence"]] %in% old_TF_sgRNA_plates_df[["sgRNA_sequence"]])
are_obsolete <- !(old_TF_sgRNA_plates_df[["sgRNA_sequence"]] %in% TF_sgRNA_plates_df[["sgRNA_sequence"]])
DfToTSV(TF_sgRNA_plates_df[are_new, ],
        file.path(TF_folder_name, "CRISPRko_TF_ordered_by_well__new_guides"),
        add_primers = TRUE
        )
DfToTSV(old_TF_sgRNA_plates_df[are_obsolete, ],
        file.path(TF_folder_name, "CRISPRko_TF_ordered_by_well__obsolete_guides"),
        add_primers = TRUE
        )




# Save data ---------------------------------------------------------------

save(list = "TF_v1_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "17) Export libraries to .tsv files - TF_v1_CRISPRko_df.RData")
     )



