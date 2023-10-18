### 25th July 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR_4sgRNA"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory   <- file.path(RData_directory, "2) CRISPRa")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "GuideScan")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData"))
load(file.path(CRISPRa_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))






# Convert TSS data to GuideScan input for the candidate genes -------------

candidates_TSS_df <- combined_TSS_CRISPRa_df[combined_TSS_CRISPRa_df[["Combined_ID"]] %in% candidates_CRISPRa_df[["Combined_ID"]], ]

TSS_ranges_candidates_df <- TSSRangesForGuideScan(candidates_TSS_df)

candidates_GuideScan_input_vec <- unique(TSSStringForGuideScan(candidates_TSS_df))





# Convert TSS data to GuideScan input for all genes -----------------------

TSS_ranges_df <- TSSRangesForGuideScan(combined_TSS_CRISPRa_df)
combined_TSS_CRISPRa_df[["GuideScan_input"]] <- TSSStringForGuideScan(combined_TSS_CRISPRa_df)
combined_TSS_CRISPRa_df[["Region_length"]] <- TSS_ranges_df[["End"]] - TSS_ranges_df[["Start"]]

unique_GuideScan_input_vec <- unique(combined_TSS_CRISPRa_df[["GuideScan_input"]])






# Write GuideScan input files to disk -------------------------------------

write.table(candidates_GuideScan_input_vec,
            file = file.path(GuideScan_files_directory, "Input_for_GuideScan_CRISPRa_candidate_gene_TSSs.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )
write.table(unique_GuideScan_input_vec,
            file = file.path(GuideScan_files_directory, "Input_for_GuideScan_CRISPRa_all_TSSs.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )







