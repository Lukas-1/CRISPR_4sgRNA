### 8th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR_4sgRNA"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory   <- file.path(RData_directory, "4) CRISPRi")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRi", "GuideScan")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "01) Compile predefined CRISPRi libraries - CRISPRi_df.RData"))
load(file.path(CRISPRi_RData_directory, "03) Map CRISPR libraries to TSS data.RData"))





# Convert TSS data to GuideScan input -------------------------------------

TSS_ranges_df <- TSSRangesForGuideScan(combined_TSS_CRISPRi_df)
combined_TSS_CRISPRi_df[["GuideScan_input"]] <- TSSStringForGuideScan(combined_TSS_CRISPRi_df)
combined_TSS_CRISPRi_df[["Region_length"]] <- TSS_ranges_df[["End"]] - TSS_ranges_df[["Start"]]

unique_GuideScan_input_vec <- unique(combined_TSS_CRISPRi_df[["GuideScan_input"]])






# Write GuideScan input files to disk -------------------------------------

write.table(unique_GuideScan_input_vec,
            file = file.path(GuideScan_files_directory, "Input_for_GuideScan_CRISPRi_all_TSSs.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )







