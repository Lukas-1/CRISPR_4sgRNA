### 30th October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory  <- file.path(RData_directory, "3) CRISPRko")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "GuideScan")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData"))





# Prepare the input to GuideScan ------------------------------------------

submit_df <- extended_CRISPRko_df[!(is.na(extended_CRISPRko_df[, "Start"])), ]

submit_df[, "GuideScan_input_sgRNA"] <- sgRNAStringForGuideScan(submit_df)




# Check for duplicated chromosomal positions ------------------------------
# (Some of these duplications are not actually duplications, instead, one sgRNA is on the + strand, and the other is on the - strand.)

num_occurrences <- table(submit_df[, "GuideScan_input_sgRNA"])[submit_df[, "GuideScan_input_sgRNA"]]

multiplicates_df <- submit_df[num_occurrences > 1, ]

multiplicates_df <- multiplicates_df[order(match(multiplicates_df[, "GuideScan_input_sgRNA"], multiplicates_df[, "GuideScan_input_sgRNA"])), ]
rownames(multiplicates_df) <- NULL





# Define all unique chromosomal positions ---------------------------------

GuideScan_input_vec <- unique(submit_df[, "GuideScan_input_sgRNA"])





# Write GuideScan input files to disk -------------------------------------

write.table(GuideScan_input_vec,
            file = file.path(GuideScan_files_directory, "Input_for_GuideScan_individual_CRISPRko_sgRNAs.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )



