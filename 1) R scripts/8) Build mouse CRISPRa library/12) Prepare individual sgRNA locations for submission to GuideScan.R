### 9th October 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR_4sgRNA"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory   <- file.path(RData_directory, "7) Mouse - CRISPRa")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Mouse - CRISPRa", "GuideScan")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "11) Refine the genomic locations of sgRNA sequences after fixing 5'G substitutions.RData"))





# Prepare the input to GuideScan ------------------------------------------

were_not_searched <- !(merged_replaced_CRISPRa_df[["TSS_searched_by_GuideScan"]] %in% "Yes")
were_mapped <- !(is.na(merged_replaced_CRISPRa_df[["Start"]]))

submit_df <- merged_replaced_CRISPRa_df[were_not_searched & were_mapped, ]

submit_df[["GuideScan_input_sgRNA"]] <- sgRNAStringForGuideScan(submit_df)





# Check for duplicated chromosomal positions ------------------------------
# (Some of these duplications are not actually duplications, but rather,
# one sgRNA is on the + strand, and the other is on the - strand.)

num_occurrences <- table(submit_df[["GuideScan_input_sgRNA"]])[submit_df[["GuideScan_input_sgRNA"]]]

multiplicates_df <- submit_df[num_occurrences > 1, ]

choose_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "Source", "mCRISPRa_v2_transcript",
                    "sgRNA_sequence", "PAM", "Caprano_rank", "GPP_rank", "mCRISPRa_v2_rank",
                    "Chromosome", "Strand", "Start", "End",
                    "Distance_from_TSS", "TSS_searched_by_GuideScan",
                    "GuideScan_efficiency", "GuideScan_specificity",
                    "Num_0MM", "Num_5G_MM", "Num_1MM", "GuideScan_Num_2MM", "GuideScan_Num_3MM", "GuideScan_offtarget_category",
                    "GuideScan_input_sgRNA"
                    )

multiplicates_df <- multiplicates_df[order(match(multiplicates_df[["GuideScan_input_sgRNA"]], multiplicates_df[["GuideScan_input_sgRNA"]])), choose_columns]
row.names(multiplicates_df) <- NULL




# Define all unique chromosomal positions ---------------------------------

GuideScan_input_vec <- unique(submit_df[["GuideScan_input_sgRNA"]])




# Write GuideScan input files to disk -------------------------------------

write.table(GuideScan_input_vec,
            file = file.path(GuideScan_files_directory, "Input_for_GuideScan_CRISPRa_individual_locations.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )



