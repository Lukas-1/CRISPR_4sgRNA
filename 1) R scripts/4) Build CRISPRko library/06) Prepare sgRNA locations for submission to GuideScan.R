### 30th October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "04) Using GuideScan.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory  <- file.path(RData_directory, "3) CRISPRko")
GuideScan_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "GuideScan")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData"))





# Prepare the input to GuideScan ------------------------------------------

submit_df <- extended_CRISPRko_df[!(is.na(extended_CRISPRko_df[["Start"]])), ]





# Check for duplicated chromosomal positions ------------------------------
# (Some of these duplications are not actually duplications, instead, one sgRNA is on the + strand, and the other is on the - strand.)

submit_df[["GuideScan_input_sgRNA"]] <- sgRNAStringForGuideScan(submit_df)

num_occurrences <- table(submit_df[["GuideScan_input_sgRNA"]])[submit_df[["GuideScan_input_sgRNA"]]]

multiplicates_df <- submit_df[num_occurrences > 1, ]

multiplicates_df <- multiplicates_df[order(match(multiplicates_df[["GuideScan_input_sgRNA"]], multiplicates_df[["GuideScan_input_sgRNA"]])), ]
row.names(multiplicates_df) <- NULL





# Retain only unique chromosomal positions --------------------------------

submit_df <- submit_df[!(duplicated(submit_df[["GuideScan_input_sgRNA"]])), ]




# Filter for already present data -----------------------------------------

previous_guidescan_sgRNAs_df <- GetCRISPRkoGuideScanOutput()

are_already_present <- submit_df[["GuideScan_input_sgRNA"]] %in% previous_guidescan_sgRNAs_df[["Region"]]



# Split the input into chunks ---------------------------------------------

chunks_list <- AppendIDsWithoutCanonicalEntrezs(entrez_chunks_list, submit_df)

chunks_df <- data.frame(
  "Combined_ID" = unlist(chunks_list, use.names = FALSE),
  "Chunk_ID" = rep(names(chunks_list), lengths(chunks_list)),
  stringsAsFactors = FALSE
)

chunk_matches_vec <- match(submit_df[["Combined_ID"]], chunks_df[["Combined_ID"]])
submit_df[["Chunk"]] <- chunks_df[["Chunk_ID"]][chunk_matches_vec]

GuideScan_input_list <- split(submit_df[["GuideScan_input_sgRNA"]], submit_df[["Chunk"]])






# Write GuideScan input files to disk -------------------------------------

for (i in seq_along(GuideScan_input_list)) {
  file_name <- paste0("Input_for_GuideScan_CRISPRko__chunk_", names(GuideScan_input_list)[[i]], ".txt")
  write.table(GuideScan_input_list[[i]],
              file = file.path(GuideScan_files_directory, file_name),
              quote = FALSE, row.names = FALSE, col.names = FALSE
              )
}


file_name <- "Input_for_GuideScan_CRISPRko_all_guides.txt"
write.table(submit_df[["GuideScan_input_sgRNA"]],
            file = file.path(GuideScan_files_directory, file_name),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )

file_name <- "Input_for_GuideScan_CRISPRko_all_guides_filtered.txt"
write.table(submit_df[["GuideScan_input_sgRNA"]][!(are_already_present)],
            file = file.path(GuideScan_files_directory, file_name),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )






