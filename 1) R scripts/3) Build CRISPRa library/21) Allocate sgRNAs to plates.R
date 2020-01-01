### 8 December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # for FormatFixedWidthInteger
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRa_RData_directory, "20) Summarize the human transcription factor sub-library.RData"))








# Define the sublibrary ---------------------------------------------------

replaced_TF_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]
rownames(replaced_TF_CRISPRa_df) <- NULL





# Exclude incomplete transcripts, or those with shared subsequences -------

are_complete_transcripts <- AreCompleteTranscripts(replaced_TF_CRISPRa_df)
have_zero_spacing        <- replaced_TF_CRISPRa_df[, "Spacing"] %in% 0 # "Zero" spacing means that there might be homologies / long shared subsequences among the chosen 4 guides!
are_top4                 <- replaced_TF_CRISPRa_df[, "Rank"] %in% 1:4
are_valid_top4           <- are_top4 & !(have_zero_spacing) & (are_complete_transcripts %in% TRUE)

stopifnot(all(replaced_TF_CRISPRa_df[are_top4 & !(are_valid_top4), "Num_TSSs"] >= 2)) # If a transcript is excluded, make sure that there is at least one other valid transcript for this gene!

# Examine the excluded sgRNAs
invalid_combo_show_columns <- c(
  "Gene_symbol", "Source", "Chromosome", "Cut_location",
  "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS", "Num_TSSs",
  "Rank", "Original_rank",
  "Num_overlaps",  "Overlaps_tolerance", "Spacing",
  "Best_combination_rank",
  "sgRNA_sequence", "PAM"
)
replaced_TF_CRISPRa_df[are_top4 & have_zero_spacing & (are_complete_transcripts %in% TRUE),  invalid_combo_show_columns]
replaced_TF_CRISPRa_df[are_top4 & have_zero_spacing & (are_complete_transcripts %in% FALSE), invalid_combo_show_columns]





# Identify problematic genes that do not meet all criteria ----------------

are_overlapping_genes      <- TF_summary_df[, "Num_overlapping_transcripts"] > 0
do_not_meet_criteria_genes <- TF_summary_df[, "Num_top4_outside_criteria"] > 0

genes_that_overlap       <- TF_summary_df[are_overlapping_genes %in% TRUE, "Entrez_ID"]
genes_that_fail_criteria <- TF_summary_df[do_not_meet_criteria_genes %in% TRUE, "Entrez_ID"]

are_problematic_genes <- are_overlapping_genes | do_not_meet_criteria_genes

problematic_genes   <- TF_summary_df[are_problematic_genes %in% TRUE,  "Entrez_ID"]
unproblematic_genes <- TF_summary_df[are_problematic_genes %in% FALSE, "Entrez_ID"]

stopifnot(all(problematic_genes %in% replaced_TF_CRISPRa_df[, "Combined_ID"]))
stopifnot(all(unproblematic_genes %in% replaced_TF_CRISPRa_df[, "Combined_ID"]))

table(are_problematic_genes, useNA = "ifany")





# Count the number of genes chosen from each of the libraries -------------

num_gene_wells <- sum(are_valid_top4) / 4

are_Calabrese_TF   <- grepl("Calabrese",   replaced_TF_CRISPRa_df[, "Source"], fixed = TRUE)
are_hCRISPRa_v2_TF <- grepl("hCRISPRa-v2", replaced_TF_CRISPRa_df[, "Source"], fixed = TRUE)
are_GPP_TFF        <- grepl("GPP",         replaced_TF_CRISPRa_df[, "Source"], fixed = TRUE)

sum(are_valid_top4 & are_Calabrese_TF)
sum(are_valid_top4 & are_hCRISPRa_v2_TF)
sum(are_valid_top4 & are_GPP_TFF)




# Count problematic genes -------------------------------------------------

are_problematic_sgRNAs   <- replaced_TF_CRISPRa_df[, "Entrez_ID"] %in% problematic_genes
are_unproblematic_sgRNAs <- replaced_TF_CRISPRa_df[, "Entrez_ID"] %in% unproblematic_genes

num_problematic   <- sum(are_valid_top4 & are_problematic_sgRNAs)   / 4
num_unproblematic <- sum(are_valid_top4 & are_unproblematic_sgRNAs) / 4

sum(are_valid_top4 & (replaced_TF_CRISPRa_df[, "Combined_ID"] %in% genes_that_overlap))       / 4
sum(are_valid_top4 & (replaced_TF_CRISPRa_df[, "Combined_ID"] %in% genes_that_fail_criteria)) / 4





# Set the number of control wells -----------------------------------------

num_control_wells <- 20





# Set the seed ------------------------------------------------------------


####################
### Set the seed ###
####################
set.seed(1)
####################
####################
####################





# Generate a pool of 4sg controls -----------------------------------------

are_controls <- merged_replaced_CRISPRa_df[, "Is_control"] == "Yes"
if (any(duplicated(toupper(merged_replaced_CRISPRa_df[are_controls, "sgRNA_sequence"])))) {
  stop("Error: Duplicated control sgRNA sequences found!")
}

are_Calabrese <- grepl("Calabrese",   merged_replaced_CRISPRa_df[, "Source"], fixed = TRUE)
are_hCRISPRa  <- grepl("hCRISPRa-v2", merged_replaced_CRISPRa_df[, "Source"], fixed = TRUE)
have_no_issues <- (merged_replaced_CRISPRa_df[, "Num_0MM"] %in% 0) &
                  (merged_replaced_CRISPRa_df[, "Num_1MM"] %in% 0) &
                  !(grepl("TTTT", merged_replaced_CRISPRa_df[, "sgRNA_sequence"], ignore.case = TRUE))

controls_Calabrese   <- merged_replaced_CRISPRa_df[are_controls & have_no_issues & are_Calabrese, "sgRNA_sequence"]
controls_hCRISPRa_v2 <- merged_replaced_CRISPRa_df[are_controls & have_no_issues & are_hCRISPRa,  "sgRNA_sequence"]

controls_hCRISPRa_v2_selected <- sample(controls_hCRISPRa_v2, length(controls_Calabrese))

guides_pool <- toupper(c(controls_Calabrese, controls_hCRISPRa_v2_selected))

indices_list <- RandomizeAllIndices(n_total = length(controls_Calabrese) * 2, n_per_plate = 4)
indices_list <- indices_list[lengths(indices_list) == 4]

guides_pool_list <- lapply(indices_list, function(x) guides_pool[x])

num_homologies_vec <- vapply(guides_pool_list, NumHomologousPairs, integer(1))

guides_pool_list <- guides_pool_list[num_homologies_vec == 0]




# Select a subset of the possible pool of 4sg controls --------------------

guides_pool_selected_indices <- sample(seq_along(guides_pool_list), num_control_wells)
guides_pool_selected <- guides_pool_list[guides_pool_selected_indices]

guides_pool_vec <- unlist(guides_pool_selected)
guides_pool_well_number <- rep(seq_along(guides_pool_selected), each = 4)
guides_pool_rep_number <-rep(1:4, times = length(guides_pool_selected))

control_sgRNAs_df <- merged_replaced_CRISPRa_df[are_controls & have_no_issues, ]

guides_pool_matches <- match(guides_pool_vec, toupper(control_sgRNAs_df[, "sgRNA_sequence"]))
control_sgRNAs_df <- control_sgRNAs_df[guides_pool_matches, ]
rownames(control_sgRNAs_df) <- NULL

control_sgRNAs_df[, "Combined_ID"] <- paste0("Control_", FormatFixedWidthInteger(guides_pool_well_number))
control_sgRNAs_df[, "Rank"] <- guides_pool_rep_number





# Shuffle the control sgRNAs ----------------------------------------------

controls_indices_vec <- RandomizeAllIndices(n_per_plate_vec = 20)[[1]]

controls_df <- RetrieveIndices(control_sgRNAs_df, controls_indices_vec, "Combined_ID")





# Build the data frames for targeting sgRNAs ------------------------------

problematic_TF_df   <- replaced_TF_CRISPRa_df[are_problematic_sgRNAs   & are_valid_top4, ]
unproblematic_TF_df <- replaced_TF_CRISPRa_df[are_unproblematic_sgRNAs & are_valid_top4, ]





# Select the targeting sgRNAs ---------------------------------------------

problematic_indices_list   <- RandomizeAllIndices(n_total = num_problematic)
unproblematic_indices_list <- RandomizeAllIndices(n_total = num_unproblematic)

problematic_df_list   <- lapply(problematic_indices_list,   function(x) RetrieveIndices(problematic_TF_df,   x, "AltTSS_ID"))
unproblematic_df_list <- lapply(unproblematic_indices_list, function(x) RetrieveIndices(unproblematic_TF_df, x, "AltTSS_ID"))





# Combine the data frames for each plate ----------------------------------

combined_df_list <- c(unproblematic_df_list, problematic_df_list)





# Randomly shuffle each plate ---------------------------------------------

shuffled_indices_list <- lapply(combined_df_list, function(x) sample(seq_len(nrow(x) / 4)))

combined_df_shuffled_list <- lapply(seq_along(shuffled_indices_list), function(x) {
  my_df <- combined_df_list[[x]]
  IDs_vec <- ifelse(my_df[, "Is_control"] == "Yes", my_df[, "Combined_ID"], my_df[, "AltTSS_ID"])
  unique_IDs_vec <- unique(IDs_vec)
  my_indices <- shuffled_indices_list[[x]]
  stopifnot(length(unique_IDs_vec) == length(my_indices))
  IDs_vec_shuffled <- unique_IDs_vec[my_indices]
  my_df <- my_df[order(match(IDs_vec, IDs_vec_shuffled)), ]
  return(my_df)
})





# Add the control sgRNAs --------------------------------------------------

use_index <- length(combined_df_shuffled_list)
combined_df_shuffled_list[[use_index]] <- rbind.data.frame(combined_df_shuffled_list[[use_index]], controls_df,
                                                           stringsAsFactors = FALSE, make.row.names = FALSE
                                                           )




# Add plate and well numbers to the data frame ----------------------------

combined_df_shuffled_list <- lapply(seq_along(combined_df_shuffled_list), function(x) {
  my_df <- combined_df_shuffled_list[[x]]
  results_df <- data.frame(
    "Plate_number" = x,
    "Well_number"  = rep(seq_len(nrow(my_df) / 4), each = 4),
    my_df,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
})




# Combine all the data frames ---------------------------------------------

TF_sgRNA_plates_df <- do.call(rbind.data.frame, c(combined_df_shuffled_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))






# Re-number the control wells ---------------------------------------------

are_controls <- TF_sgRNA_plates_df[, "Is_control"] == "Yes"

TF_sgRNA_plates_df[are_controls, "Combined_ID"] <- paste0("Control_", rep(seq_len(sum(are_controls) / 4), each = 4))






# Save data ---------------------------------------------------------------

save(list = "TF_sgRNA_plates_df",
     file = file.path(CRISPRa_RData_directory, "21) Allocate sgRNAs to plates.RData")
     )






















