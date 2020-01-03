### 1st January 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # for AreCompleteTranscripts
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library.RData"))






# Define the sublibrary ---------------------------------------------------

merged_TF_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]





# Check for invalid 4sg combinations --------------------------------------

stopifnot(!(any(merged_TF_CRISPRko_df[, "Spacing"] %in% 0)))

are_complete_genes <- sapply(unique(merged_TF_CRISPRko_df[, "Combined_ID"]), function(x) {
  all(1:4 %in% merged_TF_CRISPRko_df[merged_TF_CRISPRko_df[, "Combined_ID"] == x, "Rank"])
})
stopifnot(all(are_complete_genes))





# Define the final selection of sgRNAs ------------------------------------

top4_df <- merged_TF_CRISPRko_df[merged_TF_CRISPRko_df[, "Rank"] %in% 1:4, ]
rownames(top4_df) <- NULL





# Identify problematic genes that do not meet all criteria ----------------

entrez_sets <- ReturnProblematicGenes(TF_summary_df)
lengths(entrez_sets)






# Count problematic genes -------------------------------------------------

are_problematic_sgRNAs   <- top4_df[, "Entrez_ID"] %in% entrez_sets[["problematic"]]
are_unproblematic_sgRNAs <- top4_df[, "Entrez_ID"] %in% entrez_sets[["unproblematic"]]

num_problematic   <- sum(are_problematic_sgRNAs)   / 4
num_unproblematic <- sum(are_unproblematic_sgRNAs) / 4

sum(top4_df[, "Combined_ID"] %in% entrez_sets[["overlap"]])       / 4
sum(top4_df[, "Combined_ID"] %in% entrez_sets[["fail_criteria"]]) / 4







# Count the number of genes chosen from each of the libraries -------------

num_gene_wells <- nrow(top4_df) / 4

are_Brunello_top4 <- grepl("Brunello", top4_df[, "Source"], fixed = TRUE)
are_TKOv3_top4    <- grepl("TKOv3",    top4_df[, "Source"], fixed = TRUE)
are_GPP_top4      <- grepl("GPP",      top4_df[, "Source"], fixed = TRUE)

table(are_Brunello_top4)
table(are_TKOv3_top4)
table(are_GPP_top4)







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

are_Brunello <- grepl("Brunello", merged_CRISPRko_df[, "Source"], fixed = TRUE)
are_good_controls <- AreGoodControls(merged_CRISPRko_df)

# The TKOv3 controls target EGFP, LacZ, luciferase, etc., which we may not want
control_sequences_pool <- toupper(merged_CRISPRko_df[are_Brunello & are_good_controls, "sgRNA_sequence"])

control_sequences_pool_list <- Make4sgControlsList(control_sequences_pool)





# Randomly select control guides, and construct a data frame --------------

selected_indices <- sample(seq_along(control_sequences_pool_list), num_control_wells)
selected_controls_list <- control_sequences_pool_list[selected_indices]

controls_df <- MakeControlGuidesDf(selected_controls_list, merged_CRISPRko_df)






# Select the targeting sgRNAs ---------------------------------------------

problematic_top4_df   <- top4_df[are_problematic_sgRNAs, ]
unproblematic_top4_df <- top4_df[are_unproblematic_sgRNAs, ]

problematic_indices_list   <- RandomizeAllIndices(n_total = num_problematic)
unproblematic_indices_list <- RandomizeAllIndices(n_total = num_unproblematic)

problematic_df_list   <- lapply(problematic_indices_list,   function(x) RetrieveIndices(problematic_top4_df,   x, "Combined_ID"))
unproblematic_df_list <- lapply(unproblematic_indices_list, function(x) RetrieveIndices(unproblematic_top4_df, x, "Combined_ID"))

combined_df_list <- c(unproblematic_df_list, problematic_df_list)






# Randomly shuffle each plate ---------------------------------------------

combined_df_shuffled_list <- RandomlyShufflePlates(combined_df_list)





# Add the control sgRNAs --------------------------------------------------

use_index <- length(combined_df_shuffled_list)
combined_df_shuffled_list[[use_index]] <- rbind.data.frame(combined_df_shuffled_list[[use_index]], controls_df,
                                                           stringsAsFactors = FALSE, make.row.names = FALSE
                                                           )





# Combine into a data frame -----------------------------------------------

TF_sgRNA_plates_df <- CombinePlateDfList(combined_df_shuffled_list)






# Save data ---------------------------------------------------------------

save(list = "TF_sgRNA_plates_df",
     file = file.path(CRISPRko_RData_directory, "14) Allocate sgRNAs to plates.RData")
     )











