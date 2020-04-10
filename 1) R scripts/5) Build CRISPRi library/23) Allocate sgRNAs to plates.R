### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R")) # for FormatFixedWidthInteger
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRi_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRi_RData_directory, "21) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRi_RData_directory, "22) Summarize the human secretome sub-library.RData"))





# Exclude genes that are not protein-coding -------------------------------

TF_entrezs <- TF_overview_df[["Entrez_ID"]][!(is.na(TF_overview_df[["Num_total"]]))]
secretome_entrezs <- secretome_overview_df[["Entrez_ID"]][!(is.na(secretome_overview_df[["Num_total"]]))]
all_entrezs <- unique(c(collected_entrez_IDs, TF_entrezs, secretome_entrezs))

all_CRISPRi_df <- merged_replaced_CRISPRi_df[merged_replaced_CRISPRi_df[["Entrez_ID"]] %in% all_entrezs, ]




# Exclude incomplete transcripts, or those with shared subsequences -------

are_top4_mat <- CRISPRaAreTop4Mat(all_CRISPRi_df)





# Examine problematic genes -----------------------------------------------

invalid_combo_show_columns <- c(
  "Gene_symbol", "Source", "Chromosome", "Cut_location",
  "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS", "Num_TSSs",
  "Rank", "Original_rank",
  "Num_overlaps",  "Overlaps_tolerance", "Spacing",
  "Best_combination_rank",
  "sgRNA_sequence", "PAM"
)

all_CRISPRi_df[are_top4_mat[, "Are_valid_or_only_top4"] & !(are_top4_mat[, "Have_complete_guides"]), invalid_combo_show_columns]





# Define the sublibrary ---------------------------------------------------

are_TF <- all_CRISPRi_df[["Combined_ID"]] %in% TF_overview_df[["Combined_ID"]]

replaced_TF_CRISPRi_df <- all_CRISPRi_df[are_TF, ]
TF_are_top4_mat <- are_top4_mat[are_TF, ]





# Exclude unspaced/incomplete TF sgRNAs -----------------------------------

are_valid_top4 <- TF_are_top4_mat[, "Are_top4"] & TF_are_top4_mat[, "Have_valid_guides"]
are_invalid_top4 <- TF_are_top4_mat[, "Are_top4"] & !(are_valid_top4)


# Examine the excluded TF sgRNAs
replaced_TF_CRISPRi_df[are_invalid_top4, invalid_combo_show_columns]






# Define the final selection of sgRNAs ------------------------------------

top4_df <- replaced_TF_CRISPRi_df[are_valid_top4, ]
row.names(top4_df) <- NULL

table(table(top4_df[["AltTSS_ID"]]))





# Identify problematic genes that do not meet all criteria ----------------

entrez_sets <- ReturnProblematicGenes(TF_overview_df)
lengths(entrez_sets)





# Count problematic genes -------------------------------------------------

are_problematic_sgRNAs   <- top4_df[["Entrez_ID"]] %in% entrez_sets[["problematic"]]
are_unproblematic_sgRNAs <- top4_df[["Entrez_ID"]] %in% entrez_sets[["unproblematic"]]

num_problematic   <- sum(are_problematic_sgRNAs)   / 4
num_unproblematic <- sum(are_unproblematic_sgRNAs) / 4

sum(top4_df[["Combined_ID"]] %in% entrez_sets[["overlap"]])       / 4
sum(top4_df[["Combined_ID"]] %in% entrez_sets[["fail_criteria"]]) / 4





# Count the number of genes chosen from each of the libraries -------------

num_gene_wells <- nrow(top4_df) / 4

are_Dolcetto_top4    <- grepl("Dolcetto",    top4_df[["Source"]], fixed = TRUE)
are_hCRISPRi_v2_top4 <- grepl("hCRISPRi-v2", top4_df[["Source"]], fixed = TRUE)
are_GPP_top4         <- grepl("GPP",         top4_df[["Source"]], fixed = TRUE)

table(are_Dolcetto_top4)
table(are_hCRISPRi_v2_top4)
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

are_Dolcetto <- grepl("Dolcetto",    merged_replaced_CRISPRi_df[["Source"]], fixed = TRUE)
are_hCRISPRi <- grepl("hCRISPRi-v2", merged_replaced_CRISPRi_df[["Source"]], fixed = TRUE)

are_good_controls <- AreGoodControls(merged_replaced_CRISPRi_df)

controls_Dolcetto    <- merged_replaced_CRISPRi_df[["sgRNA_sequence"]][are_good_controls & are_Dolcetto]
controls_hCRISPRi_v2 <- merged_replaced_CRISPRi_df[["sgRNA_sequence"]][are_good_controls & are_hCRISPRi]

controls_hCRISPRi_v2_selected <- sample(controls_hCRISPRi_v2, length(controls_Dolcetto))

control_sequences_pool <- toupper(c(controls_Dolcetto, controls_hCRISPRi_v2_selected))

control_sequences_pool_list <- Make4sgControlsList(control_sequences_pool)





# Randomly select control guides, and construct a data frame --------------

selected_indices <- sample(seq_along(control_sequences_pool_list), num_control_wells)
selected_controls_list <- control_sequences_pool_list[selected_indices]

controls_df <- MakeControlGuidesDf(selected_controls_list, merged_replaced_CRISPRi_df)





# Select the targeting sgRNAs ---------------------------------------------

problematic_top4_df   <- top4_df[are_problematic_sgRNAs, ]
unproblematic_top4_df <- top4_df[are_unproblematic_sgRNAs, ]

problematic_indices_list   <- RandomizeAllIndices(n_total = num_problematic)
unproblematic_indices_list <- RandomizeAllIndices(n_total = num_unproblematic)

problematic_df_list   <- lapply(problematic_indices_list,   function(x) RetrieveIndices(problematic_top4_df,   x, "AltTSS_ID"))
unproblematic_df_list <- lapply(unproblematic_indices_list, function(x) RetrieveIndices(unproblematic_top4_df, x, "AltTSS_ID"))

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
     file = file.path(CRISPRi_RData_directory, "23) Allocate sgRNAs to plates.RData")
     )









