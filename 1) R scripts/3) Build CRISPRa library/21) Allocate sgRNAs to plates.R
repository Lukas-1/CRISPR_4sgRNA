### 8 December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "15) Separate sgRNAs for genes with multiple relevant TSSs.RData"))
load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRa_RData_directory, "20) Summarize the human transcription factor sub-library.RData"))





# Define functions --------------------------------------------------------


RandomizeAllIndices <- function(n_total = NULL, n_per_plate_vec = NULL, n_per_plate = (384 - 38)) {

  if (is.null(n_total) && is.null(n_per_plate_vec)) {
    stop("Either n_total or n_per_plate_vec must be specified!")
  }

  if (is.null(n_per_plate_vec)) {
    num_full_plates <- floor(n_total / n_per_plate)
    last_plate_n <- n_total - (num_full_plates * n_per_plate)
    n_per_plate_vec <- rep(n_per_plate, num_full_plates)
    if (last_plate_n > 0) {
      n_per_plate_vec <- c(n_per_plate_vec, last_plate_n)
    }
  } else {
    n_total <- sum(n_per_plate_vec)
  }

  indices_pool <- seq_len(n_total)
  num_plates <- length(n_per_plate_vec)
  results_list <- vector(mode = "list", length = num_plates)

  for (plate_number in seq_len(num_plates - 1L)) {
    this_sample <- sample(indices_pool, n_per_plate_vec[[plate_number]])
    results_list[[plate_number]] <- this_sample
    indices_pool <- indices_pool[!(indices_pool %in% this_sample)]
  }
  results_list[[num_plates]] <- sample(indices_pool, n_per_plate_vec[[num_plates]])
  return(results_list)
}


RetrieveIndices <- function(CRISPR_df, indices_vec, ID_column) {
  unique_combined_IDs <- unique(CRISPR_df[, ID_column])
  combined_ID_matches <- match(CRISPR_df[, ID_column], unique_combined_IDs)
  results_df <- CRISPR_df[combined_ID_matches %in% indices_vec, ]
  return(results_df)
}






# Define the sublibrary ---------------------------------------------------

replaced_TF_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[, "Combined_ID"] %in% TF_summary_df[, "Combined_ID"], ]






# Identify genes with < 4 non-overlapping guides --------------------------

are_problematic_genes <- TF_summary_df[, "Num_overlapping_transcripts"] > 0

problematic_genes    <- TF_summary_df[are_problematic_genes %in% TRUE, "Entrez_ID"]
nonoverlapping_genes <- TF_summary_df[are_problematic_genes %in% FALSE, "Entrez_ID"]

all(problematic_genes %in% replaced_TF_CRISPRa_df[, "Combined_ID"])
all(nonoverlapping_genes %in% replaced_TF_CRISPRa_df[, "Combined_ID"])

table(are_problematic_genes %in% FALSE)





# Count the number of genes chosen from each of the libraries -------------

are_complete_transcripts <- AreCompleteTranscripts(replaced_TF_CRISPRa_df)
are_top4 <- (replaced_TF_CRISPRa_df[, "Rank"] %in% 1:4) & are_complete_transcripts
table(are_top4)
sum(are_top4) / 4

are_Calabrese_TF   <- grepl("Calabrese",   replaced_TF_CRISPRa_df[, "Source"], fixed = TRUE)
are_hCRISPRa_v2_TF <- grepl("hCRISPRa-v2", replaced_TF_CRISPRa_df[, "Source"], fixed = TRUE)
are_GPP_TFF        <- grepl("GPP",         replaced_TF_CRISPRa_df[, "Source"], fixed = TRUE)

sum(are_top4 & are_Calabrese_TF)   / 4
sum(are_top4 & are_hCRISPRa_v2_TF) / 4
sum(are_top4 & are_GPP_TFF)         / 4




# Calculate the number of sgRNAs per plate --------------------------------

controls_fraction <- 0.015

num_controls_standard_plate <- floor(384 * controls_fraction)

are_problematic_sgRNAs    <- replaced_TF_CRISPRa_df[, "Entrez_ID"] %in% problematic_genes
are_nonoverlapping_sgRNAs <- replaced_TF_CRISPRa_df[, "Entrez_ID"] %in% nonoverlapping_genes

num_problematic <- sum(are_top4 & are_problematic_sgRNAs) / 4
num_nonoverlapping <- sum(are_top4 & are_nonoverlapping_sgRNAs) / 4

num_plate_5 <- num_nonoverlapping - (4 * (384 - 38))

num_controls_plate_5 <- (num_plate_5 * controls_fraction)
num_controls_problematic_plate <- (num_problematic * controls_fraction)

num_controls_plate_5 <- floor(num_controls_plate_5)
num_controls_problematic_plate <- floor(num_controls_problematic_plate)

(sum(are_top4 & are_nonoverlapping_sgRNAs) / 4) / (384 - 38)








# Set the seed ------------------------------------------------------------


####################
### Set the seed ###
####################
set.seed(1)
####################
####################
####################





# Generate a pool of 4sg controls -----------------------------------------

num_controls_vec <- c(rep(num_controls_standard_plate, 4), num_controls_plate_5, num_controls_problematic_plate)

are_controls <- merged_replaced_CRISPRa_df[, "Is_control"] == "Yes"
if (any(duplicated(toupper(merged_replaced_CRISPRa_df[are_controls, "sgRNA_sequence"])))) {
  stop("Error: Duplicatd control sgRNA sequences found!")
}

are_Calabrese <- grepl("Calabrese",   merged_replaced_CRISPRa_df[, "Source"], fixed = TRUE)
are_hCRISPRa  <- grepl("hCRISPRa-v2", merged_replaced_CRISPRa_df[, "Source"], fixed = TRUE)
have_no_issues <- (merged_replaced_CRISPRa_df[, "Num_0MM"] %in% 0) &
                  (merged_replaced_CRISPRa_df[, "Num_1MM"] %in% 0) &
                  !(grepl("TTTT", merged_replaced_CRISPRa_df[, "sgRNA_sequence"], ignore.case = TRUE))

controls_Calabrese   <- merged_replaced_CRISPRa_df[are_controls & have_no_issues & are_Calabrese, "sgRNA_sequence"]
controls_hCRISPRa_v2 <- merged_replaced_CRISPRa_df[are_controls & have_no_issues & are_hCRISPRa, "sgRNA_sequence"]

controls_hCRISPRa_v2_selected <- sample(controls_hCRISPRa_v2, length(controls_Calabrese))

guides_pool <- toupper(c(controls_Calabrese, controls_hCRISPRa_v2_selected))

indices_list <- RandomizeAllIndices(n_total = length(controls_Calabrese) * 2, n_per_plate = 4)
indices_list <- indices_list[lengths(indices_list) == 4]

guides_pool_list <- lapply(indices_list, function(x) guides_pool[x])

num_homologies_vec <- vapply(guides_pool_list, NumHomologousPairs, integer(1))

guides_pool_list <- guides_pool_list[num_homologies_vec == 0]




# Select a subset of the possible pool of 4sg controls --------------------

guides_pool_selected_indices <- sample(seq_along(guides_pool_list), sum(num_controls_vec))
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





# Add new columns to the control sgRNA data frame -------------------------

for (column_name in c("Best_combination_rank", "Spacing", "Overlaps_tolerance", "Num_overlaps", "Original_rank")) {
  control_sgRNAs_df[, column_name] <- NA_integer_
}
control_sgRNAs_df <- control_sgRNAs_df[, colnames(replaced_TF_CRISPRa_df)]




# Select the control sgRNAs for each plate --------------------------------

controls_indices_list <- RandomizeAllIndices(n_per_plate_vec = num_controls_vec)

controls_df_list <- lapply(controls_indices_list, function(x) RetrieveIndices(control_sgRNAs_df, x, "Combined_ID"))





# Build the data frames for targeting sgRNAs ------------------------------

problematic_TF_df    <- replaced_TF_CRISPRa_df[are_problematic_sgRNAs & are_top4, ]
nonoverlapping_TF_df <- replaced_TF_CRISPRa_df[are_nonoverlapping_sgRNAs & are_top4, ]





# Select the targeting sgRNAs for each plate ------------------------------

nonoverlapping_indices_list <- RandomizeAllIndices(n_total = num_nonoverlapping)

nonoverlapping_df_list <- lapply(nonoverlapping_indices_list, function(x) RetrieveIndices(nonoverlapping_TF_df, x, "AltTSS_ID"))

targeting_df_list <- c(nonoverlapping_df_list[1:2], list(problematic_TF_df), nonoverlapping_df_list[3:5])

controls_df_list <- controls_df_list[c(1:2, 6, 3:5)]




# Combine the data frames for each plate ----------------------------------

combined_df_list <- lapply(seq_along(targeting_df_list),
                           function(x) rbind.data.frame(targeting_df_list[[x]], controls_df_list[[x]], stringsAsFactors = FALSE, make.row.names = FALSE)
                           )





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

sgRNA_plates_df <- do.call(rbind.data.frame, c(combined_df_shuffled_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))






# Re-number the control wells ---------------------------------------------

are_controls <- sgRNA_plates_df[, "Is_control"] == "Yes"

sgRNA_plates_df[are_controls, "Combined_ID"] <- paste0("Control_", rep(seq_len(sum(are_controls) / 4), each = 4))






# Save data ---------------------------------------------------------------

save(list = "sgRNA_plates_df",
     file = file.path(CRISPRa_RData_directory, "21) Allocate sgRNAs to plates.RData")
     )






















