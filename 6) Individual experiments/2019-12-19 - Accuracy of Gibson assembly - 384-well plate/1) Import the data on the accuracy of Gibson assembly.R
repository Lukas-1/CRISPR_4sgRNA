### 19th December 2019 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR_4sgRNA"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
file_input_directory             <- file.path(file_directory, "1) Input")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")






# Read in data ------------------------------------------------------------

accuracies_df <- read.table(file.path(file_input_directory, "CCS5", "CCS5_identical.perc.txt"),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                            check.names = FALSE
                            )

guides_df <- data.frame(read_excel(file.path(file_input_directory, "1-384 oligo sequence.xlsx"),
                                            sheet = "1-384 primers", n_max = 384, col_names = FALSE
                                            ),
                                 stringsAsFactors = FALSE, check.names = FALSE
                                 )






# Process guides_df -------------------------------------------------------

are_empty_columns <- vapply(seq_len(ncol(guides_df)), function(x) all(is.na(guides_df[[x]])), logical(1))
guides_df <- guides_df[, !(are_empty_columns)]

symbols_mat <- do.call(cbind, unique(lapply(seq(1, ncol(guides_df), by = 2),
                                            function(x) sapply(strsplit(guides_df[[x]], "[_ ]"), "[[", 1)
                                            )
                                     )
                       )

are_identical <- apply(symbols_mat, 1, function(x) length(unique(x)) == 1)
symbols_mat[!(are_identical), ] # MIR20A is very likely an error introduced by Excel




# Compile data ------------------------------------------------------------

sequences_mat <- as.matrix(guides_df[, seq(2, ncol(guides_df), by = 2)])
colnames(sequences_mat) <- paste0("sg_", 1:4)

accuracy_colnames <- paste0("sg", 1:4, "_cr", 1:4, "_identical")

accuracies_1sg_mat <- do.call(cbind, lapply(accuracy_colnames, function(x) accuracies_df[[x]] / accuracies_df[["tot"]]))
colnames(accuracies_1sg_mat) <- paste0("Fraction_correct_sg", 1:4)

all_4_sums_vec <- rowSums(as.matrix(accuracies_df[, accuracy_colnames]))

assembly_df <- data.frame(
  "Gene_symbol"           = symbols_mat[, 1],
  "Num_reads"             = accuracies_df[["tot"]],
  "Fraction_correct_4sg"  = all_4_sums_vec / (accuracies_df[["tot"]] * 4),
  "Fraction_correct_3sg"  = rowSums(as.matrix(accuracies_df[, paste0("sg", 1:3, "_cr", 1:3, "_identical")])) /
                            (accuracies_df[["tot"]] * 3),
  accuracies_1sg_mat,
  "Longest_subsequence"   = NA,
  "Are_accuracy_outliers" = all_4_sums_vec == 0,
  sequences_mat,
  stringsAsFactors        = FALSE,
  row.names               = NULL
)

guides_list <- lapply(1:384, function(x) unname(sequences_mat[x, ]))

assembly_df[["Longest_subsequence"]] <- vapply(guides_list, LongestSharedSubsequence, integer(1))





# Re-order assembly_df ----------------------------------------------------

assembly_df <- assembly_df[order(assembly_df[["Fraction_correct_4sg"]], decreasing = FALSE), ]

row.names(assembly_df) <- NULL






# Save data ---------------------------------------------------------------

save(list = c("assembly_df", "guides_df"),
     file = file.path(intermediate_R_objects_directory,
                      "1) Import the data on the accuracy of Gibson assembly.RData"
                      )
     )









