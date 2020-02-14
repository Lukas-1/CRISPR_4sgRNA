### 7th February 2020 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-02-07 - add primer sequences")
file_input_directory   <- file.path(file_directory, "1) Input")
file_output_directory  <- file.path(file_directory, "2) Output")





# Read in data ------------------------------------------------------------

input_df_list <- sapply(list.files(file_input_directory),
                        function(x) {
                          results_df <- as.data.frame(read_excel(file.path(file_input_directory, x), col_names = FALSE), stringsAsFactors = FALSE)[, 1:2]
                          colnames(results_df) <- c("Gene_symbol", "sgRNA_sequence")
                          results_df <- results_df[order(match(results_df[["Gene_symbol"]], results_df[["Gene_symbol"]])), ]
                          rownames(results_df) <- NULL
                          results_df <- data.frame(
                            results_df["Gene_symbol"],
                            "Rank" = unlist(lapply(unique(results_df[["Gene_symbol"]]),
                                                   function(x) seq_len(sum(x == results_df[["Gene_symbol"]]))
                                                   )
                                            ),
                            "Longest_shared_subsequence" = unlist(tapply(toupper(results_df[["sgRNA_sequence"]]),
                                                                         factor(results_df[["Gene_symbol"]], levels = unique(results_df[["Gene_symbol"]])),
                                                                         VectorizedLongestSubsequenceSingle
                                                                         ), use.names = FALSE
                                                                  ),
                            "Contains_polyT" = grepl("TTTT", results_df[["sgRNA_sequence"]], ignore.case = TRUE),
                            results_df["sgRNA_sequence"],
                            stringsAsFactors = FALSE,
                            row.names = NULL
                          )
                          results_df[["Sequence_with_primers"]] <- AddPrimers(results_df)
                          results_df[["Color"]] <- OnesAndZeros(results_df[["Gene_symbol"]]) + 1L
                          if (any(results_df[["Contains_polyT"]])) {
                            message(paste0(sum(results_df[["Contains_polyT"]]), " sgRNAs with TTTT sequences were found!"))
                            results_df[["Color"]][results_df[["Contains_polyT"]]] <- 3L
                          }
                          results_df[["Contains_polyT"]] <- ifelse(results_df[["Contains_polyT"]], "Yes", "No")
                          colnames(results_df)[colnames(results_df) == "Rank"] <- "Number"
                          return(results_df)
                        }, simplify = FALSE)





# Export data -------------------------------------------------------------

file_names <- sub(".xlsx", "_primers_added.tsv", names(input_df_list), fixed = TRUE)
for (i in seq_along(input_df_list)) {
  write.table(input_df_list[[i]][, colnames(input_df_list[[i]]) != "Number"],
              file = file.path(file_output_directory, file_names[[i]]),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )
}










