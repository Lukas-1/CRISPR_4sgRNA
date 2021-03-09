### 27th January 2021 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

file_directory        <- "~/CRISPR/6) Individual experiments/2021-01-27 - calculate gene-wise p values for a pooled CRISPR screen"
file_input_directory  <- file.path(file_directory, "1) Input", "Sent by Davide")
file_output_directory <- file.path(file_directory, "2) Output")
RData_directory       <- file.path(file_directory, "3) R objects")





# Load data ---------------------------------------------------------------

load(file.path(RData_directory, "1) Calculate gene-wise p values for a pooled CRISPR screen.RData"))




# Read in data ------------------------------------------------------------

tsv_file_names <- list.files(file.path(file_input_directory, "1) Before TMM"))
before_df_list <- lapply(tsv_file_names, function(x) {
  read.table(file.path(file_input_directory, "1) Before TMM", x),
             header = TRUE,
             check.names = FALSE,
             stringsAsFactors = FALSE,
             sep = "\t"
             )
})
names(before_df_list) <- sub(".txt", "", tsv_file_names, fixed = TRUE)

common_columns <- c("ID", "TargetID", "Sequence", "GeneSymbol", "isControl")
stopifnot(length(unique(lapply(before_df_list, function(x) x[, common_columns]))) == 1)

before_mat <- do.call(cbind, lapply(before_df_list, function(x) x[["Count"]]))




# Define functions --------------------------------------------------------

GetTopGenes <- function(use_list) {
  results_list <- lapply(use_list, function(x) {
    are_top_genes <- x[["summary_df"]][["Top2_q_value"]] <= 0.05
    top_genes <- x[["summary_df"]][["Gene_symbol"]][are_top_genes]
    ### CHANGE THIS LATER ###
    ### CHANGE THIS LATER ###
    stopifnot(!(anyNA(top_genes)))
    ### CHANGE THIS LATER ###
    ### CHANGE THIS LATER ###
    return(top_genes)
  })
  return(results_list)
}
top_genes <- GetTopGenes(by_transcript_list)






