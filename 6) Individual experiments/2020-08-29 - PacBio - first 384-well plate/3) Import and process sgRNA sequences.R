### 10 September 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory  <- file.path(file_directory, "2) Input")
R_objects_directory   <- file.path(file_directory, "3) R objects")

metadata_directory    <- file.path(file_input_directory, "Metadata")
sg_sequences_file     <- file.path(metadata_directory, "1-384 sgRNA summary.xlsm")






# Read in data ------------------------------------------------------------

gRNAs_excel_df <- as.data.frame(read_excel(sg_sequences_file, sheet = 2),
                                stringsAsFactors = FALSE, check.names = FALSE
                                )





# Process guides ----------------------------------------------------------

column_names <- paste0("name of sg", 1:4)

column_splits <- lapply(gRNAs_excel_df[, column_names],
                        function(x) strsplit(x, "[_ ]")
                        )

gene_names_list <- lapply(column_splits, function(x) sapply(x, "[[", 1))
modality_list <- lapply(column_splits, function(x) substr(sapply(x, "[[", 2), 1, 1))
stopifnot(length(unique(modality_list)) == 1)

guide_seq_mat <- as.matrix(gRNAs_excel_df[, paste0("sg", 1:4, " sequence")])
colnames(guide_seq_mat) <- paste0("Sequence_sg", 1:4)

sg_sequences_df <- data.frame(
  "Well_number"         = seq_len(384),
  "Target_gene"         = gene_names_list[[1]],
  "Modality"            = paste0("CRISPR", modality_list[[1]]),
  guide_seq_mat,
  "Longest_subsequence" = apply(guide_seq_mat, 1, LongestSharedSubsequence),
  "Empty_well"          = seq_len(384) %in% c(51, 74),
  stringsAsFactors      = FALSE,
  row.names             = NULL
)





# Save data ---------------------------------------------------------------

save(list = "sg_sequences_df",
     file = file.path(R_objects_directory, "3) Import and process sgRNA sequences.RData")
     )




