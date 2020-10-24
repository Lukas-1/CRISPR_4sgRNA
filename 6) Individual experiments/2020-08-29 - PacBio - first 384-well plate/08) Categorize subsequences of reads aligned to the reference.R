### 15th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "07) Perform pairwise alignments with the reference sequence.RData"))





# Define all important features, and generate per-well coordinates --------

features_mat <- do.call(rbind, features_list)
mode(features_mat) <- "integer"
rownames(features_mat) <- NULL
colnames(features_mat) <- c("Start", "End")
features_df <- data.frame("Feature" = names(features_list),
                          features_mat,
                          stringsAsFactors = FALSE
                          )

features_mat_list <- lapply(seq_len(384), function(x) {
  sg_vec <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  AdjustForsgRNALength(features_df, sg_vec)
})

features_indices_list <- lapply(features_mat_list, function(x) {
  lapply(seq_len(nrow(x)), function(y) seq(from = x[y, "Start"], to = x[y, "End"]))
})


barcoded_plasmids <- paste0(column_bc_vec, plasmids_vec, row_bc_vec)

features_templates_list <- lapply(seq_len(384), function(x) {
  use_mat <- features_mat_list[[x]]
  vapply(seq_len(nrow(use_mat)),
         function(y) substr(barcoded_plasmids[[x]], use_mat[y, "Start"], use_mat[y, "End"]),
         ""
         )
})



# Extract sequences -------------------------------------------------------

sl7_extracted_df <- ExtractAlignedSequences(use_sl7 = TRUE)
sl9_extracted_df <- ExtractAlignedSequences(use_sl7 = FALSE)





# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_extracted_df"),
     file = file.path(R_objects_directory, "09) Categorize subsequences of reads aligned to the reference.RData")
     )








