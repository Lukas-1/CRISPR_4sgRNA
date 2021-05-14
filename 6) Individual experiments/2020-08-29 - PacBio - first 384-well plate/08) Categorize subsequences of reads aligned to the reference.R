### 15th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))





# Define all important features, and generate per-well coordinates --------

features_df <- FeaturesListToDf(features_list)

features_mat_list <- lapply(seq_len(384), function(x) {
  sg_vec <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  AdjustForsgRNALength(features_df, sg_vec)
})

features_indices_list <- lapply(features_mat_list, function(x) {
  lapply(seq_len(nrow(x)), function(y) seq(from = x[y, "Start"], to = x[y, "End"]))
})

features_templates_list <- lapply(seq_len(384), function(x) {
  use_mat <- features_mat_list[[x]]
  vapply(seq_len(nrow(use_mat)),
         function(y) substr(sg_sequences_df[x, "Barcoded_plasmid"], use_mat[y, "Start"], use_mat[y, "End"]),
         ""
         )
})



# Extract sequences -------------------------------------------------------

sl7_extracted_df <- ExtractAlignedSequences(sl7_ccs_df, sl7_alignments_df)
sl9_extracted_df <- ExtractAlignedSequences(sl9_ccs_df, sl9_alignments_df)





# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_extracted_df"),
     file = file.path(R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData")
     )








