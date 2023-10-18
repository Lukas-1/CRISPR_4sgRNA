### 20th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data.RData"))
load(file.path(p2_R_objects_directory, "04) Perform pairwise alignments with the reference sequence.RData"))




# Define all important features, and generate per-well coordinates --------

features_df <- FeaturesListToDf(features_list)

empty_list <-  vector(mode = "list", length = 384)
features_mat_list <- empty_list
features_indices_list <- empty_list
features_templates_list <- empty_list

wells_vec <- sg_sequences_df[["Well_number"]]

features_mat_list[wells_vec] <- lapply(seq_along(wells_vec), function(x) {
  sg_vec <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  AdjustForsgRNALength(features_df, sg_vec)
})

features_indices_list[wells_vec] <- lapply(features_mat_list[wells_vec], function(x) {
  lapply(seq_len(nrow(x)), function(y) seq(from = x[y, "Start"], to = x[y, "End"]))
})

features_templates_list[wells_vec] <- lapply(seq_along(wells_vec), function(x) {
  use_mat <- features_mat_list[[wells_vec[[x]]]]
  vapply(seq_len(nrow(use_mat)),
         function(y) substr(sg_sequences_df[x, "Barcoded_plasmid"], use_mat[y, "Start"], use_mat[y, "End"]),
         ""
         )
})




# Extract sequences -------------------------------------------------------

sl7_extracted_df <- ExtractAlignedSequences(sl7_ccs_df, sl7_alignments_df, unique_IDs = sg_sequences_df[["Well_number"]])
sl9_extracted_df <- ExtractAlignedSequences(sl9_ccs_df, sl9_alignments_df, unique_IDs = sg_sequences_df[["Well_number"]])





# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_extracted_df"),
     file = file.path(p2_R_objects_directory, "06) Categorize subsequences of reads aligned to the reference.RData")
     )







