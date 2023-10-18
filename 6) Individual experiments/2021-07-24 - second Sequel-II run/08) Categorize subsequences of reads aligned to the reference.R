### 30th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r2_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(s2r2_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))





# Define all important features, and generate per-well coordinates --------

features_df <- FeaturesListToDf(features_list)
for (column in c("Start", "End")) {
  features_df[[column]] <- features_df[[column]] + 17L
}

stopifnot(all(nchar(as.matrix(sg_sequences_df[, paste0("Sequence_sg", 1:4)])) == 20))

all_indices <- seq_len(nrow(sg_sequences_df))

features_mat_list <- lapply(all_indices, function(x) {
  AdjustForsgRNALength(features_df, as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)]))
})
features_indices_list <- lapply(features_mat_list, function(x) {
  lapply(seq_len(nrow(x)), function(y) seq(from = x[y, "Start"], to = x[y, "End"]))
})
features_templates_list <- lapply(all_indices, function(x) {
  use_mat <- features_mat_list[[x]]
  vapply(seq_len(nrow(use_mat)),
         function(y) substr(sg_sequences_df[x, "Barcoded_plasmid"], use_mat[y, "Start"], use_mat[y, "End"]),
         ""
         )
})




# Extract sequences -------------------------------------------------------

extracted_df <- ExtractAlignedSequences(ccs_df,
                                        alignments_df,
                                        ID_column = "Combined_ID",
                                        unique_IDs = sg_sequences_df[["Combined_ID"]]
                                        )



# Save data ---------------------------------------------------------------

save(list = c("extracted_df",
              "features_df", "features_indices_list", "features_templates_list"
              ),
     file = file.path(s2r2_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData")
     )




