## 2022-02-03


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

project_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")

source(file.path(project_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_combine_aligned_reads.RData"))




# Prepare data ------------------------------------------------------------

alignments_df[, "Mean_quality"] <- GetMeanQuality(alignments_df[, "Read_sequence"])

features_df <- FeaturesListToDf(features_list)

for (column in c("Start", "End")) {
  features_df[[column]] <- features_df[[column]] - 10L
}

features_df <- features_df[features_df[, "Feature"] %in% paste0("sg", 1:4), ]

features_indices_list <- lapply(seq_len(nrow(features_df)),
                                function(y) seq(from = features_df[y, "Start"], to = features_df[y, "End"])
                                )




# Extract aligned sequences -----------------------------------------------

num_reads <- nrow(alignments_df)
reads_per_chunk <- 10000
num_chunks <- ceiling(num_reads / reads_per_chunk)
chunks_vec <- rep(seq_len(num_chunks), each = reads_per_chunk)[seq_len(num_reads)]
chunks_list <- vector(mode = "list", length = num_chunks)
first_vec <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[1]]))
last_vec  <- format(tapply(seq_len(num_reads), chunks_vec, function(x) x[[length(x)]]))
chunk_numbers <- format(seq_len(num_chunks))

for (i in seq_len(num_chunks)) {
  are_this_chunk <- chunks_vec == i
  message("Processing chunk #", chunk_numbers[[i]], " of ",
          chunk_numbers[[length(chunk_numbers)]],  " (extracting reads ",
          first_vec[[i]], " to ", last_vec[[i]], ")..."
          )
  sub_df <- ExtractAlignedSequences(alignments_df[are_this_chunk, ])
  chunks_list[[i]] <- sub_df
}

extracted_df <- do.call(rbind.data.frame,
                        c(chunks_list,
                          stringsAsFactors = FALSE,
                          make.row.names = FALSE
                          )
                        )
extracted_df[, "Sg_number"] <- as.integer(sub("sg", "", extracted_df[, "Feature"], fixed = TRUE))

extracted_df_list <- split(extracted_df[, !(names(extracted_df) %in% c("Feature", "Sg_number"))],
                           extracted_df[, "Sg_number"]
                           )

for (i in 1:4) {
  names(extracted_df_list[[i]]) <- paste0(names(extracted_df_list[[i]]), "_sg", i)
}

extracted_df <- data.frame(
  alignments_df[, c("Orientation_fwd", "Score_fwd", "Score_rev", "Mean_quality")],
  extracted_df_list[[1]], extracted_df_list[[2]],
  extracted_df_list[[3]], extracted_df_list[[4]]
)




# Save data ---------------------------------------------------------------

save(extracted_df, file = file.path(rdata_dir, "05_extract_aligned_sgRNAs.RData"))




