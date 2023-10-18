## 2023-10-07


# Define paths ------------------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")
rdata_dir <- file.path(project_dir, "03_R_objects")




# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "05_obtain_read_counts.RData"))
for (i in 1:8) {
  load(file.path(rdata_dir, paste0("03_look_up_sgRNAs_chunk", i, ".RData")))
}



# Create a combined matched_df --------------------------------------------

matched_df_chunks <- paste0("chunk", 1:8, "_matched_df")
matched_df <- do.call(rbind.data.frame,
                      c(lapply(matched_df_chunks, get),
                        stringsAsFactors = FALSE, make.row.names = FALSE
                        ))
rm(list = matched_df_chunks)




# Check whether sgRNAs from the library are found in the data -------------

all_sg2_vec <- unique(ifelse(matched_df[, "Num_MM_sg2"] == 0,
                             matched_df[, "Sequence_sg2"],
                             matched_df[, "Correct_sgRNA_sg2"]
                             ))
all_sg3_vec <- unique(ifelse(matched_df[, "Num_MM_sg3"] == 0,
                             matched_df[, "Sequence_sg3"],
                             matched_df[, "Correct_sgRNA_sg3"]
                             ))

found_sg2 <- toupper(sg_sequences_df[, "Sequence_sg2"]) %in% all_sg2_vec
found_sg3 <- toupper(sg_sequences_df[, "Sequence_sg3"]) %in% all_sg3_vec





# Extend counts_df --------------------------------------------------------

counts_df <- data.frame(
  counts_df[, 1:3],
  "Data_contains_sg2" = found_sg2,
  "Data_contains_sg3" = found_sg3,
  counts_df[, 4:ncol(counts_df)],
  stringsAsFactors = FALSE
)



# Save data ---------------------------------------------------------------

save(list = "counts_df",
     file = file.path(rdata_dir, "06_extend_read_counts.RData")
     )


