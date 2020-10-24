### 20th October 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
file_input_directory   <- file.path(plate2_directory, "1) Input")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")
file_output_directory  <- file.path(plate2_directory, "3) Output")

metadata_directory     <- file.path(file_input_directory, "Metadata")
sg_sequences_file      <- file.path(metadata_directory, "sgRNA reference.xlsx")




# Read in data ------------------------------------------------------------

gRNAs_excel_df <- as.data.frame(read_excel(sg_sequences_file),
                                stringsAsFactors = FALSE, check.names = FALSE
                                )



# Define maps -------------------------------------------------------------

modalities_map <- c(
  "a1" = "CRISPRa",
  "Ca" = "CRISPRa",
  "Ci" = "CRISPRi",
  "i1" = "CRISPRi",
  "ko" = "CRISPRo",
  "o1" = "CRISPRo"
)



# Define well numbers -----------------------------------------------------

wells_df <- expand.grid(1:24, LETTERS[1:16])[2:1]
wells_map <- seq_len(384)
names(wells_map) <- paste0(wells_df[[1]], wells_df[[2]])




# Process guides ----------------------------------------------------------

is_empty_row <- apply(gRNAs_excel_df, 1, function(x) all(is.na(x)))

block_vec <- rep(1L, nrow(gRNAs_excel_df))
block_vec[seq_along(block_vec) > which(is_empty_row)] <- 2L

block_vec <- block_vec[!(is_empty_row)]
gRNAs_excel_df <- gRNAs_excel_df[!(is_empty_row), ]

plasmid_names_splits <- strsplit(gRNAs_excel_df[["Plasmid ID"]], " ", fixed = TRUE)
have_space <- lengths(plasmid_names_splits) > 1
plasmid_names_splits[have_space] <- lapply(plasmid_names_splits[have_space],
                                           function(x) {
                                             if (length(x) == 3) {
                                               c(paste0(x[[1]], " ", x[[2]]), x[[3]])
                                             } else {
                                               x
                                             }
                                           })

plasmid_names_splits[!(have_space)] <- strsplit(gRNAs_excel_df[["Plasmid ID"]][!(have_space)],
                                                "_",
                                                fixed = TRUE
                                                )

plasmid_names <- sapply(plasmid_names_splits, "[[", 1)
raw_modalities <- sapply(plasmid_names_splits, "[[", 2)

guide_seq_mat <- as.matrix(gRNAs_excel_df[, paste0("sg", 1:4, " (5'-3')")])
colnames(guide_seq_mat) <- paste0("Sequence_sg", 1:4)

sg_sequences_df <- data.frame(
  "Well_number"         = wells_map[gRNAs_excel_df[["Wells"]]],
  "Well_code"           = gRNAs_excel_df[["Wells"]],
  "Block"               = block_vec,
  "Target_gene"         = plasmid_names,
  "Modality"            = modalities_map[raw_modalities],
  guide_seq_mat,
  "Longest_subsequence" = apply(guide_seq_mat, 1, LongestSharedSubsequence),
  stringsAsFactors      = FALSE,
  row.names             = NULL
)



# Check for plasmids with guides <20 bp -----------------------------------

have_short_guides <- apply(guide_seq_mat, 1, function(x) any(nchar(x) != 20))
t(apply(guide_seq_mat[have_short_guides, ], 1, nchar))
sg_sequences_df[have_short_guides, ]




# Export the sgRNA sequences data frame (including well number) -----------

write.table(sg_sequences_df,
            file      = file.path(file_output_directory, paste0("sgRNAs_and_wells.tsv")),
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE
            )



# Save data ---------------------------------------------------------------

save(list = "sg_sequences_df",
     file = file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData")
     )




