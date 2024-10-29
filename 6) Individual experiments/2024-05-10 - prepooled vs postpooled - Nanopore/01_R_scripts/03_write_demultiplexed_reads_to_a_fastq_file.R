## 2024-06-01


# Load packages and source code -------------------------------------------

library("Biostrings")
library("data.table")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
rdata_dir             <- file.path(project_dir, "03_R_objects")
export_dir            <- file.path(project_dir, "04_intermediate_files", "1_output_from_R")




# Define functions --------------------------------------------------------

DivideIntoParts <- function(file_number_vec, max_per_part = 6) {
  stopifnot(length(unique(file_number_vec)) == 1)
  if (length(file_number_vec) <= max_per_part) {
    paste0("file_", file_number_vec)
  } else {
    part_number_vec <- rep(seq_len(ceiling(length(file_number_vec) / max_per_part)),
                           each = max_per_part
                           )[seq_along(file_number_vec)]
    paste0("file_", file_number_vec, "_part_", part_number_vec)
  }
}



# Identify demultiplexed chunks of data -----------------------------------

chunk_files <- grep("^02_demultiplex_reads__file_[0-9]+_chunk_[0-9]+\\.RData$", list.files(rdata_dir), value = TRUE)

stripped_files <- sub("02_demultiplex_reads__", "", chunk_files, fixed = TRUE)
stripped_files <- sub(".RData", "", stripped_files, fixed = TRUE)

file_splits <- strsplit(stripped_files, "_", fixed = TRUE)
file_numbers <- as.integer(sapply(file_splits, "[[", 2))
chunk_numbers <- as.integer(sapply(file_splits, "[[", 4))
chunks_df <- data.frame(
  "File_number"  = as.integer(sapply(file_splits, "[[", 2)),
  "Chunk_number" = chunk_numbers,
  "Chunk_name"   = stripped_files,
  "Chunk_path"   = file.path(rdata_dir, chunk_files)
)
new_order <- order(chunks_df[, "File_number"], chunks_df[, "Chunk_number"])
chunks_df <- chunks_df[new_order, ]
row.names(chunks_df) <- NULL
chunks_df[, "Part_name"] <- unlist(tapply(chunks_df[, "File_number"], chunks_df[, "File_number"], DivideIntoParts), use.names = FALSE)
chunks_df[, "Part_number"] <- match(chunks_df[, "Part_name"], unique(chunks_df[, "Part_name"]))



# Export demultiplex reads as FASTQ files ---------------------------------

for (part_name in unique(chunks_df[, "Part_name"])) {

  part_index <- match(part_name, chunks_df[, "Part_name"])
  file_number <- chunks_df[, "File_number"][[part_index]]
  part_number <- chunks_df[, "Part_number"][[part_index]]

  message("Loading the reads for ", part_name, "...")

  load_paths <- chunks_df[, "Chunk_path"][chunks_df[, "Part_name"] %in% part_name]
  object_names <- vapply(load_paths, load, "", envir = globalenv(), USE.NAMES = FALSE)


  message("Preparing data...")

  demuxed_df <- data.table::rbindlist(lapply(object_names, get))
  data.table::setDF(demuxed_df)
  rm(list = object_names)

  demuxed_df[, "File_number"] <- file_number
  demuxed_df[, "Part_number"] <- part_number
  object_name <- paste0(part_name, "_meta_df")
  retain_columns <- c("File_number", "Part_number", "Replicate", "Condition", "Is_fwd", "Mean_quality")
  assign(object_name, demuxed_df[, retain_columns], envir = globalenv())

  fastq_file_name <- paste0(part_name, ".fastq.gz")

  export_seq <- DNAStringSet(demuxed_df[, "Sequence"])
  names(export_seq) <- paste0("p", part_number, "_", seq_len(nrow(demuxed_df)))


  message("Creating the QualityScaledBStringSet object for ", fastq_file_name, "...")

  fastq_object <- QualityScaledBStringSet(export_seq, PhredQuality(demuxed_df[, "Quality"]))


  message("Writing the file to disk: ", fastq_file_name, "...")

  writeQualityScaledXStringSet(fastq_object,
                               filepath = file.path(export_dir, fastq_file_name),
                               compress = TRUE
                               )
}



# Save metadata -----------------------------------------------------------

meta_df_names <- grep("^file_[0-9]+_(part_[0-9]+_)?meta_df$", ls(), value = TRUE)
demuxed_meta_df <- data.table::rbindlist(lapply(meta_df_names, get))
data.table::setDF(demuxed_meta_df)



# Save data ---------------------------------------------------------------

save(demuxed_meta_df,
     file = file.path(rdata_dir, "03_demuxed_meta_df.RData")
     )


