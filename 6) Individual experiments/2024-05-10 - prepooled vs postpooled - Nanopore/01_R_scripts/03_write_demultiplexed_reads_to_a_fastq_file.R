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
  "Name"         = stripped_files,
  "Path"         = file.path(rdata_dir, chunk_files)
)



# Export demultiplex reads as FASTQ files ---------------------------------

for (file_number in unique(chunks_df[, "File_number"])) {

  message("Loading the reads for file ", file_number, "...")

  load_paths <- chunks_df[, "Path"][chunks_df[, "File_number"] %in% file_number]
  object_names <- vapply(load_paths, load, "", envir = globalenv(), USE.NAMES = FALSE)


  message("Preparing data...")

  demuxed_df <- data.table::rbindlist(lapply(object_names, get))
  data.table::setDF(demuxed_df)
  rm(list = object_names)

  demuxed_df[, "File_number"] <- file_number
  object_name <- paste0("file_", file_number, "_meta_df")
  retain_columns <- c("File_number", "Replicate", "Condition", "Is_fwd", "Mean_quality")
  assign(object_name, demuxed_df[, retain_columns], envir = globalenv())

  fastq_file_name <- paste0("file_", file_number, ".fastq.gz")

  export_seq <- DNAStringSet(demuxed_df[, "Sequence"])
  names(export_seq) <- paste0("read_", seq_len(nrow(demuxed_df)))


  message("Creating the QualityScaledBStringSet object for file ", file_number, "...")

  fastq_object <- QualityScaledBStringSet(export_seq, PhredQuality(demuxed_df[, "Quality"]))


  message("Writing the file to disk: ", fastq_file_name, "...")

  writeQualityScaledXStringSet(fastq_object,
                               filepath = file.path(export_dir, fastq_file_name),
                               compress = TRUE
                               )

}



# Save metadata -----------------------------------------------------------

meta_df_names <- grep("^file_[0-9]+_meta_df$", ls(), value = TRUE)
demuxed_meta_df <- data.table::rbindlist(lapply(meta_df_names, get))
data.table::setDF(demuxed_meta_df)



# Save data ---------------------------------------------------------------

save(demuxed_meta_df,
     file = file.path(rdata_dir, "03_demuxed_meta_df.RData")
     )


