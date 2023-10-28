## 2023-10-18


# Load packages and source code -------------------------------------------

library("Biostrings")
library("data.table")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2023-10-05 - prepooled vs postpooled - Nanopore")
rdata_dir             <- file.path(project_dir, "03_R_objects")
export_dir            <- file.path(project_dir, "04_intermediate_files", "1_output_from_R")



# Load data ---------------------------------------------------------------

rdata_files <- file.path(rdata_dir, paste0("01_demultiplex_reads__chunk_", 1:6, ".RData"))
object_names <- vapply(rdata_files, load, "", envir = globalenv(), USE.NAMES = FALSE)



# Export demultiplexed FASTQ reads ----------------------------------------

demuxed_df <- data.table::rbindlist(lapply(object_names, get))
data.table::setDF(demuxed_df)

export_seq <- DNAStringSet(demuxed_df[, "Sequence"])
names(export_seq) <- paste0("read_", seq_len(nrow(demuxed_df)))

fastq_object <- QualityScaledBStringSet(export_seq, PhredQuality(demuxed_df[, "Quality"]))

writeQualityScaledXStringSet(fastq_object,
                             filepath = file.path(export_dir, "All_demuxed_reads.fastq.gz"),
                             compress = TRUE
                             )



# Save metadata -----------------------------------------------------------

demuxed_meta_df <- demuxed_df[, c("Replicate", "Condition", "Is_fwd", "Mean_quality")]



# Save data ---------------------------------------------------------------

save(demuxed_meta_df,
     file = file.path(rdata_dir, "02_demuxed_meta_df.RData")
     )


