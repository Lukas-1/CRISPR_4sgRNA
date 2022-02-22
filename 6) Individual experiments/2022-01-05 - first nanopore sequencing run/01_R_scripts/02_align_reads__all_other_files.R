## 2022-01-04


# Load packages and source code -------------------------------------------

project_dir <- file.path("~", "NP1_4sg")
source(file.path(project_dir, "01_R_scripts", "1_R_functions", "01_aligning_reads.R"))



# Read in data ------------------------------------------------------------

fastq_files <- setdiff(list.files(reads_dir), "4sg_1.fastq")
fastq_reads <- ShortRead::readFastq(file.path(reads_dir, fastq_files))



# Align reads -------------------------------------------------------------

alignments_other_df <- ParallelAlignInChunks(fastq_reads)



# Save data ---------------------------------------------------------------

save(alignments_other_df, file = file.path(rdata_dir, "02_align_reads__all_other_files.RData"))



