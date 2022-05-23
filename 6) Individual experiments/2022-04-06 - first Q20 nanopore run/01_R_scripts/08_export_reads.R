## 2022-04-17



# Load packages and source code -------------------------------------------

library("Biostrings")
CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "07_exporting_reads.R"))



# Define paths ------------------------------------------------------------

project_dir   <- file.path(experiments_directory, "2022-04-06 - first Q20 nanopore run")
input_dir     <- file.path(project_dir, "02_input_data")
rdata_dir     <- file.path(project_dir, "03_R_objects")
output_dir    <- file.path(project_dir, "04_output_data")
figures_dir   <- file.path(output_dir, "Figures")
fastq_dir     <- file.path(output_dir, "Fastq")
reference_dir <- file.path(output_dir, "Reference sequences")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids.RData"))



# Read in data ------------------------------------------------------------

fastq_file_name <- "20220330_1138_X3_FAS18220_5478982b.pass.fastq.gz"
fastq_reads <- ShortRead::readFastq(file.path(input_dir, fastq_file_name))



# Create a data frame of raw reads ----------------------------------------

raw_df <- data.frame(
  "Sequence" = as.character(ShortRead::sread(fastq_reads)),
  "Quality"  = as.character(as(Biostrings::quality(fastq_reads), "PhredQuality")),
  stringsAsFactors = FALSE
)



# Create short descriptors for individual reads ---------------------------

raw_df[, "Read_string"] <- CreateReadDescriptors(raw_df, nano_df)




# Create short descriptors for plasmids -----------------------------------

CRISPRa_df <- TidyCRISPRaDf(sg_sequences_df)




# Insert gRNA sequences into the reference 4sg cassette -------------------

library_seq <- seq_len(nrow(CRISPRa_df))
plasmid_lines_list <- lapply(library_seq, function(x) {
  use_lines <- plasmid_gbk[["sequence"]]
  sg_seqs <- as.character(CRISPRa_df[x, paste0("Sequence_sg", 1:4)])
  ReplaceNNNLines(use_lines, sg_seqs)
})




# Define genes of interest ------------------------------------------------

set.seed(1)
have_one_TSS <- CRISPRa_df[, "Num_plasmids_for_gene"] %in% 1
random_genes <- unique(CRISPRa_df[, "Gene_symbol"][sample(which(have_one_TSS), 10)])[1:5]
other_genes <- c("PRNP")
all_genes <- c(random_genes, other_genes)




# Export reference sequences for genes of interest ------------------------

for (each_gene in all_genes) {
  ExportReferences(each_gene)
}



# Export reads ------------------------------------------------------------

set.seed(1)
chunk_size <- 50L
num_chunks <- 5L
random_indices <- sample(seq_len(nrow(raw_df)), size = chunk_size * num_chunks)
random_indices_list <- split(random_indices, rep(seq_len(num_chunks), each = chunk_size))

for (i in seq_len(num_chunks)) {
  these_reads <- SelectReads(raw_df, random_indices_list[[i]])
  ExportReads(these_reads, fastq_dir, file_name = paste0("50_random_reads_set", i))
}

for (each_gene in all_genes) {
  print(each_gene)
  ExportReadsForGene(each_gene, raw_df, nano_df)
}




