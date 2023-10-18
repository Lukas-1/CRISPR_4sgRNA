### 2022-06-19



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "07_exporting_reads.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
rdata_dir   <- file.path(project_dir, "03_R_objects")
fastq_dir   <- file.path(project_dir, "04_output_data", "Fastq")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))

load(file.path(rdata_dir, "01_read_in_data_run1.RData"))
load(file.path(rdata_dir, "01_read_in_data_run2_chunk1.RData"))
load(file.path(rdata_dir, "01_read_in_data_run2_chunk2.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run1.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk1.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk2.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__lumi_df.RData"))



# Create a combined matched_df --------------------------------------------

matched_df <- rbind.data.frame(data.frame("Run" = 1L, run1_matched_df,
                                          run1_fastq_df[, c("Quality_sg1", "Quality_sg2")]
                                          ),
                               data.frame("Run" = 2L, run2_chunk1_matched_df,
                                          run2_chunk1_fastq_df[, c("Quality_sg1", "Quality_sg2")]
                                          ),
                               data.frame("Run" = 2L, run2_chunk2_matched_df,
                                          run2_chunk2_fastq_df[, c("Quality_sg1", "Quality_sg2")]
                                          ),
                               stringsAsFactors = FALSE,
                               make.row.names = FALSE
                               )
matched_df <- data.frame("Read_number" = seq_len(nrow(matched_df)),
                         matched_df, stringsAsFactors = FALSE
                         )
rm(list = c("run1_matched_df", "run2_chunk1_matched_df", "run2_chunk2_matched_df"))
rm(list = c("run1_fastq_df", "run2_chunk1_fastq_df", "run2_chunk2_fastq_df"))
gc()




# Do stuff ----------------------------------------------------------------

only_has_sg1 <- (counts_df[, "Data_contains_sg1"] & (!(counts_df[, "Data_contains_sg2"])))
only_has_sg2 <- (counts_df[, "Data_contains_sg2"] & (!(counts_df[, "Data_contains_sg1"])))

only_sg1_plasmids <- counts_df[, "Plasmid_ID"][only_has_sg1]
only_sg2_plasmids <- counts_df[, "Plasmid_ID"][only_has_sg2]

set.seed(1)
five_plasmids <- sample(only_sg1_plasmids, 5)



AllReadsForPlasmid <- function(plasmid_ID, sg_number, exclude_template_switch = FALSE) {
  stopifnot(all(c("lumi_df", "matched_df") %in% ls(envir = globalenv())))
  use_regex <- paste0(plasmid_ID, "(, |$)")
  are_to_select <- grepl(use_regex, lumi_df[, paste0("Plasmid_sg", sg_number)])
  if (exclude_template_switch) {
    are_to_select[are_to_select] <- !(lumi_df[, "Has_template_switch"][are_to_select])
  }
  matches_vec <- match(lumi_df[, "Read_number"][are_to_select],
                       matched_df[, "Read_number"]
                       )
  results_df <- data.frame(lumi_df[are_to_select, ],
                           matched_df[matches_vec, c("Quality_sg1", "Quality_sg2")],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )
  print(nrow(results_df))
  return(results_df)
}



ReadLevelFASTQ <- function(input_df) {
  read_separator <- paste0(rep("N", 10), collapse = "")
  quality_separator <- paste0(rep("I", 10), collapse = "")
  sequence_vec <- paste0(input_df[, "Aligned_read_sg1"], read_separator, input_df[, "Aligned_read_sg2"])
  quality_vec <- paste0(input_df[, "Quality_sg1"], quality_separator, input_df[, "Quality_sg2"])
  quality_object <- PhredQuality(quality_vec)
  fastq_object <- QualityScaledBStringSet(sequence_vec, quality_object)
  sample_names_vec <- sub("(_off)?_2sg_", "_", input_df[, "Sample_name"])
  read_names_vec <- paste0(sample_names_vec, "_run", input_df[, "Run"],
                           "-", input_df[, "Read_number"]
                           )
  names(fastq_object) <- read_names_vec
  return(fastq_object)
}



UnmatchedSecondGuide <- function(plasmid_ID, unmatched_sg) {
  stopifnot(unmatched_sg %in% 1:2)
  stopifnot("lumi_df" %in% ls(envir = globalenv()))
  if (unmatched_sg == 2) {
    matched_sg <- 1L
  } else {
    matched_sg <- 2L
  }
  use_regex <- paste0(plasmid_ID, "(, |$)")
  are_to_select <- grepl(use_regex, lumi_df[, paste0("Plasmid_sg", matched_sg)])
  are_to_select[are_to_select] <- !(lumi_df[, "Has_template_switch"][are_to_select])
  unmatched_table <- table(lumi_df[, paste0("Aligned_read_sg", unmatched_sg)][are_to_select])
  new_order <- order(as.integer(unmatched_table), decreasing = TRUE)
  results_df <- data.frame(
    "Plasmid_ID"      = plasmid_ID,
    "Unmatched_sg"    = unmatched_sg,
    "Sequence"        = names(unmatched_table)[new_order],
    "Num_occurrences" = as.integer(unmatched_table)[new_order],
    stringsAsFactors  = FALSE
  )
  return(results_df)
}



AlignWithLibrary <- function(use_sequence, sg_number, top_n = NULL) {
  stopifnot(unmatched_sg %in% 1:2)
  stopifnot("CRISPRoff_df" %in% ls(envir = globalenv()))
  sg_column <- paste0("protospacer_", c("A", "B")[[sg_number]])
  sg_vec <- toupper(substr(CRISPRoff_df[, sg_column], 2, 20))
  results_object <- pairwiseAlignment(sg_vec, use_sequence)
  new_order <- order(score(results_object), decreasing = TRUE)
  if (!(is.null(top_n))) {
    new_order <- new_order[seq_len(top_n)]
  }
  results_df <- data.frame(
    "Plasmid_ID"         = CRISPRoff_df[, "Plasmid_ID"][new_order],
    "Sg_number"          = sg_number,
    "Library_reference"  = sg_vec[new_order],
    "Aligned_library"    = as.character(alignedPattern(results_object[new_order])),
    "Aligned_sequence"   = as.character(alignedSubject(results_object[new_order])),
    "Score"              = score(results_object)[new_order],
    stringsAsFactors     = FALSE
  )
  return(results_df)
}






SummaryFasta <- function(input_df, sg_number = 2) {

}



try_df <- AllReadsForPlasmid(five_plasmids[[1]], 1, exclude_template_switch = TRUE)
try2_object <- PrepareFASTQ(try_df)


try3_df <- AllReadsForPlasmid("ACTR6", 1, exclude_template_switch = TRUE)
try3_object <- PrepareFASTQ(try3_df)


try4_df <- UnmatchedSecondGuide(five_plasmids[[1]], 2)
try5_df <- AlignWithLibrary("TGCGCCTGCGCGTCCCCGG", 2, top_n = 5)



try4_df <- UnmatchedSecondGuide(five_plasmids[[2]], 2)
try5_df <- AlignWithLibrary("GAGGTGAATCCGGTCAGAA", 2, top_n = 5)






ExportReferenceForPlasmid <- function(plasmid_ID) {
  stopifnot(all(c("fastq_dir", "CRISPRoff_df") %in% ls(envir = globalenv())))
  plasmid_index <- match(plasmid_ID, CRISPRoff_df[, "Plasmid_ID"])
  stopifnot(!(is.na(plasmid_index)))
  use_separator <- paste0(rep("N", 10), collapse = "")
  combined_sequence <- paste0(substr(CRISPRoff_df[, "protospacer_A"][[plasmid_index]], 2, 20),
                              use_separator,
                              substr(CRISPRoff_df[, "protospacer_B"][[plasmid_index]], 2, 20),
                              )
  fasta_vec <- c(paste(">", plasmid_ID), combined_sequence)
  file_name <- paste0(plasmid_ID, "__reference.fasta")
  write.table(fasta_vec, file = file.path(fastq_dir, file_name),
              col.names = FALSE, row.names = FALSE, quote = FALSE
              )

}

ExportReferenceForPlasmid(five_plasmids[[1]])

writeQualityScaledXStringSet(try2_object,
                             filepath = file.path(fastq_dir, paste0(plasmid_ID, "__reads.fastq"))
                             )


ExportReferenceForPlasmid("ACTR6")
writeQualityScaledXStringSet(try3_object,
                             filepath = file.path(fastq_dir, paste0("ACTR6", "__reads.fastq"))
                             )




