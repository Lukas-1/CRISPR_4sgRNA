### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
CRISPOR_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRa", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRa_RData_directory, "11) Refine the genomic locations of sgRNA sequences after fixing 5'G substitutions.RData")) # The output of files 12-15) is not necessary for this step





# Include sgRNAs with missing Entrez IDs ----------------------------------

chunks_list <- AppendIDsWithoutEntrezs(entrez_chunks_list, merged_replaced_CRISPRa_df)




# Prepare data frames that can be exported to .bed files ------------------

bed_df_list <- lapply(chunks_list, function(x) MakeBedDf(merged_replaced_CRISPRa_df, x))




# Prepare objects that can be exported to FASTA files ---------------------

FASTA_df_list <- lapply(chunks_list, function(x) MakeFASTADf(merged_replaced_CRISPRa_df, x))
FASTA_df_list <- CombineDfChunks(FASTA_df_list)
FASTA_vec_list <- lapply(FASTA_df_list, MakeFASTAvec)




# Write input files for CRISPOR to disk -----------------------------------

WriteCRISPORInputFiles(bed_df_list, file_ending = "__CRISPRa.bed",
                       CRISPOR_input_directory = CRISPOR_files_directory
                       )

WriteCRISPORInputFiles(FASTA_vec_list, file_ending = "__CRISPRa.fa",
                       CRISPOR_input_directory = CRISPOR_files_directory
                       )



# Write input files (filtered by already present data) --------------------

filtered_bed_df_list    <- FilterBedDfList(bed_df_list)
filtered_FASTA_df_list  <- FilterFASTADfList(FASTA_df_list)
filtered_FASTA_vec_list <- lapply(filtered_FASTA_df_list, MakeFASTAvec)

filtered_input_directory <- file.path(CRISPOR_files_directory, "Input - filtered by already processed")

WriteCRISPORInputFiles(filtered_bed_df_list, file_ending = "__CRISPRa.bed",
                       CRISPOR_input_directory = filtered_input_directory
                       )

WriteCRISPORInputFiles(filtered_FASTA_vec_list, file_ending = "__CRISPRa.fa",
                       CRISPOR_input_directory = filtered_input_directory
                       )













