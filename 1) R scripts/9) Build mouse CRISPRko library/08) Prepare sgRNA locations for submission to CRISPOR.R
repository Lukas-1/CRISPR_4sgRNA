### 21st February 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "6) Mouse - General")
CRISPRko_RData_directory <- file.path(RData_directory, "8) Mouse - CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "Mouse - CRISPRko", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "04) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData")) # The output of "07) Integrate the output from GuideScan" is not yet needed





# Include sgRNAs with missing Entrez IDs ----------------------------------

chunks_list <- AppendIDsWithoutCanonicalEntrezs(entrez_chunks_list, extended_CRISPRko_df)




# Prepare data frames that can be exported to .bed files ------------------

bed_df_list <- BreakIntoChunks(MakeBedDf, extended_CRISPRko_df, chunks_list)





# Prepare objects that can be exported to FASTA files ---------------------

FASTA_df_list <- BreakIntoChunks(MakeFASTADf, extended_CRISPRko_df, chunks_list)
FASTA_df_list <- CombineDfChunks(FASTA_df_list)
FASTA_vec_list <- lapply(FASTA_df_list, MakeFASTAvec)





# Write input files for CRISPOR to disk -----------------------------------

WriteCRISPORInputFiles(bed_df_list,
                       file_ending = "__CRISPRko.bed",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_bed")
                       )

WriteCRISPORInputFiles(FASTA_vec_list, file_ending = "__CRISPRko.fa",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_FASTA")
                       )




# Write input files (filtered by already present data) --------------------

filtered_bed_df_list <- FilterBedDfList(bed_df_list)
filtered_FASTA_df_list <- FilterFASTADfList(FASTA_df_list)
filtered_FASTA_vec_list <- lapply(filtered_FASTA_df_list, MakeFASTAvec)

WriteCRISPORInputFiles(filtered_bed_df_list, file_ending = "__CRISPRko.bed",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_bed_filtered")
                       )

WriteCRISPORInputFiles(filtered_FASTA_vec_list, file_ending = "__CRISPRko.fa",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_FASTA_filtered")
                       )






