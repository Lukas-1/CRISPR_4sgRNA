### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")
CRISPOR_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRi", "CRISPOR")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRi_RData_directory, "11) Refine the genomic locations of sgRNA sequences after fixing 5'G substitutions.RData")) # The output of files 12-15) is not necessary for this step





# Include sgRNAs with missing Entrez IDs ----------------------------------

chunks_list <- AppendIDsWithoutCanonicalEntrezs(entrez_chunks_list, merged_replaced_CRISPRi_df)





# Prepare data frames that can be exported to .bed files ------------------

bed_df_list <- BreakIntoChunks(MakeBedDf, merged_replaced_CRISPRi_df, chunks_list)
bed_df_list <- c(BreakIntoManageableChunks(bed_df_list[c(1:2, 4:7)],
                                           80050L,
                                           MakeBedDf,
                                           merged_replaced_CRISPRi_df,
                                           chunks_list[c(1:2, 4:7)]
                                           ),
                 BreakIntoManageableChunks(bed_df_list[c(3, 8:13)],
                                           70000L,
                                           MakeBedDf,
                                           merged_replaced_CRISPRi_df,
                                           chunks_list[c(3, 8:13)]
                                           )
                 )





# Prepare objects that can be exported to FASTA files ---------------------

FASTA_df_list <- BreakIntoChunks(MakeFASTADf, merged_replaced_CRISPRi_df, chunks_list)
FASTA_df_list <- CombineDfChunks(FASTA_df_list)
FASTA_vec_list <- lapply(FASTA_df_list, MakeFASTAvec)





# Write input files for CRISPOR to disk -----------------------------------

WriteCRISPORInputFiles(bed_df_list, file_ending = "__CRISPRi.bed",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_bed")
                       )

WriteCRISPORInputFiles(FASTA_vec_list, file_ending = "__CRISPRi.fa",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_FASTA")
                       )





# Write input files (filtered by already present data) --------------------

filtered_bed_df_list <- FilterBedDfList(bed_df_list)
filtered_FASTA_df_list <- FilterFASTADfList(FASTA_df_list)
filtered_FASTA_vec_list <- lapply(filtered_FASTA_df_list, MakeFASTAvec)

WriteCRISPORInputFiles(filtered_bed_df_list, file_ending = "__CRISPRi.bed",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_bed_filtered")
                       )

WriteCRISPORInputFiles(filtered_FASTA_vec_list, file_ending = "__CRISPRi.fa",
                       CRISPOR_input_directory = file.path(CRISPOR_files_directory, "Input_FASTA_filtered")
                       )




