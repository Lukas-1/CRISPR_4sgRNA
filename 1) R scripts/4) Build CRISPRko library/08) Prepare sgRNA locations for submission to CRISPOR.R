### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData")) # The output of "07) Integrate the output from GuideScan" is not yet needed





# Prepare data frames that can be exported to .bed files ------------------

TF_bed_df_list <- lapply(entrez_chunks_list, function(x) MakeBedDf(extended_CRISPRko_df, x))





# Prepare objects that can be exported to FASTA files ---------------------

FASTA_df_list <- lapply(entrez_chunks_list, function(x) MakeFASTADf(extended_CRISPRko_df, x))
FASTA_vec_list <- lapply(FASTA_df_list, MakeFASTAvec)






# Write input files for CRISPOR to disk -----------------------------------

for (chunk_ID in names(entrez_chunks_list)) {
  for (i in 1:2) {
    file_name <- paste0("Input_for_CRISPOR__chunk_", chunk_ID, "__CRISPRko")
    write.table(get(c("TF_bed_df_list", "FASTA_vec_list")[[i]])[[chunk_ID]],
                file = file.path(CRISPOR_files_directory, paste0(file_name, c(".bed", ".fa")[[i]])),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
                )
  }
}



