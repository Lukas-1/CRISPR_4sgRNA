### 23rd December 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "19) Using CRISPOR.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
CRISPOR_files_directory  <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "CRISPOR")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "07) Integrate the output from GuideScan.RData"))






# Select subsets of genes for submission to CRISPOR -----------------------

TF_combined_IDs <- intersect(all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Combined_ID"], merged_CRISPRko_df[, "Combined_ID"])




# Prepare data frames that can be exported to .bed files ------------------

TF_bed_df <- MakeBedDf(merged_CRISPRko_df, TF_combined_IDs)




# Prepare objects that can be exported to FASTA files ---------------------

FASTA_df <- MakeFASTADf(merged_CRISPRko_df, TF_combined_IDs)
FASTA_vec <- MakeFASTAvec(FASTA_df)





# Write input files for CRISPOR to disk -----------------------------------

write.table(TF_bed_df,
            file = file.path(CRISPOR_files_directory, "Input_for_CRISPOR_CRISPRko_TFs.bed"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )

write.table(FASTA_vec,
            file = file.path(CRISPOR_files_directory, "Input_for_CRISPOR_CRISPRko_TFs.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )














