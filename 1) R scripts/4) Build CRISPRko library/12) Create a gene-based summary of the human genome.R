### 30th October 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))






# Collect all Entrez IDs from various sources -----------------------------

CRISPRko_entrez_IDs <- unique(merged_CRISPRko_df[, "Entrez_ID"])
CRISPRko_entrez_IDs <- CRISPRko_entrez_IDs[!(is.na(CRISPRko_entrez_IDs))]
CRISPRko_entrez_IDs <- unique(unlist(strsplit(CRISPRko_entrez_IDs, ", ", fixed = TRUE)))

unique_entrez_IDs <- c(collected_entrez_IDs, CRISPRko_entrez_IDs[!(CRISPRko_entrez_IDs %in% collected_entrez_IDs)])





# Create an sgRNA overview data frame, with one row per gene --------------

sgRNAs_summary_df <- SummarizeCRISPRDf(merged_CRISPRko_df)
sgRNAs_summary_df[!(is.na(sgRNAs_summary_df[, "Entrez_ID"])) & !(sgRNAs_summary_df[, "Entrez_ID"] %in% collected_entrez_IDs), ]

sgRNAs_all_genes_df <- ReorganizeSummaryDf(sgRNAs_summary_df, unique_entrez_IDs)
sgRNAs_all_genes_df[, "Entrez_ID"] <- sgRNAs_all_genes_df[, "Combined_ID"]
sgRNAs_all_genes_df <- sgRNAs_all_genes_df[, colnames(sgRNAs_all_genes_df) != "Combined_ID"]

sgRNAs_overview_df <- FixSymbolsForSummaryDf(sgRNAs_all_genes_df)
sgRNAs_overview_df[sgRNAs_overview_df[, "Original_entrez"] != "", ]





# Count the number of genes without a full complement of sgRNAs -----------

table(sgRNAs_overview_df[, "Num_total"] < 4)
table(sgRNAs_overview_df[, "Num_meeting_criteria"] < 4)





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  "Entrez_ID", "Gene_symbol", "Original_symbol", "Original_entrez",
  "Gene_present",
  "Num_overlaps", "Spacing", "Longest_subsequence", "GuideScan_specificity", "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity",
  "Num_total", "Num_meeting_criteria", "Num_overlapping_with_SNP"
)

sgRNAs_overview_excel_df <- sgRNAs_overview_df[, columns_for_excel]

for (i in seq_len(ncol(sgRNAs_overview_excel_df))) {
  sgRNAs_overview_excel_df[, i] <- ifelse(is.na(sgRNAs_overview_excel_df[, i]), "", as.character(sgRNAs_overview_excel_df[, i]))
}

write.table(sgRNAs_overview_excel_df,
            file = file.path(file_output_directory, "Overview_CRISPRko_all_genes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
            )
















