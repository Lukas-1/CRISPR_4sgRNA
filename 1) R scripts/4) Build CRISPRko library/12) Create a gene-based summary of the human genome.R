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

CRISPRko_entrez_IDs <- unique(merged_CRISPRko_df[["Entrez_ID"]])
CRISPRko_entrez_IDs <- CRISPRko_entrez_IDs[!(is.na(CRISPRko_entrez_IDs))]
CRISPRko_entrez_IDs <- unique(unlist(strsplit(CRISPRko_entrez_IDs, ", ", fixed = TRUE)))

unique_entrez_IDs <- union(collected_entrez_IDs, CRISPRko_entrez_IDs)





# Create an sgRNA overview data frame, with one row per gene --------------

sgRNAs_summary_df <- SummarizeCRISPRDf(merged_CRISPRko_df)
sgRNAs_summary_df[!(is.na(sgRNAs_summary_df[["Entrez_ID"]])) & !(sgRNAs_summary_df[["Entrez_ID"]] %in% collected_entrez_IDs), ]

sgRNAs_all_genes_df <- ReorganizeSummaryDf(sgRNAs_summary_df, unique_entrez_IDs)
sgRNAs_all_genes_df[["Entrez_ID"]] <- sgRNAs_all_genes_df[["Combined_ID"]]
sgRNAs_all_genes_df <- sgRNAs_all_genes_df[, names(sgRNAs_all_genes_df) != "Combined_ID"]

sgRNAs_overview_df <- FixSymbolsForSummaryDf(sgRNAs_all_genes_df)
sgRNAs_overview_df[sgRNAs_overview_df[["Original_entrez"]] != "", ]





# Count the number of genes without a full complement of sgRNAs -----------

table(sgRNAs_overview_df[["Num_total"]] < 4)
table(sgRNAs_overview_df[["Num_meeting_criteria"]] < 4)





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  all_genes_annotation_columns,
  "Gene_present",
  selected_metrics,
  "Num_overlapping_with_SNP"
)

columns_for_excel_with_comments <- c(setdiff(columns_for_excel, c("Gene_present", "Num_overlapping_with_SNP")), "Annotation")

WriteOverviewDfToDisk(sgRNAs_overview_df[, columns_for_excel], file_name = "Overview_CRISPRko_all_genes")
WriteOverviewDfToDisk(sgRNAs_overview_df[, columns_for_excel_with_comments], file_name = "Overview_CRISPRko_all_genes_with_comments")






# Save data ---------------------------------------------------------------

save(list = "sgRNAs_overview_df",
     file = file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData")
     )



