### 13th February 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Determine the number of available sgRNAs for secretome genes ------------

CRISPRko_secretome_sgRNAs_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% secretome_df[["Combined_ID"]], ]

secretome_sgRNAs_summary_df <- SummarizeCRISPRDf(CRISPRko_secretome_sgRNAs_df)

secretome_df <- secretome_df[!(secretome_df[["Ensembl_gene_ID"]] %in% "ENSG00000284779"), ]

reorganized_df <- ReorganizeSummaryDf(secretome_sgRNAs_summary_df, secretome_df[["Combined_ID"]])

secretome_overview_df <- data.frame(secretome_df[match(reorganized_df[["Combined_ID"]], secretome_df[["Combined_ID"]]), ],
                                    reorganized_df[, !(names(reorganized_df) %in% names(secretome_df))],
                                    stringsAsFactors = FALSE,
                                    row.names = NULL
                                    )






# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  setdiff(TF_annotation_columns, "Original_entrez"),
  selected_metrics,
  "UniProt_accession", "Annotated_category"
)

WriteOverviewDfToDisk(secretome_overview_df[, columns_for_excel], file_name = "Overview_CRISPRko_secretome")





# Save data ---------------------------------------------------------------

save(list = "secretome_overview_df",
     file = file.path(CRISPRko_RData_directory, "14) Summarize the human secretome sub-library.RData")
     )









