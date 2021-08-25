### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRi")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRi_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))





# Determine the number of available sgRNAs for secretome genes ------------

CRISPRi_secretome_sgRNAs_df <- merged_replaced_CRISPRi_df[merged_replaced_CRISPRi_df[["Combined_ID"]] %in% secretome_df[["Combined_ID"]], ]

secretome_sgRNAs_summary_df <- SummarizeCRISPRDf(CRISPRi_secretome_sgRNAs_df,
                                                 sublibraries_all_entrezs_list
                                                 )

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
  TSS_columns,
  setdiff(selected_metrics, "Num_overlaps"),
  "UniProt_accession", "Annotated_category"
)
WriteOverviewDfToDisk(secretome_overview_df[, columns_for_excel], file_name = "Overview_CRISPRi_secretome")





# Save data ---------------------------------------------------------------

save(list = "secretome_overview_df",
     file = file.path(CRISPRi_RData_directory, "22) Summarize the human secretome sub-library.RData")
     )









