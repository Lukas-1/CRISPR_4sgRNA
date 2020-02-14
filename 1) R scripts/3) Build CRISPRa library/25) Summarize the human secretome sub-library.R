### 13th February 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))





# Determine the number of available sgRNAs for secretome genes ------------

CRISPRa_secretome_sgRNAs_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% secretome_df[["Combined_ID"]], ]

secretome_sgRNAs_summary_df <- SummarizeCRISPRDf(CRISPRa_secretome_sgRNAs_df)

secretome_df <- secretome_df[!(secretome_df[["Ensembl_gene_ID"]] %in% "ENSG00000284779"), ]

reorganized_df <- ReorganizeSummaryDf(secretome_sgRNAs_summary_df, secretome_df[["Combined_ID"]])

secretome_overview_df <- data.frame(secretome_df[match(reorganized_df[["Combined_ID"]], secretome_df[["Combined_ID"]]), ],
                                    reorganized_df[, !(names(reorganized_df) %in% names(secretome_df))],
                                    stringsAsFactors = FALSE,
                                    row.names = NULL
                                    )




# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  annotation_columns[annotation_columns != "Original_Entrez_ID"],
  "Num_hCRISPRa_v2_transcripts", "Num_transcripts", "Num_overlapping_transcripts", "Num_incomplete_transcripts",
  selected_metrics[selected_metrics != "Num_overlaps"],
  "Annotated_category", "UniProt_accession"
)

secretome_overview_excel_df <- secretome_overview_df[, columns_for_excel]

secretome_overview_excel_df[["Num_total"]] <- ifelse(is.na(secretome_overview_excel_df[["Num_total"]]), 0L, secretome_overview_excel_df[["Num_total"]])

secretome_overview_excel_df <- FormatOverviewDfForExport(secretome_overview_excel_df)

write.table(secretome_overview_excel_df,
            file = file.path(file_output_directory, "Overview_CRISPRa_secretome.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
            )






# Save data ---------------------------------------------------------------

save(list = "secretome_overview_df",
     file = file.path(CRISPRa_RData_directory, "25) Summarize the human secretome sub-library.RData")
     )









