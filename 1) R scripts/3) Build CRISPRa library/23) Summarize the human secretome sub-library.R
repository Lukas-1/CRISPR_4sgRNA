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

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) Pick the top 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"))
load(file.path(CRISPRa_RData_directory, "20) Integrate the guide choices using relaxed and strict locations.RData"))





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




# Does the use of relaxed locations change the choice of guides? ----------

are_different <- DifferUsingRelaxedLocations(secretome_overview_df[["Entrez_ID"]], merged_replaced_CRISPRa_df, lax_CRISPRa_df)
secretome_overview_df[["Lax_locations_differ"]] <- ifelse(are_different, "Yes", "No")





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  setdiff(TF_annotation_columns, "Original_entrez"),
  "Num_hCRISPRa_v2_transcripts", "Num_transcripts", "Num_overlapping_transcripts", "Num_incomplete_transcripts",
  setdiff(selected_metrics, "Num_overlaps"),
  "UniProt_accession", "Annotated_category"
)
WriteOverviewDfToDisk(secretome_overview_df[, columns_for_excel], file_name = "Overview_CRISPRa_secretome")





# Save data ---------------------------------------------------------------

save(list = "secretome_overview_df",
     file = file.path(CRISPRa_RData_directory, "23) Summarize the human secretome sub-library.RData")
     )









