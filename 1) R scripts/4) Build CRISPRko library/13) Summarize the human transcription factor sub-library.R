### 26th September 2019 ###



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
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Determine the number of available sgRNAs for TF genes -------------------

CRISPRko_TF_sgRNAs_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Combined_ID"]] %in% all_TF_df[["Combined_ID"]], ]

TF_sgRNAs_summary_df <- SummarizeCRISPRDf(CRISPRko_TF_sgRNAs_df,
                                          sublibraries_all_entrezs_list
                                          )

reorganized_df <- ReorganizeSummaryDf(TF_sgRNAs_summary_df, all_TF_df[["Combined_ID"]])
all_TF_summary_df <- data.frame(all_TF_df[match(reorganized_df[["Combined_ID"]], all_TF_df[["Combined_ID"]]), ],
                                reorganized_df[, !(names(reorganized_df) %in% names(all_TF_df))],
                                stringsAsFactors = FALSE,
                                row.names = NULL
                                )




# Filter TF_overview_df for only bona fide transcription factors ----------

TF_combined_IDs <- all_TF_df[["Combined_ID"]][all_TF_df[["Is_TF"]] == "Yes"]

TF_overview_df <- all_TF_summary_df[all_TF_summary_df[["Combined_ID"]] %in% TF_combined_IDs, names(all_TF_summary_df) != "Is_TF"]
row.names(TF_overview_df) <- NULL

table(TF_overview_df[["Num_total"]] == 0)
table(TF_overview_df[["Num_total"]] < 4)
table(TF_overview_df[["Num_meeting_criteria"]] < 4)





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  TF_annotation_columns,
  selected_metrics,
  "DNA_binding_domain", "TF_assessment", "Binding_mode",
  "Is_TF_CisBP", "Is_TF_TFClass", "Is_TF_GO", "Is_C2H2_ZF"
)

WriteOverviewDfToDisk(TF_overview_df[, columns_for_excel], file_name = "Overview_CRISPRko_transcription_factors")





# Save data ---------------------------------------------------------------

save(list = "TF_overview_df",
     file = file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library - TF_overview_df.RData")
     )





