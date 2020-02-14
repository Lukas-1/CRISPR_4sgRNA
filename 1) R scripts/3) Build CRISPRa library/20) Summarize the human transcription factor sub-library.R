### 26th September 2019 ###



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

load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))





# Determine the number of available sgRNAs for TF genes -------------------

CRISPRa_TF_sgRNAs_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% all_TF_df[["Combined_ID"]], ]

TF_sgRNAs_summary_df <- SummarizeCRISPRDf(CRISPRa_TF_sgRNAs_df)

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
  annotation_columns,
  "Num_hCRISPRa_v2_transcripts", "Num_transcripts", "Num_overlapping_transcripts", "Num_incomplete_transcripts",
  selected_metrics[selected_metrics != "Num_overlaps"],
  "DNA_binding_domain", "TF_assessment", "Binding_mode",
  "Is_TF_CisBP", "Is_TF_TFClass", "Is_TF_GO", "Is_C2H2_ZF"
)

TF_summary_excel_df <- TF_overview_df[, columns_for_excel]
TF_summary_excel_df[["Num_total"]] <- ifelse(is.na(TF_summary_excel_df[["Num_total"]]), 0L, TF_summary_excel_df[["Num_total"]])

TF_summary_excel_df <- FormatOverviewDfForExport(TF_summary_excel_df)

write.table(TF_summary_excel_df,
            file = file.path(file_output_directory, "Overview_CRISPRa_transcription_factors.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
            )





# Save data ---------------------------------------------------------------

save(list = "TF_overview_df",
     file = file.path(CRISPRa_RData_directory, "20) Summarize the human transcription factor sub-library - TF_overview_df.RData")
     )











