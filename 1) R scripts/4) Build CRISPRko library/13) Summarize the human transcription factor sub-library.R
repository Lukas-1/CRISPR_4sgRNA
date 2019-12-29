### 26th September 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))






# Determine the number of available sgRNAs for TF genes -------------------

CRISPRko_TF_sgRNAs_df <- merged_CRISPRko_df[merged_CRISPRko_df[, "Combined_ID"] %in% all_TF_df[, "Combined_ID"], ]

TF_sgRNAs_summary_df <- SummarizeCRISPRDf(CRISPRko_TF_sgRNAs_df)

reorganized_df <- ReorganizeSummaryDf(TF_sgRNAs_summary_df, all_TF_df[, "Combined_ID"])
all_TF_summary_df <- data.frame(all_TF_df[match(reorganized_df[, "Combined_ID"], all_TF_df[, "Combined_ID"]), ],
                                reorganized_df[, !(colnames(reorganized_df) %in% colnames(all_TF_df))],
                                stringsAsFactors = FALSE,
                                row.names = NULL
                                )




# Filter TF_summary_df for only bona fide transcription factors -----------

TF_combined_IDs <- all_TF_df[all_TF_df[, "Is_TF"] == "Yes", "Combined_ID"]

TF_summary_df <- all_TF_summary_df[all_TF_summary_df[, "Combined_ID"] %in% TF_combined_IDs, colnames(all_TF_summary_df) != "Is_TF"]
rownames(TF_summary_df) <- NULL

table(TF_summary_df[, "Num_total"] == 0)
table(TF_summary_df[, "Num_total"] < 4)
table(TF_summary_df[, "Num_meeting_criteria"] < 4)

# table(TF_summary_df[(TF_summary_df[, "Num_total"] > 0) & (TF_summary_df[, "Num_meeting_criteria"] < 4), "Submitted_to_GuideScan"])
# TF_summary_df[(TF_summary_df[, "Num_total"] > 0) & (TF_summary_df[, "Num_meeting_criteria"] < 4) & (TF_summary_df[, "Submitted_to_GuideScan"] == "No"), ]





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  "Gene_symbol", "ENSEMBL_gene_ID", "Entrez_ID", "Original_symbol", "Original_Entrez_ID",
  "Num_overlaps", "Spacing", "Longest_subsequence", "GuideScan_specificity", "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity",
  "Num_top4_outside_criteria", "Num_total",
  "DNA_binding_domain", "TF_assessment", "Binding_mode",
  "Is_TF_CisBP", "Is_TF_TFClass", "Is_TF_GO", "Is_C2H2_ZF"
)

TF_summary_excel_df <- TF_summary_df[, columns_for_excel]

TF_summary_excel_df[, "Num_total"] <- ifelse(is.na(TF_summary_excel_df[, "Num_total"]), 0L, TF_summary_excel_df[, "Num_total"])

for (i in seq_len(ncol(TF_summary_excel_df))) {
  TF_summary_excel_df[, i] <- ifelse(is.na(TF_summary_excel_df[, i]), "", as.character(TF_summary_excel_df[, i]))
}

write.table(TF_summary_excel_df,
            file = file.path(file_output_directory, "Overview_CRISPRko_transcription_factors.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
            )





# Save data ---------------------------------------------------------------

save(list = "TF_summary_df",
     file = file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library.RData")
     )











