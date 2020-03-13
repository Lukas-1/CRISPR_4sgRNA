### 14th October 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRa_RData_directory, "19) Pick the top 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"))
load(file.path(CRISPRa_RData_directory, "20) Integrate the guide choices using relaxed and strict locations.RData"))






# Create an overview data frame -------------------------------------------

sgRNAs_overview_df <- ProduceGenomeOverviewDf(merged_replaced_CRISPRa_df, lax_CRISPRa_df)





# Create an overview data frame for relaxed locations ---------------------

sgRNAs_lax_overview_df <- ProduceGenomeOverviewDf(merged_replaced_CRISPRa_df, lax_CRISPRa_df, use_lax_df = TRUE)






# Count the number of genes without a full complement of sgRNAs -----------

table(sgRNAs_overview_df[["Num_total"]] < 4)
table(sgRNAs_overview_df[["Num_meeting_criteria"]] < 4)






# Count the number of sgRNAs that overlap with a SNP ----------------------

were_mapped <- !(is.na(merged_replaced_CRISPRa_df[["Cut_location"]]))
table((merged_replaced_CRISPRa_df[[preferred_AF_max_column]][were_mapped] > SNP_frequency_cutoff) %in% TRUE, useNA = "ifany")







# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  all_genes_annotation_columns,
  CRISPRa_columns,
  setdiff(selected_metrics, "Num_overlaps")
)

untargetable_annotations <- c("Not protein-coding", "Only annotated on alternate loci", "Not in current annotation release")

are_targetable <- !(sgRNAs_overview_df[["Gene_annotation_status"]] %in% untargetable_annotations)

WriteOverviewDfToDisk(sgRNAs_overview_df[, columns_for_excel],
                      file_name = "Overview_CRISPRa_all_genes"
                      )
WriteOverviewDfToDisk(sgRNAs_overview_df[are_targetable, columns_for_excel],
                      file_name = "Overview_CRISPRa_all_targetable_genes"
                      )




# Write the summary to disk for relaxed locations -------------------------

are_targetable <- !(sgRNAs_lax_overview_df[["Gene_annotation_status"]] %in% untargetable_annotations)

WriteOverviewDfToDisk(sgRNAs_lax_overview_df[are_targetable, columns_for_excel],
                      file_name = "Overview_CRISPRa_relaxed_all_targetable_genes"
                      )





# Save data ---------------------------------------------------------------

save(list = c("sgRNAs_overview_df", "sgRNAs_lax_overview_df"),
     file = file.path(CRISPRa_RData_directory, "21) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData")
     )











