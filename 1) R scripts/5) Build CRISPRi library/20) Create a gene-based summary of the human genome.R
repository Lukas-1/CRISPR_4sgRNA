### 9th April 2020 ###




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
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRi")

gene_lists_directory    <- file.path(CRISPR_root_directory, "2) Input data", "Gene lists", "Vacuolation")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRi_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))




# Read in the list of genes that are "vacuolation hits" -------------------

vacuolation_genes <- read.table(file.path(gene_lists_directory, "hit_gene_IDs_1.5allhitnoderatio.csv"),
                                header = FALSE, stringsAsFactors = FALSE
                                )[[1]]
LSM_genes         <- read.table(file.path(gene_lists_directory, "Additional_LSM_genes.txt"),
                                header = FALSE, stringsAsFactors = FALSE
                                )[[1]]

vacuolation_entrezs <- as.character(c(vacuolation_genes, LSM_genes))





# Create an overview data frame -------------------------------------------

sgRNAs_overview_df <- ProduceGenomeOverviewDf(merged_replaced_CRISPRi_df)







# Count the number of genes without a full complement of sgRNAs -----------

table(sgRNAs_overview_df[["Num_total"]] < 4)
table(sgRNAs_overview_df[["Num_meeting_criteria"]] < 4)






# Count the number of sgRNAs that overlap with a SNP ----------------------

were_mapped <- !(is.na(merged_replaced_CRISPRi_df[["Cut_location"]]))
table((merged_replaced_CRISPRi_df[[preferred_AF_max_column]][were_mapped] > SNP_frequency_cutoff) %in% TRUE, useNA = "ifany")






# Check for hCRISPRi-v2 TSSs that are not being targeted separately -------

sgRNAs_overview_df[which(sgRNAs_overview_df[["Num_transcripts"]] < sgRNAs_overview_df[["Num_hCRISPR_v2_transcripts"]]), ]





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  all_genes_annotation_columns,
  TSS_columns,
  setdiff(selected_metrics, "Num_overlaps")
)

untargetable_annotations <- c("Not protein-coding", "Only annotated on alternate loci", "Not in current annotation release")

are_targetable <- !(sgRNAs_overview_df[["Gene_annotation_status"]] %in% untargetable_annotations)

WriteOverviewDfToDisk(sgRNAs_overview_df[, columns_for_excel],
                      file_name = "Overview_CRISPRi_all_genes"
                      )
WriteOverviewDfToDisk(sgRNAs_overview_df[are_targetable, columns_for_excel],
                      file_name = "Overview_CRISPRi_all_targetable_genes"
                      )




# Write an overview of the vacuolation candidate genes to disk ------------

vacuolation_summary_df <- SummarizeCRISPRDf(merged_replaced_CRISPRi_df[merged_replaced_CRISPRi_df[["Combined_ID"]] %in% vacuolation_entrezs, ])

vacuolation_overview_df <- ReorganizeSummaryDf(vacuolation_summary_df, vacuolation_entrezs)
vacuolation_overview_df[["Entrez_ID"]] <- vacuolation_overview_df[["Combined_ID"]]
vacuolation_overview_df <- vacuolation_overview_df[, names(vacuolation_overview_df) != "Combined_ID"]
vacuolation_overview_df <- FixSymbolsForSummaryDf(vacuolation_overview_df)

WriteOverviewDfToDisk(vacuolation_overview_df[, columns_for_excel],
                      file_name = "Overview_CRISPRi_vacuolation"
                      )






# Save data ---------------------------------------------------------------

save(list = "sgRNAs_overview_df",
     file = file.path(CRISPRi_RData_directory, "20) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData")
     )

save(list = "vacuolation_entrezs",
     file = file.path(CRISPRi_RData_directory, "20) Create a gene-based summary of the human genome - vacuolation_entrezs.RData")
     )












