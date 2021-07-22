### 23rd November 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "10) Rat - General")
CRISPRko_RData_directory <- file.path(RData_directory, "12) Rat - CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Rat - CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Rn.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Create an empty SNP column ----------------------------------------------

merged_CRISPRko_df[[preferred_AF_max_column]] <- NA_real_





# Create an overview data frame -------------------------------------------

sgRNAs_overview_df <- ProduceGenomeOverviewDf(merged_CRISPRko_df,
                                              sublibraries_entrezs_list = list(collected_entrez_IDs),
                                              is_rat = TRUE
                                              )





# Count the number of genes without a full complement of sgRNAs -----------

table(sgRNAs_overview_df[["Num_total"]] < 4)
table(sgRNAs_overview_df[["Num_meeting_criteria"]] < 4)





# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  all_genes_annotation_columns,
  selected_metrics
)

untargetable_annotations <- c("Not protein-coding", "Only annotated on alternate loci", "Not in current annotation release")

are_targetable <- !(sgRNAs_overview_df[["Gene_annotation_status"]] %in% untargetable_annotations)

WriteOverviewDfToDisk(sgRNAs_overview_df[, columns_for_excel],
                      file_name = "Overview_CRISPRko_all_genes"
                      )
WriteOverviewDfToDisk(sgRNAs_overview_df[are_targetable, columns_for_excel],
                      file_name = "Overview_CRISPRko_all_targetable_genes"
                      )




# Save data ---------------------------------------------------------------

save("sgRNAs_overview_df",
     file = file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the rat genome - sgRNAs_overview_df.RData")
     )





