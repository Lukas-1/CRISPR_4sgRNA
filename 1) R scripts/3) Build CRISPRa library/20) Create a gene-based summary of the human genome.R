### 14th October 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")
file_output_directory   <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - num_genes_in_library.RData"))





# Create an overview data frame -------------------------------------------

sgRNAs_overview_df <- ProduceGenomeOverviewDf(merged_replaced_CRISPRa_df,
                                              sublibraries_all_entrezs_list
                                              )





# Count the number of genes without a full complement of sgRNAs -----------

table(sgRNAs_overview_df[["Num_total"]] < 4)
table(sgRNAs_overview_df[["Num_meeting_criteria"]] < 4)






# Count the number of sgRNAs that overlap with a SNP ----------------------

were_mapped <- !(is.na(merged_replaced_CRISPRa_df[["Cut_location"]]))
table((merged_replaced_CRISPRa_df[[preferred_AF_max_column]][were_mapped] > SNP_frequency_cutoff) %in% TRUE, useNA = "ifany")




# Count the number of genes for which multiple TSSs are targeted ----------

table(sgRNAs_overview_df[["Num_hCRISPR_v2_transcripts"]], useNA = "ifany")
table(sgRNAs_overview_df[["Num_hCRISPR_v2_transcripts"]] > 1, useNA = "ifany")






# Check for hCRISPRa-v2 TSSs that are not being targeted separately -------

sgRNAs_overview_df[which(sgRNAs_overview_df[["Num_transcripts"]] < sgRNAs_overview_df[["Num_hCRISPRa_v2_transcripts"]]), ]

show_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol",
                  "Nearest_symbols",
                  "Nearest_gene_distance",
                  "Calabrese_rank", "hCRISPRa_v2_rank",
                  "sgRNA_sequence",

                  "Num_TSSs", "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS",
                  "hCRISPRa_v2_transcript",
                  "Rank", "Original_rank",

                  "Best_TSS", "Entrez_chromosome", "Chromosome", "Strand", "Distance_from_TSS",
                  "CRISPOR_3MM_specificity",
                  "CRISPOR_Num_0MM", "CRISPOR_Num_1MM",
                  "CRISPOR_Num_2MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM"
                  )




# Add the number of genes in the 4sg library ------------------------------

num_genes_in_library <- c(
  num_genes_in_library,
  "4sg" = sum(sgRNAs_overview_df[["In_4sg_library"]] %in% "Yes")
)






# Write the summary data frame to disk ------------------------------------

columns_for_excel <- c(
  all_genes_annotation_columns,
  TSS_columns,
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





# Save data ---------------------------------------------------------------

save(list = c("sgRNAs_overview_df", "num_genes_in_library", "num_TSSs_hCRISPR"),
     file = file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the human genome.RData")
     )






