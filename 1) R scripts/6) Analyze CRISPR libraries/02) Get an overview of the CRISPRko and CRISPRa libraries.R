### 15th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R")) # For CRISPRaAreTop4Mat



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
annotation_intermediate_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Annotation")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))

load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRko_RData_directory, "14) Summarize the human secretome sub-library.RData"))
CRISPRko_sgRNAs_overview_df    <- sgRNAs_overview_df
CRISPRko_TF_overview_df        <- TF_overview_df
CRISPRko_secretome_overview_df <- secretome_overview_df

load(file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRa_RData_directory, "21) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRa_RData_directory, "22) Summarize the human secretome sub-library.RData"))
CRISPRa_sgRNAs_overview_df    <- sgRNAs_overview_df
CRISPRa_TF_overview_df        <- TF_overview_df
CRISPRa_secretome_overview_df <- secretome_overview_df

rm(sgRNAs_overview_df)


load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Define functions --------------------------------------------------------

ArePresentInOverviewDf <- function(entrez_IDs, overview_df) {
  matches_vec <- match(entrez_IDs, overview_df[["Entrez_ID"]])
  are_non_homologous <- !(is.na(overview_df[["Spacing"]])) & !(overview_df[["Spacing"]] %in% "None")
  results_mat <- cbind(
    "Are_present"  = (overview_df[["Num_total"]][matches_vec] >= 1) %in% TRUE,
    "Are_complete" = (overview_df[["Num_total"]][matches_vec] >= 4) %in% TRUE,
    "Are_non_homologous" = are_non_homologous[matches_vec] %in% TRUE
  )
  rownames(results_mat) <- entrez_IDs
  return(results_mat)
}


ArePresentInCRISPRDf <- function(entrez_IDs, CRISPR_df) {
  is_CRISPRa <- "Calabrese_rank" %in% colnames(CRISPR_df)
  if (is_CRISPRa) {
    are_top4_mat <- CRISPRaAreTop4Mat(CRISPR_df)
  } else {
    are_top4_mat <- CRISPRkoAreTop4Mat(CRISPR_df)
  }


  present_entrezs        <- unique(CRISPR_df[["Entrez_ID"]])
  complete_entrezs       <- unique(CRISPR_df[["Entrez_ID"]][are_top4_mat[, "Have_complete_guides"]])
  non_homologous_entrezs <- unique(CRISPR_df[["Entrez_ID"]][are_top4_mat[, "Have_non_homologous_guides"]])

  present_entrezs        <- present_entrezs[!(is.na(present_entrezs))]
  complete_entrezs       <- complete_entrezs[!(is.na(complete_entrezs))]
  non_homologous_entrezs <- non_homologous_entrezs[!(is.na(non_homologous_entrezs))]

  results_mat <- cbind(
    "Are_present"        = entrez_IDs %in% present_entrezs,
    "Are_complete"       = entrez_IDs %in% complete_entrezs,
    "Are_non_homologous" = entrez_IDs %in% non_homologous_entrezs
  )
  rownames(results_mat) <- entrez_IDs
  return(results_mat)
}




# Define the set of included genes ----------------------------------------

TF_entrezs <- unique(c(CRISPRa_TF_overview_df[["Entrez_ID"]][!(is.na(CRISPRa_TF_overview_df[["Num_total"]]))],
                       CRISPRko_TF_overview_df[["Entrez_ID"]][!(is.na(CRISPRko_TF_overview_df[["Num_total"]]))]
                       )
                     )
secretome_entrezs <- unique(c(CRISPRa_secretome_overview_df[["Entrez_ID"]][!(is.na(CRISPRa_secretome_overview_df[["Num_total"]]))],
                              CRISPRko_secretome_overview_df[["Entrez_ID"]][!(is.na(CRISPRko_secretome_overview_df[["Num_total"]]))]
                              )
                            )
all_entrezs <- unique(c(collected_entrez_IDs, TF_entrezs, secretome_entrezs))




# Construct a data frame summarizing which genes are present --------------

CRISPRa_are_present_mat <- ArePresentInOverviewDf(all_entrezs, CRISPRa_sgRNAs_overview_df)
stopifnot(identical(CRISPRa_are_present_mat, ArePresentInCRISPRDf(all_entrezs, merged_replaced_CRISPRa_df)))

CRISPRko_are_present_mat <- ArePresentInOverviewDf(all_entrezs, CRISPRko_sgRNAs_overview_df)
stopifnot(identical(CRISPRko_are_present_mat, ArePresentInCRISPRDf(all_entrezs, merged_CRISPRko_df)))

colnames(CRISPRa_are_present_mat)  <- paste0(colnames(CRISPRa_are_present_mat),  "_CRISPRa")
colnames(CRISPRko_are_present_mat) <- paste0(colnames(CRISPRko_are_present_mat), "_CRISPRko")

untargetable_annotations <- c("Only annotated on alternate loci", "Not in current annotation release")

targetable_entrezs <- collected_entrezs_df[["Entrez_ID"]][!(collected_entrezs_df[["Category"]] %in% untargetable_annotations)]
categories_vec <- collected_entrezs_df[["Category"]][match(all_entrezs, collected_entrezs_df[["Entrez_ID"]])]

present_genes_df <- data.frame(
  "Entrez_ID"               = all_entrezs,
  MapToEntrezs(all_entrezs)["Gene_symbol"],
  "Category"                = ifelse(is.na(categories_vec), "Not protein-coding", categories_vec),
  "Are_protein_coding"      = all_entrezs %in% collected_entrez_IDs,
  "Are_on_reference_genome" = all_entrezs %in% targetable_entrezs,
  CRISPRa_are_present_mat,
  CRISPRko_are_present_mat,
  stringsAsFactors = FALSE
)



# Define subsets of genes -------------------------------------------------

are_present_in_both    <- present_genes_df[["Are_present_CRISPRa"]] & present_genes_df[["Are_present_CRISPRko"]]
are_present_in_neither <- !(present_genes_df[["Are_present_CRISPRa"]]) & !(present_genes_df[["Are_present_CRISPRko"]])
are_present_in_either  <- present_genes_df[["Are_present_CRISPRa"]] | present_genes_df[["Are_present_CRISPRko"]]

are_non_homologous_in_both <- present_genes_df[["Are_non_homologous_CRISPRa"]] & present_genes_df[["Are_non_homologous_CRISPRko"]]
are_non_homologous_in_neither <- !(present_genes_df[["Are_non_homologous_CRISPRa"]]) & !(present_genes_df[["Are_non_homologous_CRISPRko"]])

present_in_neither_and_targetable <- are_present_in_neither & present_genes_df[["Are_on_reference_genome"]] & present_genes_df[["Are_protein_coding"]]
present_in_neither_and_untargetable <- are_present_in_neither & !(present_genes_df[["Are_on_reference_genome"]]) & present_genes_df[["Are_protein_coding"]]


are_exclusive_to_CRISPRa  <- present_genes_df[["Are_present_CRISPRa"]] &  !(present_genes_df[["Are_present_CRISPRko"]])
are_exclusive_to_CRISPRko <- present_genes_df[["Are_present_CRISPRko"]] & !(present_genes_df[["Are_present_CRISPRa"]])

are_non_homologous_in_CRISPRa_only  <- present_genes_df[["Are_non_homologous_CRISPRa"]] &  !(present_genes_df[["Are_non_homologous_CRISPRko"]])
are_non_homologous_in_CRISPRko_only <- present_genes_df[["Are_non_homologous_CRISPRko"]] & !(present_genes_df[["Are_non_homologous_CRISPRa"]])

are_unspaced_CRISPRa <- present_genes_df[["Are_present_CRISPRa"]] & !(present_genes_df[["Are_non_homologous_CRISPRa"]]) & present_genes_df[["Are_on_reference_genome"]]
are_unspaced_CRISPRko <- present_genes_df[["Are_present_CRISPRko"]] & !(present_genes_df[["Are_non_homologous_CRISPRko"]]) & present_genes_df[["Are_on_reference_genome"]]

are_present_but_not_on_genome <- (present_genes_df[["Are_present_CRISPRa"]] |
                                  present_genes_df[["Are_present_CRISPRko"]]
                                  ) & present_genes_df[["Are_protein_coding"]] &
                                  !(present_genes_df[["Are_on_reference_genome"]])




# Count the number of available genes -------------------------------------

sum(are_present_in_both)
sum(present_genes_df[["Are_present_CRISPRa"]])
sum(present_genes_df[["Are_present_CRISPRko"]])

sum(are_non_homologous_in_both)
sum(present_genes_df[["Are_non_homologous_CRISPRa"]])
sum(present_genes_df[["Are_non_homologous_CRISPRko"]])




# Examine subsets of genes ------------------------------------------------

present_genes_df[present_in_neither_and_targetable, ]
present_genes_df[present_in_neither_and_untargetable, ]

table(present_genes_df[present_in_neither_and_untargetable, "Category"])

present_genes_df[are_exclusive_to_CRISPRko, ]

present_genes_df[are_non_homologous_in_CRISPRa_only, ]
present_genes_df[are_non_homologous_in_CRISPRko_only, ]

present_genes_df[are_present_but_not_on_genome, ]

present_genes_df[are_unspaced_CRISPRa, ]
present_genes_df[are_unspaced_CRISPRko, ]




# Examine problematic sgRNAs for CRISPRa ----------------------------------

all_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% all_entrezs, ]

CRISPRa_are_top4_mat <- CRISPRaAreTop4Mat(all_CRISPRa_df)

show_columns_CRISPRa <- c(
  "Gene_symbol", "Source", "Chromosome", "Cut_location",
  "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS", "Num_TSSs",
  "Rank", "Original_rank",
  "Num_overlaps",  "Overlaps_tolerance", "Spacing",
  "Best_combination_rank",
  "sgRNA_sequence", "PAM"
)

all_CRISPRa_df[CRISPRa_are_top4_mat[, "Are_chosen_4sg"] & !(CRISPRa_are_top4_mat[, "Have_complete_guides"]), show_columns_CRISPRa]





# Examine problematic sgRNAs for CRISPRko ----------------------------------

all_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Entrez_ID"]] %in% all_entrezs, ]

CRISPRko_are_top4_mat <- CRISPRkoAreTop4Mat(all_CRISPRko_df)

show_columns_CRISPRko <- grep("TSS", show_columns_CRISPRa, value = TRUE, invert = TRUE)

all_CRISPRko_df[CRISPRko_are_top4_mat[, "Are_top4"] & !(CRISPRko_are_top4_mat[, "Have_complete_guides"]), show_columns_CRISPRko]





# Count the number of wells -----------------------------------------------

unique_IDs_CRISPRa     <- unique(all_CRISPRa_df[["Combined_ID"]][CRISPRa_are_top4_mat[, "Are_chosen_4sg"]])
unique_TSS_IDs_CRISPRa <- unique(all_CRISPRa_df[["AltTSS_ID"]][CRISPRa_are_top4_mat[, "Are_chosen_4sg"]])
unique_IDs_CRISPRko    <- unique(all_CRISPRko_df[["Combined_ID"]][CRISPRko_are_top4_mat[, "Are_chosen_4sg"]])
length(unique_IDs_CRISPRa)
length(unique_TSS_IDs_CRISPRa)
length(unique_IDs_CRISPRko)






# Filter the sublibraries for available genes -----------------------------

sublibrary_filtered_df <- sublibrary_df[sublibrary_df[["Entrez_ID"]] %in% present_genes_df[["Entrez_ID"]][are_present_in_either], ]
rownames(sublibrary_filtered_df) <- NULL

table(sublibrary_filtered_df[["Sublibrary"]])






# Count the number of plasmids with gRNAs derived from GPP ----------------

load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRa_by_well_df <- full_4sg_by_well_df
load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_by_well_df <- full_4sg_by_well_df
rm(full_4sg_by_gene_df)
rm(full_4sg_by_well_df)


are_selected <- !((CRISPRa_by_well_df[["Is_obsolete"]] %in% "Yes") |
                  (CRISPRa_by_well_df[["Is_control"]] %in% "Yes")
                  )
CRISPRa_use_df <- CRISPRa_by_well_df[are_selected, c("Combined_ID", "AltTSS_ID", "Source")]

have_GPP_CRISPRa <- tapply(CRISPRa_use_df[["Source"]],
                           factor(CRISPRa_use_df[["AltTSS_ID"]], levels = unique(CRISPRa_use_df[["AltTSS_ID"]])),
                           function(x) any(x == "GPP")
                           )
table(have_GPP_CRISPRa)


are_selected <- !((CRISPRko_by_well_df[["Is_obsolete"]] %in% "Yes") |
                  (CRISPRko_by_well_df[["Is_control"]] %in% "Yes")
                  )
CRISPRko_use_df <- CRISPRko_by_well_df[are_selected, c("Combined_ID", "Source")]

have_GPP_CRISPRko <- tapply(CRISPRko_use_df[["Source"]],
                            factor(CRISPRko_use_df[["Combined_ID"]], levels = unique(CRISPRko_use_df[["Combined_ID"]])),
                            function(x) any(x == "GPP")
                            )






