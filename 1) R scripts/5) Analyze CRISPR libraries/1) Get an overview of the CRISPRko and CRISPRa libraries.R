### 15th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
annotation_intermediate_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Annotation")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))

load(file.path(CRISPRko_RData_directory, "13) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))
load(file.path(CRISPRko_RData_directory, "14) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRko_RData_directory, "15) Summarize the human secretome sub-library.RData"))
CRISPRko_sgRNAs_overview_df    <- sgRNAs_overview_df
CRISPRko_TF_overview_df        <- TF_overview_df
CRISPRko_secretome_overview_df <- secretome_overview_df

load(file.path(CRISPRa_RData_directory, "21) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))
load(file.path(CRISPRa_RData_directory, "22) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRa_RData_directory, "23) Summarize the human secretome sub-library.RData"))
CRISPRa_sgRNAs_overview_df    <- sgRNAs_overview_df
CRISPRa_TF_overview_df        <- TF_overview_df
CRISPRa_secretome_overview_df <- secretome_overview_df

rm(sgRNAs_overview_df)


load(file.path(CRISPRa_RData_directory, "20) For problematic genes, pick 4 guides without reference to the TSS - merged_replaced_CRISPRa_df.RData"))
load(file.path(CRISPRa_RData_directory, "20) For problematic genes, pick 4 guides without reference to the TSS - lax_CRISPRa_df.RData"))

load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))
load(file.path(CRISPRko_RData_directory, "12) Pick 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"))






# Define functions --------------------------------------------------------

ArePresentInOverviewDf <- function(entrez_IDs, overview_df) {
  matches_vec <- match(entrez_IDs, overview_df[["Entrez_ID"]])
  are_spaced <- !(is.na(overview_df[["Spacing"]])) & !(overview_df[["Spacing"]] %in% "None")
  results_mat <- cbind(
    "Are_present"  = (overview_df[["Num_total"]][matches_vec] >= 1) %in% TRUE,
    "Are_complete" = (overview_df[["Num_total"]][matches_vec] >= 4) %in% TRUE,
    "Are_spaced"   = are_spaced[matches_vec] %in% TRUE
  )
  rownames(results_mat) <- entrez_IDs
  return(results_mat)
}


ArePresentInCRISPRDf <- function(entrez_IDs, CRISPR_df) {
  is_CRISPRa <- "Calabrese_rank" %in% colnames(CRISPR_df)
  if (is_CRISPRa) {
    are_top4_mat <- CRISPRaAreTop4Mat(CRISPR_df)
    are_chosen <- are_top4_mat[, "Are_valid_or_only_top4"]
  } else {
    are_top4_mat <- CRISPRkoAreTop4Mat(CRISPR_df)
    are_chosen <- are_top4_mat[, "Are_top4"]
  }

  present_entrezs  <- unique(CRISPR_df[["Entrez_ID"]])
  complete_entrezs <- unique(CRISPR_df[["Entrez_ID"]][are_top4_mat[, "Have_complete_guides"]])
  spaced_entrezs   <- unique(CRISPR_df[["Entrez_ID"]][are_top4_mat[, "Have_spaced_guides"]])

  present_entrezs  <- present_entrezs[!(is.na(present_entrezs))]
  complete_entrezs <- complete_entrezs[!(is.na(complete_entrezs))]
  spaced_entrezs   <- spaced_entrezs[!(is.na(spaced_entrezs))]

  results_mat <- cbind(
    "Are_present"  = entrez_IDs %in% present_entrezs,
    "Are_complete" = entrez_IDs %in% complete_entrezs,
    "Are_spaced"   = entrez_IDs %in% spaced_entrezs
  )
  rownames(results_mat) <- entrez_IDs
  return(results_mat)
}





# Check for genes without sgRNAs from GPP ---------------------------------

are_present_CRISPRko <- collected_entrez_IDs %in% merged_CRISPRko_df[["Entrez_ID"]]
are_present_CRISPRko_GPP <- collected_entrez_IDs %in% merged_CRISPRko_df[["Entrez_ID"]][grepl("GPP", merged_CRISPRko_df[["Source"]], fixed = TRUE)]

are_present_CRISPRa <- collected_entrez_IDs %in% merged_replaced_CRISPRa_df[["Entrez_ID"]]
are_present_CRISPRa_GPP <- collected_entrez_IDs %in% merged_replaced_CRISPRa_df[["Entrez_ID"]][grepl("GPP", merged_replaced_CRISPRa_df[["Source"]], fixed = TRUE)]

only_CRISPRa_GPP <- are_present_CRISPRa_GPP & !(are_present_CRISPRko_GPP)
only_CRISPRko_GPP <- are_present_CRISPRko_GPP & !(are_present_CRISPRa_GPP)





# Examine missing genes ---------------------------------------------------

missing_only_for_CRISPRa_df <- collected_entrezs_df[only_CRISPRko_GPP, c("Entrez_ID", "Gene_symbol", "Category")]
rownames(missing_only_for_CRISPRa_df) <- NULL

missing_only_for_CRISPRko_df <- collected_entrezs_df[only_CRISPRa_GPP, c("Entrez_ID", "Gene_symbol", "Category")]
rownames(missing_only_for_CRISPRko_df) <- NULL

# All genes that are not in current annotation release, or only annotated on alternate loci, are missing
table(collected_entrezs_df[["Category"]], !(are_present_CRISPRa_GPP))
table(collected_entrezs_df[["Category"]], !(are_present_CRISPRko_GPP))
table(collected_entrezs_df[["Category"]], !(are_present_CRISPRa_GPP & are_present_CRISPRko_GPP))




# Check for targetable genes ----------------------------------------------

not_targetable_categories <- c(
  "Not in current annotation release",
  "Only annotated on alternate loci"
)

are_targetable <- !(collected_entrezs_df[["Category"]] %in% not_targetable_categories)

CRISPRa_redo_entrezs_vec  <- collected_entrez_IDs[are_targetable & !(are_present_CRISPRa_GPP)]
CRISPRko_redo_entrezs_vec <- collected_entrez_IDs[are_targetable & !(are_present_CRISPRko_GPP)]





# Export targetable genes -------------------------------------------------

write.table(CRISPRa_redo_entrezs_vec,
            file = file.path(annotation_intermediate_files_directory, "Entrez_IDs_not_covered_by_GPP_CRISPRa.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )
write.table(CRISPRko_redo_entrezs_vec,
            file = file.path(annotation_intermediate_files_directory, "Entrez_IDs_not_covered_by_GPP_CRISPRko.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )






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

are_spaced_in_both   <- present_genes_df[["Are_spaced_CRISPRa"]] & present_genes_df[["Are_spaced_CRISPRko"]]
are_spacd_in_neither <- !(present_genes_df[["Are_spaced_CRISPRa"]]) & !(present_genes_df[["Are_spaced_CRISPRko"]])

present_in_neither_and_targetable <- are_present_in_neither & present_genes_df[["Are_on_reference_genome"]] & present_genes_df[["Are_protein_coding"]]
present_in_neither_and_untargetable <- are_present_in_neither & !(present_genes_df[["Are_on_reference_genome"]]) & present_genes_df[["Are_protein_coding"]]


are_exclusive_to_CRISPRa  <- present_genes_df[["Are_present_CRISPRa"]] &  !(present_genes_df[["Are_present_CRISPRko"]])
are_exclusive_to_CRISPRko <- present_genes_df[["Are_present_CRISPRko"]] & !(present_genes_df[["Are_present_CRISPRa"]])

are_spaced_in_CRISPRa_only  <- present_genes_df[["Are_spaced_CRISPRa"]] &  !(present_genes_df[["Are_spaced_CRISPRko"]])
are_spaced_in_CRISPRko_only <- present_genes_df[["Are_spaced_CRISPRko"]] & !(present_genes_df[["Are_spaced_CRISPRa"]])

are_unspaced_CRISPRa <- present_genes_df[["Are_present_CRISPRa"]] & !(present_genes_df[["Are_spaced_CRISPRa"]]) & present_genes_df[["Are_on_reference_genome"]]
are_unspaced_CRISPRko <- present_genes_df[["Are_present_CRISPRko"]] & !(present_genes_df[["Are_spaced_CRISPRko"]]) & present_genes_df[["Are_on_reference_genome"]]

are_present_but_not_on_genome <- (present_genes_df[["Are_present_CRISPRa"]] |
                                  present_genes_df[["Are_present_CRISPRko"]]
                                  ) & present_genes_df[["Are_protein_coding"]] &
                                  !(present_genes_df[["Are_on_reference_genome"]])





# Count the number of available genes -------------------------------------

sum(are_present_in_both)
sum(present_genes_df[["Are_present_CRISPRa"]])
sum(present_genes_df[["Are_present_CRISPRko"]])

sum(are_spaced_in_both)
sum(present_genes_df[["Are_spaced_CRISPRa"]])
sum(present_genes_df[["Are_spaced_CRISPRko"]])





# Examine subsets of genes ------------------------------------------------

present_genes_df[present_in_neither_and_targetable, ]
present_genes_df[present_in_neither_and_untargetable, ]

table(present_genes_df[present_in_neither_and_untargetable, "Category"])

present_genes_df[are_exclusive_to_CRISPRko, ]

present_genes_df[are_spaced_in_CRISPRa_only, ]
present_genes_df[are_spaced_in_CRISPRko_only, ]

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

all_CRISPRa_df[CRISPRa_are_top4_mat[, "Are_valid_or_only_top4"] & !(CRISPRa_are_top4_mat[, "Have_complete_guides"]), show_columns_CRISPRa]






# Examine problematic sgRNAs for CRISPRko ----------------------------------

all_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Entrez_ID"]] %in% all_entrezs, ]

CRISPRko_are_top4_mat <- CRISPRkoAreTop4Mat(all_CRISPRko_df)

show_columns_CRISPRko <- grep("TSS", show_columns_CRISPRa, value = TRUE, invert = TRUE)

all_CRISPRko_df[CRISPRko_are_top4_mat[, "Are_top4"] & !(CRISPRko_are_top4_mat[, "Have_complete_guides"]), show_columns_CRISPRko]






# Count the number of wells -----------------------------------------------

## Guides picked using strict locations
unique_IDs_CRISPRa_strict     <- unique(all_CRISPRa_df[["Combined_ID"]][CRISPRa_are_top4_mat[, "Are_valid_or_only_top4"]])
unique_TSS_IDs_CRISPRa_strict <- unique(all_CRISPRa_df[["AltTSS_ID"]][CRISPRa_are_top4_mat[, "Are_valid_or_only_top4"]])
unique_IDs_CRISPRko_strict    <- unique(all_CRISPRko_df[["Combined_ID"]][CRISPRko_are_top4_mat[, "Are_top4"]])
length(unique_IDs_CRISPRa_strict)
length(unique_TSS_IDs_CRISPRa_strict)
length(unique_IDs_CRISPRko_strict)


## Guides picked using relaxed locations
lax_CRISPRa_df <- lax_CRISPRa_df[lax_CRISPRa_df[["Entrez_ID"]] %in% all_entrezs, ]
lax_CRISPRko_df <- lax_CRISPRko_df[lax_CRISPRko_df[["Entrez_ID"]] %in% all_entrezs, ]
lax_CRISPRa_are_top4_mat <- CRISPRaAreTop4Mat(lax_CRISPRa_df)
lax_CRISPRko_are_top4_mat <- CRISPRkoAreTop4Mat(lax_CRISPRko_df)
unique_IDs_CRISPRa_lax     <- unique(lax_CRISPRa_df[["Combined_ID"]][lax_CRISPRa_are_top4_mat[, "Are_valid_or_only_top4"]])
unique_TSS_IDs_CRISPRa_lax <- unique(lax_CRISPRa_df[["AltTSS_ID"]][lax_CRISPRa_are_top4_mat[, "Are_valid_or_only_top4"]])
unique_IDs_CRISPRko_lax    <- unique(lax_CRISPRko_df[["Combined_ID"]][lax_CRISPRko_are_top4_mat[, "Are_top4"]])
length(unique_IDs_CRISPRa_lax)
length(unique_TSS_IDs_CRISPRa_lax)
length(unique_IDs_CRISPRko_lax)


identical(unique_IDs_CRISPRko_lax, unique_IDs_CRISPRko_strict)
setdiff(unique_IDs_CRISPRa_lax, unique_IDs_CRISPRa_strict)
setdiff(unique_TSS_IDs_CRISPRa_lax, unique_TSS_IDs_CRISPRa_strict)








# Filter the sublibraries for available genes -----------------------------

sublibrary_filtered_df <- sublibrary_df[sublibrary_df[["Entrez_ID"]] %in% present_genes_df[["Entrez_ID"]][are_present_in_either], ]
rownames(sublibrary_filtered_df) <- NULL

table(sublibrary_filtered_df[["Sublibrary"]])






