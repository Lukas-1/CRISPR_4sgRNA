### 15th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))




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

load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))
CRISPRko_sgRNAs_overview_df <- sgRNAs_overview_df

load(file.path(CRISPRa_RData_directory, "19) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))
CRISPRa_sgRNAs_overview_df <- sgRNAs_overview_df

rm(sgRNAs_overview_df)

load(file.path(CRISPRko_RData_directory, "11) Re-order the library to prioritize non-overlapping sgRNAs.RData"))
load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))






# Check for missing genes -------------------------------------------------

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









