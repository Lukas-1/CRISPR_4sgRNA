### 11th January 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Define functions --------------------------------------------------------




# Do stuff ----------------------------------------------------------------

example_entrezs <- unique(merged_CRISPRko_df[["Entrez_ID"]])[1:1000]
# example_entrezs <- unlist(sublibraries_all_entrezs_list, use.names = FALSE)
are_example_genes <- merged_CRISPRko_df[["Entrez_ID"]] %in% example_entrezs


CRISPR_df <- merged_CRISPRko_df[are_example_genes, ]


stopifnot(!(anyNA(CRISPR_df[["Entrez_ID"]])))
stopifnot(!(any(grepl(",", CRISPR_df[["Entrez_ID"]], fixed = TRUE))))


genes_df <- PrepareGenesDf(exon_locations_df)


# are_same <- mapply(identical,
#                    lengths(strsplit(genes_df[["Entrez_IDs"]], ", ", fixed = TRUE)),
#                    lengths(strsplit(genes_df[["Gene_symbols"]], ", ", fixed = TRUE))
#                    )


rename_gene_columns <- c(
  "Gene_IDs" = "Affected_gene_IDs",
  "Source"   = "Gene_source"
)

for (column_name in names(rename_gene_columns)) {
  names(genes_df)[names(genes_df) == column_name] <- rename_gene_columns[[column_name]]
}

split_results <- Process0MMLoci(CRISPR_df)
unique_loci_df <- split_results[["unique_df"]]
expanded_0MM_df <- split_results[["expanded_df"]]
rm(split_results)

combined_df <- FindOverlappingHits(unique_loci_df, genes_df)
full_combined_df <- AlignHits(expanded_0MM_df, combined_df)



identical(lengths(strsplit(full_combined_df[["Affected_Entrez_IDs"]], ", ", fixed = TRUE)),
          lengths(strsplit(full_combined_df[["Affected_gene_symbols"]], ", ", fixed = TRUE))
          )


combined_columns <- c(
  "Locus_0MM",
  "Index", "Is_primary_location", "Intended_Entrez_ID", "Affected_Entrez_IDs",
  "Number_of_gene_IDs", "Number_of_Entrez_IDs",
  "Intended_gene_symbol", "Affected_gene_symbols", "Num_loci",
  "Guide_locus", "Gene_locus", "Chromosome",
  "Affected_gene_IDs", "Ensembl_gene_IDs", "Ensembl_transcript_ID", "Exon_ID", "Strand",
  "Gene_types", "Gene_source"
)

full_combined_df <- full_combined_df[, combined_columns]

stopifnot(identical(unique(full_combined_df[["Index"]]),
                    seq_len(nrow(CRISPR_df))
                    )
          )


summary_df <- SummarizeFullDf(full_combined_df)

stopifnot(nrow(summary_df) == nrow(CRISPR_df))


summary_df <- SummarizeSummaryDf(summary_df)









