### 11th January 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R")) # For Are4sg
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R")) # For CheckThatFactorIsInOrder




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Find overlapping genes for 4sg combinations (predicted deletions) -------

are_4sg <- Are4sg(merged_CRISPRko_df, sublibraries_all_entrezs_list)

deletions_CDS_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                CDS_or_exon_locations_df
                                                )
deletions_exon_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                 exon_locations_df
                                                 )




# Determine which sgRNAs can be analyzed ----------------------------------

have_single_entrez <- !(is.na(merged_CRISPRko_df[["Entrez_ID"]])) &
                      !(grepl(",", merged_CRISPRko_df[["Entrez_ID"]], fixed = TRUE))

indices_vec <- rep(NA, nrow(merged_CRISPRko_df))
indices_vec[have_single_entrez] <- seq_len(sum(have_single_entrez))




# Find overlapping genes for sgRNAs ---------------------------------------

guides_CDS_list <- FindOverlappingGenes(merged_CRISPRko_df[have_single_entrez, ],
                                        CDS_or_exon_locations_df
                                        )
guides_exon_list <- FindOverlappingGenes(merged_CRISPRko_df[have_single_entrez, ],
                                         exon_locations_df
                                         )



# Prepare for saving summary data frames ----------------------------------

deletions_CDS_df <- deletions_CDS_list[["summary_df"]]
deletions_exon_df <- deletions_exon_list[["summary_df"]]

guides_CDS_df <- guides_CDS_list[["summary_df"]][have_single_entrez, ]
row.names(guides_CDS_df) <- NULL

guides_exon_df <- guides_CDS_list[["summary_df"]][have_single_entrez, ]
row.names(guides_exon_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = c("deletions_CDS_df", "deletions_exon_df",
              "guides_CDS_df", "guides_exon_df"
              ),
     file = file.path(CRISPRa_RData_directory, "24) Find all genes targeted by each sgRNA.RData")
     )




# gene_ID_df <- delete_gene_ID_df
# symbol_df <- delete_symbol_df
# shared_columns <- delete_shared_columns
#
# for (column in shared_columns) {
#   print(paste0(column, "  ", identical(gene_ID_df[[column]], symbol_df[[column]])))
# }
#
# assign("delete_gene_ID_df", gene_ID_df, envir = globalenv())
# assign("delete_symbol_df", symbol_df, envir = globalenv())
# assign("delete_shared_columns", shared_columns, envir = globalenv())
#
#
#
# are_same <- mapply(identical,
#                    gene_ID_df[["Num_affected_genes"]],
#                    symbol_df[["Num_affected_genes"]]
#                    )
#
# gene_ID_df[!(are_same), ]
# symbol_df[!(are_same), ]
#
#
#
# head(delete_symbol_df)








