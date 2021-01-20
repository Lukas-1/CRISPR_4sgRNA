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
deletions_exons_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                  exon_locations_df
                                                  )




# Find overlapping genes for sgRNAs ---------------------------------------

guides_CDS_list <- FindOverlappingGenes(merged_CRISPRko_df[are_example_genes, ],
                                        CDS_or_exon_locations_df
                                        )

guides_exons_list <- FindOverlappingGenes(merged_CRISPRko_df[are_example_genes, ],
                                          exon_locations_df
                                          )

table(overlap_list[["summary_df"]][["Affects_intended_gene"]])
table(overlap_list_2[["summary_df"]][["Affects_intended_gene"]])

summary_df <- overlap_list[["summary_df"]]

all_entrezs <- merged_CRISPRko_df[["Entrez_ID"]][are_4sg]






comb_blocks <- Get4sgProjectedDeletions(merged_CRISPRko_df[are_4sg, ])
comb_no_blocks <- Get4sgProjectedDeletions(merged_CRISPRko_df[are_4sg, ],
                                           split_large_deletions = FALSE
                                           )
comb_only_primary_blocks <- Get4sgProjectedDeletions(merged_CRISPRko_df[are_4sg, ],
                                                     only_primary_location = TRUE
                                                     )
comb_only_primary_no_blocks <- Get4sgProjectedDeletions(merged_CRISPRko_df[are_4sg, ],
                                                        only_primary_location = TRUE,
                                                        split_large_deletions = FALSE
                                                        )

reduced_df <- comb_blocks[["reduced_df"]]

reduced_df[["Factor"]] <- reduced_df[["Span"]] / 10^6

reduced_df <- reduced_df[order(reduced_df[["Span"]], decreasing = TRUE), ]

goo






reduced_df <- comb_no_blocks[["reduced_df"]]

reduced_df[["Factor"]] <- reduced_df[["Span"]] / 10^6

reduced_df <- reduced_df[order(reduced_df[["Span"]], decreasing = TRUE), ]
head(reduced_df)


reduced_df[reduced_df[["Span"]] > 10^6, ]


gaps <- reduced_df[["End"]] - reduced_df[["Start"]]


reduced_df[reduced_df[["Span"]] > 10^6, ]



expanded_df <- comb_blocks[["expanded_df"]]

expanded_df[expanded_df[["List_index"]] %in% 19471, ]


reduced_df[((gaps > 0) %in% TRUE) & ((gaps < 3) %in% TRUE), ]

are_chosen <- are_first | are_last

expanded_df <- expanded_df[are_chosen, ]

index_chromosomes <- index_chromosomes[are_chosen]

expanded_df[expanded_df[["List_index"]] %in% 1763, ]





ex_sub_df <- ex_CRISPR_df[ex_CRISPR_df[["Entrez_ID"]] %in% "276", ]



# summary_df <- overlap_list[["summary_df"]]
#
#
#
# shared_columns <- intersect(names(summary_df),
#                             names(new_summary_df)
#                             )
# shared_columns <- setdiff(shared_columns, c("Affected_gene_IDs", "Affected_gene_symbols"))
#
#
# identical(summary_df[, shared_columns], new_summary_df[, shared_columns])
#
# are_same <- lengths(affected_symbols_results[["list"]]) == lengths(affected_genes_results[["list"]])
#
# lengths(affected_symbols_results[["list"]])[!(are_same)]
# lengths(affected_genes_results[["list"]])[!(are_same)]
#
# summary_df[!(are_same), ]
#
#
# stopifnot(identical(lengths(affected_symbols_results[["list"]]),
#                     lengths(affected_genes_results[["list"]])
#                     ))
#
#
#
# try_df <- full_df[full_df[["Index"]] %in% c(469005, 469008, 469014, 469029 ), ]
#
#
# try_df <- full_df[full_df[["Index"]] %in% c(469005), ]
#
# try_df <- unique(try_df[, c("Index", "Affected_gene_symbols", "Affected_Entrez_IDs", "Affected_gene_IDs")])
# num_occurrences_vec <- table(try_df[["Affected_gene_symbols"]])[try_df[["Affected_gene_symbols"]]]
#
# try_df[num_occurrences_vec > 1, ]
#
#
#
#
# CDS_or_exon_locations_df[CDS_or_exon_locations_df[["Gene_symbols"]] %in% "LINC00904", ]
#
#

#
# symbol_df <- delete_symbol_df
# gene_ID_df <- delete_gene_ID_df
#
#
# ## DELETE THIS::
# #
shared_columns <- intersect(names(symbol_df), names(gene_ID_df ))

identical(symbol_df[, shared_columns], gene_ID_df[, shared_columns])

for (column in shared_columns) {
  print(paste0(column, ":", identical(symbol_df[[column]], gene_ID_df[[column]])))
}

are_same <- mapply(identical, symbol_df[["Num_affected_genes"]], gene_ID_df[["Num_affected_genes"]])
symbol_df[!(are_same), ]
gene_ID_df[!(are_same), ]

symbol_df[!(are_same), "Num_affected_genes"]
gene_ID_df[!(are_same), "Num_affected_genes"]

# symbol_df <- symbol_df[, shared_columns]
# gene_ID_df <- gene_ID_df[, shared_columns]
# identical(symbol_df, gene_ID_df)
#
#
# are_same <- mapply(identical,
#                    gene_ID_df[["Loci_targeting_intended_gene"]],
#                    symbol_df[["Loci_targeting_intended_gene"]]
#                    )
#
# table(are_same)
#
# head(gene_ID_df[!(are_same), ])
# head(symbol_df[!(are_same), ])



are_4sg <- Are4sg(merged_CRISPRko_df, sublibraries_all_entrezs_list)

CRISPR_df <- merged_CRISPRko_df[are_4sg, ]
input_genes_df <- CDS_or_exon_locations_df

only_protein_coding  <- FALSE
exclude_pseudogenes  <- TRUE
only_known_gene_type <- TRUE
require_Entrez_ID    <- FALSE



shared_columns <- intersect(names(symbol_df), names(gene_ID_df))
symbol_df[!(are_same), "Num_affected_genes"]
gene_ID_df[!(are_same), "Num_affected_genes"]












