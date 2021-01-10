### 18th October 2020 ###




# Import packages and source code -----------------------------------------

library("Matrix")
library("igraph")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R")) # For MakeLocationStrings
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R")) # For GetCutLocations
# source(file.path(general_functions_directory, "29) Determining gene types.R"))

source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))





# Define functions --------------------------------------------------------

MergeCommonElements <- function(input_list) {
  # from https://stackoverflow.com/a/47328161
  i <- rep(seq_along(input_list), lengths(input_list))
  j <- factor(unlist(input_list))
  tab <- sparseMatrix(i = i, j = as.integer(j), x = TRUE, dimnames = list(NULL, levels(j)))
  connects <- tcrossprod(tab, boolArith = TRUE)
  group <- clusters(graph_from_adjacency_matrix(as(connects, "lsCMatrix")))[["membership"]]
  results_list <- tapply(input_list,
                         group,
                         function(x) sort(unique(unlist(x))),
                         simplify = FALSE
                         )
  return(results_list)
}




Process0MMLoci <- function(CRISPR_df) {
  loci_splits <- strsplit(CRISPR_df[, "Locations_0MM"], "; ", fixed = TRUE)

  primary_locus_vec <- MakeLocationStrings(CRISPR_df)
  expanded_df <- ExpandList(loci_splits)

  names(expanded_df) <- c("Locus_0MM", "Index")

  index_vec <- expanded_df[["Index"]]
  are_primary <- primary_locus_vec[index_vec] == expanded_df[["Locus_0MM"]]
  expanded_df[["Is_primary_location"]] <- are_primary

  expanded_df[["Intended_Entrez_IDs"]] <- CRISPR_df[["Entrez_ID"]][index_vec]
  expanded_df[["Intended_gene_symbols"]] <- CRISPR_df[["Gene_symbol"]][index_vec]
  expanded_df[["Num_loci"]] <- CRISPR_df[["Num_0MM"]][index_vec]

  are_not_NA <- !(is.na(expanded_df[["Locus_0MM"]]))
  unique_loci_vec <- unique(expanded_df[["Locus_0MM"]][are_not_NA])
  unique_loci_df <- data.frame("Locus_0MM" = unique_loci_vec,
                               LocationStringToDf(unique_loci_vec),
                               stringsAsFactors = FALSE
                               )
  unique_loci_df[["Cut_location"]] <- GetCutLocations(unique_loci_df)
  results_list <- list(
    "expanded_df" = expanded_df,
    "unique_df"   = unique_loci_df
  )
  return(results_list)
}



StandardLocationString <- function(ranges_df) {
  paste0(ranges_df[, "Chromosome"],
         ":",
         ranges_df[, "Start"],
         "-",
         ranges_df[, "End"],
         ":",
         ranges_df[, "Strand"]
         )
}




# Do stuff ----------------------------------------------------------------


if (FALSE) {
  CRISPRa_splits <- strsplit(merged_replaced_CRISPRa_df[["Locations_0MM"]], "; ", fixed = TRUE)
  CRISPRa_splits[is.na(CRISPRa_splits)] <- list(integer(0))
  stopifnot(identical(merged_replaced_CRISPRa_df[["Num_0MM"]], unname(lengths(CRISPRa_splits))))
}




example_entrezs <- unique(merged_replaced_CRISPRa_df[["Entrez_ID"]])[1:1000]
example_entrezs <- unlist(sublibraries_all_entrezs_list, use.names = FALSE)
are_example_genes <- merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% example_entrezs

CRISPR_df <- merged_replaced_CRISPRa_df[are_example_genes, ]


stopifnot(!(anyNA(CRISPR_df[["Entrez_ID"]])))
stopifnot(!(any(grepl(",", CRISPR_df[["Entrez_ID"]], fixed = TRUE))))


split_results <- Process0MMLoci(CRISPR_df)
unique_loci_df <- split_results[["unique_df"]]
expanded_0MM_df <- split_results[["expanded_df"]]

sgRNA_GPos_object <- LocationsDfToGPosObject(unique_loci_df)

unique_loci_df[["Guide_locus"]] <- StandardLocationString(unique_loci_df)



TSS_df <- PrepareTSSDf(all_TSS_df)


TSS_GRanges_object <- RangesDfToGRangesObject(TSS_df)
TSS_df[["Gene_locus"]] <- StandardLocationString(TSS_df)


rename_TSS_columns <- c(
  "Strand"       = "TSS_strand",
  "Entrez_IDs"   = "Affected_Entrez_IDs",
  "Gene_symbols" = "Affected_gene_symbols",
  "TSS"          = "TSS_position",
  "Score"        = "FANTOM5_score",
  "Source"       = "TSS_source"
)

for (column_name in names(rename_TSS_columns)) {
  names(TSS_df)[names(TSS_df) == column_name] <- rename_TSS_columns[[column_name]]
}

hits_object <- findOverlaps(sgRNA_GPos_object,
                            TSS_GRanges_object,
                            ignore.strand = TRUE,
                            select = "all"
                            )


anyNA(TSS_df[["Chromosome"]])

TSS_columns <- c(
  "Gene_locus", "Affected_Entrez_IDs", "Number_of_Entrez_IDs",
  "Affected_gene_symbols", "Is_main_TSS", "TSS_source", "FANTOM5_score",
  "Chromosome", "TSS_strand", "TSS_position"
)

unique_loci_aligned_df <- unique_loci_df[queryHits(hits_object), ]
TSS_aligned_df <- TSS_df[subjectHits(hits_object), ]

are_identical <- mapply(identical,
                        unique_loci_aligned_df[["Chromosome"]],
                        TSS_aligned_df[["Chromosome"]]
                        )
identical(unique_loci_aligned_df[["Chromosome"]],
          TSS_aligned_df[["Chromosome"]]
          )




combined_df <- data.frame(
  unique_loci_df[queryHits(hits_object), c("Locus_0MM", "Guide_locus")],
  TSS_df[subjectHits(hits_object), TSS_columns],
  stringsAsFactors = FALSE,
  row.names = NULL
)


# matches_vec <- match(expanded_0MM_df[["Locus_0MM"]], combined_df[["Locus_0MM"]])
#
# combined_0MM_df <- data.frame(
#   expanded_0MM_df,
#   combined_df[matches_vec, names(combined_df) != "Locus_0MM"],
#   stringsAsFactors = FALSE,
#   row.names = NULL
# )

# combined_0MM_df <- combined_0MM_df[, combined_columns]

are_to_retain <- (expanded_0MM_df[["Locus_0MM"]] %in% combined_df[["Locus_0MM"]]) |
                 (!(duplicated(expanded_0MM_df[["Index"]])))


indices_list <- split(
  seq_len(nrow(combined_df)),
  factor(combined_df[["Locus_0MM"]],
         levels = unique(expanded_0MM_df[["Locus_0MM"]][are_to_retain]),
         exclude = c()
         )
)

indices_list[lengths(indices_list) == 0] <- list(NA_integer_)

matches_vec <- match(expanded_0MM_df[["Locus_0MM"]], names(indices_list))

num_reps_vec <- lengths(indices_list)[matches_vec]

expanded_indices <- rep(seq_len(nrow(expanded_0MM_df)),
                        times = num_reps_vec
                        )
combined_indices <- unlist(indices_list[matches_vec], use.names = FALSE)

full_df <- data.frame(expanded_0MM_df[expanded_indices, ],
                      combined_df[combined_indices, names(combined_df) != "Locus_0MM"],
                      stringsAsFactors = FALSE,
                      row.names = NULL
                      )

combined_columns <- c(
  "Locus_0MM",
  "Index", "Is_primary_location", "Intended_Entrez_IDs", "Affected_Entrez_IDs", "Number_of_Entrez_IDs",
  "Intended_gene_symbols", "Affected_gene_symbols", "Num_loci",
  "Guide_locus", "Gene_locus", "Chromosome",
  "Is_main_TSS", "TSS_source", "TSS_strand", "TSS_position"
)

full_df <- full_df[, combined_columns]



unique_entrezs <- unique(full_df[["Affected_Entrez_IDs"]])
min_unique_entrezs <- GetMinEntrez(unique_entrezs)
min_entrezs_vec <- min_unique_entrezs[match(full_df[["Affected_Entrez_IDs"]], unique_entrezs)]


new_order <- order(
  full_df[["Index"]],
  full_df[["Intended_Entrez_IDs"]] != full_df[["Affected_Entrez_IDs"]],
  min_entrezs_vec,
  match(full_df[["Locus_0MM"]], full_df[["Locus_0MM"]])
)



reordered_by_gene_df <- full_df
reordered_by_gene_df[["Original_order"]] <- seq_len(nrow(reordered_by_gene_df))
reordered_by_gene_df <- full_df[new_order, ]
row.names(reordered_by_gene_df) <- NULL


SplitAndRejoin <- function(long_df, guides_df, use_column) {

  splits_list <- strsplit(long_df[, use_column], ", ", fixed = TRUE)

  joined_list <- tapply(splits_list,
                        factor(long_df[["Index"]], levels = seq_len(nrow(guides_df))),
                        function(x) {
                          results_vec <- unique(unlist(x, use.names = FALSE))
                          results_vec <- results_vec[!(is.na(results_vec))]
                          return(results_vec)
                        },
                        simplify = FALSE
                        )

  joined_vec <- vapply(joined_list,
                       function(x) {
                         if (length(x) == 0) {
                           NA_character_
                         } else {
                           paste0(x, collapse = ", ")
                         }
                       }, "")
  results_list <- list(
    "list"   = joined_list,
    "vector" = joined_vec
  )
  return(results_list)
}


affected_symbols_results <- SplitAndRejoin(reordered_by_gene_df, CRISPR_df, "Affected_gene_symbols")
affected_entrezs_results <- SplitAndRejoin(reordered_by_gene_df, CRISPR_df, "Affected_Entrez_IDs")

intended_symbols_results <- SplitAndRejoin(reordered_by_gene_df, CRISPR_df, "Intended_gene_symbols")
intended_entrezs_results <- SplitAndRejoin(reordered_by_gene_df, CRISPR_df, "Intended_Entrez_IDs")

first_indices <- match(unique(full_df[["Index"]]),
                       full_df[["Index"]]
                       )

intended_entrezs <- full_df[["Intended_Entrez_IDs"]][first_indices]
intended_symbols <- full_df[["Intended_gene_symbols"]][first_indices]
total_loci       <- full_df[["Num_loci"]][first_indices]

identical(intended_entrezs, unname(intended_entrezs_results[["vector"]]))


identical(lengths(affected_symbols_results[["list"]]),
          lengths(affected_entrezs_results[["list"]])
          )


summary_df <- data.frame(
  "Intended_Entrez_ID"    = intended_entrezs,
  "Affected_Entrez_IDs"   = affected_entrezs_results[["vector"]],
  "Intended_gene_symbol"  = intended_symbols,
  "Affected_gene_symbols" = affected_symbols_results[["vector"]],
  "Num_affected_genes"    = lengths(affected_entrezs_results[["list"]]),
  "Total_loci"            = total_loci,
  stringsAsFactors        = FALSE
)


index_loci_vec <- paste0(full_df[["Index"]], "__", full_df[["Locus_0MM"]])

are_duplicated <- duplicated(index_loci_vec)

locus_fac <- factor(index_loci_vec, levels = index_loci_vec[!(are_duplicated)])

are_NA_entrezs <- is.na(full_df[["Affected_Entrez_IDs"]])


entrezs_locus_splits <- split(full_df[["Affected_Entrez_IDs"]][!(are_NA_entrezs)],
                              locus_fac[!(are_NA_entrezs)]
                              )

index_vec <- full_df[["Index"]][!(duplicated(index_loci_vec))]

second_splits <- split(entrezs_locus_splits, index_vec)
second_splits <- lapply(second_splits, function(x) x[lengths(x) > 0])


have_targets <- tapply(!(are_NA_entrezs), locus_fac, any)
summary_df[["Loci_with_targets"]] <- as.vector(tapply(have_targets, index_vec, sum))
stopifnot(identical(unname(lengths(second_splits)), summary_df[["Loci_with_targets"]]))


summary_df[["Loci_targeting_intended_gene"]] <- vapply(seq_len(nrow(CRISPR_df)),
                                                       function(x) {
                                                         intended_entrez <- summary_df[["Intended_Entrez_ID"]][[x]]
                                                         sum(vapply(second_splits[[x]], function(y) intended_entrez %in% y, logical(1)))
                                                       }, integer(1))



distinct_splits <- lapply(second_splits, unique)

have_multiple_loci <- summary_df[["Loci_with_targets"]] > 1

distinct_splits[have_multiple_loci] <- lapply(distinct_splits[have_multiple_loci],
                                              MergeCommonElements
                                              )

split_distinct_splits <- lapply(distinct_splits,
                                function(x) lapply(x, function(y) unlist(strsplit(y, ", ", fixed = TRUE)))
                                )


summary_df[["Distinct_loci"]] <- lengths(distinct_splits, use.names = FALSE)

summary_df[["Distinct_loci_with_multiple_targets"]] <- vapply(split_distinct_splits,
                                                              function(x) sum(lengths(x) > 1),
                                                              integer(1)
                                                              )






TSS_strands_result <- SplitAndRejoin(reordered_by_gene_df, CRISPR_df, "TSS_strand")

summary_df[["Affected_genes_strand"]] <- ifelse(TSS_strands_result[["vector"]] == "+",
                                                "+",
                                                ifelse(TSS_strands_result[["vector"]] == "-",
                                                       "-",
                                                       "Both"
                                                       )
                                                )




head(summary_df[summary_df[["Distinct_loci"]] < summary_df[["Loci_with_targets"]], ])


gooo


# symbols_targeted_list <- tapply(symbol_splits,
#                                 factor(full_df[["Index"]], levels = seq_len(nrow(CRISPR_df))),
#                                 function(x) {
#                                   results_vec <- unique(unlist(x, use.names = FALSE))
#                                   results_vec <- results_vec[!(is.na(results_vec))]
#                                   return(results_vec)
#                                 },
#                                 simplify = FALSE
#                                 )
#
# symbols_targeted_vec <- vapply(symbols_targeted_list,
#                                function(x) {
#                                  if (length(x) == 0) {
#                                    NA_character_
#                                  } else {
#                                    paste0(x, collapse = ", ")
#                                  }
#                                }, "")
#
#
# strands_targeted <- tapply(symbol_splits,
#                            factor(full_df[["Index"]], levels = seq_len(nrow(CRISPR_df))),
#                            function(x) {
#                               results_vec <- unique(unlist(x, use.names = FALSE))
#                               results_vec <- results_vec[!(is.na(results_vec))]
#                               return(results_vec)
#                            })


table(seq_len(nrow(CRISPR_df)) %in% full_df[["Index"]])







hits_entrezs_vec  <- TSS_df[["Entrez_ID"]][subjectHits(hits_object)]
hits_symbols_vec  <- TSS_df[["Gene_symbol"]][subjectHits(hits_object)]
hits_strands_vec <- TSS_df[["Strand"]][subjectHits(hits_object)]
hits_TSS_vec <- TSS_df[["TSS"]][subjectHits(hits_object)]

new_order <- order(queryHits(hits_object),
                   GetMinEntrez(hits_entrezs_vec)
                   )

hits_entrezs_vec <- hits_entrezs_vec[new_order]
hits_symbols_vec  <- hits_symbols_vec[new_order]
hits_strands_vec <- hits_strands_vec[new_order]
hits_TSS_vec <- hits_TSS_vec[new_order]

query_vec <- queryHits(hits_object)[new_order]
query_fac <- factor(query_vec, levels = seq_len(nrow(unique_locations_df)))

sub_df <- data.frame(
  "Gene_IDs" = hits_entrezs_vec,
  "Symbols"  = hits_symbols_vec,
  "Strands"  = hits_strands_vec,
  "Group"    = query_fac,
  stringsAsFactors = FALSE
)

are_duplicated <- duplicated(sub_df)


entrezs_vec <- tapply(hits_entrezs_vec[!(are_duplicated)],
                      query_fac[!(are_duplicated)],
                      paste0,
                      collapse = ", "
                      )
symbols_vec <- tapply(hits_symbols_vec[!(are_duplicated)],
                      query_fac[!(are_duplicated)],
                      paste0,
                      collapse = ", "
                      )
strands_vec <- tapply(hits_strands_vec[!(are_duplicated)],
                  query_fac[!(are_duplicated)],
                  paste0,
                  collapse = ", "
                  )
TSS_vec <- tapply(hits_TSS_vec[!(are_duplicated)],
                  query_fac[!(are_duplicated)],
                  paste0,
                  collapse = ", "
                  )

results_df <- data.frame(
  "Gene_IDs"       = entrezs_vec,
  "Gene_symbols"   = symbols_vec,
  "Strands"        = strands_vec,
  "TSS"            = TSS_vec,
  "Num_genes"      = as.integer(table(query_fac[!(are_duplicated)])),
  stringsAsFactors = FALSE
)




matches_vec <- match(expanded_df[["Value"]], expanded_df[["Value"]])


final_results_df <- results_df[matches_vec, ]







