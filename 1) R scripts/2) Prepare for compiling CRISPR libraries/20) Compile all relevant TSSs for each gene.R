### 31st March 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder
source(file.path(general_functions_directory, "28) Merging the FANTOM5 and BioMart TSS data.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(general_RData_directory, "19) Compile the information on gene type.RData"))





# Check for inconsistent chromosome data ----------------------------------

head(FANTOM5_df)

FANTOM5_entrez_matches <- match(FANTOM5_df[["Entrez_ID"]],
                                entrez_to_symbol_df[["Entrez_ID"]]
                                )


# FANTOM5_df[["Symbol_from_entrez"]] <- MapToEntrezs(entrez_IDs_vec = FANTOM5_df[["Entrez_ID"]])[["Gene_symbol"]]

# grep(", ", FANTOM5_df[["Entrez_ID"]], fixed = TRUE, value = TRUE)
# grep(", ", FANTOM5_df[["Gene_symbol"]], fixed = TRUE, value = TRUE)





# Use the chromosome to disambiguate FANTOM5 Entrez IDs -------------------

have_entrezs <- !(is.na(FANTOM5_df[["Entrez_ID"]]))
unique_entrez_ID_strings <- unique(FANTOM5_df[["Entrez_ID"]][have_entrezs])

ambiguous_entrezs_vec <- unique(grep(", ", FANTOM5_df[["Entrez_ID"]], fixed = TRUE, value = TRUE))
ambiguous_entrezs_list <- strsplit(ambiguous_entrezs_vec, ", ", fixed = TRUE)

expanded_entrezs_df <- ExpandList(ambiguous_entrezs_list)
colnames(expanded_entrezs_df)[[1]] <- "Entrez_ID"

entrez_matches <- match(expanded_entrezs_df[["Entrez_ID"]],
                        entrez_to_symbol_df[["Entrez_ID"]]
                        )
expanded_entrezs_df[["Entrez_chromosome"]] <- entrez_to_symbol_df[["Chromosome"]][entrez_matches]
expanded_entrezs_df[["Entrez_chromosome"]] <- ifelse(expanded_entrezs_df[["Entrez_chromosome"]] == "chrX, chrY",
                                              "chrX",
                                               expanded_entrezs_df[["Entrez_chromosome"]]
                                              )

matches_vec <- match(ambiguous_entrezs_vec, FANTOM5_df[["Entrez_ID"]])
original_chromosomes_vec <- FANTOM5_df[["Chromosome"]][matches_vec]

expanded_entrezs_df[["Original_chromosome"]] <- original_chromosomes_vec[expanded_entrezs_df[["List_index"]]]
expanded_entrezs_df[["Are_matching"]] <- expanded_entrezs_df[["Original_chromosome"]] == expanded_entrezs_df[["Entrez_chromosome"]]

are_matching_list <- split(expanded_entrezs_df[["Are_matching"]], expanded_entrezs_df[["List_index"]])
entrez_IDs_list <- split(expanded_entrezs_df[["Entrez_ID"]], expanded_entrezs_df[["List_index"]])
disambiguated_entrezs_list <- lapply(seq_along(ambiguous_entrezs_vec),
                                     function(x) entrez_IDs_list[[x]][which(are_matching_list[[x]])]
                                     )

disambiguated_entrezs_vec <- vapply(disambiguated_entrezs_list,
                                    function(x) paste0(x, collapse = ", "),
                                    ""
                                    )

FANTOM5_df[["Entrez_ID"]][matches_vec] <- disambiguated_entrezs_vec







# Merge the FANTOM5 and BioMart data frames -------------------------------

FANTOM5_have_entrez <- !(is.na(FANTOM5_df[["Entrez_ID"]]))
BioMart_have_entrez <- !(is.na(BioMart_filtered_df[["Entrez_ID"]]))

all_TSS_df <- MergeTSSData(FANTOM5_df[FANTOM5_have_entrez, ],
                           BioMart_filtered_df[BioMart_have_entrez, ],
                           common_columns = c("Entrez_ID",
                                              "Chromosome",
                                              "Strand"
                                              ),
                           use_legacy_TSS_order = FALSE
                           )





# Add the original gene symbols -------------------------------------------

MakeTSSIDs <- function(TSS_df, TSS_column = "TSS") {
  paste0(TSS_df[["Entrez_ID"]], " | ",
         TSS_df[["Chromosome"]], " | ",
         TSS_df[["Strand"]], " | ",
         TSS_df[[TSS_column]]
         )
}

FANTOM5_TSS_IDs <- MakeTSSIDs(FANTOM5_df, TSS_column = "TSS_start")[FANTOM5_have_entrez]
BioMart_TSS_IDs <- MakeTSSIDs(BioMart_filtered_df)[BioMart_have_entrez]

long_TSS_IDs <- c(FANTOM5_TSS_IDs, BioMart_TSS_IDs)

all_original_symbols <- c(FANTOM5_df[["Gene_symbol"]][FANTOM5_have_entrez],
                          BioMart_filtered_df[["Gene_symbol"]][BioMart_have_entrez]
                          )

all_original_symbols_list <- split(all_original_symbols, long_TSS_IDs)
all_original_symbols_list <- lapply(all_original_symbols_list, unique)
all_original_symbols_vec <- vapply(all_original_symbols_list, paste0, collapse = ", ", "")
all_original_symbols_list <- strsplit(all_original_symbols_vec, ", ", fixed = TRUE)
all_original_symbols_list <- lapply(all_original_symbols_list, unique)
all_original_symbols_vec <- vapply(all_original_symbols_list, paste0, collapse = ", ", "")

combined_TSS_IDs <- MakeTSSIDs(all_TSS_df)

TSS_ID_matches <- match(combined_TSS_IDs, names(all_original_symbols_vec))
all_TSS_df[["Original_symbol"]] <- all_original_symbols_vec[TSS_ID_matches]






# Collect information for Entrez IDs --------------------------------------

entrezs_list <- strsplit(all_TSS_df[["Entrez_ID"]], ", ", fixed = TRUE)

entrezs_vec <- unique(unlist(entrezs_list, use.names = FALSE))
entrezs_vec <- entrezs_vec[order(as.integer(entrezs_vec))]

head(entrez_to_symbol_df)

entrez_matches <- match(entrezs_vec, entrez_to_symbol_df[["Entrez_ID"]])
symbols_vec <- ifelse(is.na(entrez_to_symbol_df[["Symbol_Org_Hs_eg_db"]]),
                      entrez_to_symbol_df[["Symbol_NCBI_Hs_info"]],
                      entrez_to_symbol_df[["Symbol_Org_Hs_eg_db"]]
                      )[entrez_matches]
gene_type_matches <- match(entrezs_vec, entrez_to_gene_type_df[["Entrez_ID"]])

entrez_info_df <- data.frame(
  "Entrez_gene_ID"    = entrezs_vec,
  "Gene_symbol"       = symbols_vec,
  "Entrez_chromosome" = entrez_to_symbol_df[["Chromosome"]][entrez_matches],
  "Gene_type"         = entrez_to_gene_type_df[["Consensus_type"]][gene_type_matches],
  stringsAsFactors    = FALSE
)

entrez_info_df[["Entrez_chromosome"]] <- ifelse(entrez_info_df[["Entrez_chromosome"]] == "chrX, chrY",
                                                "chrX",
                                                entrez_info_df[["Entrez_chromosome"]]
                                                )



# Add the data linked to the Entrez ID to all_TSS_df ----------------------

expanded_entrezs_df <- ExpandList(entrezs_list)
expanded_matches <- match(expanded_entrezs_df[["Value"]],
                          entrez_info_df[["Entrez_gene_ID"]]
                          )

expanded_info_df <- entrez_info_df[expanded_matches, ]

expanded_groups <- factor(expanded_entrezs_df[["List_index"]])

info_list <- sapply(names(expanded_info_df), function(x) {
  tapply(expanded_info_df[[x]],
         expanded_groups,
         function(y) {
           if (all(is.na(y))) {
             NA_character_
           } else {
             y_vec <- unique(y[!(is.na(y))])
             paste0(unique(y_vec), collapse = ", ")
           }
         }
         )
}, simplify = FALSE)

all_TSS_df <- data.frame(
  all_TSS_df,
  do.call(cbind, info_list),
  stringsAsFactors = FALSE,
  row.names = NULL
)

stopifnot(identical(all_TSS_df[["Entrez_ID"]],
                    all_TSS_df[["Entrez_gene_ID"]]
                    )
          )
all_TSS_df <- all_TSS_df[, colnames(all_TSS_df) != "Entrez_gene_ID"]

all_TSS_columns <- c(
  "Entrez_ID", "Gene_symbol", "Original_symbol",
  "Entrez_chromosome", "Chromosome",
  "Gene_type",
  "Source" ,"Strand", "TSS", "Score"
)
all_TSS_df <- all_TSS_df[, all_TSS_columns]





# Drop some TSS entries that seem to be on the wrong chromosome -----------

are_identical <- mapply(identical,
                        all_TSS_df[["Chromosome"]],
                        all_TSS_df[["Entrez_chromosome"]]
                        )
are_not_NA <- !(is.na(all_TSS_df[["Entrez_chromosome"]]))

problematic_entrezs <- all_TSS_df[["Entrez_ID"]][!(are_identical) & are_not_NA]

are_problematic_entrezs <- all_TSS_df[["Entrez_ID"]] %in% problematic_entrezs

problematic_entrezs_vec <- all_TSS_df[["Entrez_ID"]][are_problematic_entrezs]
problematic_entrezs_fac <- factor(problematic_entrezs_vec,
                                  levels = unique(problematic_entrezs_vec)
                                  )
CheckThatFactorIsInOrder(problematic_entrezs_fac)

are_discordant_list <- tapply(seq_along(problematic_entrezs_fac),
                              problematic_entrezs_fac,
                              function(x) {
                                use_indices <- which(are_problematic_entrezs)[x]
                                all_TSS_df[["Chromosome"]][use_indices] != all_TSS_df[["Entrez_chromosome"]][use_indices]
                              }, simplify = FALSE)

are_discordant_chromosome <- unlist(are_discordant_list, use.names = FALSE)

are_all_discordant <- unlist(lapply(are_discordant_list,
                                    function(x) rep(all(x), length(x))
                                    ),
                             use.names = FALSE
                             )

long_are_all_discordant <- rep(FALSE, nrow(all_TSS_df))
long_are_all_discordant[are_problematic_entrezs] <- are_all_discordant

nrow(all_TSS_df)
all_TSS_df <- all_TSS_df[!(are_problematic_entrezs & !(long_are_all_discordant)), ]
nrow(all_TSS_df)
row.names(all_TSS_df) <- NULL

CheckThatFactorIsInOrder(factor(all_TSS_df[["Entrez_ID"]]))




# Identify Entrez IDs for which no unambiguous TSS is available -----------

have_multiple_entrezs <- grepl(", ", all_TSS_df[["Entrez_ID"]], fixed = TRUE)

single_entrezs <- all_TSS_df[["Entrez_ID"]][!(have_multiple_entrezs)]
entrez_splits <- strsplit(all_TSS_df[["Entrez_ID"]], ", ", fixed = TRUE)
all_entrezs <- unique(unlist(entrez_splits, use.names = FALSE))
only_multiple_entrezs <- setdiff(all_entrezs, single_entrezs)

have_only_multiple_entrezs <- vapply(entrez_splits, function(x) {
  all(x %in% only_multiple_entrezs)
}, logical(1))

are_eligible <- !(have_multiple_entrezs) | have_only_multiple_entrezs




# Identify the TSS with the highest FANTOM5 score -------------------------

entrezs_fac <- factor(all_TSS_df[["Entrez_ID"]],
                      levels = unique(all_TSS_df[["Entrez_ID"]])
                      )

all_TSS_df[["Is_main_TSS"]] <- unlist(tapply(
  seq_len(nrow(all_TSS_df)),
  entrezs_fac,
  function(x) {
    score_vec <- all_TSS_df[["Score"]][x]
    results_vec <- rep(FALSE, length(x))
    are_available <- are_eligible[x] & !(is.na(score_vec))
    if (!(any(are_available))) {
      return(results_vec)
    } else {
      max_score <- max(score_vec, na.rm = TRUE)
      are_max <- are_available & (score_vec %in% max_score)
      results_vec[which(are_max)[[1]]] <- TRUE
      if (sum(are_max) > 1) {
        message("\nMore than one maximum score found:")
        print(all_TSS_df[x, ])
        message(paste0("The first TSS was chosen: ",
                       all_TSS_df[["Chromosome"]][x][results_vec], ":",
                       all_TSS_df[["TSS"]][x][results_vec], ":",
                       all_TSS_df[["Strand"]][x][results_vec],
                       "\n"
                       )
                )
      }
      return(results_vec)
    }
  }, simplify = FALSE),
  use.names = FALSE
)





# Choose a TSS for each Entrez ID -----------------------------------------

all_TSS_df[["Is_chosen_TSS"]] <- unlist(tapply(
  seq_len(nrow(all_TSS_df)),
  entrezs_fac,
  function(x) {
    results_vec <- rep(FALSE, length(x))
    if (!(any(are_eligible[x]))) {
      return(results_vec)
    } else {
      are_main <- all_TSS_df[["Is_main_TSS"]][x]
      if (any(are_main)) {
        return(are_main)
      } else {
        results_vec[[1]] <- TRUE # This is assuming the order is already correct, i.e. highest or lowest position, depending on the strand
        return(results_vec)
      }
    }
  }, simplify = FALSE),
  use.names = FALSE
)




# Add a column indicating the consistency of chromosomal locations --------

stopifnot(!(anyNA(all_TSS_df[["Chromosome"]])))
all_TSS_df[["Has_consistent_chromosome"]] <- mapply(identical,
                                                    all_TSS_df[["Chromosome"]],
                                                    all_TSS_df[["Entrez_chromosome"]]
                                                    )




# Save data ---------------------------------------------------------------

save(list = "all_TSS_df",
     file = file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData")
     )





