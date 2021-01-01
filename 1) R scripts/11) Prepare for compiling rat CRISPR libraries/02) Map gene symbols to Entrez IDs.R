### 18 November 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files", "10) Rat - General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Rn.eg.db Bioconductor database.RData"))





# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz
# on 17 November 2020

NCBI_Rn_info_df <- read.table(file.path(CRISPR_input_directory, "Rat genome", "NCBI", "Rattus_norvegicus.gene_info"),
                              sep = "\t", quote = "", header = TRUE, row.names = NULL,
                              fill = TRUE, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE
                              )




# Define functions --------------------------------------------------------

ReverseList <- function(my_list) {
  my_entries <- unlist(my_list, recursive = FALSE, use.names = FALSE)
  my_names <- rep(names(my_list), times = lengths(my_list, use.names = FALSE))
  names(my_names) <- my_entries
  return(my_names)
}

ReverseAndRelist <- function(my_list) {
  reversed_list <- ReverseList(my_list)
  results_list <- sapply(unique(names(reversed_list)), function(x) unique(unlist(reversed_list[names(reversed_list) == x])), simplify = FALSE)
  return(results_list)
}

OrderEntrezsList <- function(my_list) {
  results_list <- lapply(my_list, function(x) x[order(as.integer(x))])
  names(results_list) <- names(my_list)
  return(results_list)
}

SynonymsToEntrez <- function(my_symbol) {
  if (my_symbol %in% names(entrez_to_symbol_list)) {
    my_entrez <- entrez_to_symbol_list[[my_entrez]]
  }
}

ConsolidateList <- function(input_list) {
  sapply(unique(names(input_list)), function(x) unlist(input_list[names(input_list) == x], use.names = FALSE), simplify = FALSE)
}


NamesToAllCaps <- function(input_list) {
  names(input_list) <- toupper(names(input_list))
  return(input_list)
}







# Convert alias_ and symbol_to_entrez_list to all capital letters ---------

symbol_to_entrez_list <- ConsolidateList(NamesToAllCaps(symbol_to_entrez_list))
alias_to_entrez_list  <- OrderEntrezsList(ConsolidateList(NamesToAllCaps(alias_to_entrez_list)))







# Perform checks on alias_ and symbol_to_entrez_list ----------------------

symbol_contains_alias <- vapply(names(symbol_to_entrez_list), function(x) all(symbol_to_entrez_list[[x]] %in% alias_to_entrez_list[[x]]), logical(1))
stopifnot(all(symbol_contains_alias))

stopifnot(identical(symbol_to_entrez_list, OrderEntrezsList(symbol_to_entrez_list)))






# Check for conflicting mappings of Entrez IDs to gene symbols ------------

any(duplicated(names(entrez_to_symbol_list)))
any(duplicated(names(NCBI_Rn_info_df[["GeneID"]])))

unique_entrez_IDs <- unique(c(unlist(unname(symbol_to_entrez_list)), names(entrez_to_symbol_list), as.character(NCBI_Rn_info_df[["GeneID"]])))
unique_entrez_IDs <- unique_entrez_IDs[order(as.integer(unique_entrez_IDs))]

chromosome_vec <- vapply(unique_entrez_IDs, function(x) {
  if (x %in% names(entrez_to_chromosome_list)) {
    chromosome_vec <- ifelse(entrez_to_chromosome_list[[x]] == "MT",
                             "chrM",
                             ifelse(entrez_to_chromosome_list[[x]] == "Un", NA_character_, paste0("chr", entrez_to_chromosome_list[[x]]))
                             )
    if (all(is.na(chromosome_vec))) {
      return(NA_character_)
    } else {
      return(paste0(chromosome_vec, collapse = ", "))
    }
  } else {
    return(NA_character_)
  }
}, "")

Entrez_map_df <- data.frame(
  "Entrez_ID"           = unique_entrez_IDs,
  "Symbol_Org_Rn_eg_db" = EntrezIDsToSymbols(unique_entrez_IDs),
  "Symbol_NCBI_Rn_info" = NCBI_Rn_info_df[["Symbol"]][match(unique_entrez_IDs, NCBI_Rn_info_df[["GeneID"]])],
  "Chromosome"          = chromosome_vec,
  stringsAsFactors = FALSE
)

are_identical <- mapply(identical, Entrez_map_df[["Symbol_Org_Rn_eg_db"]], Entrez_map_df[["Symbol_NCBI_Rn_info"]])

Entrez_map_df[!(are_identical), ]

Entrez_map_df[is.na(Entrez_map_df[["Symbol_Org_Rn_eg_db"]]) & !(is.na(Entrez_map_df[["Symbol_NCBI_Rn_info"]])), ]

Entrez_map_df[!(are_identical) & !(is.na(Entrez_map_df[["Symbol_Org_Rn_eg_db"]])) & !(is.na(Entrez_map_df[["Symbol_NCBI_Rn_info"]])), ]







# Create a list of gene synonyms from NCBI's Rn.gene_info -----------------

are_empty <- NCBI_Rn_info_df[["Synonyms"]] == "-"

synonyms_list <- strsplit(NCBI_Rn_info_df[["Synonyms"]][!(are_empty)], "|", fixed = TRUE)

entrez_to_synonyms_list <- synonyms_list
names(entrez_to_synonyms_list) <- NCBI_Rn_info_df[["GeneID"]][!(are_empty)]

synonyms_to_entrez_list <- ReverseAndRelist(entrez_to_synonyms_list)

# symbol_to_synonyms_list <- synonyms_list
# names(symbol_to_synonyms_list) <- NCBI_Rn_info_df[["Symbol"]][!(are_empty)]
# synonyms_to_symbol_list <- ReverseAndRelist(symbol_to_synonyms_list)




# Convert the lists from NCBI's Rn.info file to all capital letters -------

synonyms_to_entrez_list <- OrderEntrezsList(NamesToAllCaps(ConsolidateList(synonyms_to_entrez_list)))
stopifnot(identical(synonyms_to_entrez_list, OrderEntrezsList(synonyms_to_entrez_list)))





# Collect all known symbols for which mappings exist ----------------------

collected_symbols <- sort(unique(c(names(symbol_to_entrez_list),
                                   names(alias_to_entrez_list),
                                   toupper(unlist(entrez_to_symbol_list, recursive = FALSE, use.names = FALSE)),
                                   toupper(NCBI_Rn_info_df[["Symbol"]]),
                                   names(synonyms_to_entrez_list)
                                   )
                                 )
                          )





# Build a data frame for assigning gene symbols to Entrez IDs -------------

combined_list <- sapply(collected_symbols, function(x) unique(c(alias_to_entrez_list[[x]], synonyms_to_entrez_list[[x]])))

symbols_df <- data.frame("Symbol"            = collected_symbols,
                         "Chosen_entrez"     = NA_character_,
                         "Chosen_via"        = NA_character_,
                         "Number_OrgRn"      = vapply(collected_symbols, function(x) if (x %in% names(symbol_to_entrez_list)) length(symbol_to_entrez_list[[x]]) else NA_integer_, integer(1)),
                         "Number_all"        = lengths(combined_list),
                         "Entrez_IDs_OrgRn"  = vapply(collected_symbols, function(x) if (x %in% names(symbol_to_entrez_list)) paste0(symbol_to_entrez_list[[x]], collapse = ", ") else NA_character_, ""),
                         "Entrez_IDs_all"    = vapply(combined_list, paste0, collapse = ", ", ""),
                         stringsAsFactors = FALSE,
                         row.names = NULL
                         )


unique_OrgRn_symbols_list <- symbol_to_entrez_list[lengths(symbol_to_entrez_list) == 1]
unique_OrgRn_alias_list   <- alias_to_entrez_list[lengths(alias_to_entrez_list) == 1]
unique_NCBI_synonyms_list <- synonyms_to_entrez_list[lengths(synonyms_to_entrez_list) == 1]

intersect_OrgRn_NCBI_list <- sapply(intersect(names(alias_to_entrez_list), names(synonyms_to_entrez_list)),
                                    function(x) intersect(alias_to_entrez_list[[x]], synonyms_to_entrez_list[[x]])
                                    , simplify = FALSE
                                    )
unique_intersect_OrgRn_NCBI_list <- intersect_OrgRn_NCBI_list[lengths(intersect_OrgRn_NCBI_list) == 1]

symbols_df[["Unique_in_OrgRn_symbol"]] <- vapply(collected_symbols, function(x) if (x %in% names(unique_OrgRn_symbols_list)) unique_OrgRn_symbols_list[[x]] else NA_character_, "")
symbols_df[["Present_in_NCBI_symbol"]] <- as.character(NCBI_Rn_info_df[["GeneID"]][match(collected_symbols, toupper(NCBI_Rn_info_df[["Symbol"]]))])
symbols_df[["Unambiguous_alias"]] <- vapply(collected_symbols, function(x) {
  unique_OrgRn_alias_lookup <- unique_OrgRn_alias_list[[x]]
  unique_NCBI_synonyms_lookup <- unique_NCBI_synonyms_list[[x]]
  if (is.null(unique_OrgRn_alias_lookup) || is.null(unique_NCBI_synonyms_lookup)) {
    return(NA_character_)
  } else {
    unique_symbols <- unique(unlist(c(unique_OrgRn_alias_lookup, unique_NCBI_synonyms_lookup), use.names = FALSE))
    if (length(unique_symbols) == 1) {
      return(unique_symbols)
    } else {
      return(NA_character_)
    }
  }
}, "")

symbols_df[["Unique_in_overlap"]]      <- vapply(collected_symbols, function(x) if (x %in% names(unique_intersect_OrgRn_NCBI_list)) unique_intersect_OrgRn_NCBI_list[[x]] else NA_character_, "")
symbols_df[["Unique_in_OrgRn_alias"]]  <- vapply(collected_symbols, function(x) if (x %in% names(unique_OrgRn_alias_list)) unique_OrgRn_alias_list[[x]] else NA_character_, "")
symbols_df[["Unique_in_NCBI_synonym"]] <- vapply(collected_symbols, function(x) if (x %in% names(unique_NCBI_synonyms_list)) unique_NCBI_synonyms_list[[x]] else NA_character_, "")

symbols_df[["Min_OrgRn_symbol"]]  <- vapply(collected_symbols, function(x) {
  my_lookup <- symbol_to_entrez_list[[x]]
  if (is.null(my_lookup)) {
    return(NA_character_)
  } else {
    return(my_lookup[[1]])
  }
}, "")

symbols_df[["Min_OrgRn_alias"]]  <- vapply(collected_symbols, function(x) {
  my_lookup <- alias_to_entrez_list[[x]]
  if (is.null(my_lookup)) {
    return(NA_character_)
  } else {
    return(my_lookup[[1]])
  }
}, "")

symbols_df[["Min_NCBI_synonym"]] <- vapply(collected_symbols, function(x) {
  my_lookup <- synonyms_to_entrez_list[[x]]
  if (is.null(my_lookup)) {
    return(NA_character_)
  } else {
    return(my_lookup[[1]])
  }
}, "")


symbols_assignment_map <- c(
  "Unique_in_OrgRn_symbol" = "Only one mapping in Org.Rn.egSYMBOL",
  "Present_in_NCBI_symbol" = "No unique mapping in Org.Rn.egSYMBOL, but found in the Symbol column of Rattus_norvegicus.gene_info",
  "Unambiguous_alias"      = "One and only one mapping in both Org.Rn.egALIAS and Rattus_norvegicus.gene_info",
  "Unique_in_overlap"      = "Only one mapping that is coRnon to both Org.Rn.egALIAS and Rattus_norvegicus.gene_info",
  "Unique_in_OrgRn_alias"  = "Only present in Org.Rn.egALIAS, and only one mapping found",
  "Unique_in_NCBI_synonym" = "Only present in Rattus_norvegicus.gene_info, and only one mapping found",
  "Min_OrgRn_symbol"       = "Maps to multiple IDs in Org.Rn.egSYMBOL; the first (lowest) ID was chosen",
  "Min_OrgRn_alias"        = "Maps to multiple IDs in Org.Rn.egALIAS; the first (lowest) ID was chosen",
  "Min_NCBI_synonym"       = "Not found in Org.Rn.egALIAS, but maps to multiple IDs in Rattus_norvegicus.gene_info; the first (lowest) ID was chosen"
)

column_indices <- match(names(symbols_assignment_map), names(symbols_df))

are_available_mat <- !(apply(symbols_df[, column_indices], 2, function(x) is.na(x)))
which_first <- apply(are_available_mat, 1, function(x) if (any(x)) which(x)[[1]] else NA_integer_)

symbols_df[["Chosen_via"]] <- factor(unname(symbols_assignment_map[which_first]), levels = unname(symbols_assignment_map))
symbols_df[["Chosen_entrez"]] <- vapply(seq_along(which_first), function(x) symbols_df[x, column_indices[[which_first[[x]]]]], "")

symbol_to_entrez_df <- symbols_df
entrez_to_symbol_df <- Entrez_map_df





# Save data ---------------------------------------------------------------

save(list = c("symbol_to_entrez_df", "entrez_to_symbol_df"),
     file = file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData")
     )



















