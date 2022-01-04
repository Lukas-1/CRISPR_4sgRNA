### 19th April 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")
gene_lists_directory    <- file.path(CRISPR_root_directory, "2) Input data", "Mouse gene lists")

membrane_hits_path      <- file.path(gene_lists_directory,
                                     "Membrane heterogeneity",
                                     "Hits list for membrane heterogeneity_SAM library_Mouse_GT1-7.xlsx"
                                     )

synthetic_lethal_path   <- file.path(gene_lists_directory,
                                     "Prions synthetic lethality",
                                     "Hits list.xlsx"
                                     )



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))



# Read in data ------------------------------------------------------------

membrane_df <- data.frame(read_excel(membrane_hits_path),
                          stringsAsFactors = FALSE,
                          check.names = FALSE
                          )

lethality_df <- data.frame(read_excel(synthetic_lethal_path),
                           stringsAsFactors = FALSE,
                           check.names = FALSE
                           )




# Collect Entrez IDs from the membrane screen -----------------------------

membrane_symbols_vec <- c(membrane_df[["Enriched Hits"]], membrane_df[["Decreased Hits"]])
membrane_symbols_vec <- membrane_symbols_vec[!(is.na(membrane_symbols_vec))]

membrane_mapped_df <- MapToEntrezs(symbols_vec = membrane_symbols_vec, is_mouse = TRUE)

stopifnot(all(membrane_mapped_df[["Original_symbol"]] == ""))


are_enriched <- membrane_mapped_df[["Gene_symbol"]] %in% membrane_df[["Enriched Hits"]]
are_depleted <- membrane_mapped_df[["Gene_symbol"]] %in% membrane_df

stopifnot(!(any(are_enriched & are_depleted)))

membrane_het_df <- data.frame(
  membrane_mapped_df,
  "Category" = ifelse(are_enriched, "Enriched", "Depleted"),
  "Num_occurrences" = table(membrane_mapped_df[["Entrez_ID"]])[membrane_mapped_df[["Entrez_ID"]]],
  stringsAsFactors = FALSE
)

membrane_het_df <- membrane_het_df[!(duplicated(membrane_het_df)), ]
row.names(membrane_het_df) <- NULL




# Collect Entrez IDs from the prion synthetic lethality screen ------------

lethal_symbols_vec <- lethality_df[, "Hits in total_decreased_P2 P4 P8 and P2 P4 and P4 P8 overlapped preference"]
lethal_mapped_df <- MapToEntrezs(symbols_vec = lethal_symbols_vec, is_mouse = TRUE)
prions_lethality_df <- lethal_mapped_df[, 1:2]

prions_lethality_df[["Gene_rank"]] <- seq_len(nrow(prions_lethality_df))

prions_category_vec <- vapply(lethal_mapped_df[, "Gene_symbol"], function(x) {
  col_names <- names(lethality_df)[4:2]
  is_in_list <- vapply(col_names, function(y) x %in% lethality_df[, y], logical(1))
  col_names[[which(is_in_list)[[1]]]]
}, "")

prions_category_fac <- factor(prions_category_vec, names(lethality_df)[c(4, 2, 3)])

prions_lethality_df[["Category"]] <- prions_category_fac




# Save data ---------------------------------------------------------------

save(list = c("membrane_het_df", "prions_lethality_df"),
     file = file.path(general_RData_directory, "06) Read in gene lists.RData")
     )






