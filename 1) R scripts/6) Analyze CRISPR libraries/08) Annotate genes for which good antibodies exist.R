### 21st April 2021 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Analysis", "Alexa", "Antibodies")
gene_lists_directory     <- file.path(CRISPR_root_directory, "2) Input data", "Gene lists")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "22) Compile data on essential genes - essential_df.RData"))
load(file.path(general_RData_directory, "23) Compile data on protein localization.RData"))
load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRa_by_gene_df <- full_4sg_by_gene_df
load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_by_gene_df <- full_4sg_by_gene_df
rm(full_4sg_by_gene_df)




# Read in data ------------------------------------------------------------

original_4i_df <- data.frame(read_excel(file.path(gene_lists_directory, "4i_Antibodies_PelkmansLab.xlsx")),
                             stringsAsFactors = FALSE, check.names = FALSE
                             )




# Define functions --------------------------------------------------------

ReplaceNAs <- function(input_df) {
  NA_columns <- c("Achilles_common", "CRISPR_common", "Hart_3_or_more_lines",
                  "Hart_HeLa", "Blomen_HAP1_KBM7_intersect"
                  )
  no_entrez <- is.na(input_df[["Entrez_ID"]])
  for (column_name in NA_columns) {
    input_df[, column_name] <- ifelse(is.na(input_df[, column_name]),
                                      ifelse(no_entrez, "", "N/A"),
                                      as.character(input_df[, column_name])
                                      )
  }

  for (column_name in setdiff(names(input_df), NA_columns)) {
    input_df[, column_name] <- ifelse(is.na(input_df[, column_name]),
                                      "",
                                      as.character(input_df[, column_name])
                                      )
  }
  return(input_df)
}




# Translate 4i antibody names ---------------------------------------------

are_all_NA <- apply(original_4i_df[, 2:ncol(original_4i_df)], 1, function(x) all(is.na(x)))
original_4i_df <- original_4i_df[!(are_all_NA), ]

initial_entrezs_4i <- SymbolsToEntrezIDs(toupper(original_4i_df[["Name"]]))

intermediate_4i_df <- data.frame(
  original_4i_df[1],
  "Entrez_ID"        = initial_entrezs_4i,
  "Gene_symbol"      = EntrezIDsToSymbols(initial_entrezs_4i),
  "Gene_description" = EntrezIDsToGeneNames(initial_entrezs_4i),
  original_4i_df[, 2:ncol(original_4i_df)],
  stringsAsFactors = FALSE
)




# Export intermediate 4i antibody gene IDs --------------------------------

write.table(intermediate_4i_df,
            file = file.path(file_output_directory, "Antibodies_4i_intermediate_gene_IDs.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )




# Process 4i table --------------------------------------------------------

antibodies_4i_df <- data.frame(read_excel(file.path(file_output_directory, "4i_antibodies_manually_annotated.xlsx")),
                               stringsAsFactors = FALSE, check.names = FALSE
                               )

are_all_NA <- apply(antibodies_4i_df[, 2:ncol(antibodies_4i_df)], 1, function(x) all(is.na(x)))
antibodies_4i_df <- antibodies_4i_df[!(are_all_NA), ]


rename_vec <- c(
  "Entrez ID"    = "Entrez_ID",
  "Gene symbol"  = "Gene_symbol",
  "Category"     = "Category",
  "Antibody #"   = "Antibody_number_4i",
  "Name"         = "Antibody_name_4i",
  "Multimer?"    = "Is_multimer",
  ">1 gene?"     = "Encoded_by_multiple_genes",
  "PTM/isoform?" = "Specific_PTM_or_isoform",
  "Comment"      = "Comment_Lukas",
  "Species"      = "Host_species",
  "Clonality"    = "Clonality",
  "Manufacturer" = "Manufacturer",
  "Catalogue #"  = "Catalog_number",
  "Clone ID"     = "Clone_ID",
  "Remarks"      = "Remarks",
  "Conc. (1/x)"  = "Dilution_factor",
  "Cycle Used"   = "Cycle_used"
)

new_4i_df <- antibodies_4i_df
new_4i_df[["Category"]] <- NA

for (column_name in names(rename_vec)) {
  names(new_4i_df)[names(new_4i_df) == column_name] <- rename_vec[[column_name]]
}

new_4i_df <- new_4i_df[, rename_vec]

entrezs_4i <- SymbolsToEntrezIDs(new_4i_df[["Gene_symbol"]])
symbols_4i <- EntrezIDsToSymbols(entrezs_4i)

stopifnot(!(any(is.na(entrezs_4i) & !(is.na(new_4i_df[["Gene_symbol"]])))))
stopifnot(!(any(!(is.na(entrezs_4i)) & is.na(symbols_4i))))
stopifnot(identical(symbols_4i, MapToEntrezs(entrez_IDs_vec = entrezs_4i)[["Gene_symbol"]]))

new_4i_df[["Entrez_ID"]] <- entrezs_4i
new_4i_df[["Gene_symbol"]] <- symbols_4i




# Define gene lists -------------------------------------------------------

Alexa_list <- list(
  "High essentiality" = c(
    "DDX6", "SON", "SRRM1", "SRRM2", "Sec13"
  ),
  "Low essentiality" = c(
    "DCP1A", "G3BP1", "NONO", "SP100"
  ),
  "FISH probes available" = c(
    "EEA1", "TFRC", "LAMP1", "VPS35", "ACTB"
  ),
  "Antibodies to validate" = c(
    "PPP1CA", "PPP1CB", "PPP1CC", "PPP2CA", "DYRK1A", "DYRK2", "DYRK3"
  )
)

Alexa_df <- data.frame(
  "Category"    = rep(names(Alexa_list), lengths(Alexa_list)),
  "Gene_symbol" = toupper(unlist(Alexa_list, use.names = FALSE)),
  stringsAsFactors  = FALSE
)

Alexa_df[["Entrez_ID"]] <- SymbolsToEntrezIDs(Alexa_df[["Gene_symbol"]])

stopifnot(!(anyNA(Alexa_df[["Entrez_ID"]])))
stopifnot(identical(Alexa_df[["Gene_symbol"]], EntrezIDsToSymbols(Alexa_df[["Entrez_ID"]])))




# Combine the "Alexa" and "4i" gene lists ---------------------------------

Alexa_matches_list <- lapply(Alexa_df[["Entrez_ID"]], function(x) {
  are_this_gene <- new_4i_df[["Entrez_ID"]] %in% x
  if (!(any(are_this_gene))) {
    return(NA)
  } else {
    return(which(are_this_gene))
  }
})

Alexa_matches <- unlist(Alexa_matches_list)

Alexa_4i_df <- new_4i_df[Alexa_matches, ]

original_indices <- rep(seq_len(nrow(Alexa_df)), lengths(Alexa_matches_list))

for (column_name in c("Category", "Entrez_ID", "Gene_symbol")) {
  Alexa_4i_df[[column_name]] <- Alexa_df[[column_name]][original_indices]
}

non_Alexa_indices <- setdiff(seq_len(nrow(new_4i_df)), Alexa_matches)

non_Alexa_4i_df <- new_4i_df[non_Alexa_indices, ]

non_Alexa_4i_df[["Category"]] <- "Other 4i"

combined_4i_df <- rbind.data.frame(
  Alexa_4i_df,
  non_Alexa_4i_df,
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)





# Create a combined data frame of protein-coding genes --------------------

all_entrezs <- unique(c(CRISPRa_by_gene_df[["Entrez_ID"]],
                        CRISPRko_by_gene_df[["Entrez_ID"]]
                      ))
all_entrezs <- sort(as.integer(all_entrezs))

sublibrary_entrezs <- unlist(sublibraries_all_entrezs_list)
sublibrary_vec <- rep(names(sublibraries_all_entrezs_list),
                      lengths(sublibraries_all_entrezs_list)
                      )
sublibrary_matches <- match(all_entrezs, sublibrary_entrezs)

all_4sg_df <- data.frame(
  "Entrez_ID"    = all_entrezs,
  "Gene_symbol"  = MapToEntrezs(entrez_IDs_vec = as.character(all_entrezs))[["Gene_symbol"]],
  "CRISPRa_4sg"  = ifelse(all_entrezs %in% CRISPRa_by_gene_df[["Entrez_ID"]],
                          "Yes", "No"
                          ),
  "CRISPRko_4sg" = ifelse(all_entrezs %in% CRISPRko_by_gene_df[["Entrez_ID"]],
                          "Yes", "No"
                          ),
  "Sublibrary"   = sublibrary_vec[sublibrary_matches],
  stringsAsFactors = FALSE
)

stopifnot(!(anyNA(all_4sg_df[["Sublibrary"]])))





# Add additional data -----------------------------------------------------

surfaceome_matches <- match(all_4sg_df[["Entrez_ID"]], surfaceome_df[["Entrez_ID"]])
CD_matches         <- match(all_4sg_df[["Entrez_ID"]], names(entrez_to_CD_map))
HPA_matches        <- match(all_4sg_df[["Entrez_ID"]], HPA_df[["Entrez_ID"]])
essential_matches  <- match(all_4sg_df[["Entrez_ID"]], essential_df[["Entrez_ID"]])

HPA_columns <- setdiff(c("Subcellular_location", names(HPA_df)),
                       c("Ensembl_gene_ID", "Original_symbol", "Antibody_IDs",
                         names(all_4sg_df)
                       ))

essential_columns <- setdiff(names(essential_df),
                             c("Achilles_num_essential", "Achilles_num_cell_lines",
                               "Achilles_mean_probability",
                               "DEMETER2_mean_probability", "CRISPR_mean_probability",
                               "BlomenHart_intersect", "BlomenHart_intersect_DepMap",
                               "Three_categories", "Four_categories",
                               names(all_4sg_df)
                               )
                             )

all_genes_df <- data.frame(all_4sg_df,
                           "CD" = entrez_to_CD_map[CD_matches],
                           surfaceome_df[surfaceome_matches, !(names(surfaceome_df) %in% names(all_4sg_df))],
                           HPA_df[HPA_matches, HPA_columns],
                           essential_df[essential_matches, essential_columns],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )







# Select surface proteins that are not expressed for CRISPRa FACS ---------

are_CRISPRa <- all_genes_df[["CRISPRa_4sg"]] %in% "Yes"
are_surfaceome <- all_genes_df[, "SURFY_surfaceome_2018"] %in% "Yes"

use_surface_columns <- c("CSC_HeLa_2018", "CSC_HeLa_2015", "CSC_HEK_2015", "CSC_LN229_2015")
are_not_expressed <- rowSums((as.matrix(all_genes_df[, use_surface_columns]) == "Yes"), na.rm = TRUE) == 0

are_plasma_membrane <- all_genes_df[["Subcellular_location"]] %in% "Plasma membrane"

are_eligible <- are_CRISPRa & are_surfaceome & are_not_expressed & are_plasma_membrane

FACS_df <- all_genes_df[are_eligible, ]
row.names(FACS_df) <- NULL

set.seed(1)
FACS_df[["Randomized_rank"]] <- sample(seq_len(nrow(FACS_df)))
shuffled_FACS_df <- FACS_df[order(FACS_df[["Randomized_rank"]]), ]
row.names(shuffled_FACS_df) <- NULL

export_FACS_df <- ReplaceNAs(shuffled_FACS_df)

write.table(export_FACS_df,
            file = file.path(file_output_directory, "CRISPRa_for_FACS.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )




# Integrate the 4i data ---------------------------------------------------

matches_4i_4sg <- match(combined_4i_df[["Entrez_ID"]], all_genes_df[["Entrez_ID"]])

genes_4i_4sg_df <- all_genes_df[matches_4i_4sg, ]

are_other_4sg <- !(all_genes_df[["Entrez_ID"]] %in% genes_4i_4sg_df[["Entrez_ID"]])

expanded_4i_indices <- c(seq_len(nrow(combined_4i_df)),
                         rep(NA, sum(are_other_4sg))
                         )

expanded_4i_df <- combined_4i_df[expanded_4i_indices, ]
for (column_name in intersect(names(all_genes_df), names(combined_4i_df))) {
  expanded_4i_df[[column_name]] <- c(combined_4i_df[[column_name]],
                                     all_genes_df[[column_name]][are_other_4sg]
                                     )
}

expanded_4i_matches <- match(expanded_4i_df[["Entrez_ID"]],
                             all_genes_df[["Entrez_ID"]]
                             )

all_4i_4sg_df <- data.frame(
  expanded_4i_df[, 1:5],
  all_genes_df[expanded_4i_matches, !(names(all_genes_df) %in% names(expanded_4i_df))],
  expanded_4i_df[, 6:ncol(expanded_4i_df)],
  stringsAsFactors = FALSE
)




# Replace NAs with empty strings ------------------------------------------

export_4i_4sg_df <- ReplaceNAs(all_4i_4sg_df)
export_all_genes_df <- ReplaceNAs(all_genes_df)





# Export the integrated data frame ----------------------------------------

write.table(export_4i_4sg_df,
            file = file.path(file_output_directory, "All_4i_4sg.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )


write.table(export_all_genes_df,
            file = file.path(file_output_directory, "All_4sg.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )

