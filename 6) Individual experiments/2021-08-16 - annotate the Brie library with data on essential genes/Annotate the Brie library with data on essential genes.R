### 7th March 2021 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))
source(file.path(general_functions_directory, "32) Compiling data on essential genes.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
CRISPR_input_directory    <- file.path(CRISPR_root_directory, "2) Input data")
essential_genes_directory <- file.path(CRISPR_input_directory, "Human genome", "Essential genes")
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRko_Brie_path        <- file.path(CRISPR_input_directory,
                                       "Mouse CRISPR libraries",
                                       "CRISPRko",
                                       "Brie",
                                       "2016 - Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9 - Table S22 - Brie.xlsx"
                                       )

file_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-08-16 - annotate the Brie library with data on essential genes")
file_input_directory      <- file.path(file_directory, "1) Input")
file_output_directory     <- file.path(file_directory, "3) Output")
R_objects_directory       <- file.path(file_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "22) Compile data on essential genes - dependencies data.RData"))




# Read in data ------------------------------------------------------------

Brie_df <- data.frame(read_excel(CRISPRko_Brie_path),
                      stringsAsFactors = FALSE, check.names = FALSE
                      )


homologs_df <- read.delim(file.path(file_input_directory, "HOM_MouseHumanSequence.rpt.txt"),
                          stringsAsFactors = FALSE, check.names = FALSE
                          )





# Map from mouse Entrez gene IDs to human gene IDs ------------------------

Brie_df[, "Target Gene ID"] <- as.integer(Brie_df[, "Target Gene ID"])
mouse_entrezs <- sort(unique(Brie_df[, "Target Gene ID"]))

are_mouse <- homologs_df[, "Common Organism Name"] == "mouse, laboratory"
mouse_to_key_list <- lapply(mouse_entrezs, function(x) {
  are_this_entrez <- are_mouse & (homologs_df[, "EntrezGene ID"] == x)
  results_vec <- homologs_df[are_this_entrez, "DB Class Key"]
  return(results_vec)
})

homologs_df[homologs_df[["DB Class Key"]] == 38833823, ]

mouse_to_key_list <- lapply(mouse_to_key_list, unique)

are_present <- lengths(mouse_to_key_list) != 0

mappings_mat <- cbind(
  "Mouse_Entrez_gene_ID" = mouse_entrezs[are_present],
  "DB_key" = unlist(mouse_to_key_list[are_present])
)

key_to_human_list <- lapply(mappings_mat[, "DB_key"], function(x) {
  are_this_key <- !(are_mouse) & (homologs_df[, "DB Class Key"] == x)
  stopifnot(!(anyNA(are_this_key)))
  results_vec <- homologs_df[["EntrezGene ID"]][are_this_key]
  return(results_vec)
})

indices_vec <- rep(seq_len(nrow(mappings_mat)), lengths(key_to_human_list))
mappings_mat <- cbind(mappings_mat[indices_vec, ],
                      "Human_Entrez_gene_ID" = unlist(key_to_human_list)
                      )




# Process genes from the Brie library -------------------------------------

all_entrezs <- sort(unique(mappings_mat[, "Human_Entrez_gene_ID"]))

essential_df <- GetGeneEssentiality(all_entrezs, essential_datasets_list)

stopifnot(identical(all_entrezs, essential_df[["Entrez_ID"]]))






# Annotate mouse genes with human essentiality data -----------------------

stopifnot(!(any(duplicated(essential_df[, "Entrez_ID"]))))

symbols_list <- lapply(mouse_entrezs, function(x) {
  are_this_gene <- Brie_df[, "Target Gene ID"] == x
  unique(Brie_df[["Target Gene Symbol"]][are_this_gene])
})
stopifnot(all(lengths(symbols_list) == 1))

df_list <- lapply(mouse_entrezs, function(x) {
  are_this_gene <- mappings_mat[, "Mouse_Entrez_gene_ID"] == x
  if (any(are_this_gene)) {
    human_entrezs <- unique(mappings_mat[are_this_gene, "Human_Entrez_gene_ID"])
    matches_vec <- match(human_entrezs, essential_df[, "Entrez_ID"])
  } else {
    matches_vec <- NA_integer_
  }
  results_df <- essential_df[matches_vec, ]
  return(results_df)
})
reordered_df <- do.call(rbind.data.frame,
                        c(df_list,
                          stringsAsFactors = FALSE,
                          make.row.names = FALSE
                          )
                        )

num_entrezs <- vapply(df_list, nrow, integer(1))
mouse_essential_df <- data.frame(
  "Mouse_Entrez_ID" = rep(mouse_entrezs, num_entrezs),
  "Mouse_symbol"    = rep(unlist(symbols_list), num_entrezs),
  "Mapping"         = NA,
  "Human_Entrez_ID" = reordered_df[, "Entrez_ID"],
  "Human_symbol"    = reordered_df[, "Gene_symbol"],
  reordered_df[, 3:ncol(reordered_df)],
  stringsAsFactors = FALSE
)

mapping_types <- vapply(mouse_essential_df[, "Mouse_Entrez_ID"], function(x) {
  are_this_gene <- mouse_essential_df[, "Mouse_Entrez_ID"] == x
  human_entrezs <- mouse_essential_df[["Human_Entrez_ID"]][are_this_gene]
  if (all(is.na(human_entrezs))) {
    return("No mapping")
  } else {
    if (length(human_entrezs) > 1) {
      return("One-to-many")
    } else {
      are_these <- mouse_essential_df[["Human_Entrez_ID"]] %in% human_entrezs
      mouse_human_mouse <- unique(mouse_essential_df[["Mouse_Entrez_ID"]][are_these])
      if (length(mouse_human_mouse) > 1) {
        return("Many-to-one")
      } else {
        return("One-to-one")
      }
    }

  }
}, "")

mouse_essential_df[["Mapping"]] <- mapping_types




# Get an overview over the mappings ---------------------------------------

mappings_vec <- vapply(mouse_entrezs, function(x) {
  are_this_gene <- mouse_essential_df[, "Mouse_Entrez_ID"] == x
  mappings_vec <- unique(mouse_essential_df[["Mapping"]][are_this_gene])
  return(mappings_vec)
}, "")


data_columns <- names(mouse_essential_df)[6:ncol(mouse_essential_df)]
are_NA_mat <- do.call(cbind, lapply(mouse_essential_df[, data_columns], is.na))
are_all_NA <- apply(are_NA_mat, 1, all)

have_no_data <- vapply(mouse_entrezs, function(x) {
  are_this_gene <- mouse_essential_df[, "Mouse_Entrez_ID"] == x
  results_vec <- all(are_all_NA[are_this_gene])
  return(results_vec)
}, logical(1))

new_mappings_vec <- ifelse(mappings_vec == "None",
                           "None",
                           paste0(mappings_vec,
                                  ifelse(have_no_data, ", but there is no data", "")
                                  )

                           )



# Export data -------------------------------------------------------------

essential_export_df <- mouse_essential_df
NA_empty_columns <- c("CRISPR_mean_probability",   "CRISPR_num_essential",   "CRISPR_num_cell_lines",
                      "CRISPR_mean_effect",
                      "Achilles_mean_probability", "Achilles_num_essential", "Achilles_num_cell_lines",
                      "DEMETER2_mean_probability", "DEMETER2_num_essential", "DEMETER2_num_cell_lines",
                      "Human_Entrez_ID", "Human_symbol"
                      )
for (column_name in NA_empty_columns) {
  essential_export_df[, column_name] <- ifelse(is.na(essential_export_df[, column_name]),
                                               "",
                                               essential_export_df[, column_name]
                                               )
}

essential_export_df <- essential_export_df[, !(names(essential_export_df) %in% c("Achilles_mean_probability", "Achilles_num_essential", "Achilles_num_cell_lines"))]

write.table(essential_export_df,
            file = file.path(file_output_directory, "Mouse_to_human_essential_genes.tsv"),
            sep = "\t", na = "N/A", quote = FALSE, row.names = FALSE
            )





# Save data ---------------------------------------------------------------

save(list = "mouse_essential_df",
     file = file.path(R_objects_directory, "Annotate the Brie library with data on essential genes.RData")
     )




