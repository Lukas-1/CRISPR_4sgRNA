### 7th March 2021 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
CRISPR_input_directory    <- file.path(CRISPR_root_directory, "2) Input data")
essential_genes_directory <- file.path(CRISPR_input_directory, "Human genome", "Essential genes")
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
file_output_directory     <- file.path(CRISPR_root_directory, "5) Output", "Analysis", "Essential genes")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))




# Define functions --------------------------------------------------------

ProcessAchillesGenesDf <- function(genes_df) {
  # genes_df is a single column in the format A1BG (1) {1 being the Entrez ID}

  gene_splits <- strsplit(genes_df[, "gene"], " ", fixed = TRUE)
  entrezs_vec <- sapply(gene_splits, "[[", 2)
  stopifnot(all(substr(entrezs_vec, 1, 1) == "("))
  stopifnot(all(substr(entrezs_vec, nchar(entrezs_vec), nchar(entrezs_vec)) == ")"))

  stripped_entrezs <- as.integer(substr(entrezs_vec, 2, nchar(entrezs_vec) - 1))
  stripped_entrezs <- as.integer(stripped_entrezs)
  original_symbols <- sapply(gene_splits, "[[", 1)
  results_df <- data.frame(
    "Entrez_ID"       = stripped_entrezs,
    "Original_symbol" = original_symbols,
    stringsAsFactors = FALSE
  )
  CheckEntrezs(stripped_entrezs, results_df)
  return(results_df)
}



CheckEntrezs <- function(entrezs_vec, results_df = NULL) {
  are_present <- as.character(entrezs_vec) %in% names(entrez_to_symbol_list)
  if (!(all(are_present))) {
    message(paste0("The following Entrez IDs were not found in ",
                   "org.Hs.egSYMBOL: ",
                   paste0(entrezs_vec[!(are_present)], collapse = ", ")
                   ))
    if (!(is.null(results_df))) {
      print(results_df[!(are_present), ])
    }
    message("")
  }
  return(invisible(NULL))
}


ProcessAchillesDataDf <- function(data_df, cell_lines_df) {

  data_mat <- t(as.matrix(data_df[, 2:ncol(data_df)]))
  cell_line_matches <- match(data_df[, "DepMap_ID"],
                             cell_lines_df[["DepMap_ID"]]
                             )
  colnames(data_mat) <- cell_lines_df[["stripped_cell_line_name"]][cell_line_matches]
  genes_df <- data.frame("gene" = row.names(data_mat),
                         stringsAsFactors = FALSE
                         )
  genes_df <- ProcessAchillesGenesDf(genes_df)

  results_df <- data.frame(
    genes_df,
    "Mean" = rowMeans(data_mat, na.rm = TRUE),
    data_mat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



ProcessDEMETERDataDf <- function(data_df) {
  stopifnot(all(c("depmap_samples_df", "DEMETER2_samples_df") %in% ls(envir = globalenv())))

  new_order <- order(as.integer(sub("ACH-", "", data_df[[1]])))
  data_df <- data_df[new_order, ]
  row.names(data_df) <- NULL

  data_mat <- t(as.matrix(data_df[, 2:ncol(data_df)]))
  cell_line_matches <- match(data_df[[1]], depmap_samples_df[["DepMap_ID"]])

  DEMETER_matches <- match(data_df[[1]], DEMETER2_samples_df[[1]])
  DEMETER_names <- DEMETER2_samples_df[["Marcotte_name"]][DEMETER_matches]
  have_no_match <- is.na(cell_line_matches)
  replaced_names <- DEMETER_names[have_no_match]
  stopifnot(!(anyNA(replaced_names)))
  replaced_names <- toupper(replaced_names)
  stopifnot(!(any(replaced_names %in% depmap_samples_df[["stripped_cell_line_name"]])))

  colnames(data_mat) <- ifelse(have_no_match,
                               replaced_names,
                               depmap_samples_df[["stripped_cell_line_name"]][cell_line_matches]
                               )

  gene_splits <- strsplit(rownames(data_mat), " ", fixed = TRUE)
  gene_symbols <- sapply(gene_splits, "[[", 1)
  entrez_IDs <- sapply(gene_splits, "[[", 2)
  entrez_IDs <- substr(entrez_IDs, 2, nchar(entrez_IDs) - 1)
  symbol_splits <- strsplit(gene_symbols, "&", fixed = TRUE)
  entrez_splits <- strsplit(entrez_IDs, "&", fixed = TRUE)
  stopifnot(identical(lengths(symbol_splits), lengths(entrez_splits)))
  stopifnot(!(any(unlist(entrez_IDs[lengths(entrez_IDs) == 1]) %in% unlist(entrez_IDs[lengths(entrez_IDs) > 1]))))

  symbols_vec <- unlist(symbol_splits, use.names = FALSE)
  entrezs_vec <- as.integer(unlist(entrez_splits, use.names = FALSE))

  CheckEntrezs(entrezs_vec)
  expanded_indices <- rep(seq_along(entrez_IDs), lengths(entrez_splits))

  expanded_df <- data.frame(
    "Entrez_ID" = entrezs_vec,
    "Gene_symbol" = symbols_vec,
    "Num_targeted_genes" = lengths(entrez_splits)[expanded_indices],
    "Mean" = rowMeans(data_mat[expanded_indices, ], na.rm = TRUE),
    data_mat[expanded_indices, ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  expanded_df <- expanded_df[order(expanded_df[["Entrez_ID"]]), ]
  row.names(expanded_df) <- NULL

  return(expanded_df)
}




GetGeneEssentiality <- function(entrezs_vec, datasets_list) {

  required_objects <- "achilles_depend_df"
  stopifnot(all(required_objects %in% ls(envir = globalenv())))

  stopifnot(identical("integer", typeof(entrezs_vec)))

  AreEssential <- function(essential_genes, all_genes) {
    stopifnot(all(essential_genes %in% all_genes))
    ifelse(!(entrezs_vec %in% all_genes),
           NA,
           entrezs_vec %in% essential_genes
           )
  }

  symbols_vec <- MapToEntrezs(entrez_IDs_vec = as.character(all_entrezs))[["Gene_symbol"]]

  Achilles_matches <- match(entrezs_vec, achilles_depend_df[["Entrez_ID"]])
  Achilles_depend_mat <- as.matrix(achilles_depend_df[, 4:ncol(achilles_depend_df)], na.rm = TRUE)
  Achilles_depend_mat <- Achilles_depend_mat[Achilles_matches, ]

  DEMETER_matches <-  match(entrezs_vec, DEMETER2_combined_depend_df[["Entrez_ID"]])
  DEMETER_depend_mat <- as.matrix(DEMETER2_combined_depend_df[, 4:ncol(DEMETER2_combined_depend_df)], na.rm = TRUE)
  DEMETER_depend_mat <- DEMETER_depend_mat[DEMETER_matches, ]

  genes_list <- list(
    "Entrez_ID"                 = entrezs_vec,
    "Gene_symbol"               = symbols_vec,
    "Achilles_mean_probability" = achilles_depend_df[["Mean"]][Achilles_matches],
    "Achilles_num_essential"    = ifelse(is.na(Achilles_matches),
                                         NA,
                                         rowSums(Achilles_depend_mat > 0.5, na.rm = TRUE)
                                         ),
    "Achilles_num_cell_lines"   = ifelse(is.na(Achilles_matches),
                                         NA,
                                         rowSums(!(is.na(Achilles_depend_mat)))
                                         ),
    "DEMETER2_mean_probability" = DEMETER2_combined_depend_df[["Mean"]][DEMETER_matches],
    "DEMETER2_num_essential"    = ifelse(is.na(DEMETER_matches),
                                         NA,
                                         rowSums(DEMETER_depend_mat > 0.5, na.rm = TRUE)
                                         ),
    "DEMETER2_num_cell_lines"   = ifelse(is.na(DEMETER_matches),
                                         NA,
                                         rowSums(!(is.na(DEMETER_depend_mat)))
                                         )
  )

  genes_list <- c(
    genes_list,
    lapply(datasets_list, function(x) AreEssential(x[["essential"]], x[["all"]]))
  )

  genes_df <- do.call(data.frame, c(genes_list, stringsAsFactors = FALSE))

  for (i in seq_along(genes_df)) {
    if (is.logical(genes_df[[i]])) {
      genes_df[[i]] <- ifelse(genes_df[[i]], "Essential", "Non-essential")
    }
  }

  return(genes_df)
}





# Read in lists of essential genes from original publications -------------

blomen_file_prefix <- "2015 - Gene essentiality and synthetic lethality in haploid human cells - Table S"
blomen_df <- data.frame(read_excel(file.path(essential_genes_directory,
                                             "Publications", paste0(blomen_file_prefix, "3.xlsx")
                                             ),
                                   skip = 1, na = "NA", .name_repair = "minimal"
                                   ),
                        stringsAsFactors = FALSE, check.names = FALSE
                        )
names(blomen_df)[3:8] <- paste0("KBM7_", names(blomen_df)[3:8])
names(blomen_df)[11:16] <- paste0("HAP1_", names(blomen_df)[11:16])
blomen_df <- blomen_df[, !(duplicated(as.list(blomen_df)))]




blomen_KBM7_all_genes_df <- data.frame(read_excel(file.path(essential_genes_directory,
                                                            "Publications", paste0(blomen_file_prefix, "1.xlsx")
                                                            ),
                                                  skip = 1, na = "NA"
                                                  ),
                                       stringsAsFactors = FALSE, check.names = FALSE
                                       )

blomen_HAP1_all_genes_df <- data.frame(read_excel(file.path(essential_genes_directory,
                                                            "Publications", paste0(blomen_file_prefix, "2.xlsx")
                                                            ),
                                                  skip = 1, na = "NA"
                                                  ),
                                       stringsAsFactors = FALSE, check.names = FALSE
                                       )



hart_file_name <- "2015 - High-Resolution CRISPR Screens Reveal Fitness Genes - Table S2.xlsx"
hart_df <- data.frame(read_excel(file.path(essential_genes_directory,
                                           "Publications", hart_file_name
                                           )
                                 ),
                      stringsAsFactors = FALSE
                      )





# Read in lists of essential genes from DepMap / Project Achilles ---------
## Download link: https://depmap.org/portal/download/all/

ReadDepMapFile <- function(file_name, sub_folder = "2021_Q1") {
  read.csv(file.path(essential_genes_directory,
                     "DepMap", sub_folder,
                     paste0(file_name, ".csv")
                     ),
           stringsAsFactors = FALSE, check.names = FALSE
           )
}

depmap_hart_blomen_df        <- ReadDepMapFile("common_essentials_21Q1")

achilles_common_df           <- ReadDepMapFile("Achilles_common_essentials_21Q1")
achilles_original_effects_df <- ReadDepMapFile("Achilles_gene_effect_21Q1")
achilles_original_depend_df  <- ReadDepMapFile("Achilles_gene_dependency_21Q1")

CRISPR_common_df             <- ReadDepMapFile("CRISPR_common_essentials_21Q1")
CRISPR_original_effects_df   <- ReadDepMapFile("CRISPR_gene_effect_21Q1")
CRISPR_original_depend_df    <- ReadDepMapFile("CRISPR_gene_dependency_21Q1")

depmap_samples_df            <- ReadDepMapFile("sample_info_21Q1")

DEMETER2_combined_original_depend_df <- ReadDepMapFile("gene_dependency", sub_folder = "DEMETER 2 Combined RNAi")
DEMETER2_combined_original_effects_df <- ReadDepMapFile("gene_effect", sub_folder = "DEMETER 2 Combined RNAi")
DEMETER2_samples_df <- ReadDepMapFile("sample_info", sub_folder = "DEMETER 2 Combined RNAi")





# Create Entrez IDs for data of Blomen et al. (2015) ----------------------


### Essential genes ###

blomen_essential_genes_df <- data.frame("Ensembl_gene_ID" = blomen_df[, "ENSEMBL_ID"],
                                        "Gene_symbol"     = blomen_df[, "GENE_SYMBOL"],
                                        stringsAsFactors = FALSE
                                        )
blomen_essential_mappings_df <- MapEnsemblIDs(blomen_essential_genes_df)

blomen_df[["Entrez_ID"]] <- blomen_essential_mappings_df[["Consensus_entrez"]]



### All genes ###

stopifnot(identical(sort(blomen_KBM7_all_genes_df[["ENSEMBL_ID"]]),
                    sort(blomen_HAP1_all_genes_df[["ENSEMBL_ID"]])
                    ))

blomen_all_genes_df <- data.frame("Ensembl_gene_ID" = blomen_KBM7_all_genes_df[, "ENSEMBL_ID"],
                                  "Gene_symbol"     = blomen_KBM7_all_genes_df[, "GENE_SYMBOL"],
                                  stringsAsFactors = FALSE
                                  )
blomen_all_mappings_df <- MapEnsemblIDs(blomen_all_genes_df)
blomen_all_entrezs <- blomen_all_mappings_df[["Consensus_entrez"]]





# Create Entrez IDs for data of Hart et al. (2015) ------------------------

### Essential genes ###

hart_df[, "numTKOHits"] <- as.integer(hart_df[, "numTKOHits"])

hart_symbols_list <- lapply(hart_df[, "Gene"], function(x) {
  if (x %in% names(symbol_to_entrez_list)) {
    results_vec <- symbol_to_entrez_list[[x]]
  } else if (x %in% names(alias_to_entrez_list)) {
    results_vec <- alias_to_entrez_list[[x]]
  } else {
    results_vec <- NA_character_
  }
  return(results_vec)
})

hart_df[["Entrez_IDs"]] <- vapply(hart_symbols_list, function(x) {
  if (is.null(x)) {
    return(NA_character_)
  } else {
    paste0(x, collapse = ", ")
  }
}, "")

hart_df[["Num_Entrez_ID_mappings"]] <- ifelse(is.na(hart_symbols_list),
                                              0L,
                                              lengths(hart_symbols_list)
                                              )

expanded_hart_indices <- rep(seq_len(nrow(hart_df)), lengths(hart_symbols_list))
expanded_hart_df <- hart_df[expanded_hart_indices, ]
expanded_hart_df[["Entrez_ID"]] <- as.integer(unlist(hart_symbols_list, use.names = FALSE))





# Intersect the Hart and Blomen datasets (for comparison...) --------------

blomen_are_selected <- (blomen_df[, "KBM7_selected"] == "YES") &
                       (blomen_df[, "HAP1_selected"] == "YES")

hart_2_or_more <- expanded_hart_df[, "numTKOHits"] >= 2

blomen_hart_2_intersect <- intersect(blomen_df[["Entrez_ID"]][blomen_are_selected],
                                     expanded_hart_df[["Entrez_ID"]][hart_2_or_more]
                                     )
blomen_hart_2_intersect <- blomen_hart_2_intersect[!(is.na(blomen_hart_2_intersect))]





# Process data frames from the Achilles project ---------------------------


### Split Entrez IDs and gene symbols for lists of common genes ###

depmap_hart_blomen_df <- ProcessAchillesGenesDf(depmap_hart_blomen_df)
achilles_common_df    <- ProcessAchillesGenesDf(achilles_common_df)
CRISPR_common_df      <- ProcessAchillesGenesDf(CRISPR_common_df)


### Achilles gene effects ###

achilles_effects_df <- ProcessAchillesDataDf(achilles_original_effects_df,
                                             depmap_samples_df
                                             )
achilles_depend_df  <- ProcessAchillesDataDf(achilles_original_depend_df,
                                             depmap_samples_df
                                             )
CRISPR_effects_df   <- ProcessAchillesDataDf(CRISPR_original_effects_df,
                                             depmap_samples_df
                                             )
CRISPR_depend_df    <- ProcessAchillesDataDf(CRISPR_original_depend_df,
                                             depmap_samples_df
                                             )

stopifnot(identical(names(achilles_depend_df), names(achilles_effects_df)))
stopifnot(identical(achilles_effects_df[, "Entrez_ID"], achilles_depend_df[, "Entrez_ID"]))

stopifnot(identical(names(CRISPR_effects_df), names(CRISPR_depend_df)))
stopifnot(identical(CRISPR_effects_df[, "Entrez_ID"], CRISPR_depend_df[, "Entrez_ID"]))



# Process the combined RNAi dataset ---------------------------------------

DEMETER2_combined_depend_df  <- ProcessDEMETERDataDf(DEMETER2_combined_original_depend_df)
DEMETER2_combined_effects_df <- ProcessDEMETERDataDf(DEMETER2_combined_original_effects_df)

stopifnot(identical(names(DEMETER2_combined_depend_df), names(DEMETER2_combined_effects_df)))
stopifnot(identical(DEMETER2_combined_depend_df[, "Entrez_ID"], DEMETER2_combined_effects_df[, "Entrez_ID"]))





# Define "essential genes" and "all genes" for each dataset ---------------

hart_are_selected <- expanded_hart_df[, "numTKOHits"] >= 3

essential_blomen_hart_intersect <- intersect(blomen_df[["Entrez_ID"]][blomen_are_selected],
                                             expanded_hart_df[["Entrez_ID"]][hart_are_selected]
                                             )
essential_blomen_hart_intersect <- essential_blomen_hart_intersect[!(is.na(essential_blomen_hart_intersect))]

all_blomen_hart_intersect <- intersect(expanded_hart_df[["Entrez_ID"]],
                                       blomen_all_entrezs
                                       )

essential_datasets_list <- list(

    "Achilles_common"             = list("essential" = achilles_common_df[["Entrez_ID"]],
                                         "all"       = union(achilles_common_df[["Entrez_ID"]],
                                                             achilles_depend_df[["Entrez_ID"]]
                                                             )
                                         ),

    "CRISPR_common"               = list("essential" = CRISPR_common_df[["Entrez_ID"]],
                                         "all"       = union(achilles_common_df[["Entrez_ID"]],
                                                             achilles_depend_df[["Entrez_ID"]]
                                                             )
                                         ),

    "BlomenHart_intersect"        = list("essential" = essential_blomen_hart_intersect,
                                         "all"       = all_blomen_hart_intersect
                                         ),

    "BlomenHart_intersect_DepMap" = list("essential" = depmap_hart_blomen_df[["Entrez_ID"]],
                                         "all"       = union(depmap_hart_blomen_df[["Entrez_ID"]],
                                                             intersect(all_blomen_hart_intersect,
                                                                       achilles_depend_df[["Entrez_ID"]]
                                                                       )
                                                             )
                                         ),

    "Hart_3_or_more_lines"        = list("essential" = expanded_hart_df[["Entrez_ID"]][hart_are_selected],
                                         "all"       = expanded_hart_df[["Entrez_ID"]]
                                         ),

    "Hart_HeLa"                   = list("essential" = expanded_hart_df[["Entrez_ID"]][expanded_hart_df[, "BF_hela"] > 15.47],
                                         "all"       = expanded_hart_df[["Entrez_ID"]]
                                         ),

    "Blomen_HAP1_KBM7_intersect"  = list("essential" = blomen_df[["Entrez_ID"]][blomen_are_selected],
                                         "all"       = blomen_all_entrezs
                                         )

)



# Process protein-coding genes --------------------------------------------

all_entrezs <- union(collected_entrez_IDs,
                     unlist(sublibraries_all_entrezs_list, use.names = FALSE)
                     )
all_entrezs <- sort(as.integer(all_entrezs))

essential_df <- GetGeneEssentiality(all_entrezs, essential_datasets_list)

stopifnot(identical(all_entrezs, essential_df[["Entrez_ID"]]))





# Export data -------------------------------------------------------------

essential_export_df <- essential_df
for (column_name in c("Achilles_mean_probability", "Achilles_num_essential", "Achilles_num_cell_lines")) {
  essential_export_df[, column_name] <- ifelse(is.na(essential_export_df[, column_name]),
                                               "",
                                               essential_export_df[, column_name]
                                               )
}

write.table(essential_export_df,
            file = file.path(file_output_directory, "Essential_genes.tsv"),
            sep = "\t", na = "N/A", quote = FALSE, row.names = FALSE
            )






# Categorize into three groups by essentiality ----------------------------

essential_fraction <- essential_df[["Achilles_num_essential"]] /
                      essential_df[["Achilles_num_cell_lines"]]
essential_vec <- ifelse(essential_fraction < 0.1,
                        "Non-essential",
                        ifelse(essential_fraction > 0.9,
                               "Essential",
                               "Intermediate"
                               )
                        )
essential_fac <- factor(essential_vec,
                        levels = c("Essential", "Intermediate", "Non-essential"),
                        ordered = TRUE
                        )

table(essential_fac, essential_df[["Achilles_common"]])

essential_df[["Three_categories"]] <- essential_fac







# Save data ---------------------------------------------------------------

save(list = "essential_df",
     file = file.path(general_RData_directory, "22) Compile data on essential genes.RData")
     )









