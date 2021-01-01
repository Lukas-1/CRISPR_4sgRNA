### 22nd December 2020 ###




# Import packages and source code -----------------------------------------

library("readxl")
library("org.Hs.eg.db")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
file_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-12-22 - collect data on essential genes")
file_input_directory      <- file.path(file_directory, "1) Input")
essential_genes_directory <- file.path(file_input_directory, "Essential gene lists")
vacuolation_directory     <- file.path(file_input_directory, "Vacuolation genes")
file_output_directory     <- file.path(file_directory, "2) Output")




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

depmap_hart_blomen_df <- read.csv(file.path(essential_genes_directory,
                                            "DepMap",
                                            "common_essentials_20Q4.csv"
                                            ),
                                  stringsAsFactors = FALSE
                                  )


achilles_common_df <- read.csv(file.path(essential_genes_directory,
                                         "DepMap",
                                         "Achilles_common_essentials_20Q4_v2.csv"
                                         ),
                               stringsAsFactors = FALSE
                               )

achilles_original_effects_df <- read.csv(file.path(essential_genes_directory,
                                                   "DepMap",
                                                   "Achilles_gene_effect_20Q4_v2.csv"
                                                   ),
                                         stringsAsFactors = FALSE, check.names = FALSE
                                         )

achilles_original_depend_df <- read.csv(file.path(essential_genes_directory,
                                                  "DepMap",
                                                  "Achilles_gene_dependency_20Q4_v2.csv"
                                                  ),
                                        stringsAsFactors = FALSE, check.names = FALSE
                                        )

depmap_samples_df <- read.csv(file.path(essential_genes_directory,
                                        "DepMap",
                                        "sample_info_20Q4.csv"
                                        ),
                              stringsAsFactors = FALSE
                              )



# Read in the hit list of vacuolation-related genes -----------------------

vacuolation_genes <- read.table(file.path(vacuolation_directory, "hit_gene_IDs_1.5allhitnoderatio.csv"),
                                header = FALSE, stringsAsFactors = FALSE
                                )[[1]]
LSm_genes         <- read.table(file.path(vacuolation_directory, "Additional_LSm_genes.txt"),
                                header = FALSE, stringsAsFactors = FALSE
                                )[[1]]

vac_CRISPRi_df    <- data.frame(read_excel(file.path(vacuolation_directory,
                                                     "CRISPRi - 4sg vacuolation - re-ordered.xlsx"
                                                     ),
                                           guess_max = 10000
                                           ),
                                check.names = FALSE, stringsAsFactors = FALSE
                                )




# Define gene annotation lookup lists -------------------------------------

ensembl_to_entrez_list <- as.list(org.Hs.egENSEMBL2EG[mappedkeys(org.Hs.egENSEMBL2EG)])
symbol_to_entrez_list  <- as.list(org.Hs.egSYMBOL2EG[mappedkeys(org.Hs.egSYMBOL2EG)])
alias_to_entrez_list   <- as.list(org.Hs.egALIAS2EG[mappedkeys(org.Hs.egALIAS2EG)])
entrez_to_symbol_list  <- as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])




# Define functions --------------------------------------------------------

ENSEMBLandHUGOtoEntrez <- function(ensg_vec, symbols_vec) {

  stopifnot(length(ensg_vec) == length(symbols_vec))

  from_ensembl_list <- ensembl_to_entrez_list[ensg_vec]
  from_symbol_list <- symbol_to_entrez_list[symbols_vec]

  results_list <- lapply(seq_along(symbols_vec), function(x) {
    if (length(from_ensembl_list[[x]]) == 1) {
      result <- from_ensembl_list[[x]]
    } else if (is.null(from_ensembl_list[[x]])) {
      if (is.null(from_symbol_list[[x]])) {
        result <- NA_character_
      } else {
        result <- from_symbol_list[[x]]
      }
    } else {
      if (any(from_ensembl_list[[x]] %in% from_symbol_list[[x]])) {
        result <- intersect(from_ensembl_list[[x]], from_symbol_list[[x]])
      } else {
        result <- from_ensembl_list[[x]]
      }
    }
    return(as.integer(result))
  })
  return(results_list)
}



ProcessAchillesGenesDf <- function(genes_df) {
  # genes_df is a single column in the format A1BG (1) {1 being the Entrez ID}

  gene_splits <- strsplit(genes_df[, "gene"], " ", fixed = TRUE)
  entrezs_vec <- sapply(gene_splits, "[[", 2)
  stopifnot(all(substr(entrezs_vec, 1, 1) == "("))
  stopifnot(all(substr(entrezs_vec, nchar(entrezs_vec), nchar(entrezs_vec)) == ")"))

  stripped_entrezs <- as.integer(substr(entrezs_vec, 2, nchar(entrezs_vec) - 1))
  are_present <- stripped_entrezs %in% names(entrez_to_symbol_list)
  if (!(all(are_present))) {
    message(paste0("The following Entrez IDs were not found in ",
                   "org.Hs.egSYMBOL: ",
                   paste0(stripped_entrezs[!(are_present)], collapse = ", ")
                   ))
  }

  stripped_entrezs <- as.integer(stripped_entrezs)
  original_symbols <- sapply(gene_splits, "[[", 1)
  results_df <- data.frame(
    "Entrez_ID"       = stripped_entrezs,
    "Original_symbol" = original_symbols,
    stringsAsFactors = FALSE
  )
  return(results_df)
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
    "Mean" = rowMeans(data_mat),
    data_mat,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
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


  TwoByTwoTable <- function(essential_genes, all_genes, dataset_label = NULL) {
    assign("delete_essential_genes", essential_genes, envir = globalenv())
    assign("delete_all_genes", all_genes, envir = globalenv())
    assign("delete_dataset_label", dataset_label, envir = globalenv())
    stopifnot(all(essential_genes %in% all_genes))
    are_candidates <- all_genes %in% entrezs_vec
    are_essential <- all_genes %in% essential_genes
    two_by_two_table <- table(are_candidates, are_essential)
    fisher_results <- fisher.test(two_by_two_table,
                                  alternative = "two.sided"
                                  )
    results_list <- list(
      "Total_num_genes"                 = length(all_genes),
      "Num_essential_genes"             = sum(are_essential),
      "Num_available_candidates"        = sum(entrezs_vec %in% all_genes),
      "Neither_candidate_nor_essential" = two_by_two_table[1, 1],
      "Essential_but_not_candidate"     = two_by_two_table[1, 2],
      "Candidate_but_not_essential"     = two_by_two_table[2, 1],
      "Both_candidate_and_essential"    = two_by_two_table[2, 2],
      "Odds_ratio"                      = fisher_results[["estimate"]][[1]],
      "OR_lower_bound"                  = fisher_results[["conf.int"]][[1]],
      "OR_upper_bound"                  = fisher_results[["conf.int"]][[2]],
      "p_value"                         = fisher_results[["p.value"]]
    )
    if (!(is.null(dataset_label))) {
      results_list <- c(
        list("Dataset" = dataset_label),
        results_list
      )
    }
    return(results_list)
  }



  symbols_vec <- vapply(entrez_to_symbol_list[as.character(entrezs_vec)],
                        function(x) if (is.null(x)) NA_character_ else paste0(x, collapse = ", "),
                        ""
                        )



  genes_list <- list(
    "Entrez_ID"                 = entrezs_vec,
    "Gene_symbol"               = symbols_vec,
    "Achilles_mean_probability" = achilles_depend_df[["Mean"]][match(entrezs_vec, achilles_depend_df[["Entrez_ID"]])]
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

  fisher_list <- lapply(seq_along(datasets_list),
                        function(x) TwoByTwoTable(datasets_list[[x]][["essential"]],
                                                  datasets_list[[x]][["all"]],
                                                  names(datasets_list)[[x]]
                                                  )
                        )

  fisher_df <- do.call(rbind.data.frame, c(fisher_list, stringsAsFactors = FALSE, make.row.names = FALSE))

  results_list <- list("genes_df" = genes_df, "fisher_df" = fisher_df)
  return(results_list)
}




# Process vacuolation-related genes ---------------------------------------

vac_df <- data.frame("Category"            = c(rep("LSm genes", length(LSm_genes)),
                                               rep("Vacuolation hits", length(vacuolation_genes))
                                               ),
                     "Entrez_ID"            = c(LSm_genes, vacuolation_genes),
                     stringsAsFactors = FALSE
                     )
vac_df[["Gene_symbol"]] <- vapply(as.character(vac_df[, "Entrez_ID"]), function(x) {
  gene_symbol <- entrez_to_symbol_list[[x]]
  if (is.null(gene_symbol)) {
    return(NA_character_)
  } else {
    return(gene_symbol)
  }
}, "")

vac_df[["Included_for_CRISPRi"]] <- ifelse(vac_df[["Entrez_ID"]] %in% vac_CRISPRi_df[, "Entrez ID"],
                                           "Yes", "No"
                                           )

new_order <- order(match(vac_df[["Entrez_ID"]], vac_CRISPRi_df[, "Entrez ID"]))
vac_df <- vac_df[new_order, ]
row.names(vac_df) <- NULL






# Create Entrez IDs for data of Blomen et al. (2015) ----------------------



### Essential genes ###

blomen_list <- ENSEMBLandHUGOtoEntrez(blomen_df[, "ENSEMBL_ID"],
                                      blomen_df[, "GENE_SYMBOL"]
                                      )
stopifnot(all(lengths(blomen_list) == 1))

blomen_df[["Entrez_ID"]] <- unlist(blomen_list)


### All genes ###

stopifnot(identical(sort(blomen_KBM7_all_genes_df[["ENSEMBL_ID"]]),
                    sort(blomen_HAP1_all_genes_df[["ENSEMBL_ID"]])
                    ))

blomen_all_entrezs_list <- ENSEMBLandHUGOtoEntrez(blomen_KBM7_all_genes_df[["ENSEMBL_ID"]],
                                                  blomen_KBM7_all_genes_df[["GENE_SYMBOL"]]
                                                  )
blomen_all_entrezs <- unlist(blomen_all_entrezs_list)





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
achilles_common_df <- ProcessAchillesGenesDf(achilles_common_df)



### Achilles gene effects ###

achilles_effects_df <- ProcessAchillesDataDf(achilles_original_effects_df,
                                             depmap_samples_df
                                             )
achilles_depend_df  <- ProcessAchillesDataDf(achilles_original_depend_df,
                                             depmap_samples_df
                                             )

stopifnot(identical(names(achilles_depend_df), names(achilles_depend_df)))
stopifnot(identical(achilles_effects_df[, "Entrez_ID"], achilles_depend_df[, "Entrez_ID"]))




### Check for correspondence between the overlap of Hart and Blomen provided
### by DepMap, and my own overlap of the same datasets

all_intersect_genes <- union(blomen_hart_2_intersect, depmap_hart_blomen_df[["Entrez_ID"]])
in_mine <- all_intersect_genes %in% blomen_hart_2_intersect
in_depmap <- all_intersect_genes %in% depmap_hart_blomen_df[["Entrez_ID"]]

table(in_mine, in_depmap)
table(all_intersect_genes[!(in_depmap)] %in% achilles_depend_df[["Entrez_ID"]])






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





# Annotate vacuolation hits with gene essentiality data -------------------

vac_essential_list <- GetGeneEssentiality(vac_df[["Entrez_ID"]],
                                          essential_datasets_list
                                          )

stopifnot(identical(vac_df[["Entrez_ID"]],
                    vac_essential_list[["genes_df"]][["Entrez_ID"]]
                    )
          )
stopifnot(identical(vac_df[["Gene_symbol"]],
                    vac_essential_list[["genes_df"]][["Gene_symbol"]]
                    )
          )
vacuolation_df <- data.frame(
  vac_df[, c("Category", "Included_for_CRISPRi")],
  vac_essential_list[["genes_df"]],
  stringsAsFactors = FALSE
)




# Draw plots --------------------------------------------------------------

hist(as.matrix(achilles_depend_df[, 4:ncol(achilles_depend_df)]), breaks = 200)

are_common_essentials <- achilles_effects_df[["Entrez_ID"]] %in% achilles_common_df[["Entrez_ID"]]
boxplot(achilles_effects_df[["Mean"]] ~ are_common_essentials)

boxplot(achilles_depend_df[["Mean"]] ~ are_common_essentials)

plot(achilles_effects_df[["Mean"]], achilles_depend_df[["Mean"]])


are_hart_essentials <- hart_df[["numTKOHits"]] != 0
boxplot(hart_df[["BF_hela"]] ~ are_hart_essentials)

range(hart_df[["BF_hela"]], na.rm = TRUE)
range(hart_df[["BF_hela"]][are_hart_essentials], na.rm = TRUE)
range(hart_df[["BF_hela"]][!(are_hart_essentials)], na.rm = TRUE)


table(hart_df[["BF_hela"]] > 15.47)


hist(hart_df[["BF_hela"]], breaks = 200)
hist(hart_df[["BF_hela"]][are_hart_essentials], breaks = 200)
hist(hart_df[["BF_hela"]][!(are_hart_essentials)], breaks = 200)





# Export data -------------------------------------------------------------

vac_export_df <- vacuolation_df
for (column_name in "Achilles_mean_probability") {
  vac_export_df[, column_name] <- ifelse(is.na(vac_export_df[, column_name]),
                                         "",
                                         vac_export_df[, column_name]
                                         )
}

write.table(vac_export_df,
            file = file.path(file_output_directory, "Vacuolation_essential_genes.tsv"),
            sep = "\t", na = "N/A", quote = FALSE, row.names = FALSE
            )

write.table(vac_essential_list[["fisher_df"]],
            file = file.path(file_output_directory, "Vacuolation_overrepresentation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE
            )





