### 2022-04-22 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))
source(file.path(general_functions_directory, "32) Compiling data on essential genes.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

project_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 4sg - first trial")
input_dir <- file.path(project_dir, "02_input_data")
rdata_dir <- file.path(project_dir, "03_R_objects")
output_dir <- file.path(project_dir, "04_output_data", "Essential genes")

essential_directory <- file.path(input_dir, "Essential genes")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")

DepMap2020Q2_dir <- file.path(input_dir, "Essential genes", "Used for Nunez et al")



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))
load(file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))




# Read in lists of essential genes from original publications -------------

blomen_file_prefix <- "2015 - Gene essentiality and synthetic lethality in haploid human cells - Table S"
blomen_df <- data.frame(read_excel(file.path(essential_directory,
                                             "Publications", paste0(blomen_file_prefix, "3.xlsx")
                                             ),
                                   skip = 1, na = "NA", .name_repair = "minimal"
                                   ),
                        stringsAsFactors = FALSE, check.names = FALSE
                        )
names(blomen_df)[3:8] <- paste0("KBM7_", names(blomen_df)[3:8])
names(blomen_df)[11:16] <- paste0("HAP1_", names(blomen_df)[11:16])
blomen_df <- blomen_df[, !(duplicated(as.list(blomen_df)))]




blomen_KBM7_all_genes_df <- data.frame(read_excel(file.path(essential_directory,
                                                            "Publications", paste0(blomen_file_prefix, "1.xlsx")
                                                            ),
                                                  skip = 1, na = "NA"
                                                  ),
                                       stringsAsFactors = FALSE, check.names = FALSE
                                       )

blomen_HAP1_all_genes_df <- data.frame(read_excel(file.path(essential_directory,
                                                            "Publications", paste0(blomen_file_prefix, "2.xlsx")
                                                            ),
                                                  skip = 1, na = "NA"
                                                  ),
                                       stringsAsFactors = FALSE, check.names = FALSE
                                       )



hart_file_name <- "2015 - High-Resolution CRISPR Screens Reveal Fitness Genes - Table S2.xlsx"
hart_df <- data.frame(read_excel(file.path(essential_directory,
                                           "Publications", hart_file_name
                                           )
                                 ),
                      stringsAsFactors = FALSE
                      )



# Read in lists of essential genes from DepMap / Project Achilles ---------
## Download link: https://depmap.org/portal/download/all/

ReadDepMapFile <- function(file_name, sub_folder = "2022_Q1") {
  read.csv(file.path(essential_directory,
                     "DepMap", sub_folder,
                     paste0(file_name, ".csv")
                     ),
           stringsAsFactors = FALSE, check.names = FALSE
           )
}

depmap_hart_blomen_df        <- ReadDepMapFile("common_essentials")

achilles_common_df           <- ReadDepMapFile("Achilles_common_essentials")
achilles_original_effects_df <- ReadDepMapFile("Achilles_gene_effect")
achilles_original_depend_df  <- ReadDepMapFile("Achilles_gene_dependency")

CRISPR_common_df             <- ReadDepMapFile("CRISPR_common_essentials")
CRISPR_original_effects_df   <- ReadDepMapFile("CRISPR_gene_effect")
CRISPR_original_depend_df    <- ReadDepMapFile("CRISPR_gene_dependency")

depmap_samples_df            <- ReadDepMapFile("sample_info")

DEMETER2_combined_original_depend_df <- ReadDepMapFile("gene_dependency", sub_folder = "DEMETER 2 Combined RNAi")
DEMETER2_combined_original_effects_df <- ReadDepMapFile("gene_effect", sub_folder = "DEMETER 2 Combined RNAi")
DEMETER2_samples_df <- ReadDepMapFile("sample_info", sub_folder = "DEMETER 2 Combined RNAi")



# Read in the data used for Nunez et al., 2021 (Cell) ---------------------

essentials_2020Q2_df     <- read.csv(file.path(DepMap2020Q2_dir, "DepMap_20Q2__common_essentials.csv"), stringsAsFactors = FALSE)
non_essentials_2020Q2_df <- read.csv(file.path(DepMap2020Q2_dir, "DepMap_20Q2__nonessentials.csv"), stringsAsFactors = FALSE)



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

DEMETER2_combined_depend_df  <- ProcessDEMETERDataDf(DEMETER2_combined_original_depend_df, check_replacements = FALSE)
DEMETER2_combined_effects_df <- ProcessDEMETERDataDf(DEMETER2_combined_original_effects_df, check_replacements = FALSE)

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

    "CRISPR_common"               = list("essential" = CRISPR_common_df[["Entrez_ID"]],
                                         "all"       = union(achilles_common_df[["Entrez_ID"]],
                                                             achilles_depend_df[["Entrez_ID"]]
                                                             )
                                         ),

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



# Process CRISPRoff 2sg library genes -------------------------------------

all_entrezs <- CRISPRoff_df[, "Entrez_ID"]
all_entrezs <- all_entrezs[!(is.na(all_entrezs))]
all_entrezs <- as.integer(unlist(strsplit(all_entrezs, ", ", fixed = TRUE)))

essential_df <- GetGeneEssentiality(all_entrezs, essential_datasets_list)

stopifnot(identical(all_entrezs, essential_df[["Entrez_ID"]]))





# Export data -------------------------------------------------------------

essential_export_df <- essential_df
NA_empty_columns <- c("CRISPR_mean_probability",   "CRISPR_num_essential",   "CRISPR_num_cell_lines",
                      "CRISPR_mean_effect",
                      "Achilles_mean_probability", "Achilles_num_essential", "Achilles_num_cell_lines",
                      "DEMETER2_mean_probability", "DEMETER2_num_essential", "DEMETER2_num_cell_lines"
                      )
for (column_name in NA_empty_columns) {
  essential_export_df[, column_name] <- ifelse(is.na(essential_export_df[, column_name]),
                                               "",
                                               essential_export_df[, column_name]
                                               )
}

essential_export_df <- essential_export_df[, !(names(essential_export_df) %in% c("Achilles_mean_probability", "Achilles_num_essential", "Achilles_num_cell_lines"))]

write.table(essential_export_df,
            file = file.path(output_dir, "Essential_genes_CRISPRoff_2sg.tsv"),
            sep = "\t", na = "N/A", quote = FALSE, row.names = FALSE
            )








# Categorize into four categories by essentiality -------------------------

essential_fraction <- essential_df[["CRISPR_num_essential"]] /
                      essential_df[["CRISPR_num_cell_lines"]]
essential_vec <- ifelse(essential_fraction == 0,
                        "Non-essential",
                        ifelse(essential_fraction > 0.95,
                               "Essential",
                               ifelse((essential_fraction > 0.1) & (essential_fraction < 0.15),
                                      "Intermediate",
                                      "Other"
                                      )
                               )
                        )
essential_fac <- factor(essential_vec,
                        levels = c("Non-essential", "Intermediate", "Essential", "Other")
                        )

essential_df[["Four_categories"]] <- essential_fac




# Draw histograms ---------------------------------------------------------

categ_mat <- sapply(levels(essential_df[["Four_categories"]]), function(x) {
  are_this_category <- essential_df[["Four_categories"]] %in% x
  CRISPR_effects_df[["Entrez_ID"]] %in% essential_df[["Entrez_ID"]][are_this_category]
})
DrawEssentialityHistograms()



# Tidy the data used for Nunez et al., 2021 (Cell) ------------------------

TidyDepMapGeneList <- function(input_df) {
  splits_list <- strsplit(input_df[, 1], " (", fixed = TRUE)
  results_df <- data.frame(
    "Gene_symbol" = sapply(splits_list, "[[", 1),
    "Entrez_ID"   = sub(")", "", sapply(splits_list, "[[", 2), fixed = TRUE),
    stringsAsFactors = FALSE
  )
  results_df[, "Entrez_ID"] <- as.integer(results_df[, "Entrez_ID"])
  return(results_df)
}

essentials_2020Q2_df <- TidyDepMapGeneList(essentials_2020Q2_df)
non_essentials_2020Q2_df <- TidyDepMapGeneList(non_essentials_2020Q2_df)




# Save data ---------------------------------------------------------------

save(list = c("essentials_2020Q2_df", "non_essentials_2020Q2_df"),
     file = file.path(rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData")
     )

save(list = "essential_df",
     file = file.path(rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData")
     )

save(list = c("essential_datasets_list", "achilles_depend_df",
              "CRISPR_depend_df", "DEMETER2_combined_depend_df",
              "CRISPR_effects_df"
              ),
     file = file.path(rdata_dir, "05_compile_data_on_essential_genes__datasets.RData")
     )



