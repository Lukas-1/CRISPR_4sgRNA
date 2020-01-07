### 26th September 2019 ###





# Import packages and source code -----------------------------------------

library("biomaRt")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
human_TFs_directory     <- file.path(CRISPR_input_directory, "Human genome", "Human TFs - Lambert and Jolma")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))





# Read in data ------------------------------------------------------------

# Downloaded from http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
# on 3 October 2019

tf_new_df <- read.csv(file.path(human_TFs_directory, "DatabaseExtract_v_1.01.csv"), stringsAsFactors = FALSE, check.names = FALSE)[, -1]





# Rename the columns ------------------------------------------------------

names(tf_new_df) <- c(
  "ENSEMBL_gene_ID",
  "Gene_symbol",
  "DNA_binding_domain",
  "Is_TF",
  "TF_assessment",
  "Binding_mode",
  "Motif_status",
  "Final_notes",
  "Final_comments",
  "Interpro_ID",
  "Original_Entrez_ID",
  "Entrez_description",
  "PBD_ID",
  "Tested_by_HT_SELEX",
  "Tested_by_PBM",
  "Conditional_binding_requirements",
  "Original_comments",

  "Vaquerizas_2009_TF_classification",
  "Is_TF_CisBP",
  "TFCAT_classification",
  "Is_TF_GO",

  "Initial_assessment",
  "Curator1_name",
  "Curator2_name",

  "Is_TF_TFClass",

  "GO_evidence",
  "Pfam_domains",

  "Is_C2H2_ZF"
)






# Explore the TF databases ------------------------------------------------

table(tf_new_df[, "Original_Entrez_ID"])[table(tf_new_df[, "Original_Entrez_ID"]) > 1]

tf_new_df[tf_new_df[, "Original_Entrez_ID"] == "107984053", ]
tf_new_df[tf_new_df[, "Original_Entrez_ID"] == "159119", ]
tf_new_df[tf_new_df[, "Original_Entrez_ID"] == "728689", ]
tf_new_df[tf_new_df[, "Original_Entrez_ID"] == "Has no Entrez identifier", ]
tf_new_df[tf_new_df[, "Original_Entrez_ID"] == "None", ]

table(tf_new_df[, "Binding_mode"], useNA = "ifany")





# Map TF gene symbols to Entrez IDs ---------------------------------------

symbol_mappings_df <- MapToEntrezs(symbols_vec = tf_new_df[, "Gene_symbol"])
symbol_mappings_df <- symbol_mappings_df[, names(symbol_mappings_df) != "Original_entrez"]
problematic_symbols_df <- symbol_mappings_df[!(symbol_mappings_df[, "Entrez_source"] %in% 1), ]





# Create an ENSEMBL-Entrez mapping using BioMart --------------------------

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

ensembl_to_entrez_df <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = tf_new_df[, "ENSEMBL_gene_ID"],
  mart = mart
)





# Map TF ENSEMBL IDs to Entrez IDs for the new file -----------------------

ensembl_mappings_df <- data.frame(
  tf_new_df[, c("ENSEMBL_gene_ID", "Gene_symbol", "Original_Entrez_ID", "Is_TF")],
  "OrgHs_entrez"         = vapply(ensembl_to_entrez_list[tf_new_df[, "ENSEMBL_gene_ID"]], function(x) if (length(x) == 0) NA_character_ else paste0(x, collapse = ", "), ""),
  "OrgHs_num_mappings"   = lengths(ensembl_to_entrez_list[tf_new_df[, "ENSEMBL_gene_ID"]]),
  "BioMart_entrez"       = as.character(ensembl_to_entrez_df[match(tf_new_df[, "ENSEMBL_gene_ID"], ensembl_to_entrez_df[, "ensembl_gene_id"]), "entrezgene_id"]),
  "Symbol_to_entrez"     = symbol_mappings_df[, "Entrez_ID"],
  "Symbol_entrez_symbol" = symbol_mappings_df[, "Gene_symbol"],
  "Symbol_entrez_source" = symbol_mappings_df[, "Entrez_source"],
  stringsAsFactors       = FALSE,
  row.names              = NULL
)

are_not_identical <- !(mapply(identical, ensembl_mappings_df[, "OrgHs_entrez"], ensembl_mappings_df[, "BioMart_entrez"], USE.NAMES = FALSE))
are_NA <- is.na(ensembl_mappings_df[, "BioMart_entrez"])

ensembl_mappings_df[, "Consensus_entrez"] <- vapply(seq_len(nrow(ensembl_mappings_df)), function(x) {
  this_vec <- unlist(ensembl_mappings_df[x, c("BioMart_entrez", "OrgHs_entrez", "Symbol_to_entrez")])
  if (!(are_NA[[x]]) && !(are_not_identical[[x]])) {
    this_vec[["BioMart_entrez"]]
  } else if (!(is.na(this_vec[["BioMart_entrez"]]))) {
    if ((!(is.na(this_vec[["OrgHs_entrez"]]))) &&
        identical(this_vec[["OrgHs_entrez"]], this_vec[["Symbol_to_entrez"]])
        ) { # majority vote
      this_vec[["OrgHs_entrez"]]
    } else if ((!(is.na(this_vec[["OrgHs_entrez"]]))) &&
                any(strsplit(this_vec[["OrgHs_entrez"]], ", ", fixed = TRUE)[[1]] %in% this_vec[["Symbol_to_entrez"]])
               ) { # this had to be introduced to keep the number of transcription factors at 2765
      this_vec[["Symbol_to_entrez"]]
    } else {
      this_vec[["BioMart_entrez"]]
    }
  } else if (!(is.na(this_vec[["OrgHs_entrez"]]))) {
    if (grepl(", ", this_vec[["OrgHs_entrez"]], fixed = TRUE)) {
      if (!(grepl(", ", this_vec[["Symbol_to_entrez"]], fixed = TRUE)) && grepl(this_vec[["Symbol_to_entrez"]], this_vec[["OrgHs_entrez"]], fixed = TRUE)) {
        this_vec[["Symbol_to_entrez"]]
      } else {
        NA_character_
      }
    } else {
      this_vec[["OrgHs_entrez"]]
    }
  } else if (!(is.na(this_vec[["Symbol_to_entrez"]])) && !(grepl(", ", this_vec[["Symbol_to_entrez"]], fixed = TRUE))) {
    this_vec[["Symbol_to_entrez"]]
  } else {
    NA_character_
  }
}, "")



if (any(grepl(", ", ensembl_mappings_df[, "Consensus_entrez"]))) {
  stop("Ambiguous Entrez ID mappings were found!")
}

# Check for duplicate mappings that would lead to fewer than 2765 entries being found in the final database!
num_occurrences_vec <- table(ensembl_mappings_df[, "Consensus_entrez"])[ensembl_mappings_df[, "Consensus_entrez"]]
if (any(num_occurrences_vec > 1, na.rm = TRUE)) {
  message("\nProblematic genes:")
  print(ensembl_mappings_df[(num_occurrences_vec > 1) %in% TRUE, ])
  stop("Multiple ensembl IDs seemed to map to same gene or Entrez ID!")
}



final_symbols_df <- MapToEntrezs(ensembl_mappings_df[, "Consensus_entrez"], ensembl_mappings_df[, "Gene_symbol"])

ensembl_mappings_df[, "Consensus_symbol"] <- ifelse(is.na(ensembl_mappings_df[, "Consensus_entrez"]), NA_character_, final_symbols_df[, "Gene_symbol"])
ensembl_mappings_df[, "Original_symbol"] <- final_symbols_df[, "Original_symbol"]
ensembl_mappings_df[, "Combined_ID"] <- ifelse(is.na(ensembl_mappings_df[, "Consensus_entrez"]),
                                               toupper(ensembl_mappings_df[, "Original_symbol"]),
                                               ensembl_mappings_df[, "Consensus_entrez"]
                                               )






# Examine transcription factors with problematic mappings -----------------

ensembl_mappings_df[are_NA | are_not_identical, ]

identical_to_original_entrez <- mapply(identical, ensembl_mappings_df[, "Consensus_entrez"], ensembl_mappings_df[, "Original_Entrez_ID"], USE.NAMES = FALSE)
original_entrez_available <- !(ensembl_mappings_df[, "Original_Entrez_ID"] %in% c("Has no Entrez identifier", "None"))

ensembl_mappings_df[is.na(ensembl_mappings_df[, "Consensus_entrez"]), ]
ensembl_mappings_df[(!(identical_to_original_entrez)) & original_entrez_available, ]






# Select columns for a summary TF data frame ------------------------------

ensembl_mappings_df[identical_to_original_entrez, "Original_Entrez_ID"] <- ""

all_TF_df <- data.frame(
  "Combined_ID"        = ensembl_mappings_df[, "Combined_ID", drop = FALSE],
  "Entrez_ID"          = ensembl_mappings_df[, "Consensus_entrez"],
  "Gene_symbol"        = ensembl_mappings_df[, "Consensus_symbol"],
  ensembl_mappings_df[, "Original_symbol", drop = FALSE],
  "Original_Entrez_ID" = ifelse(identical_to_original_entrez, "", ensembl_mappings_df[, "Original_Entrez_ID"]),
  tf_new_df[, c("ENSEMBL_gene_ID", "Is_TF", "DNA_binding_domain", "TF_assessment", "Binding_mode", "Is_TF_CisBP", "Is_TF_TFClass", "Is_TF_GO")],
  "Is_C2H2_ZF"         = c("False" = "No", "True" = "Yes")[tf_new_df[, "Is_C2H2_ZF"]],
  stringsAsFactors     = FALSE,
  row.names            = NULL
)






# Save data ---------------------------------------------------------------

save(list = "all_TF_df",
     file = file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData")
     )




