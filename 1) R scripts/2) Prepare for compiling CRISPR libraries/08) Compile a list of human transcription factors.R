### 26th September 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
human_TFs_directory     <- file.path(CRISPR_input_directory, "Sublibraries", "Human TFs - Lambert and Jolma")
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
  "Ensembl_gene_ID",
  "Gene_symbol",
  "DNA_binding_domain",
  "Is_TF",
  "TF_assessment",
  "Binding_mode",
  "Motif_status",
  "Final_notes",
  "Final_comments",
  "Interpro_ID",
  "Original_entrez",
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

table(tf_new_df[["Original_entrez"]])[table(tf_new_df[["Original_entrez"]]) > 1]

tf_new_df[tf_new_df[["Original_entrez"]] == "107984053", ]
tf_new_df[tf_new_df[["Original_entrez"]] == "159119", ]
tf_new_df[tf_new_df[["Original_entrez"]] == "728689", ]
tf_new_df[tf_new_df[["Original_entrez"]] == "Has no Entrez identifier", ]
tf_new_df[tf_new_df[["Original_entrez"]] == "None", ]

table(tf_new_df[["Binding_mode"]], useNA = "ifany")







# Map transcription factors to Entrez IDs ---------------------------------

ensembl_mappings_df <- MapEnsemblIDs(tf_new_df)







# Examine transcription factors with problematic mappings -----------------

ensembl_mappings_df[ensembl_mappings_df[["Are_problematic"]], ]

identical_to_original_entrez <- mapply(identical, ensembl_mappings_df[["Consensus_entrez"]], tf_new_df[["Original_entrez"]], USE.NAMES = FALSE)
original_entrez_available <- !(ensembl_mappings_df[["Original_entrez"]] %in% c("Has no Entrez identifier", "None"))

ensembl_mappings_df[is.na(ensembl_mappings_df[["Consensus_entrez"]]), ]
ensembl_mappings_df[(!(identical_to_original_entrez)) & original_entrez_available, ]






# Select columns for the TF data frame ------------------------------------

ensembl_mappings_df[["Original_entrez"]][identical_to_original_entrez] <- ""

all_TF_df <- data.frame(
  ensembl_mappings_df["Combined_ID"],
  "Entrez_ID"       = ensembl_mappings_df[["Consensus_entrez"]],
  "Gene_symbol"     = ensembl_mappings_df[["Consensus_symbol"]],
  ensembl_mappings_df["Original_symbol"],
  "Original_entrez" = ifelse(identical_to_original_entrez, "", tf_new_df[["Original_entrez"]]),
  tf_new_df[, c("Ensembl_gene_ID", "Is_TF", "DNA_binding_domain", "TF_assessment", "Binding_mode", "Is_TF_CisBP", "Is_TF_TFClass", "Is_TF_GO")],
  "Is_C2H2_ZF"      = c("False" = "No", "True" = "Yes")[tf_new_df[["Is_C2H2_ZF"]]],
  stringsAsFactors  = FALSE,
  row.names         = NULL
)






# Save data ---------------------------------------------------------------

save(list = "all_TF_df",
     file = file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData")
     )




