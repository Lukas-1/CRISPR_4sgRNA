### 21st July 2019 ###



# Legacy mode -------------------------------------------------------------

legacy_mode <- TRUE





# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
CRISPR_input_directory           <- file.path(CRISPR_root_directory, "2) Input data")
CRISPR_libraries_directory       <- file.path(CRISPR_input_directory, "CRISPR libraries")
gene_lists_directory             <- file.path(CRISPR_input_directory, "Gene lists")
RData_directory                  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory          <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory          <- file.path(RData_directory, "2) CRISPRa")

CRISPR_libraries_directory       <- file.path(CRISPR_input_directory, "CRISPR libraries")
CRISPRa_datasets_directory       <- file.path(CRISPR_libraries_directory, "CRISPRa")

CRISPRa_Horlbeck2016_path        <- file.path(CRISPRa_datasets_directory, "Horlbeck, Kampmann, Weissman - eLife 2016")
CRISPRa_Horlbeck2016_sgRNAs_path <- file.path(CRISPRa_Horlbeck2016_path, "2016 - Compact and highly active next-generation libraries - Table S5.xlsx")
CRISPRa_Horlbeck2016_TSSs_path   <- file.path(CRISPRa_Horlbeck2016_path, "2016 - Compact and highly active next-generation libraries - Table S2.xlsx")
CRISPRa_Doench2018_path          <- file.path(CRISPRa_datasets_directory, "Sanson, Doench - Nat Comm 2018", "2018 - Optimized libraries for CRISPR-Cas9 genetic screens - Data S5.xlsx")
hand_picked_CRISPRa_path         <- file.path(CRISPRa_datasets_directory, "Manually curated CRISPRa gRNAs.xlsx")
GPP_CRISPRa_path                 <- file.path(CRISPR_root_directory, "4) Intermediate files/CRISPRa/GPP sgRNA designer/2) Output files")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Read in data ------------------------------------------------------------

hCRISPRa_v2_df <- data.frame(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, skip = 7)[-1, ], stringsAsFactors = FALSE)
names(hCRISPRa_v2_df) <- names(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, n_max = 1))

hCRISPRa_v2_TSS_df <- data.frame(read_excel(CRISPRa_Horlbeck2016_TSSs_path), stringsAsFactors = FALSE, check.names = FALSE)

CalabreseA_df  <- data.frame(read_excel(CRISPRa_Doench2018_path, sheet = "SetA sgRNA Annotations"), stringsAsFactors = FALSE, check.names = FALSE)
CalabreseB_df  <- data.frame(read_excel(CRISPRa_Doench2018_path, sheet = "SetB sgRNA Annotations"), stringsAsFactors = FALSE, check.names = FALSE)

### Read in custom selections of genes
hand_picked_df <- data.frame(read_excel(hand_picked_CRISPRa_path), stringsAsFactors = FALSE, check.names = FALSE)
first_trial_df <- read.table(file.path(gene_lists_directory, "Trial_genes.txt"), stringsAsFactors = FALSE, row.names = NULL, quote = "")


### Read in the output from the GPP sgRNA designer tool
GPP_CRISPRa_file_names <- list.files(GPP_CRISPRa_path)
GPP_CRISPRa_all_df <- ReadGPPOutputFiles(GPP_CRISPRa_file_names, GPP_CRISPRa_path)



# Annotate the hCRISPRa_v2 data -------------------------------------------

hCRISPRa_v2_sublibrary_map <- c(
  "h1" = "Kinases/Phosphatases/Drug Targets",
  "h2" = "Cancer/Apoptosis",
  "h3" = "Stress/Proteostasis",
  "h4" = "Mitochondria/Trafficking/Motility",
  "h5" = "Gene Expression",
  "h6" = "Membrane Proteins",
  "h7" = "Unassigned"
)

hCRISPRa_v2_df[, "Sublibrary"] <- hCRISPRa_v2_sublibrary_map[hCRISPRa_v2_df[, "Sublibrary"]]

hCRISPRa_v2_df_transcript_IDs     <- paste0(hCRISPRa_v2_df[, "gene"], "__", hCRISPRa_v2_df[, "transcript"])
hCRISPRa_v2_TSS_df_transcript_IDs <- paste0(hCRISPRa_v2_TSS_df[, "gene"], "__", hCRISPRa_v2_TSS_df[, "transcript"])
TSS_source_matches <- match(hCRISPRa_v2_df_transcript_IDs, hCRISPRa_v2_TSS_df_transcript_IDs)
hCRISPRa_v2_df[, "TSS_source"] <- hCRISPRa_v2_TSS_df[TSS_source_matches, "TSS source"]



# Combine the Calabrese data ----------------------------------------------

Calabrese_df <- rbind.data.frame(
  data.frame(CalabreseA_df, "Calabrese_rank" = "1/2/3", stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(CalabreseB_df, "Calabrese_rank" = "4/5/6", stringsAsFactors = FALSE, check.names = FALSE),
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)




# Process the data from the Broad Institute's GPP portal ------------------

GPP_CRISPRa_df <- TidyGPPOutputDf(GPP_CRISPRa_all_df, CRISPRa_GPP_output_columns)
# GPP_CRISPRa_df <- GPP_CRISPRa_df[GPP_CRISPRa_df[, "Pick Order"] %in% 1:100, ]

table(GPP_CRISPRa_df[, "Pick Order"])
sort(table(GPP_CRISPRa_df[, "Target Gene Symbol"]))






# Collect all unique HUGO gene symbols ------------------------------------

are_Calabrese_controls   <- ((Calabrese_df[, "Annotated Gene Symbol"] == "CONTROL") & (Calabrese_df[, "Annotated Gene ID"] == "CONTROL")) |
                            ((Calabrese_df[, "Annotated Gene Symbol"] == "NO-TARGET") & (Calabrese_df[, "Annotated Gene ID"] == "UNKNOWN_NO-TARGET"))

are_hCRISPRa_v2_controls <- hCRISPRa_v2_df[, "gene"] == "negative_control"

hCRISPRa_v2_symbols_vec <- unique(hCRISPRa_v2_df[!(are_hCRISPRa_v2_controls), "gene"])
Calabrese_symbols_entrezs_df <- unique(Calabrese_df[!(are_Calabrese_controls), c("Annotated Gene Symbol", "Annotated Gene ID")], MARGIN = 1)

names(Calabrese_symbols_entrezs_df) <- c("Gene_symbol", "Entrez_ID")

unique_symbols_vec <- unique(c(hCRISPRa_v2_symbols_vec, Calabrese_symbols_entrezs_df[, "Gene_symbol"][!(toupper(Calabrese_symbols_entrezs_df[, "Gene_symbol"]) %in% toupper(hCRISPRa_v2_symbols_vec))]))





# Check for discrepancies in the Entrez IDs -------------------------------

### The following gene symbols are associated with multiple Entrez IDs in the Calabrese dataset
Calabrese_symbols_entrezs_df[table(Calabrese_symbols_entrezs_df[, "Gene_symbol"])[Calabrese_symbols_entrezs_df[, "Gene_symbol"]] == 2, ]

Calabrese_symbols_to_entrezs <- SymbolsToEntrezIDs(Calabrese_symbols_entrezs_df[, "Gene_symbol"])
Calabrese_entrezs_to_symbols <- EntrezIDsToSymbols(Calabrese_symbols_entrezs_df[, "Entrez_ID"])

Calabrese_symbols_entrezs_df[is.na(Calabrese_symbols_to_entrezs), ] ### The following gene symbols could not be mapped to Entrez IDs
Calabrese_symbols_entrezs_df[is.na(Calabrese_entrezs_to_symbols), ] ### The following Entrez IDs could not be translated to gene symbols

### The following genes have discrepant Entrez IDs (those that were translated from gene symbols, versus those provided by the authors)
are_different_entrezs <- !(is.na(Calabrese_symbols_to_entrezs)) & (Calabrese_symbols_to_entrezs != Calabrese_symbols_entrezs_df[, "Entrez_ID"])
data.frame(Calabrese_symbols_entrezs_df, "Symbol_to_Entrez" = Calabrese_symbols_to_entrezs, stringsAsFactors = FALSE)[are_different_entrezs, ]

### The following genes have discrepant gene symbols (those that were translated from Entrez IDs, versus those provided by the authors)
are_different_symbols <- !(is.na(Calabrese_entrezs_to_symbols)) & (Calabrese_entrezs_to_symbols != Calabrese_symbols_entrezs_df[, "Gene_symbol"])
data.frame(Calabrese_symbols_entrezs_df, "Entrez_to_Symbol" = Calabrese_entrezs_to_symbols, stringsAsFactors = FALSE)[are_different_symbols, ]





# Harmonize the data frames -----------------------------------------------

mapped_Calabrese_df   <- MapToEntrezs(Calabrese_df[, "Annotated Gene ID"], Calabrese_df[, "Annotated Gene Symbol"])
mapped_hCRISPRa_v2_df <- MapToEntrezs(symbols_vec = hCRISPRa_v2_df[, "gene"])
mapped_GPP_df         <- MapToEntrezs(as.character(GPP_CRISPRa_df[, "Target Gene ID"]), GPP_CRISPRa_df[, "Target Gene Symbol"])
mapped_manual_df      <- MapToEntrezs(symbols_vec = hand_picked_df[, "Gene symbol"])


new_Calabrese_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "Calabrese",
  "hCRISPRa_v2_transcript"    = NA_character_,
  "Is_control"                = ifelse(are_Calabrese_controls, "Yes", "No"),
  mapped_Calabrese_df[, 1:4],
  "Entrez_source_Calabrese"   = mapped_Calabrese_df[, "Entrez_source"],
  "Entrez_source_hCRISPRa_v2" = NA_integer_,
  "sgRNA_sequence"            = Calabrese_df[, "sgRNA Sequence"],
  "Original_PAM"              = NA_character_,
  "Calabrese_rank"            = Calabrese_df[, "Calabrese_rank"],
  "GPP_rank"                  = NA_integer_,
  "hCRISPRa_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "hCRISPRa_v2_ID"            = NA_character_,
  "hCRISPRa_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)

new_hCRISPRa_v2_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "hCRISPRa-v2",
  "hCRISPRa_v2_transcript"    = ifelse(hCRISPRa_v2_df[, "transcript"] == "na",
                                       NA_character_,
                                       gsub(",", ", ", hCRISPRa_v2_df[, "transcript"], fixed = TRUE)
                                       ),
  "Is_control"                = ifelse(are_hCRISPRa_v2_controls, "Yes", "No"),
  mapped_hCRISPRa_v2_df[, 1:4],
  "Entrez_source_Calabrese"   = NA_integer_,
  "Entrez_source_hCRISPRa_v2" = mapped_hCRISPRa_v2_df[, "Entrez_source"],
  "sgRNA_sequence"            = hCRISPRa_v2_df[, "protospacer sequence"],
  "Original_PAM"              = NA_character_,
  "Calabrese_rank"            = NA_integer_,
  "GPP_rank"                  = NA_integer_,
  "hCRISPRa_v2_rank"          = ifelse(is.na(hCRISPRa_v2_df[, "selection rank"]), hCRISPRa_v2_df[, "Sublibrary half"], hCRISPRa_v2_df[, "selection rank"]),
  "Predicted_score"           = hCRISPRa_v2_df[, "predicted score"],
  "Empirical_score"           = hCRISPRa_v2_df[, "empirical score"],
  "Off_target_stringency"     = hCRISPRa_v2_df[, "off-target stringency"],
  "Sublibrary"                = hCRISPRa_v2_df[, "Sublibrary"],
  "hCRISPRa_v2_ID"            = sub(".23", "", hCRISPRa_v2_df[, "sgID"], fixed = TRUE),
  "hCRISPRa_TSS_source"       = hCRISPRa_v2_df[, "TSS_source"],
  stringsAsFactors            = FALSE
)


new_GPP_CRISPRa_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "GPP",
  "hCRISPRa_v2_transcript"    = NA_character_,
  "Is_control"                = "No",
  mapped_GPP_df[, 1:4],
  "Entrez_source_Calabrese"   = NA_integer_,
  "Entrez_source_hCRISPRa_v2" = NA_integer_,
  "sgRNA_sequence"            = GPP_CRISPRa_df[, "sgRNA Sequence"],
  "Original_PAM"              = GPP_CRISPRa_df[, "PAM Sequence"],
  "Calabrese_rank"            = NA_integer_,
  "GPP_rank"                  = GPP_CRISPRa_df[, "Pick Order"],
  "hCRISPRa_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "hCRISPRa_v2_ID"            = NA_character_,
  "hCRISPRa_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)


manual_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "Curated",
  "hCRISPRa_v2_transcript"    = NA_character_,
  "Is_control"                = "No",
  mapped_manual_df[, 1:4],
  "Entrez_source_Calabrese"   = NA_integer_,
  "Entrez_source_hCRISPRa_v2" = NA_integer_,
  "sgRNA_sequence"            = hand_picked_df[, "gRNA sequence"],
  "Original_PAM"              = NA_character_,
  "Calabrese_rank"            = NA_character_,
  "GPP_rank"                  = NA_integer_,
  "hCRISPRa_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "hCRISPRa_v2_ID"            = NA_character_,
  "hCRISPRa_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)





# Check for multiple occurrences of the same sgRNA sequence ---------------

num_occurrences_hCRISPRa_v2 <- table(new_hCRISPRa_v2_df[, "sgRNA_sequence"])[new_hCRISPRa_v2_df[, "sgRNA_sequence"]]

multiple_occurrences_df <- data.frame(new_hCRISPRa_v2_df, "Num_occurrences" = as.integer(num_occurrences_hCRISPRa_v2))
multiple_occurrences_df <- multiple_occurrences_df[order(multiple_occurrences_df[, "sgRNA_sequence"]), ]
multiple_occurrences_df[multiple_occurrences_df[, "Num_occurrences"] > 1, ]

any(duplicated(new_Calabrese_df[, "sgRNA_sequence"]))




# Build a combined data frame of all genes --------------------------------

any(is.na(new_Calabrese_df[, "Entrez_ID"]) & (new_Calabrese_df[, "Original_symbol"] == ""))

combined_df <- rbind.data.frame(manual_df, new_Calabrese_df, new_hCRISPRa_v2_df, new_GPP_CRISPRa_df, stringsAsFactors = FALSE, make.row.names = FALSE)
are_controls <- combined_df[, "Is_control"] == "Yes"

control_names_map <- c("CONTROL" = "Control", "NEGATIVE_CONTROL" = "Control", "NO-TARGET" = "No-target")
combined_df[are_controls, "Original_symbol"] <- control_names_map[combined_df[are_controls, "Original_symbol"]]

combined_df[, "Combined_ID"] <- ifelse(are_controls,
                                       "Control",
                                       ifelse(is.na(combined_df[, "Entrez_ID"]),
                                              toupper(combined_df[, "Original_symbol"]),
                                              combined_df[, "Entrez_ID"]
                                              )
                                       )

if (legacy_mode) {
  CRISPRa_df <- ResolveDuplicates(combined_df, concatenate_columns = c("Sublibrary", "hCRISPRa_v2_ID"))
} else {
  CRISPRa_df <- ResolveDuplicates(combined_df, concatenate_columns = c("Sublibrary", "hCRISPRa_v2_ID", "hCRISPRa_TSS_source"))
}





# Search for duplicated sgRNAs within the same gene -----------------------

sgID_vec <- paste0(CRISPRa_df[, "Combined_ID"], "__", toupper(CRISPRa_df[, "sgRNA_sequence"]))
num_occurrences <- table(sgID_vec)[sgID_vec]
have_duplications_sgRNA_list <- split(CRISPRa_df[num_occurrences > 1, ],
                                      factor(sgID_vec[num_occurrences > 1], levels = unique(sgID_vec[num_occurrences > 1]))
                                      )







# Count genes with multiple transcripts in the hCRISPRa-v2 database -------

are_hCRISPRa_v2 <- grepl("hCRISPRa-v2", CRISPRa_df[, "Source"], fixed = TRUE)
transcripts_list <- tapply(CRISPRa_df[are_hCRISPRa_v2, "hCRISPRa_v2_transcript"],
                           factor(CRISPRa_df[are_hCRISPRa_v2, "Combined_ID"], levels = unique(CRISPRa_df[are_hCRISPRa_v2, "Combined_ID"])),
                           function(x) unique(x[!(is.na(x))])
                           )
table(lengths(transcripts_list))





# Build a data frame from the candidate genes -----------------------------

candidate_genes_vec <- c(unique(hand_picked_df[, "Gene symbol"]), first_trial_df[, 1])

candidates_CRISPRa_df <- GetGenes(candidate_genes_vec, CRISPRa_df)





# Save data ---------------------------------------------------------------

save(list = c("CRISPRa_df", "candidates_CRISPRa_df"),
     file = file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData")
     )


















