### 20 February 2021 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory        <- "~/CRISPR"
CRISPR_input_directory       <- file.path(CRISPR_root_directory, "2) Input data")
CRISPR_libraries_directory   <- file.path(CRISPR_input_directory, "Mouse CRISPR libraries")
CRISPRko_datasets_directory  <- file.path(CRISPR_libraries_directory, "CRISPRko")
CRISPRko_Brie_path           <- file.path(CRISPRko_datasets_directory, "Brie",
                                          "2016 - Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9 - Table S22 - Brie.xlsx"
                                          )
CRISPRko_Brie_controls_path  <- file.path(CRISPRko_datasets_directory, "Brie",
                                          "broadgpp-brie-library-controls.csv"
                                          )
RData_directory              <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory      <- file.path(RData_directory, "6) Mouse - General")
CRISPRko_RData_directory     <- file.path(RData_directory, "8) Mouse - CRISPRko")

GPP_CRISPRko_path            <- file.path(CRISPR_root_directory, "4) Intermediate files",
                                          "Mouse - CRISPRko", "GPP sgRNA designer",
                                          "2) Output files",
                                          "All genes"
                                          )




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "18) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData"))






# Read in data ------------------------------------------------------------

Brie_df <- data.frame(read_excel(CRISPRko_Brie_path),
                      stringsAsFactors = FALSE, check.names = FALSE
                      )

Brie_controls_df <- read.csv(CRISPRko_Brie_controls_path,
                             stringsAsFactors = FALSE, check.names = FALSE
                             )

GPP_CRISPRko_full_df <- ReadGPPOutputFiles(GPP_CRISPRko_path)






# Process the data from the Broad Institute's GPP portal ------------------

GPP_CRISPRko_df <- TidyGPPOutputDf(GPP_CRISPRko_full_df, CRISPRko_GPP_output_columns)
GPP_CRISPRko_df <- FilterGPPOutputDf(GPP_CRISPRko_df, problematic_entrezs, n_unproblematic = 10, n_problematic = 50)
# GPP_CRISPRko_df <- FilterGPPOutputDf(GPP_CRISPRko_df, c(), n_unproblematic = 50, n_problematic = 50)





# Tidy the data frame -----------------------------------------------------

mapped_Brie_df <- MapToEntrezs(entrez_IDs_vec = as.character(Brie_df[["Target Gene ID"]]),
                               symbols_vec    = Brie_df[["Target Gene Symbol"]],
                               is_mouse       = TRUE
                               )

mapped_GPP_df <- MapToEntrezs(entrez_IDs_vec = as.character(GPP_CRISPRko_df[["Target Gene ID"]]),
                              symbols_vec    = GPP_CRISPRko_df[["Target Gene Symbol"]],
                              is_mouse       = TRUE
                              )


new_Brie_df <- data.frame(
  "Combined_ID"            = ifelse(is.na(mapped_Brie_df[["Entrez_ID"]]),
                                    mapped_Brie_df[["Original_symbol"]],
                                    mapped_Brie_df[["Entrez_ID"]]
                                    ),
  "Source"                 = "Brie",
  "Is_control"             = "No",
  mapped_Brie_df[, 1:4],
  "Original_order"         = NA_integer_,
  "Entrez_source_Brie"     = mapped_Brie_df[["Entrez_source"]],
  "Transcript_ID"          = Brie_df[["Target Transcript"]],
  "Genomic_sequence_ID"    = Brie_df[["Genomic Sequence"]],
  "Exon_number_Brie"       = as.integer(Brie_df[["Exon Number"]]),
  "Exon_number_GPP"        = NA_integer_,
  "sgRNA_sequence"         = Brie_df[["sgRNA Target Sequence"]],
  "sgRNA_context_sequence" = Brie_df[["Target Context Sequence"]],
  "Original_PAM"           = Brie_df[["PAM Sequence"]],
  "Original_cut_position"  = as.integer(Brie_df[["Position of Base After Cut (1-based)"]]),
  "Original_orientation"   = Brie_df[["Strand"]],
  "GPP_rank"               = NA_integer_,
  "Rule_set_2_score"       = Brie_df[["Rule Set 2 score"]],
  stringsAsFactors         = FALSE
)


Brie_controls_df <- data.frame(
  "Combined_ID"            = paste0("Control_", Brie_controls_df[["Public ID"]]),
  "Source"                 = "Brie",
  "Is_control"             = "Yes",
  "Entrez_ID"              = NA_character_,
  "Gene_symbol"            = NA_character_,
  "Original_entrez"        = NA_character_,
  "Original_symbol"        = NA_character_,
  "Original_order"         = NA_integer_,
  "Entrez_source_Brie"     = NA_integer_,
  "Transcript_ID"          = NA_character_,
  "Genomic_sequence_ID"    = NA_character_,
  "Exon_number_Brie"       = NA_integer_,
  "Exon_number_GPP"        = NA_integer_,
  "sgRNA_sequence"         = Brie_controls_df[["Target Sequence"]],
  "sgRNA_context_sequence" = NA_character_,
  "Original_PAM"           = NA_character_,
  "Original_cut_position"  = NA_integer_,
  "Original_orientation"   = NA_character_,
  "GPP_rank"               = NA_integer_,
  "Rule_set_2_score"       = NA_real_,
  stringsAsFactors         = FALSE
)


new_GPP_CRISPRko_df <- data.frame(
  "Combined_ID"            = NA_character_,
  "Source"                 = "GPP",
  "Is_control"             = "No",
  mapped_GPP_df[, 1:4],
  "Original_order"         = NA_integer_,
  "Entrez_source_Brie"     = NA_integer_,
  "Transcript_ID"          = NA_character_, # GPP_CRISPRko_df[["Target Transcript"]]
  "Genomic_sequence_ID"    = NA_character_,
  "Exon_number_Brie"       = NA_integer_,
  "Exon_number_GPP"        = GPP_CRISPRko_df[["Exon Number"]],
  "sgRNA_sequence"         = GPP_CRISPRko_df[["sgRNA Sequence"]],
  "sgRNA_context_sequence" = GPP_CRISPRko_df[["sgRNA Context Sequence"]],
  "Original_PAM"           = GPP_CRISPRko_df[["PAM Sequence"]],
  "Original_cut_position"  = GPP_CRISPRko_df[["sgRNA Cut Position (1-based)"]],
  "Original_orientation"   = GPP_CRISPRko_df[["Orientation"]],
  "GPP_rank"               = GPP_CRISPRko_df[["Pick Order"]],
  "Rule_set_2_score"       = NA_real_, # GPP uses a different version of Rule Set 2... GPP_CRISPRko_df[["On-Target Efficacy Score"]]
  stringsAsFactors         = FALSE
)





# Count the number of genes in each library -------------------------------

num_genes_in_library <- c(
  "Brie"          = length(unique(Brie_df[["Target Gene ID"]])),
  "GPP"           = length(unique(GPP_CRISPRko_df[["Target Gene ID"]])),
  "GPP_4_or_more" = sum(table(GPP_CRISPRko_df[["Target Gene ID"]]) >= 4)
)





# Build a combined data frame ---------------------------------------------

combined_df <- rbind.data.frame(new_Brie_df,
                                Brie_controls_df,
                                new_GPP_CRISPRko_df,
                                stringsAsFactors = FALSE,
                                make.row.names = FALSE
                                )

combined_df[["Combined_ID"]] <- ifelse(!(is.na(combined_df[["Combined_ID"]])),
                                       combined_df[["Combined_ID"]],
                                       ifelse(is.na(combined_df[["Entrez_ID"]]),
                                              toupper(combined_df[["Original_symbol"]]),
                                              combined_df[["Entrez_ID"]]
                                              )
                                       )

CRISPRko_df <- ResolveDuplicates(combined_df, concatenate_columns = c("Exon_number_GPP"))





# Make sure that sgRNAs for the same gene are on consecutive rows ---------

new_CRISPRko_df <- CRISPRko_df[order(match(CRISPRko_df[["Combined_ID"]], CRISPRko_df[["Combined_ID"]])), ]
row.names(CRISPRko_df) <- NULL





# Check for multiple occurrences of the same sgRNA sequence ---------------

any(duplicated(CRISPRko_df[["sgRNA_sequence"]]))





# Save data ---------------------------------------------------------------

save(list = "CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData")
     )

save(list = "num_genes_in_library",
     file = file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries - num_genes_in_library.RData")
     )



