### 29 October 2019 ###




# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory        <- "~/CRISPR_4sgRNA"
CRISPR_input_directory       <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory              <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory      <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory     <- file.path(RData_directory, "3) CRISPRko")

CRISPR_libraries_directory   <- file.path(CRISPR_input_directory, "CRISPR libraries")
CRISPRko_datasets_directory  <- file.path(CRISPR_libraries_directory, "CRISPRko")

CRISPRko_Brunello_path       <- file.path(CRISPRko_datasets_directory, "Brunello",
                                          "2016 - Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9 - Table S21 - Brunello.xlsx"
                                          )
CRISPRko_Brunello_2018_path  <- file.path(CRISPRko_datasets_directory, "Brunello",
                                          "2018 - Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities - Data S1.xlsx"
                                          )
CRISPRko_TKOv3_path          <- file.path(CRISPRko_datasets_directory, "TKOv3",  "tkov3_guide_sequence.xlsx")

GPP_CRISPRko_path            <- file.path(CRISPR_root_directory, "4) Intermediate files/CRISPRko/GPP sgRNA designer/2) Output files")
GPP_priority_CRISPRko_path   <- file.path(GPP_CRISPRko_path, "1) High-priority")
GPP_optional_CRISPRko_path   <- file.path(GPP_CRISPRko_path, "2) Optional")
GPP_semimanual_CRISPRko_path <- file.path(GPP_CRISPRko_path, "3) Semi-manual")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "18) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData"))






# Read in data ------------------------------------------------------------

Brunello_df <- data.frame(read_excel(CRISPRko_Brunello_path),
                          stringsAsFactors = FALSE, check.names = FALSE
                          )[, -1]

Brunello_2018_df <- data.frame(read_excel(CRISPRko_Brunello_2018_path, sheet = "sgRNA annotations"),
                               stringsAsFactors = FALSE, check.names = FALSE
                               )

TKOv3_df <- data.frame(read_excel(CRISPRko_TKOv3_path), stringsAsFactors = FALSE, check.names = FALSE)

GPP_CRISPRko_full_df <- ReadGPPOutputFiles(c(GPP_priority_CRISPRko_path, GPP_optional_CRISPRko_path, GPP_semimanual_CRISPRko_path))





# Find control sgRNA sequences --------------------------------------------

are_controls_Brunello <- startsWith(Brunello_2018_df[["Annotated Gene Symbol"]], "NO_CURRENT_")





# Explore the Brunello library --------------------------------------------

head(Brunello_2018_df)
head(Brunello_df)

identical(sort(Brunello_df[["sgRNA Target Sequence"]]),
          sort(Brunello_2018_df[["sgRNA Sequence"]][!(are_controls_Brunello)])
          )

all(Brunello_df[["sgRNA Target Sequence"]] %in% Brunello_2018_df[["sgRNA Sequence"]])
length(unique(Brunello_df[["sgRNA Target Sequence"]]))
length(unique(Brunello_2018_df[["sgRNA Sequence"]]))

head(Brunello_2018_df[!(Brunello_2018_df[["sgRNA Sequence"]] %in% Brunello_df[["sgRNA Target Sequence"]]), ])






# Explore the TKOv3 library -----------------------------------------------

TKOv3_df[sub("exon_", "", TKOv3_df[["TARGET EXON"]], fixed = TRUE) == "CTRL", ]





# Take the updated gene symbols from the 2018 version of Brunello ---------

stopifnot(!(anyDuplicated(Brunello_2018_df[["sgRNA Sequence"]])))

Brunello_matches <- match(toupper(Brunello_df[["sgRNA Target Sequence"]]), toupper(Brunello_2018_df[["sgRNA Sequence"]]))
Brunello_symbols_noNA <- Brunello_df[["Target Gene Symbol"]][!(is.na(Brunello_matches))]
Brunello_2018_symbols <- Brunello_2018_df[["Annotated Gene Symbol"]][Brunello_matches]
Brunello_2018_symbols_noNA <- Brunello_2018_symbols[!(is.na(Brunello_matches))]

are_identical_symbols <- mapply(identical, Brunello_symbols_noNA, Brunello_2018_symbols_noNA)
Brunello_symbols_df <- data.frame("Older" = Brunello_symbols_noNA,
                                  "Newer" = Brunello_2018_symbols_noNA,
                                  stringsAsFactors = FALSE,
                                  row.names = NULL
                                  )
Brunello_symbols_df[!(are_identical_symbols), ]

Brunello_entrezs_noNA <- as.character(Brunello_df[["Target Gene ID"]][!(is.na(Brunello_matches))])
Brunello_2018_entrezs_noNA <- Brunello_2018_df[["Annotated Gene ID"]][Brunello_matches[!(is.na(Brunello_matches))]]
are_identical_entrezs <- mapply(identical, Brunello_entrezs_noNA, Brunello_2018_entrezs_noNA)
stopifnot(all(are_identical_entrezs))





# Process the data from the Broad Institute's GPP portal ------------------

GPP_CRISPRko_df <- TidyGPPOutputDf(GPP_CRISPRko_full_df, CRISPRko_GPP_output_columns)
GPP_CRISPRko_df <- FilterGPPOutputDf(GPP_CRISPRko_df, problematic_entrezs, n_unproblematic = 10, n_problematic = 50)






# Tidy the data frame -----------------------------------------------------

Brunello_combined_symbols <- ifelse(is.na(Brunello_2018_symbols),
                                    Brunello_df[["Target Gene Symbol"]],
                                    Brunello_2018_symbols
                                    )

mapped_Brunello_df <- MapToEntrezs(entrez_IDs_vec = as.character(Brunello_df[["Target Gene ID"]]),
                                   symbols_vec = Brunello_combined_symbols
                                   )

mapped_Brunello_df[mapped_Brunello_df[["Original_entrez"]] != "", ]

mapped_TKOv3_df <- MapToEntrezs(symbols_vec = TKOv3_df[["GENE"]])

mapped_GPP_df <- MapToEntrezs(entrez_IDs_vec = as.character(GPP_CRISPRko_df[["Target Gene ID"]]),
                              symbols_vec    = GPP_CRISPRko_df[["Target Gene Symbol"]]
                              )


new_Brunello_df <- data.frame(
  "Combined_ID"            = ifelse(is.na(mapped_Brunello_df[["Entrez_ID"]]),
                                    toupper(mapped_Brunello_df[["Original_symbol"]]),
                                    mapped_Brunello_df[["Entrez_ID"]]
                                    ),
  "Source"                 = "Brunello",
  "Is_control"             = "No",
  mapped_Brunello_df[, 1:4],
  "Original_order"         = NA_integer_,
  "Entrez_source_Brunello" = mapped_Brunello_df[["Entrez_source"]],
  "Entrez_source_TKOv3"    = NA_integer_,
  "Transcript_ID"          = Brunello_df[["Target Transcript"]],
  "Genomic_sequence_ID"    = Brunello_df[["Genomic Sequence"]],
  "Exon_number_Brunello"   = as.integer(Brunello_df[["Exon Number"]]),
  "Exon_number_TKOv3"      = NA_integer_,
  "Exon_number_GPP"        = NA_integer_,
  "sgRNA_sequence"         = Brunello_df[["sgRNA Target Sequence"]],
  "sgRNA_context_sequence" = Brunello_df[["Target Context Sequence"]],
  "Original_PAM"           = Brunello_df[["PAM Sequence"]],
  "Original_cut_position"  = Brunello_df[["Position of Base After Cut (1-based)"]],
  "Original_orientation"   = Brunello_df[["Strand"]],
  "GPP_rank"               = NA_integer_,
  "Rule_set_2_score"       = Brunello_df[["Rule Set 2 score"]],
  "TKOv3_ID"               = NA_character_,
  stringsAsFactors         = FALSE
)

new_Brunello_df[["Original_order"]] <- unlist(tapply(seq_len(nrow(new_Brunello_df)),
                                                     factor(new_Brunello_df[["Combined_ID"]], levels = unique(new_Brunello_df[["Combined_ID"]])),
                                                     function(x) seq_len(length(x)),
                                                     simplify = FALSE
                                                     )
                                              )



Brunello_controls_df <- data.frame(
  "Combined_ID"            = paste0("Control_", sub("^NO_CURRENT_", "", Brunello_2018_df[["Annotated Gene Symbol"]][are_controls_Brunello])),
  "Source"                 = "Brunello",
  "Is_control"             = "Yes",
  "Entrez_ID"              = NA_character_,
  "Gene_symbol"            = NA_character_,
  "Original_entrez"        = NA_character_,
  "Original_symbol"        = NA_character_,
  "Original_order"         = NA_integer_,
  "Entrez_source_Brunello" = NA_integer_,
  "Entrez_source_TKOv3"    = NA_integer_,
  "Transcript_ID"          = NA_character_,
  "Genomic_sequence_ID"    = NA_character_,
  "Exon_number_Brunello"   = NA_integer_,
  "Exon_number_TKOv3"      = NA_integer_,
  "Exon_number_GPP"        = NA_integer_,
  "sgRNA_sequence"         = Brunello_2018_df[["sgRNA Sequence"]][are_controls_Brunello],
  "sgRNA_context_sequence" = NA_character_,
  "Original_PAM"           = NA_character_,
  "Original_cut_position"  = NA_integer_,
  "Original_orientation"   = NA_character_,
  "GPP_rank"               = NA_integer_,
  "Rule_set_2_score"       = NA_real_,
  "TKOv3_ID"               = NA_character_,
  stringsAsFactors         = FALSE
)


TKOv3_exons <- ifelse(TKOv3_df[["TARGET EXON"]] == "CTRL",
                      NA,
                      sub("exon_", "", TKOv3_df[["TARGET EXON"]], fixed = TRUE)
                      )
TKOv3_exons <- as.integer(TKOv3_exons)

are_TKOv3_controls <- TKOv3_df[["TARGET EXON"]] == "CTRL"
new_TKOv3_df <- data.frame(
  "Combined_ID"            = ifelse(are_TKOv3_controls, "Control", NA_character_),
  "Source"                 = "TKOv3",
  "Is_control"             = ifelse(are_TKOv3_controls, "Yes", "No"),
  mapped_TKOv3_df[, 1:4],
  "Original_order"         = NA_integer_,
  "Entrez_source_Brunello" = NA_integer_,
  "Entrez_source_TKOv3"    = mapped_TKOv3_df[["Entrez_source"]],
  "Transcript_ID"          = NA_character_,
  "Genomic_sequence_ID"    = NA_character_,
  "Exon_number_Brunello"   = NA_integer_,
  "Exon_number_TKOv3"      = TKOv3_exons,
  "Exon_number_GPP"        = NA_integer_,
  "sgRNA_sequence"         = TKOv3_df[["SEQUENCE"]],
  "sgRNA_context_sequence" = NA_character_,
  "Original_PAM"           = NA_character_,
  "Original_cut_position"  = NA_integer_,
  "Original_orientation"   = NA_character_,
  "GPP_rank"               = NA_integer_,
  "Rule_set_2_score"       = NA_real_,
  "TKOv3_ID"               = TKOv3_df[["GUIDE_ID"]],
  stringsAsFactors         = FALSE
)



new_GPP_CRISPRko_df <- data.frame(
  "Combined_ID"            = NA_character_,
  "Source"                 = "GPP",
  "Is_control"             = "No",
  mapped_GPP_df[, 1:4],
  "Original_order"         = NA_integer_,
  "Entrez_source_Brunello" = NA_integer_,
  "Entrez_source_TKOv3"    = NA_integer_,
  "Transcript_ID"          = NA_character_, # GPP_CRISPRko_df[["Target Transcript"]] # most of the transcript IDs are a different version from Brunello
  "Genomic_sequence_ID"    = NA_character_,
  "Exon_number_Brunello"   = NA_integer_,
  "Exon_number_TKOv3"      = NA_integer_,
  "Exon_number_GPP"        = GPP_CRISPRko_df[["Exon Number"]],
  "sgRNA_sequence"         = GPP_CRISPRko_df[["sgRNA Sequence"]],
  "sgRNA_context_sequence" = GPP_CRISPRko_df[["sgRNA Context Sequence"]],
  "Original_PAM"           = GPP_CRISPRko_df[["PAM Sequence"]],
  "Original_cut_position"  = NA_integer_, # GPP_CRISPRko_df[["sgRNA Cut Position (1-based)"]] # the cut locations are based on hg38, not hg19 as in Brunello
  "Original_orientation"   = GPP_CRISPRko_df[["Orientation"]],
  "GPP_rank"               = GPP_CRISPRko_df[["Pick Order"]],
  "Rule_set_2_score"       = NA_real_, # GPP uses a different version of Rule Set 2... GPP_CRISPRko_df[["On-Target Efficacy Score"]],
  "TKOv3_ID"               = NA_character_,
  stringsAsFactors         = FALSE
)





# Count the number of genes in each library -------------------------------

num_genes_in_library <- c(
  "Brunello"      = length(unique(Brunello_df[["Target Gene ID"]])),
  "TKOv3"         = length(unique(TKOv3_df[["GENE"]][!(are_TKOv3_controls)])),
  "GPP"           = length(unique(GPP_CRISPRko_df[["Target Gene ID"]])),
  "GPP_4_or_more" = sum(table(GPP_CRISPRko_df[["Target Gene ID"]]) >= 4)
)




# Build a combined data frame ---------------------------------------------

combined_df <- rbind.data.frame(new_Brunello_df,
                                Brunello_controls_df,
                                new_TKOv3_df,
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

CRISPRko_df <- ResolveDuplicates(combined_df, concatenate_columns = c("TKOv3_ID", "Exon_number_GPP"))





# Make sure that sgRNAs for the same gene are on consecutive rows ---------

new_CRISPRko_df <- CRISPRko_df[order(match(CRISPRko_df[["Combined_ID"]], CRISPRko_df[["Combined_ID"]])), ]
row.names(CRISPRko_df) <- NULL





# Check for multiple occurrences of the same sgRNA sequence ---------------

anyDuplicated(CRISPRko_df[["sgRNA_sequence"]])





# Save data ---------------------------------------------------------------

save(list = "CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData")
     )

save(list = "num_genes_in_library",
     file = file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries - num_genes_in_library.RData")
     )







