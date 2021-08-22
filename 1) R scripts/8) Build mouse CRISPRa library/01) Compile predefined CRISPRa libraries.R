### 21st July 2019 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))
source(file.path(general_functions_directory, "24) Assigning genes to sublibraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
CRISPR_input_directory           <- file.path(CRISPR_root_directory, "2) Input data")
CRISPR_libraries_directory       <- file.path(CRISPR_input_directory, "Mouse CRISPR libraries")
CRISPRa_datasets_directory       <- file.path(CRISPR_libraries_directory, "CRISPRa")

gene_lists_directory             <- file.path(CRISPR_input_directory, "Gene lists")
RData_directory                  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory          <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory          <- file.path(RData_directory, "7) Mouse - CRISPRa")

CRISPRa_Horlbeck2016_path        <- file.path(CRISPRa_datasets_directory, "Horlbeck, Kampmann, Weissman - eLife 2016")
CRISPRa_Horlbeck2016_sgRNAs_path <- file.path(CRISPRa_Horlbeck2016_path, "2016 - Compact and highly active next-generation libraries - Table S6.xlsx")
Horlbeck2016_TSSs_path           <- file.path(CRISPR_input_directory, "CRISPR libraries",
                                              "TSS", "Horlbeck, Kampmann, Weissman - eLife 2016",
                                              "2016 - Compact and highly active next-generation libraries - Table S2.xlsx"
                                              )
CRISPRa_Doench2018_path          <- file.path(CRISPRa_datasets_directory, "Sanson, Doench - Nat Comm 2018")

GPP_CRISPRa_path                 <- file.path(CRISPR_root_directory, "4) Intermediate files/Mouse - CRISPRa/GPP sgRNA designer/2) Output files/All genes")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRa_RData_directory, "26) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData"))





# Read in data ------------------------------------------------------------

mCRISPRa_v2_df <- data.frame(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, skip = 7)[-1, ], stringsAsFactors = FALSE)
names(mCRISPRa_v2_df) <- names(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, n_max = 1))

Horlbeck_TSS_df <- data.frame(read_excel(Horlbeck2016_TSSs_path, sheet = "Mouse mm10"), stringsAsFactors = FALSE, check.names = FALSE)

# Downloaded from https://www.addgene.org/pooled-library/broadgpp-mouse-crispra-caprano-p65hsf/ on 6 October 2020
CapranoA_df  <- read.table(file.path(CRISPRa_Doench2018_path, "broadgpp-caprano-targets-seta.txt"),
                           sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE
                           )
CapranoB_df  <- read.table(file.path(CRISPRa_Doench2018_path, "broadgpp-caprano-targets-setb.txt"),
                           sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE
                           )

### Read in the output from the GPP sgRNA designer tool
GPP_CRISPRa_full_df <- ReadGPPOutputFiles(GPP_CRISPRa_path)





# Annotate the mCRISPRa_v2 data -------------------------------------------

mCRISPRa_v2_sublibrary_map <- hCRISPRa_v2_sublibrary_map
names(mCRISPRa_v2_sublibrary_map) <- sub("h", "m", names(mCRISPRa_v2_sublibrary_map), fixed = TRUE)

mCRISPRa_v2_df[["Sublibrary"]] <- mCRISPRa_v2_sublibrary_map[mCRISPRa_v2_df[["Sublibrary"]]]

mCRISPRa_v2_df <- AddHorlbeckTSSSource(mCRISPRa_v2_df, Horlbeck_TSS_df)





# Combine the Caprano data ------------------------------------------------

Caprano_df <- rbind.data.frame(
  data.frame(CapranoA_df, "Caprano_rank" = "1/2/3", stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(CapranoB_df, "Caprano_rank" = "4/5/6", stringsAsFactors = FALSE, check.names = FALSE),
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)




# Process the data from the Broad Institute's GPP portal ------------------

GPP_CRISPRa_df <- TidyGPPOutputDf(GPP_CRISPRa_full_df, CRISPRa_GPP_output_columns)
GPP_CRISPRa_df <- FilterGPPOutputDf(GPP_CRISPRa_df,
                                    problematic_genes = problematic_entrezs, # unique(GPP_CRISPRa_df[["Target Gene ID"]]),
                                    n_unproblematic = 10,
                                    n_problematic = NULL
                                    )




# Collect all unique HUGO gene symbols ------------------------------------

are_mCRISPRa_v2_controls <- mCRISPRa_v2_df[["gene"]] == "negative_control"

mCRISPRa_v2_symbols_vec <- unique(mCRISPRa_v2_df[["gene"]][!(are_mCRISPRa_v2_controls)])
Caprano_symbols_entrezs_df <- unique(Caprano_df[, c("Annotated Gene Symbol", "Annotated Gene ID")], MARGIN = 1)

names(Caprano_symbols_entrezs_df) <- c("Gene_symbol", "Entrez_ID")

unique_symbols_vec <- unique(c(mCRISPRa_v2_symbols_vec,
                               Caprano_symbols_entrezs_df[["Gene_symbol"]][!(toupper(Caprano_symbols_entrezs_df[["Gene_symbol"]]) %in% toupper(mCRISPRa_v2_symbols_vec))]
                               )
                             )





# Check for discrepancies in the Entrez IDs -------------------------------

### The following gene symbols are associated with multiple Entrez IDs in the Caprano dataset
Caprano_symbols_entrezs_df[table(Caprano_symbols_entrezs_df[["Gene_symbol"]])[Caprano_symbols_entrezs_df[["Gene_symbol"]]] == 2, ]

Caprano_symbols_to_entrezs <- SymbolsToEntrezIDs(Caprano_symbols_entrezs_df[["Gene_symbol"]])
Caprano_entrezs_to_symbols <- EntrezIDsToSymbols(Caprano_symbols_entrezs_df[["Entrez_ID"]])

Caprano_symbols_entrezs_df[is.na(Caprano_symbols_to_entrezs), ] ### The following gene symbols could not be mapped to Entrez IDs
Caprano_symbols_entrezs_df[is.na(Caprano_entrezs_to_symbols), ] ### The following Entrez IDs could not be translated to gene symbols

### The following genes have discrepant Entrez IDs (those that were translated from gene symbols, versus those provided by the authors)
are_different_entrezs <- !(is.na(Caprano_symbols_to_entrezs)) & (Caprano_symbols_to_entrezs != Caprano_symbols_entrezs_df[["Entrez_ID"]])
data.frame(Caprano_symbols_entrezs_df, "Symbol_to_Entrez" = Caprano_symbols_to_entrezs, stringsAsFactors = FALSE)[are_different_entrezs, ]

### The following genes have discrepant gene symbols (those that were translated from Entrez IDs, versus those provided by the authors)
are_different_symbols <- !(is.na(Caprano_entrezs_to_symbols)) & (Caprano_entrezs_to_symbols != Caprano_symbols_entrezs_df[["Gene_symbol"]])
data.frame(Caprano_symbols_entrezs_df, "Entrez_to_Symbol" = Caprano_entrezs_to_symbols, stringsAsFactors = FALSE)[are_different_symbols, ]





# Harmonize the data frames -----------------------------------------------

mapped_Caprano_df     <- MapToEntrezs(Caprano_df[["Annotated Gene ID"]], Caprano_df[["Annotated Gene Symbol"]], is_mouse = TRUE)
mapped_mCRISPRa_v2_df <- MapToEntrezs(symbols_vec = mCRISPRa_v2_df[["gene"]], is_mouse = TRUE)
mapped_GPP_df         <- MapToEntrezs(as.character(GPP_CRISPRa_df[["Target Gene ID"]]), GPP_CRISPRa_df[["Target Gene Symbol"]], is_mouse = TRUE)


new_Caprano_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "Caprano",
  "mCRISPRa_v2_transcript"    = NA_character_,
  "Is_control"                = "No",
  mapped_Caprano_df[, 1:4],
  "Entrez_source_Caprano"     = mapped_Caprano_df[["Entrez_source"]],
  "Entrez_source_mCRISPRa_v2" = NA_integer_,
  "sgRNA_sequence"            = Caprano_df[["Barcode Sequence"]],
  "Original_PAM"              = NA_character_,
  "Caprano_rank"              = Caprano_df[["Caprano_rank"]],
  "GPP_rank"                  = NA_integer_,
  "mCRISPRa_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "mCRISPRa_v2_ID"            = NA_character_,
  "mCRISPRa_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)

new_mCRISPRa_v2_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "mCRISPRa-v2",
  "mCRISPRa_v2_transcript"    = ifelse(mCRISPRa_v2_df[["transcript"]] == "na",
                                       NA_character_,
                                       gsub(",", ", ", mCRISPRa_v2_df[["transcript"]], fixed = TRUE)
                                       ),
  "Is_control"                = ifelse(are_mCRISPRa_v2_controls, "Yes", "No"),
  mapped_mCRISPRa_v2_df[, 1:4],
  "Entrez_source_Caprano"     = NA_integer_,
  "Entrez_source_mCRISPRa_v2" = mapped_mCRISPRa_v2_df[["Entrez_source"]],
  "sgRNA_sequence"            = mCRISPRa_v2_df[["protospacer sequence"]],
  "Original_PAM"              = NA_character_,
  "Caprano_rank"              = NA_integer_,
  "GPP_rank"                  = NA_integer_,
  "mCRISPRa_v2_rank"          = ifelse(is.na(mCRISPRa_v2_df[["selection rank"]]), mCRISPRa_v2_df[["Sublibrary half"]], mCRISPRa_v2_df[["selection rank"]]),
  "Predicted_score"           = mCRISPRa_v2_df[["predicted score"]],
  "Empirical_score"           = mCRISPRa_v2_df[["empirical score"]],
  "Off_target_stringency"     = mCRISPRa_v2_df[["off-target stringency"]],
  "Sublibrary"                = mCRISPRa_v2_df[["Sublibrary"]],
  "mCRISPRa_v2_ID"            = sub(".23", "", mCRISPRa_v2_df[["sgID"]], fixed = TRUE),
  "mCRISPRa_TSS_source"       = mCRISPRa_v2_df[["TSS_source"]],
  stringsAsFactors            = FALSE
)


new_GPP_CRISPRa_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "GPP",
  "mCRISPRa_v2_transcript"    = NA_character_,
  "Is_control"                = "No",
  mapped_GPP_df[, 1:4],
  "Entrez_source_Caprano"     = NA_integer_,
  "Entrez_source_mCRISPRa_v2" = NA_integer_,
  "sgRNA_sequence"            = GPP_CRISPRa_df[["sgRNA Sequence"]],
  "Original_PAM"              = GPP_CRISPRa_df[["PAM Sequence"]],
  "Caprano_rank"              = NA_integer_,
  "GPP_rank"                  = GPP_CRISPRa_df[["Pick Order"]],
  "mCRISPRa_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "mCRISPRa_v2_ID"            = NA_character_,
  "mCRISPRa_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)






# Check for multiple occurrences of the same sgRNA sequence ---------------

num_occurrences_mCRISPRa_v2 <- table(new_mCRISPRa_v2_df[["sgRNA_sequence"]])[new_mCRISPRa_v2_df[["sgRNA_sequence"]]]

multiple_occurrences_df <- data.frame(new_mCRISPRa_v2_df, "Num_occurrences" = as.integer(num_occurrences_mCRISPRa_v2))
multiple_occurrences_df <- multiple_occurrences_df[order(multiple_occurrences_df[["sgRNA_sequence"]]), ]
multiple_occurrences_df[multiple_occurrences_df[["Num_occurrences"]] > 1, ]

any(duplicated(new_Caprano_df[["sgRNA_sequence"]]))




# Build a combined data frame of all genes --------------------------------

any(is.na(new_Caprano_df[["Entrez_ID"]]) & (new_Caprano_df[["Original_symbol"]] == ""))

combined_df <- rbind.data.frame(new_Caprano_df, new_mCRISPRa_v2_df, new_GPP_CRISPRa_df, stringsAsFactors = FALSE, make.row.names = FALSE)
are_controls <- combined_df[["Is_control"]] == "Yes"

control_names_map <- c("CONTROL" = "Control", "NEGATIVE_CONTROL" = "Control", "NO-TARGET" = "No-target")
combined_df[["Original_symbol"]][are_controls] <- control_names_map[combined_df[["Original_symbol"]][are_controls]]

combined_df[["Combined_ID"]] <- ifelse(are_controls,
                                       "Control",
                                       ifelse(is.na(combined_df[["Entrez_ID"]]),
                                              toupper(combined_df[["Original_symbol"]]),
                                              combined_df[["Entrez_ID"]]
                                              )
                                       )

CRISPRa_df <- ResolveDuplicates(combined_df, concatenate_columns = c("Sublibrary", "mCRISPRa_v2_ID", "mCRISPRa_TSS_source"))





# Search for duplicated sgRNAs within the same gene -----------------------

sgID_vec <- paste0(CRISPRa_df[["Combined_ID"]], "__", toupper(CRISPRa_df[["sgRNA_sequence"]]))
num_occurrences <- table(sgID_vec)[sgID_vec]
have_duplications_sgRNA_list <- split(CRISPRa_df[num_occurrences > 1, ],
                                      factor(sgID_vec[num_occurrences > 1], levels = unique(sgID_vec[num_occurrences > 1]))
                                      )





# Count genes with multiple transcripts in the mCRISPRa-v2 database -------

are_mCRISPRa_v2 <- grepl("mCRISPRa-v2", CRISPRa_df[["Source"]], fixed = TRUE)
transcripts_list <- tapply(CRISPRa_df[["mCRISPRa_v2_transcript"]][are_mCRISPRa_v2],
                           factor(CRISPRa_df[["Combined_ID"]][are_mCRISPRa_v2], levels = unique(CRISPRa_df[["Combined_ID"]][are_mCRISPRa_v2])),
                           function(x) unique(x[!(is.na(x))])
                           )
table(lengths(transcripts_list))






# Save data ---------------------------------------------------------------

save(list = "CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData")
     )



