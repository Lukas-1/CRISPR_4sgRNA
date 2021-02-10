### 7th April 2020 ###



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
CRISPR_libraries_directory       <- file.path(CRISPR_input_directory, "CRISPR libraries")
CRISPRi_datasets_directory       <- file.path(CRISPR_libraries_directory, "CRISPRi")

RData_directory                  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory          <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory          <- file.path(RData_directory, "4) CRISPRi")

CRISPRi_Horlbeck2016_path        <- file.path(CRISPRi_datasets_directory, "Horlbeck, Kampmann, Weissman - eLife 2016")
CRISPRi_Horlbeck2016_sgRNAs_path <- file.path(CRISPRi_Horlbeck2016_path, "2016 - Compact and highly active next-generation libraries - Table S3.xlsx")
Horlbeck2016_TSSs_path           <- file.path(CRISPR_libraries_directory, "TSS", "Horlbeck, Kampmann, Weissman - eLife 2016",
                                              "2016 - Compact and highly active next-generation libraries - Table S2.xlsx"
                                              )
CRISPRi_Doench2018_path          <- file.path(CRISPRi_datasets_directory, "Sanson, Doench - Nat Comm 2018",
                                              "2018 - Optimized libraries for CRISPR-Cas9 genetic screens - Data S3.xlsx"
                                              )
GPP_CRISPRi_path                 <- file.path(CRISPR_root_directory, "4) Intermediate files/CRISPRi/GPP sgRNA designer/2) Output files")
GPP_CRISPRi_path                 <- file.path(GPP_CRISPRi_path, "All genes")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRi_RData_directory, "26) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData"))





# Read in data ------------------------------------------------------------

hCRISPRi_v2.0_df <- data.frame(read_excel(CRISPRi_Horlbeck2016_sgRNAs_path, sheet = "hCRISPRi-v2", skip = 8, col_names = FALSE), stringsAsFactors = FALSE)
names(hCRISPRi_v2.0_df) <- names(read_excel(CRISPRi_Horlbeck2016_sgRNAs_path, sheet = "hCRISPRi-v2", n_max = 1))

hCRISPRi_v2.1_df <- data.frame(read_excel(CRISPRi_Horlbeck2016_sgRNAs_path, sheet = "hCRISPRi-v2.1", skip = 2, col_names = FALSE), stringsAsFactors = FALSE)
names(hCRISPRi_v2.1_df) <- names(read_excel(CRISPRi_Horlbeck2016_sgRNAs_path, sheet = "hCRISPRi-v2.1", n_max = 1))

Horlbeck_TSS_df <- data.frame(read_excel(Horlbeck2016_TSSs_path), stringsAsFactors = FALSE, check.names = FALSE)

DolcettoA_df  <- data.frame(read_excel(CRISPRi_Doench2018_path, sheet = "SetA sgRNA annotations"), stringsAsFactors = FALSE, check.names = FALSE)
DolcettoB_df  <- data.frame(read_excel(CRISPRi_Doench2018_path, sheet = "SetB sgRNA annotations"), stringsAsFactors = FALSE, check.names = FALSE)

### Read in the output from the GPP sgRNA designer tool
GPP_CRISPRi_full_df <- ReadGPPOutputFiles(GPP_CRISPRi_path)




# Annotate the hCRISPRi_v2.0 data  ----------------------------------------

hCRISPRi_v2.0_df[["Sublibrary"]] <- hCRISPRa_v2_sublibrary_map[hCRISPRi_v2.0_df[["Sublibrary"]]]




# Annotate the hCRISPRi_v2.1 data -----------------------------------------

hCRISPRi_v2.1_df <- AddHorlbeckTSSSource(hCRISPRi_v2.1_df, Horlbeck_TSS_df)





# Combine the Dolcetto data -----------------------------------------------

Dolcetto_df <- rbind.data.frame(
  data.frame(DolcettoA_df, "Dolcetto_rank" = "1/2/3", stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(DolcettoB_df, "Dolcetto_rank" = "4/5/6", stringsAsFactors = FALSE, check.names = FALSE),
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)




# Process the data from the Broad Institute's GPP portal ------------------

GPP_CRISPRi_df <- TidyGPPOutputDf(GPP_CRISPRi_full_df, CRISPRa_GPP_output_columns)
GPP_CRISPRi_df <- FilterGPPOutputDf(GPP_CRISPRi_df, problematic_entrezs, n_unproblematic = 10, n_problematic = NULL)





# Identify control sgRNAs -------------------------------------------------

are_Dolcetto_controls <- ((Dolcetto_df[["Annotated Gene Symbol"]] == "CONTROL") &
                          (Dolcetto_df[["Annotated Gene ID"]] == "CONTROL")) |
                         ((Dolcetto_df[["Annotated Gene Symbol"]] == "NO-TARGET") &
                          (Dolcetto_df[["Annotated Gene ID"]] == "UNKNOWN_NO-TARGET"))

are_hCRISPRi_v2.0_controls <- hCRISPRi_v2.0_df[["gene"]] == "negative_control"






# Count the number of genes in each library -------------------------------

num_genes_in_library <- c(
  "Dolcetto"      = length(unique(Dolcetto_df[["Annotated Gene ID"]][!(are_Dolcetto_controls)])),
  "hCRISPRi-v2"   = length(unique(hCRISPRi_v2.1_df[["gene"]])),
  "GPP"           = length(unique(GPP_CRISPRi_df[["Target Gene ID"]])),
  "GPP_4_or_more" = sum(table(GPP_CRISPRi_df[["Target Gene ID"]]) >= 4)
)




# Harmonize the data frames -----------------------------------------------

mapped_Dolcetto_df      <- MapToEntrezs(Dolcetto_df[["Annotated Gene ID"]], Dolcetto_df[["Annotated Gene Symbol"]])
mapped_hCRISPRi_v2.1_df <- MapToEntrezs(symbols_vec = hCRISPRi_v2.1_df[["gene"]])
mapped_GPP_df           <- MapToEntrezs(as.character(GPP_CRISPRi_df[["Target Gene ID"]]), GPP_CRISPRi_df[["Target Gene Symbol"]])


new_Dolcetto_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "Dolcetto",
  "hCRISPRi_v2_transcript"    = NA_character_,
  "Is_control"                = ifelse(are_Dolcetto_controls, "Yes", "No"),
  mapped_Dolcetto_df[, 1:4],
  "Entrez_source_Dolcetto"    = mapped_Dolcetto_df[["Entrez_source"]],
  "Entrez_source_hCRISPRi_v2" = NA_integer_,
  "sgRNA_sequence"            = Dolcetto_df[["sgRNA Sequence"]],
  "Original_PAM"              = NA_character_,
  "Dolcetto_rank"             = Dolcetto_df[["Dolcetto_rank"]],
  "GPP_rank"                  = NA_integer_,
  "hCRISPRi_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "hCRISPRi_v2_ID"            = NA_character_,
  "hCRISPRi_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)

new_hCRISPRi_v2_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "hCRISPRi-v2.1",
  "hCRISPRi_v2_transcript"    = gsub(",", ", ", hCRISPRi_v2.1_df[["transcript"]], fixed = TRUE),
  "Is_control"                = "No",
  mapped_hCRISPRi_v2.1_df[, 1:4],
  "Entrez_source_Dolcetto"    = NA_integer_,
  "Entrez_source_hCRISPRi_v2" = mapped_hCRISPRi_v2.1_df[["Entrez_source"]],
  "sgRNA_sequence"            = hCRISPRi_v2.1_df[["protospacer sequence"]],
  "Original_PAM"              = NA_character_,
  "Dolcetto_rank"             = NA_integer_,
  "GPP_rank"                  = NA_integer_,
  "hCRISPRi_v2_rank"          = hCRISPRi_v2.1_df[["selection rank"]],
  "Predicted_score"           = hCRISPRi_v2.1_df[["predicted score"]],
  "Empirical_score"           = hCRISPRi_v2.1_df[["empirical score"]],
  "Off_target_stringency"     = hCRISPRi_v2.1_df[["off-target stringency"]],
  "Sublibrary"                = NA_character_,
  "hCRISPRi_v2_ID"            = sub(".23", "", hCRISPRi_v2.1_df[["sgID"]], fixed = TRUE),
  "hCRISPRi_TSS_source"       = hCRISPRi_v2.1_df[["TSS_source"]],
  stringsAsFactors            = FALSE
)


hCRISPRi_v2_controls_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "hCRISPRi-v2.0",
  "hCRISPRi_v2_transcript"    = NA_character_,
  "Is_control"                = "Yes",
  "Entrez_ID"                 = NA_character_,
  "Gene_symbol"               = NA_character_,
  "Original_entrez"           = NA_character_,
  "Original_symbol"           = NA_character_,
  "Entrez_source_Dolcetto"    = NA_integer_,
  "Entrez_source_hCRISPRi_v2" = NA_integer_,
  "sgRNA_sequence"            = hCRISPRi_v2.0_df[["protospacer sequence"]][are_hCRISPRi_v2.0_controls],
  "Original_PAM"              = NA_character_,
  "Dolcetto_rank"             = NA_integer_,
  "GPP_rank"                  = NA_integer_,
  "hCRISPRi_v2_rank"          = hCRISPRi_v2.0_df[["Sublibrary half"]][are_hCRISPRi_v2.0_controls],
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = hCRISPRi_v2.0_df[["off-target stringency"]][are_hCRISPRi_v2.0_controls],
  "Sublibrary"                = hCRISPRi_v2.0_df[["Sublibrary"]][are_hCRISPRi_v2.0_controls],
  "hCRISPRi_v2_ID"            = hCRISPRi_v2.0_df[["sgID"]][are_hCRISPRi_v2.0_controls],
  "hCRISPRi_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)


new_GPP_CRISPRi_df <- data.frame(
  "Combined_ID"               = NA_character_,
  "Source"                    = "GPP",
  "hCRISPRi_v2_transcript"    = NA_character_,
  "Is_control"                = "No",
  mapped_GPP_df[, 1:4],
  "Entrez_source_Dolcetto"    = NA_integer_,
  "Entrez_source_hCRISPRi_v2" = NA_integer_,
  "sgRNA_sequence"            = GPP_CRISPRi_df[["sgRNA Sequence"]],
  "Original_PAM"              = GPP_CRISPRi_df[["PAM Sequence"]],
  "Dolcetto_rank"             = NA_integer_,
  "GPP_rank"                  = GPP_CRISPRi_df[["Pick Order"]],
  "hCRISPRi_v2_rank"          = NA_character_,
  "Predicted_score"           = NA_real_,
  "Empirical_score"           = NA_real_,
  "Off_target_stringency"     = NA_real_,
  "Sublibrary"                = NA_character_,
  "hCRISPRi_v2_ID"            = NA_character_,
  "hCRISPRi_TSS_source"       = NA_character_,
  stringsAsFactors            = FALSE
)




# Check for multiple occurrences of the same sgRNA sequence ---------------

num_occurrences_hCRISPRi_v2 <- table(new_hCRISPRi_v2_df[["sgRNA_sequence"]])[new_hCRISPRi_v2_df[["sgRNA_sequence"]]]

multiple_occurrences_df <- data.frame(new_hCRISPRi_v2_df, "Num_occurrences" = as.integer(num_occurrences_hCRISPRi_v2))
multiple_occurrences_df <- multiple_occurrences_df[order(multiple_occurrences_df[["sgRNA_sequence"]]), ]
multiple_occurrences_df[multiple_occurrences_df[["Num_occurrences"]] > 1, ]

any(duplicated(new_Dolcetto_df[["sgRNA_sequence"]]))




# Build a combined data frame of all genes --------------------------------

any(is.na(new_Dolcetto_df[["Entrez_ID"]]) & (new_Dolcetto_df[["Original_symbol"]] == ""))

combined_df <- rbind.data.frame(new_Dolcetto_df,
                                new_hCRISPRi_v2_df,
                                hCRISPRi_v2_controls_df,
                                new_GPP_CRISPRi_df,
                                stringsAsFactors = FALSE,
                                make.row.names = FALSE
                                )
are_controls <- combined_df[["Is_control"]] == "Yes"

combined_df[["Original_symbol"]][are_controls] <- NA_character_

combined_df[["Combined_ID"]] <- ifelse(are_controls,
                                       "Control",
                                       ifelse(is.na(combined_df[["Entrez_ID"]]),
                                              toupper(combined_df[["Original_symbol"]]),
                                              combined_df[["Entrez_ID"]]
                                              )
                                       )

CRISPRi_df <- ResolveDuplicates(combined_df, concatenate_columns = c("Sublibrary", "hCRISPRi_v2_ID", "hCRISPRi_TSS_source"))




# Search for duplicated sgRNAs within the same gene -----------------------

sgID_vec <- paste0(CRISPRi_df[["Combined_ID"]], "__", toupper(CRISPRi_df[["sgRNA_sequence"]]))
num_occurrences <- table(sgID_vec)[sgID_vec]
have_duplications_sgRNA_list <- split(CRISPRi_df[num_occurrences > 1, ],
                                      factor(sgID_vec[num_occurrences > 1], levels = unique(sgID_vec[num_occurrences > 1]))
                                      )



# Count genes with multiple transcripts in the hCRISPRi-v2 database -------

are_hCRISPRi_v2 <- grepl("hCRISPRi-v2", CRISPRi_df[["Source"]], fixed = TRUE)
transcripts_list <- tapply(CRISPRi_df[["hCRISPRi_v2_transcript"]][are_hCRISPRi_v2],
                           factor(CRISPRi_df[["Combined_ID"]][are_hCRISPRi_v2],
                                  levels = unique(CRISPRi_df[["Combined_ID"]][are_hCRISPRi_v2])
                                  ),
                           function(x) unique(x[!(is.na(x))])
                           )
table(lengths(transcripts_list))




# Save data ---------------------------------------------------------------

save("CRISPRi_df",
     file = file.path(CRISPRi_RData_directory, "01) Compile predefined CRISPRi libraries - CRISPRi_df.RData")
     )

save(list = "num_genes_in_library",
     file = file.path(CRISPRi_RData_directory, "01) Compile predefined CRISPRa libraries - num_genes_in_library.RData")
     )

