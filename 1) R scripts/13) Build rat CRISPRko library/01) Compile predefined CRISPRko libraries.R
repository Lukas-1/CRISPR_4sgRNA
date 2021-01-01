### 29 October 2019 ###




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
RData_directory              <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory      <- file.path(RData_directory, "10) Rat - General")
CRISPRko_RData_directory     <- file.path(RData_directory, "12) Rat - CRISPRko")

mouse_CRISPRa_library_path   <- file.path(RData_directory, "7) Mouse - CRISPRa",
                                          "19) For problematic genes, pick 4 guides without reference to the TSS.RData"
                                          )
GPP_CRISPRko_path            <- file.path(CRISPR_root_directory, "4) Intermediate files",
                                          "Rat - CRISPRko", "GPP sgRNA designer",
                                          "2) Output files",
                                          "Pharmacology genes"
                                          )




# Load data ---------------------------------------------------------------

load(mouse_CRISPRa_library_path)
load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Rn.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Read in data ------------------------------------------------------------

GPP_CRISPRko_full_df <- ReadGPPOutputFiles(GPP_CRISPRko_path)






# Pick control sgRNA sequences (from the mouse mCRISPRa-v2 library) -------

mouse_controls_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Is_control"]] == "Yes", ]




# Process the data from the Broad Institute's GPP portal ------------------

GPP_CRISPRko_df <- TidyGPPOutputDf(GPP_CRISPRko_full_df, CRISPRko_GPP_output_columns)
GPP_CRISPRko_df <- FilterGPPOutputDf(GPP_CRISPRko_df, c(), n_unproblematic = 50)





# Tidy the data frame -----------------------------------------------------

mapped_GPP_df <- MapToEntrezs(entrez_IDs_vec = as.character(GPP_CRISPRko_df[["Target Gene ID"]]),
                              symbols_vec    = GPP_CRISPRko_df[["Target Gene Symbol"]],
                              is_rat         = TRUE
                              )

new_GPP_CRISPRko_df <- data.frame(
  "Combined_ID"            = NA_character_,
  "Source"                 = "GPP",
  "Is_control"             = "No",
  mapped_GPP_df[, 1:4],
  "Transcript_ID"          = GPP_CRISPRko_df[["Target Transcript"]],
  "Exon_number_GPP"        = GPP_CRISPRko_df[["Exon Number"]],
  "sgRNA_sequence"         = GPP_CRISPRko_df[["sgRNA Sequence"]],
  "sgRNA_context_sequence" = GPP_CRISPRko_df[["sgRNA Context Sequence"]],
  "Original_PAM"           = GPP_CRISPRko_df[["PAM Sequence"]],
  "Original_cut_position"  = GPP_CRISPRko_df[["sgRNA Cut Position (1-based)"]],
  "Original_orientation"   = GPP_CRISPRko_df[["Orientation"]],
  "GPP_rank"               = GPP_CRISPRko_df[["Pick Order"]],
  "Rule_set_2_score"       = GPP_CRISPRko_df[["On-Target Efficacy Score"]],
  stringsAsFactors         = FALSE
)







# Build a combined data frame ---------------------------------------------

for (column_name in setdiff(names(new_GPP_CRISPRko_df), names(mouse_controls_df))) {
  mouse_controls_df[[column_name]] <- NA
}

combined_df <- rbind.data.frame(new_GPP_CRISPRko_df,
                                mouse_controls_df[, names(new_GPP_CRISPRko_df)],
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






