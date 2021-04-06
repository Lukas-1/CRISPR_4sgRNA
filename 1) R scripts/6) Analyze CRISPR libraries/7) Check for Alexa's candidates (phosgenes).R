### 20 May 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Analysis", "Alexa")
gene_lists_directory     <- file.path(CRISPR_root_directory, "2) Input data", "Gene lists")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRa_by_gene_df <- full_4sg_by_gene_df
load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRo_by_gene_df <- full_4sg_by_gene_df
rm(full_4sg_by_gene_df)




# Read in data ------------------------------------------------------------

phosgenes_path <- file.path(gene_lists_directory, "crispr_phosgenes_of_interest.txt")
phosgenes_df <- read.table(phosgenes_path, stringsAsFactors = FALSE)




# Check for the availability of Alexa's phosgenes -------------------------

are_present_CRISPRa  <- phosgenes_df[[1]] %in% CRISPRa_by_gene_df[["Gene_symbol"]]
are_present_CRISPRko <- phosgenes_df[[1]] %in% CRISPRo_by_gene_df[["Gene_symbol"]]

stopifnot(identical(are_present_CRISPRa, are_present_CRISPRko))

CRISPRa_matches <- match(phosgenes_df[[1]], CRISPRa_by_gene_df[["Gene_symbol"]])
CRISPRo_matches <- match(phosgenes_df[[1]], CRISPRo_by_gene_df[["Gene_symbol"]])

CRISPRa_libraries <- CRISPRa_by_gene_df[["Sublibrary_4sg"]][CRISPRa_matches]
CRISPRo_libraries <- CRISPRo_by_gene_df[["Sublibrary_4sg"]][CRISPRo_matches]

CRISPRa_plates <- CRISPRa_by_gene_df[["Plate_number"]][CRISPRa_matches]
CRISPRo_plates <- CRISPRo_by_gene_df[["Plate_number"]][CRISPRo_matches]

table(CRISPRa_plates)
table(CRISPRo_plates)



# Select only Alexa's phosgenes -------------------------------------------

ReorderByLibrary <- function(input_df) {
  new_order <- order(input_df[["Sublibrary_4sg"]] != "Kinases/Phosphatases/Drug Targets",
                     match(input_df[["Sublibrary_4sg"]], names(sublibraries_all_entrezs_list)),
                     as.integer(input_df[["Entrez_ID"]])
                     )
  input_df <- input_df[new_order, ]
  row.names(input_df) <- NULL
  are_all_NA <- rowSums(is.na(as.matrix(input_df))) == 3
  input_df[["Gene_symbol"]][are_all_NA] <- phosgenes_df[[1]][!(are_present_CRISPRa)]
  input_df[["Entrez_ID"]][are_all_NA] <- SymbolsToEntrezIDs(phosgenes_df[[1]][!(are_present_CRISPRa)])
  return(input_df)
}

phosgenes_CRISPRa_df <- ReorderByLibrary(CRISPRa_by_gene_df[CRISPRa_matches, c("Entrez_ID", "Gene_symbol", "Sublibrary_4sg")])
phosgenes_CRISPRo_df <- ReorderByLibrary(CRISPRo_by_gene_df[CRISPRo_matches, c("Entrez_ID", "Gene_symbol", "Sublibrary_4sg")])
stopifnot(identical(phosgenes_CRISPRa_df, phosgenes_CRISPRo_df))




# Export phosgenes --------------------------------------------------------

write.table(table(CRISPRa_plates),
            file = file.path(file_output_directory, "CRISPRa_plates_overview.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )

write.table(table(CRISPRo_plates),
            file = file.path(file_output_directory, "CRISPRo_plates_overview.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )

write.table(table(CRISPRa_libraries),
            file = file.path(file_output_directory, "CRISPRa_libraries_overview.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )

write.table(table(CRISPRo_libraries),
            file = file.path(file_output_directory, "CRISPRo_libraries_overview.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )

write.table(phosgenes_CRISPRa_df,
            file = file.path(file_output_directory, "Genes_and_libraries_4sg.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )















