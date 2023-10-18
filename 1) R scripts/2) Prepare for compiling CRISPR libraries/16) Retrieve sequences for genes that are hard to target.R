### 18th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Define functions --------------------------------------------------------

GetCodingSequence <- function(entrez_ID, get_exons = FALSE) {
  if (get_exons) {
    UseFunction <- GenomicFeatures::exons
  } else {
    UseFunction <- GenomicFeatures::cds
  }
  cds_GRanges_object <- UseFunction(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    filter = list("gene_id" = entrez_ID)
                                    )
  results_df <- as.data.frame(cds_GRanges_object)
  colnames(results_df) <- c("Chromosome", "Start", "End", "Length", "Strand")
  results_df <- data.frame(
    "Entrez_ID" = entrez_ID,
    "Symbol"    = EntrezIDsToSymbols(entrez_ID),
    results_df[, c("Chromosome", "Strand", "Start", "End", "Length")],
    stringsAsFactors = FALSE
  )
  return(results_df)
}






# Find the coding sequences for select genes ------------------------------

GetCodingSequence("5350")
GetCodingSequence("100507027")

GetCodingSequence("5350")[["Start"]] - 50
GetCodingSequence("5350")[["End"]] + 50


GetCodingSequence(SymbolsToEntrezIDs("HNRNPK"))
GetCodingSequence(SymbolsToEntrezIDs("HNRNPK"), get_exons = TRUE)

GetCodingSequence(SymbolsToEntrezIDs("PRNP"))
GetCodingSequence(SymbolsToEntrezIDs("PRNP"), get_exons = TRUE)






