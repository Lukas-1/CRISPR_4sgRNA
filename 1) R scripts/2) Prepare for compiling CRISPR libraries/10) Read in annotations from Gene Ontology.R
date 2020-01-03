### 8th November 2019 ###





# Import packages and source code -----------------------------------------






# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files", "1) General")




# Read in data ------------------------------------------------------------

BioMart_GO_df <- read.table(file.path(CRISPR_input_directory, "Sublibraries", "Gene Ontology", "biomart_export_2019-08-11_Gene_Ontology.txt"),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE,
                            fill = TRUE
                            )



# Define functions --------------------------------------------------------

SearchGOTerms <- function(search_term, unique_terms = unique_GO_terms, fixed = FALSE) {
  grep(search_term, unique_terms, fixed = fixed, value = TRUE)
}





# Find GO terms matching certain keywords ---------------------------------

unique_GO_terms <- unique(BioMart_GO_df[, "GO term name"])
all_GO_terms_TFs          <- SearchGOTerms("transcription factor", fixed = TRUE)
all_GO_terms_GPCRs        <- SearchGOTerms("G[- ]protein")
all_GO_terms_kinases      <- SearchGOTerms("kinase",      fixed = TRUE)
all_GO_terms_phosphatases <- SearchGOTerms("phosphatase", fixed = TRUE)
all_GO_terms_ion_channels <- SearchGOTerms("ion channel", fixed = TRUE)
all_GO_terms_transporters <- SearchGOTerms("transporter", fixed = TRUE)
all_GO_terms_ubiquitin    <- SearchGOTerms("ubiquitin",   fixed = TRUE)





# Save data ---------------------------------------------------------------

save(list = "BioMart_GO_df",
     file = file.path(general_RData_directory, "10) Read in annotations from Gene Ontology.RData")
     )




