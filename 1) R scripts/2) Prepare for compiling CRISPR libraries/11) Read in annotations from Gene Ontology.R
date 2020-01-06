### 8th November 2019 ###





# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")




# Read in data ------------------------------------------------------------

# Downloaded from https://www.ensembl.org/biomart/martview
BioMart_GO_df <- read.table(file.path(CRISPR_input_directory, "Sublibraries", "Gene Ontology", "biomart_export_2019-08-11_Gene_Ontology.txt"),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE,
                            fill = TRUE
                            )



# Define functions --------------------------------------------------------

SearchGOTerms <- function(search_term, unique_terms = unique_GO_terms, fixed = FALSE) {
  results_vec <- grep(search_term, unique_terms, fixed = fixed, value = TRUE)
  results_vec <- results_vec[order(nchar(results_vec), results_vec)]
  return(results_vec)
}

CountNumGenes <- function(GO_terms) {
  are_these_terms <- BioMart_GO_df[, "GO term name"] %in% GO_terms
  unique_IDs <- unique(BioMart_GO_df[are_these_terms, "NCBI gene ID"])
  unique_IDs <- unique_IDs[unique_IDs != ""]
  return(length(unique_IDs))
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


all_GO_terms_receptor     <- SearchGOTerms("receptor",   fixed = TRUE)
all_GO_terms_membrane_receptor <- SearchGOTerms("membrane receptor",   fixed = TRUE)





# Count the number of genes matching GO terms -----------------------------

CountNumGenes("kinase activity")
CountNumGenes("phosphatase activity")

CountNumGenes("ion channel activity")
CountNumGenes("transporter activity")

CountNumGenes("protein ubiquitination")
CountNumGenes("protein deubiquitination")

CountNumGenes("G protein-coupled receptor activity")

CountNumGenes("receptor complex")
CountNumGenes("transmembrane receptor protein serine/threonine kinase activity")
CountNumGenes("transmembrane receptor protein tyrosine kinase activity")






# Save data ---------------------------------------------------------------

save(list = "BioMart_GO_df",
     file = file.path(general_RData_directory, "11) Read in annotations from Gene Ontology.RData")
     )




