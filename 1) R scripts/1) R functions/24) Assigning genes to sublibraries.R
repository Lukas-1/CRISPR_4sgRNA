### 13th February 2020 ###




# Define functions --------------------------------------------------------

TidyEntrezs <- function(entrez_IDs_vec) {
  entrez_IDs_vec <- entrez_IDs_vec[!(is.na(entrez_IDs_vec))]
  entrez_IDs_vec <- entrez_IDs_vec[order(as.integer(entrez_IDs_vec))]
  return(entrez_IDs_vec)
}




# Define mappings ---------------------------------------------------------

hCRISPRa_v2_sublibrary_map <- c(
  "h6" = "Membrane Proteins",
  "h1" = "Kinases/Phosphatases/Drug Targets",
  "h4" = "Mitochondria/Trafficking/Motility",
  "h3" = "Stress/Proteostasis",
  "h2" = "Cancer/Apoptosis",
  "h5" = "Gene Expression",
  "h7" = "Unassigned"
)
