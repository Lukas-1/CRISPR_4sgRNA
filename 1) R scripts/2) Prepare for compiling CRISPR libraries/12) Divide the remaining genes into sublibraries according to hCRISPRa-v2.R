### 13th February 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("eulerr")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))
source(file.path(general_functions_directory, "24) Assigning genes to sublibraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")

CRISPRa_datasets_directory       <- file.path(CRISPR_input_directory, "CRISPR libraries", "CRISPRa")
CRISPRa_Horlbeck2016_path        <- file.path(CRISPRa_datasets_directory, "Horlbeck, Kampmann, Weissman - eLife 2016")
CRISPRa_Horlbeck2016_sgRNAs_path <- file.path(CRISPRa_Horlbeck2016_path, "2016 - Compact and highly active next-generation libraries - Table S5.xlsx")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))
load(file.path(general_RData_directory, "11) Compile gene groups from HGNC.RData"))





# Read in data ------------------------------------------------------------

hCRISPRa_v2_df <- data.frame(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, skip = 7)[-1, ], stringsAsFactors = FALSE)
names(hCRISPRa_v2_df) <- names(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, n_max = 1))


# Downloaded from https://www.ensembl.org/biomart/martview
BioMart_GO_df <- read.table(file.path(CRISPR_input_directory, "Sublibraries", "Gene Ontology", "biomart_export_2020-03-25_Gene_Ontology.txt"),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE,
                            fill = TRUE
                            )
# Downloaded from https://www.proteinatlas.org/search/protein_class:Predicted+membrane+proteins
HPA_membrane_proteins_df <- read.table(file.path(CRISPR_input_directory, "Human Genome", "Human Protein Atlas", "protein_class_Predicted_Membrane_Proteins__2020_03_11.tsv"),
                                       sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE,
                                       fill = TRUE
                                       )





# Map membrane proteins from the Human Protein Atlas to Entrez IDs --------

HPA_membrane_proteins_df <- HPA_membrane_proteins_df[, c("Gene", "Ensembl", "Chromosome", "Gene description")]
colnames(HPA_membrane_proteins_df)[1:2] <- c("Gene_symbol", "Ensembl_gene_ID")
colnames(HPA_membrane_proteins_df) <- gsub(" ", "_", colnames(HPA_membrane_proteins_df), fixed = TRUE)

HPA_ensembl_membrane_df <- MapEnsemblIDs(HPA_membrane_proteins_df, warn = FALSE)





# Select gene-sublibrary mappings from the hCRISPRa-v2 library ------------

hC2_sublibrary_df <- unique(hCRISPRa_v2_df[, c("gene", "Sublibrary")])
colnames(hC2_sublibrary_df)[colnames(hC2_sublibrary_df) == "Sublibrary"] <- "Sublibrary_code"

hC2_sublibrary_df <- hC2_sublibrary_df[hC2_sublibrary_df[["gene"]] != "negative_control", ]

if (any(duplicated(hC2_sublibrary_df[["gene"]]))) {
  stop("Duplicated symbol-sublibrary combination found!")
}

hC2_sublibrary_df[["Sublibrary"]] <- hCRISPRa_v2_sublibrary_map[hC2_sublibrary_df[["Sublibrary_code"]]]





# Map gene symbols to Entrez IDs ------------------------------------------

hC2_sublibrary_df <- data.frame(MapToEntrezs(symbols_vec = hC2_sublibrary_df[["gene"]]),
                                hC2_sublibrary_df[, c("Sublibrary_code", "Sublibrary")],
                                stringsAsFactors = FALSE
                                )
hC2_sublibrary_df <- hC2_sublibrary_df[, colnames(hC2_sublibrary_df) != "Original_entrez"]





# Resolve duplicated Entrez IDs, where possible ---------------------------

num_occurrences_vec <- table(hC2_sublibrary_df[["Entrez_ID"]])[hC2_sublibrary_df[["Entrez_ID"]]]
multiplicates_df <- hC2_sublibrary_df[(num_occurrences_vec >= 2) %in% TRUE, ]
multiplicates_df <- multiplicates_df[order(GetMinEntrez(multiplicates_df[["Entrez_ID"]])), ]

unicates_df <- hC2_sublibrary_df[(num_occurrences_vec == 1) %in% TRUE, ]

resolved_df_list <- by(multiplicates_df,
                       factor(multiplicates_df[["Entrez_ID"]], levels = unique(multiplicates_df[["Entrez_ID"]])),
                       function(x) {
                         are_original_symbol <- x[["Original_symbol"]] == ""
                         if (any(are_original_symbol)) {
                           x[are_original_symbol, ]
                         } else {
                           x
                         }
                       }, simplify = FALSE
                       )

hC2_sublibrary_df <- do.call(rbind.data.frame, c(list(unicates_df), resolved_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
hC2_sublibrary_df <- hC2_sublibrary_df[order(GetMinEntrez(hC2_sublibrary_df[["Entrez_ID"]])), ]
rownames(hC2_sublibrary_df) <- NULL




# Define the Entrez IDs from prioritized sub-libraries --------------------

TF_entrez_IDs <- TidyEntrezs(all_TF_df[["Entrez_ID"]][all_TF_df[["Is_TF"]] == "Yes"])

are_secretome <- !(secretome_df[["Supercategory"]] %in% "Intracellular or membrane-bound")
secretome_entrez_IDs <- TidyEntrezs(secretome_df[["Entrez_ID"]][are_secretome])

GPCR_entrez_IDs <- TidyEntrezs(GPCR_df[["Entrez_ID"]])

membrane_protein_entrez_IDs <- TidyEntrezs(HPA_ensembl_membrane_df[["Entrez_ID"]])




# Create a combined list --------------------------------------------------

sublibrary_list <- c(
  list(
    "Transcription factors" = TF_entrez_IDs,
    "GPCRs"                 = GPCR_entrez_IDs,
    "Secretome"             = secretome_entrez_IDs
  ),
  tapply(hC2_sublibrary_df[["Entrez_ID"]],
         factor(hC2_sublibrary_df[["Sublibrary"]], levels = hCRISPRa_v2_sublibrary_map),
         function(x) TidyEntrezs(unlist(strsplit(x, ", ", fixed = TRUE))),
         simplify = FALSE
         )
)




# Assign protein-coding genes to each of the sub-libraries ----------------

are_in_library_mat <- do.call(cbind, lapply(sublibrary_list, function(x) collected_entrez_IDs %in% x))
rownames(are_in_library_mat) <- collected_entrez_IDs

assignment_vec <- apply(are_in_library_mat, 1, function(x) if (!(any(x))) "None" else colnames(are_in_library_mat)[[which(x)[[1]]]])

sublibrary_df <- data.frame("Entrez_ID"   = collected_entrez_IDs,
                            MapToEntrezs(entrez_IDs_vec = collected_entrez_IDs)["Gene_symbol"],
                            "Sublibrary"  = factor(assignment_vec, levels = c(colnames(are_in_library_mat), "None")),
                            stringsAsFactors = FALSE
                            )





# Examine genes that are found in multiple sub-libraries ------------------

sublibrary_df[are_in_library_mat[, "Transcription factors"] & are_in_library_mat[, "Secretome"], ]
sublibrary_df[are_in_library_mat[, "Secretome"] & are_in_library_mat[, "GPCRs"], ]

any(are_in_library_mat[, "Transcription factors"] & are_in_library_mat[, "GPCRs"])
any(are_in_library_mat[, "Secretome"] & are_in_library_mat[, "GPCRs"])

plot(eulerr::euler(are_in_library_mat[, c("Transcription factors", "Secretome", "GPCRs")]),
     quantities = list(font = 2, cex = 0.4)
     )




# Check the overlap between 'secretome' and 'membrane' genes --------------

## This is somewhat obsolete, since 'Intracellular or membrane-bound' genes
## are no longer included in the 'secretome' sub-library

secretome_entrezs <- secretome_df[["Entrez_ID"]][secretome_df[["Supercategory"]] %in% "Intracellular or membrane-bound"]
membrane_entrezs <- hC2_sublibrary_df[["Entrez_ID"]][hC2_sublibrary_df[["Sublibrary"]] %in% "Membrane Proteins"]

secretome_entrezs <- intersect(secretome_entrezs, collected_entrez_IDs)
membrane_entrezs <- intersect(membrane_entrezs, collected_entrez_IDs)

length(setdiff(secretome_entrezs, membrane_entrezs))
length(setdiff(membrane_entrezs, secretome_entrezs))
length(intersect(membrane_entrezs, secretome_entrezs))






# Define functions for interrogating GO terms -----------------------------

FindGOTerms <- function(search_term, unique_terms = unique_GO_terms, fixed = FALSE) {
  results_vec <- grep(search_term, unique_terms, fixed = fixed, value = TRUE)
  results_vec <- results_vec[order(nchar(results_vec), results_vec)]
  return(results_vec)
}

GenesForGoTerm <- function(GO_term) {
  # Requires 'BioMart_GO_df' and 'sublibrary_df' in the global environment
  unique_entrezs <- EntrezsForTerm(GO_term)
  results_df <- data.frame(
    "Entrez_ID"   = unique_entrezs,
    "Gene_symbol" = MapToEntrezs(entrez_IDs_vec = unique_entrezs)[["Gene_symbol"]],
    "Sublibrary"  = sublibrary_df[["Sublibrary"]][match(unique_entrezs, sublibrary_df[["Entrez_ID"]])],
    stringsAsFactors = FALSE
  )
  return(results_df)
}


LibrariesForTerm <- function(GO_term) {
  table(GenesForGoTerm(GO_term)[["Sublibrary"]])
}

EntrezsForTerm <- function(GO_term) {
  are_this_term <- BioMart_GO_df[["GO term name"]] == GO_term
  unique_entrezs <- unique(BioMart_GO_df[["NCBI gene ID"]][are_this_term])
  unique_entrezs <- as.character(unique_entrezs[!(is.na(unique_entrezs))])
  if (length(unique_entrezs) == 0) {
    stop(paste0("No Entrez IDs were found for the GO term '", GO_term, "'!"))
  }
  return(unique_entrezs)
}






# Investigate promising GO terms ------------------------------------------

unique_GO_terms <- unique(BioMart_GO_df[["GO term name"]])

## Look for GO terms that match specific keywords
all_GO_terms_kinases           <- FindGOTerms("kinase",      fixed = TRUE)
all_GO_terms_phosphatases      <- FindGOTerms("phosphatase", fixed = TRUE)
all_GO_terms_ion_channels      <- FindGOTerms("ion channel", fixed = TRUE)
all_GO_terms_transporters      <- FindGOTerms("transporter", fixed = TRUE)
all_GO_terms_receptor          <- FindGOTerms("receptor", fixed = TRUE)
all_GO_terms_membrane_receptor <- FindGOTerms("membrane receptor", fixed = TRUE)
all_GO_terms_transmembrane     <- FindGOTerms("transmembrane", fixed = TRUE)
all_GO_terms_gene_expression   <- FindGOTerms("gene expression", fixed = TRUE)



## For various GO terms, check which sublibraries were assigned to the genes that contain this term
LibrariesForTerm("kinase activity")
LibrariesForTerm("phosphatase activity")

LibrariesForTerm("membrane")
LibrariesForTerm("plasma membrane")
LibrariesForTerm("transmembrane transporter activity")
LibrariesForTerm("transmembrane receptor protein tyrosine kinase activity")
LibrariesForTerm("transmembrane receptor protein tyrosine phosphatase activity")
LibrariesForTerm("transmembrane receptor protein serine/threonine kinase activity")

LibrariesForTerm("G protein-coupled receptor activity")

LibrariesForTerm("gene expression")
LibrariesForTerm("regulation of gene expression")





## Confirm that the more general terms contain all genes present in the more specific terms
table(GenesForGoTerm("transmembrane transporter activity")[["Entrez_ID"]] %in% GenesForGoTerm("membrane")[["Entrez_ID"]])
table(c(GenesForGoTerm("transmembrane receptor protein tyrosine kinase activity")[["Entrez_ID"]],
        GenesForGoTerm("transmembrane receptor protein tyrosine phosphatase activity")[["Entrez_ID"]],
        GenesForGoTerm("transmembrane receptor protein serine/threonine kinase activity")[["Entrez_ID"]]
        ) %in% GenesForGoTerm("membrane")[["Entrez_ID"]])







# Categorize genes based on GO terms --------------------------------------

are_none_or_unassigned <- sublibrary_df[["Sublibrary"]] %in% c("None", "Unassigned")


## Assign kinases/phosphatases
kinase_entrezs <- EntrezsForTerm("kinase activity")
phosphatase_entrezs <- EntrezsForTerm("phosphatase activity")

are_new_kinases <- are_none_or_unassigned &
                   (sublibrary_df[["Entrez_ID"]] %in% c(kinase_entrezs, phosphatase_entrezs))

sublibrary_df[["Sublibrary"]][are_new_kinases] <- "Kinases/Phosphatases/Drug Targets"


## Assign membrane proteins
are_membrane_proteins <- sublibrary_df[["Entrez_ID"]] %in% HPA_ensembl_membrane_df[["Consensus_entrez"]]

receptor_activity_terms <- FindGOTerms("receptor activity$", fixed = FALSE)
are_eligible_terms <- !(grepl("^negative regulation of ", receptor_activity_terms)) &
                      !(grepl("^positive regulation of ", receptor_activity_terms)) &
                      !(grepl("^regulation of ", receptor_activity_terms))
receptor_activity_terms <- receptor_activity_terms[are_eligible_terms]

receptor_or_transporter_terms <- c(receptor_activity_terms, "transmembrane transporter activity")

receptor_entrezs <- unique(unlist(lapply(receptor_or_transporter_terms, EntrezsForTerm)))
are_receptors <- sublibrary_df[["Entrez_ID"]] %in% receptor_entrezs

are_ion_channels <- sublibrary_df[["Entrez_ID"]] %in% ion_channel_df[["Entrez_ID"]]

are_new_membrane_proteins <- are_none_or_unassigned &
                             are_membrane_proteins &
                             (are_ion_channels | are_receptors)
sublibrary_df[["Sublibrary"]][are_new_membrane_proteins] <- "Membrane Proteins"






# Define an inclusive list of Entrez IDs for sublibraries -----------------
# ... that also includes non-protein coding genes for transcription factors
# and the secretome

sublibraries_all_entrezs_list <- split(sublibrary_df[["Entrez_ID"]], sublibrary_df[["Sublibrary"]])
sublibraries_all_entrezs_list[["Transcription factors"]] <- TF_entrez_IDs
sublibraries_all_entrezs_list[["Secretome"]] <- setdiff(secretome_entrez_IDs, c(TF_entrez_IDs, GPCR_entrez_IDs))






# Save data ---------------------------------------------------------------

save(list = c("sublibrary_df", "sublibraries_all_entrezs_list"),
     file = file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData")
     )





