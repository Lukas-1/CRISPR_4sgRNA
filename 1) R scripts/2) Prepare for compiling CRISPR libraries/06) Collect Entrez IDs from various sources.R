### 30 October 2019 ###






# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "24) Assigning genes to sublibraries.R"))






# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
NCBI_directory          <- file.path(CRISPR_input_directory, "Human genome", "NCBI")
annotation_intermediate_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Annotation")




# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# on 25 March 2020
NCBI_Hs_info_df <- read.table(file.path(NCBI_directory, "Homo_sapiens.gene_info"),
                              sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                              fill = TRUE, check.names = FALSE, comment.char = ""
                              )



# Downloaded from ftp://ftp.ensembl.org/pub/release-99/tsv/homo_sapiens/Homo_sapiens.GRCh38.99.entrez.tsv.gz
# on 25 March 2020
Ensembl_Hs_entrez_df <- read.table(file.path(CRISPR_input_directory, "Human genome", "Ensembl", "Homo_sapiens.GRCh38.99.entrez.tsv"),
                                   sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                                   check.names = FALSE
                                   )


# Downloaded from https://www.ncbi.nlm.nih.gov/gene/?term=%22Homo+sapiens%22%5BOrganism%5D+AND+%22source_+genomic%22%5Bproperties%5D+AND+%22genetype+protein+coding%22%5BProperties%5D+AND+alive%5Bprop%5D
# (send to file, format ASN.1)
# on 15 February 2020
NCBI_gene_ASN1_files <- grep("NCBI_gene_protein_coding_ASN1", list.files(NCBI_directory), fixed = TRUE, value = TRUE)
NCBI_gene_ASN1_vec_list <- lapply(NCBI_gene_ASN1_files, function(x) readLines(file.path(NCBI_directory, x)))




# Collect Entrez IDs from tabulated NCBI and Ensembl data -----------------

NCBI_Hs_info_df <- NCBI_Hs_info_df[NCBI_Hs_info_df[["#tax_id"]] == "9606", ]
NCBI_entrezs_vec <- unique(as.character(NCBI_Hs_info_df[["GeneID"]][NCBI_Hs_info_df[["type_of_gene"]] == "protein-coding"]))

Ensembl_Hs_have_entrez <- !(is.na(suppressWarnings(as.integer(Ensembl_Hs_entrez_df[["xref"]]))))
Ensembl_Hs_are_protein_coding <- Ensembl_Hs_entrez_df[["protein_stable_id"]] != "-"
Ensembl_entrezs_vec <- unique(Ensembl_Hs_entrez_df[["xref"]][Ensembl_Hs_have_entrez & Ensembl_Hs_are_protein_coding])




# Parse the ASN.1 output from NCBI gene -----------------------------------

NCBI_gene_ASN1_vec <- do.call(c, NCBI_gene_ASN1_vec_list)

are_new_gene <- NCBI_gene_ASN1_vec == "Entrezgene ::= {"
counter <- 0L
number_vec <- rep(NA_integer_, length(NCBI_gene_ASN1_vec))
for (i in seq_along(NCBI_gene_ASN1_vec)) {
  if (are_new_gene[[i]]) {
    counter <- counter + 1L
  }
  number_vec[[i]] <- counter
}

ASN1_trimmed_vec <- trimws(NCBI_gene_ASN1_vec)
ASN1_list <- unname(split(ASN1_trimmed_vec, number_vec))

NCBI_gene_entrezs_vec <- vapply(ASN1_list, function(x) sub(",", "", sub("geneid ", "", x[[3]], fixed = TRUE), fixed = TRUE), "")
if (anyNA(as.integer(NCBI_gene_entrezs_vec))) {
  stop("Invalid Entrez IDs found after parsing the NCBI ASN.1 file!")
}




# Search the NCBI ASN.1 output for specific comments ----------------------

SearchNCBIGene <- function(string) {
  vapply(ASN1_list, function(x) any(grepl(string, x, fixed = TRUE)), logical(1))
}
not_in_current_release  <- SearchNCBIGene("not in current annotation release")
partial_on_assembly     <- SearchNCBIGene("partial on reference assembly")
only_on_alternate_loci  <- SearchNCBIGene("only annotated on alternate loci in reference assembly")
only_on_loci_or_patches <- SearchNCBIGene("only annotated on alternate loci and patches unit in reference")
only_on_patches         <- SearchNCBIGene("only annotated on patches unit in reference assembly")





# Integrate information on Entrez IDs -------------------------------------

collected_entrez_IDs <- TidyEntrezs(union(NCBI_entrezs_vec, Ensembl_entrezs_vec))

collected_entrezs_df <- data.frame(
  "Entrez_ID"               = collected_entrez_IDs,
  "Gene_symbol"             = EntrezIDsToSymbols(collected_entrez_IDs),
  "Is_in_NCBI_Hs_info"      = collected_entrez_IDs %in% NCBI_entrezs_vec,
  "Is_in_NCBI_gene"         = collected_entrez_IDs %in% NCBI_gene_entrezs_vec,
  "Not_on_current_assembly" = collected_entrez_IDs %in% NCBI_gene_entrezs_vec[not_in_current_release],
  "Partial_on_assembly"     = collected_entrez_IDs %in% NCBI_gene_entrezs_vec[partial_on_assembly],
  "Only_on_alternate_loci"  = collected_entrez_IDs %in% NCBI_gene_entrezs_vec[only_on_alternate_loci | only_on_loci_or_patches | only_on_patches],
  "Is_in_Ensembl"           = collected_entrez_IDs %in% Ensembl_entrezs_vec,
  stringsAsFactors = FALSE
)

collected_entrezs_df[["Category"]] <- vapply(seq_along(collected_entrez_IDs), function(x) {
  is_in_NCBI_info <- collected_entrezs_df[["Is_in_NCBI_Hs_info"]][[x]]
  is_in_Ensembl <- collected_entrezs_df[["Is_in_Ensembl"]][[x]]
  if (!(is_in_NCBI_info)) {
    return("Ensembl only (not in NCBI)")
  } else if (!(collected_entrezs_df[["Is_in_NCBI_gene"]][[x]])) {
    if (is_in_Ensembl) {
      return("Not in NCBI Gene, but in Ensembl")
    } else {
      return("Not in NCBI Gene or in Ensembl")
    }
  } else if (collected_entrezs_df[["Not_on_current_assembly"]][[x]]) {
    return("Not in current annotation release")
  } else if (collected_entrezs_df[["Partial_on_assembly"]][[x]]) {
    return("Partial on reference assembly")
  } else if (collected_entrezs_df[["Only_on_alternate_loci"]][[x]]) {
    return("Only annotated on alternate loci")
  } else {
    if (is_in_Ensembl) {
      return("Present in NCBI and Ensembl")
    } else {
      return("NCBI only (not in Ensembl)")
    }
  }
}, "")




# Export problematic Entrez IDs -------------------------------------------

if (!(all(collected_entrezs_df[["Is_in_NCBI_gene"]]))) {
  write.table(collected_entrezs_df[["Entrez_ID"]][!(collected_entrezs_df[["Is_in_NCBI_gene"]])],
              file = file.path(annotation_intermediate_files_directory, "Entrez_IDs_not_covered_by_NCBI_gene.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE
              )
}




# Save data ---------------------------------------------------------------

save(list = c("collected_entrez_IDs", "collected_entrezs_df"),
     file = file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData")
     )



