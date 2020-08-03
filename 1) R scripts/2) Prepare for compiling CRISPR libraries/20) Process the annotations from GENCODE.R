### 30th July 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))






# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
human_genome_directory <- file.path(CRISPR_root_directory, "2) Input data", "Human genome")
Gencode_directory      <- file.path(human_genome_directory, "Gencode")



# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz
# on 26 March 2020
original_gencode_df <- read.table(file.path(Gencode_directory, "gencode.v33.annotation.gff3"),
                                  sep = "\t",
                                  header = FALSE,
                                  stringsAsFactors = FALSE
                                  )[, -6]

# Downloaded from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.metadata.EntrezGene.gz
# on 31 July 2020
gencode_entrez_df <- read.table(file.path(Gencode_directory, "gencode.v33.metadata.EntrezGene"),
                                sep = "\t",
                                header = FALSE,
                                stringsAsFactors = FALSE
                                )

# Downloaded from ftp://ftp.ensembl.org/pub/release-99/tsv/homo_sapiens/Homo_sapiens.GRCh38.99.entrez.tsv.gz
# on 25 March 2020
Ensembl_Hs_entrez_df <- read.table(file.path(human_genome_directory, "Ensembl", "Homo_sapiens.GRCh38.99.entrez.tsv"),
                                   sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                                   check.names = FALSE
                                   )

BioMart_df <- read.table(file.path(human_genome_directory, "Ensembl", "BioMart_human_2020-03-25_mart_export.txt"),
                         sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                         )





# Define functions --------------------------------------------------------

Collapsify <- function(input_list) {
  vapply(input_list,
         function(x) {
           are_NA <- is.na(x)
           if (all(are_NA)) {
             return(NA_character_)
           } else {
             paste0(x[!(are_NA)], collapse = ", ")
           }
         },
         ""
         )
}

StripVersion <- function(char_vec) {
  splits_list <- strsplit(char_vec,   ".", fixed = TRUE)
  results_vec <- sapply(splits_list, "[", 1)
  return(results_vec)
}





# Assign column names -----------------------------------------------------

# See see https://m.ensembl.org/info/website/upload/gff3.html
names(original_gencode_df) <- c("seqid", "source", "type", "start", "end",
                                "strand", "phase", "attributes"
                                )
names(gencode_entrez_df) <- c("ENST_version", "Entrez_ID")





# Explore the Gencode data frame ------------------------------------------

table(original_gencode_df[["type"]])
table(original_gencode_df[["seqid"]])
table(original_gencode_df[["source"]])
table(original_gencode_df[["phase"]])






# Explore the Ensembl Homo_sapiens.GRCh38.99.entrez.tsv data frame --------

are_integer_Ensembl_Hs <- !(is.na(suppressWarnings(as.integer(Ensembl_Hs_entrez_df[["xref"]]))))

length(unique((Ensembl_Hs_entrez_df[["xref"]][!(are_integer_Ensembl_Hs)])))




# Process the Gencode-Entrez mappings -------------------------------------

gencode_entrez_df[["ENST_ID"]] <- StripVersion(gencode_entrez_df[["ENST_version"]])




# Check for duplicate Gencode-Entrez mappings -----------------------------

num_occurrences_ENST_vec <- table(gencode_entrez_df[["ENST_version"]]
                                  )[gencode_entrez_df[["ENST_version"]]]

duplicate_ENST_df <- gencode_entrez_df[num_occurrences_ENST_vec > 1, ]
order_vec <- order(duplicate_ENST_df[["ENST_ID"]],
                   duplicate_ENST_df[["ENST_version"]]
                   )
duplicate_ENST_df <- duplicate_ENST_df[order_vec, ]
row.names(duplicate_ENST_df) <- NULL

num_occurrences_entrez_vec <- table(gencode_entrez_df[["Entrez_ID"]]
                                    )[gencode_entrez_df[["Entrez_ID"]]]
table(num_occurrences_entrez_vec)




# Map Ensembl transcript IDs to Entrez IDs using the Gencode file ---------

gencode_ENST_version_to_entrez_list <- split(gencode_entrez_df[["Entrez_ID"]],
                                             gencode_entrez_df[["ENST_version"]]
                                             )
gencode_ENST_version_to_entrez_vec <- Collapsify(gencode_ENST_version_to_entrez_list)


gencode_ENST_ID_to_entrez_list <- split(gencode_entrez_df[["Entrez_ID"]],
                                        gencode_entrez_df[["ENST_ID"]]
                                        )
gencode_ENST_ID_to_entrez_vec <- Collapsify(gencode_ENST_ID_to_entrez_list)


table(lengths(gencode_ENST_version_to_entrez_list))
table(lengths(gencode_ENST_ID_to_entrez_list))



# Map Ensembl transcript IDs to Entrez IDs using the Ensembl file ---------

unique_entrezs <- unique(Ensembl_Hs_entrez_df[["transcript_stable_id"]][are_integer_Ensembl_Hs])
unique_entrezs <- unique_entrezs[order(as.integer(unique_entrezs))]

ensembl_ENST_ID_to_entrez_list <- split(Ensembl_Hs_entrez_df[["xref"]],
                                        factor(Ensembl_Hs_entrez_df[["transcript_stable_id"]],
                                               levels = unique_entrezs
                                               )
                                        )
ensembl_ENST_ID_to_entrez_vec <- Collapsify(ensembl_ENST_ID_to_entrez_list)


table(lengths(ensembl_ENST_ID_to_entrez_list))
length(unique(Ensembl_Hs_entrez_df[["xref"]][are_integer_Ensembl_Hs]))




# Map Ensembl transcript IDs to Entrez IDs using BioMart ------------------

BioMart_ENST_version_to_entrez_list <- split(BioMart_df[["NCBI gene ID"]],
                                             BioMart_df[["Transcript stable ID version"]]
                                             )
BioMart_ENST_version_to_entrez_vec <- Collapsify(BioMart_ENST_version_to_entrez_list)


BioMart_ENST_ID_to_entrez_list <- split(BioMart_df[["NCBI gene ID"]],
                                        BioMart_df[["Transcript stable ID"]]
                                        )
BioMart_ENST_ID_to_entrez_vec <- Collapsify(BioMart_ENST_ID_to_entrez_list)




# Process Gencode attributes ----------------------------------------------

ID_splits <- strsplit(original_gencode_df[["attributes"]], ";", fixed = TRUE)

ID_split_splits <- lapply(ID_splits, function(x) strsplit(x, "=", fixed = TRUE))

ID_categories <- Collapsify(ID_split_splits)
categories_table <- table(unlist(ID_categories))




# Extract specific attributes ---------------------------------------------

ExtractIDs <- function(ID_list_list, ID_tag, allow_missing = TRUE) {
  vapply(ID_list_list, function(x) {
    are_gene_ID <- vapply(x, function(y) y[[1]], "") == ID_tag
    ids_vec <- vapply(x, function(y) y[[2]], "")
    if (allow_missing && !(any(are_gene_ID))) {
      return(NA_character_)
    } else {
      return(ids_vec[[which(are_gene_ID)]])
    }
  }, "")
}

ENSG_versions_vec <- ExtractIDs(ID_split_splits, "gene_id")
ENST_versions_vec <- ExtractIDs(ID_split_splits, "transcript_id")
ENST_IDs_vec      <- StripVersion(ENST_versions_vec)
gene_symbols_vec  <- ExtractIDs(ID_split_splits, "gene_name")





# Construct a new data frame ----------------------------------------------

gencode_df <- data.frame(
  "Ensembl_gene_version"              = ENSG_versions_vec,
  "Ensembl_gene_ID"                   = StripVersion(ENSG_versions_vec),
  "Ensembl_transcript_version"        = ENST_versions_vec,
  "Ensembl_transcript_ID"             = ENST_IDs_vec,
  "Gencode_ENST_version_to_entrez_ID" = gencode_ENST_version_to_entrez_vec[ENST_versions_vec],
  "Gencode_ENST_ID_to_entrez_ID"      = gencode_ENST_ID_to_entrez_vec[ENST_IDs_vec],
  "Ensembl_Hs_ENST_ID_to_entrez_ID"   = ensembl_ENST_ID_to_entrez_vec[ENST_IDs_vec],
  "BioMart_ENST_version_to_entrez_ID" = BioMart_ENST_version_to_entrez_vec[ENST_versions_vec],
  "BioMart_ENST_ID_to_entrez_ID"      = BioMart_ENST_ID_to_entrez_vec[ENST_IDs_vec],
  "Gene_symbol"                       = gene_symbols_vec,
  "Chromosome"                        = original_gencode_df[["seqid"]],
  "Start"                             = original_gencode_df[["start"]],
  "End"                               = original_gencode_df[["end"]],
  "Strand"                            = original_gencode_df[["strand"]],
  stringsAsFactors                    = FALSE
)







