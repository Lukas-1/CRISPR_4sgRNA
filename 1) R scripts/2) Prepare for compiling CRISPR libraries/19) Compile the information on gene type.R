### 4th June 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "28) Merging the FANTOM5 and BioMart TSS data.R")) # For TidyBioMartChromosomes




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
human_genome_directory  <- file.path(CRISPR_input_directory, "Human genome")






# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(general_RData_directory, "18) Process the annotations from GENCODE.RData"))





# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# on 25 March 2020
NCBI_Hs_info_df <- read.table(file.path(human_genome_directory, "NCBI", "Homo_sapiens.gene_info"),
                              sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                              fill = TRUE, check.names = FALSE, comment.char = ""
                              )

BioMart_gene_type_df <- read.table(file.path(human_genome_directory, "Ensembl", "BioMart_human_2020-06-03_mart_export_with_gene_type.txt"),
                                   sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                                   )



# Define functions --------------------------------------------------------

CollapseFunction <- function(char_vec) {
  char_vec <- unique(char_vec)
  if (all(is.na(char_vec))) {
    return(NA_character_)
  } else {
    return(paste0(char_vec, collapse = ", "))
  }
}



# Define common gene type mappings ----------------------------------------

NCBI_gene_type_map <- c(
  "protein-coding"    = "protein-coding",
  "tRNA"              = "tRNA",
  "rRNA"              = "rRNA",
  "snoRNA"            = "snoRNA",
  "snRNA"             = "snRNA",
  "pseudo"            = "pseudogene",
  "ncRNA"             = "ncRNA",
  "scRNA"             = "other",
  "biological-region" = "other",
  "other"             = "other",
  "unknown"           = "unknown"
)


BioMart_gene_type_map <- c(

  "protein_coding"                     = "protein-coding",

  "IG_C_gene"                          = "protein-coding",
  "IG_D_gene"                          = "protein-coding",
  "IG_J_gene"                          = "protein-coding",
  "IG_V_gene"                          = "protein-coding",
  "TR_C_gene"                          = "protein-coding",
  "TR_D_gene"                          = "protein-coding",
  "TR_J_gene"                          = "protein-coding",
  "TR_V_gene"                          = "protein-coding",

  "Mt_tRNA"                            = "tRNA",
  "rRNA"                               = "rRNA",
  "Mt_rRNA"                            = "rRNA",

  "miRNA"                              = "miRNA",

  "scaRNA"                             = "scaRNA",
  "snoRNA"                             = "snoRNA",
  "snRNA"                              = "snRNA",
  "vaultRNA"                           = "vtRNA",
  "ribozyme"                           = "ribozyme",

  "lncRNA"                             = "lncRNA",

  "TR_V_pseudogene"                    = "pseudogene",
  "transcribed_processed_pseudogene"   = "pseudogene",
  "transcribed_unitary_pseudogene"     = "pseudogene",
  "transcribed_unprocessed_pseudogene" = "pseudogene",
  "translated_processed_pseudogene"    = "pseudogene",
  "translated_unprocessed_pseudogene"  = "pseudogene",
  "unitary_pseudogene"                 = "pseudogene",
  "unprocessed_pseudogene"             = "pseudogene",
  "IG_C_pseudogene"                    = "pseudogene",
  "IG_J_pseudogene"                    = "pseudogene",
  "IG_pseudogene"                      = "pseudogene",
  "IG_V_pseudogene"                    = "pseudogene",
  "polymorphic_pseudogene"             = "pseudogene",
  "processed_pseudogene"               = "pseudogene",
  "pseudogene"                         = "pseudogene",
  "rRNA_pseudogene"                    = "pseudogene",
  "TR_J_pseudogene"                    = "pseudogene",

  "scRNA"                              = "other",
  "sRNA"                               = "other",
  "misc_RNA"                           = "other",
  "non_stop_decay"                     = "other",
  "nonsense_mediated_decay"            = "other",
  "processed_transcript"               = "other",
  "retained_intron"                    = "other",

  "TEC"                                = "unknown"
)





# Explore NCBI's Homo_sapiens.gene_info data ------------------------------

anyNA(NCBI_Hs_info_df[["GeneID"]])
table(NCBI_Hs_info_df[["type_of_gene"]])

any(duplicated(NCBI_Hs_info_df[["GeneID"]]))





# Explore the BioMart export (containing gene/transcript type) ------------

table(BioMart_gene_type_df[["Gene type"]])
table(BioMart_gene_type_df[["Transcript type"]])

entrez_column <- "NCBI gene (formerly Entrezgene) ID"
are_NA_entrezs <- is.na(BioMart_gene_type_df[[entrez_column]])

BioMart_gene_type_df[are_NA_entrezs, ]

table(BioMart_gene_type_df[["Gene type"]][are_NA_entrezs])

anyNA(BioMart_gene_type_df[["Gene stable ID"]])

table(are_NA_entrezs)
table(duplicated(BioMart_gene_type_df[[entrez_column]]))





# Explore the GENCODE data frame ------------------------------------------

stopifnot(identical(sort(unique(gencode_df[["Gene_type"]])),
                    sort(unique(BioMart_gene_type_df[["Gene type"]]))
                    ))




# Reduce the BioMart table to Entrez ID and gene type ---------------------

BioMart_Entrez_type_only_df <- unique(BioMart_gene_type_df[!(are_NA_entrezs), c(entrez_column, "Gene type") ])

colnames(BioMart_Entrez_type_only_df) <- c("NCBI_gene_ID", "Gene_type")
rownames(BioMart_Entrez_type_only_df) <- NULL

any(duplicated(BioMart_Entrez_type_only_df[["NCBI_gene_ID"]]))

NCBI_IDs_char_vec <- as.character(BioMart_Entrez_type_only_df[["NCBI_gene_ID"]])

## Explore duplicates (Entrez IDs with multiple gene types)
num_occurrences <- table(NCBI_IDs_char_vec)[NCBI_IDs_char_vec]
duplicates_df <- BioMart_Entrez_type_only_df[num_occurrences > 1, ]
duplicates_df <- duplicates_df[order(duplicates_df[["NCBI_gene_ID"]]), ]





# Collect all Entrez IDs --------------------------------------------------

collected_entrezs_vec <- unique(c(NCBI_Hs_info_df[["GeneID"]]),
                                BioMart_Entrez_type_only_df[["NCBI_gene_ID"]]
                                )
collected_entrezs_vec <- sort(collected_entrezs_vec)






# Ascertain that ENSG IDs are only for genes that lack Entrez IDs ---------

NA_entrez_ENSGs_vec <- unique(BioMart_gene_type_df[["Gene stable ID"]][are_NA_entrezs])
are_NA_entrez_ENSG <- BioMart_gene_type_df[["Gene stable ID"]] %in% NA_entrez_ENSGs_vec
stopifnot(all(is.na(BioMart_gene_type_df[[entrez_column]][are_NA_entrez_ENSG])))





# Reduce the BioMart table to Ensembl gene ID and gene type ---------------

BioMart_ENSG_gene_type_df <- unique(BioMart_gene_type_df[are_NA_entrezs, c("Gene stable ID", "Gene type")])

BioMart_ENSG_gene_type_df <- BioMart_ENSG_gene_type_df[order(BioMart_ENSG_gene_type_df[["Gene stable ID"]]), ]

colnames(BioMart_ENSG_gene_type_df) <- c("Ensembl_gene_ID", "Original_type")
rownames(BioMart_ENSG_gene_type_df) <- NULL

stopifnot(!(any(duplicated(BioMart_ENSG_gene_type_df[["Ensembl_gene_ID"]]))))




# Reduce the GENCODE table to Ensembl gene ID and gene type ---------------

gencode_ENSG_gene_type_df <- unique(gencode_df[, c("Ensembl_gene_ID", "Gene_type")])

gencode_ENSG_gene_type_df <- gencode_ENSG_gene_type_df[order(gencode_ENSG_gene_type_df[["Ensembl_gene_ID"]]), ]
rownames(gencode_ENSG_gene_type_df) <- NULL
colnames(gencode_ENSG_gene_type_df)[[2]] <- "Original_type"

gencode_ENSG_matches <- match(gencode_ENSG_gene_type_df[["Original_type"]],
                              names(BioMart_gene_type_map)
                              )
gencode_ENSG_gene_type_df[, "Recoded_type"] <- BioMart_gene_type_map[gencode_ENSG_matches]
gencode_ENSG_gene_type_df <- gencode_ENSG_gene_type_df[ c(1, 3, 2)]

stopifnot(!(any(duplicated(gencode_ENSG_gene_type_df[["Ensembl_gene_ID"]]))))




# Collect all Ensembl gene IDs --------------------------------------------

collected_ENSG_vec <- sort(unique(BioMart_ENSG_gene_type_df[["Ensembl_gene_ID"]]))





# Complete the data frame for Ensembl gene IDs ----------------------------

BioMart_ENSG_matches <- match(BioMart_ENSG_gene_type_df[["Original_type"]],
                              names(BioMart_gene_type_map)
                              )
BioMart_ENSG_gene_type_df[["Recoded_type"]] <- BioMart_gene_type_map[BioMart_ENSG_matches]

ENSG_ID_matches <- match(BioMart_ENSG_gene_type_df[["Ensembl_gene_ID"]],
                         BioMart_gene_type_df[["Gene stable ID"]]
                         )

BioMart_ENSG_gene_type_df[["BioMart_gene_name"]] <- BioMart_gene_type_df[["Gene name"]][ENSG_ID_matches]

BioMart_ENSG_gene_type_df[["TSS"]] <- BioMart_gene_type_df[["Transcription start site (TSS)"]][ENSG_ID_matches]
BioMart_ENSG_gene_type_df[["Chromosome"]] <- TidyBioMartChromosomes(BioMart_gene_type_df[["Chromosome/scaffold name"]][ENSG_ID_matches])
BioMart_ENSG_gene_type_df[["Strand"]] <- ifelse(BioMart_gene_type_df[["Strand"]][ENSG_ID_matches] == -1, "-", "+")

BioMart_ENSG_columns <- c("Ensembl_gene_ID", "BioMart_gene_name",
                          "Recoded_type", "Original_type",
                          "Chromosome", "Strand", "TSS"
                          )
BioMart_ENSG_gene_type_df <- BioMart_ENSG_gene_type_df[, BioMart_ENSG_columns]





# Process NCBI's gene type information ------------------------------------

NCBI_entrez_matches <- match(collected_entrezs_vec, NCBI_Hs_info_df[["GeneID"]])
NCBI_types_vec <- NCBI_Hs_info_df[["type_of_gene"]][NCBI_entrez_matches]

NCBI_types_matches <- match(NCBI_types_vec, names(NCBI_gene_type_map))
NCBI_types_recoded_vec <- NCBI_gene_type_map[NCBI_types_matches]




# Process BioMart's gene type information for Entrez IDs ------------------

BioMart_types_list <- lapply(collected_entrezs_vec, function(x) {
  are_this_entrez <- which(BioMart_Entrez_type_only_df[["NCBI_gene_ID"]] == x)
  if (any(are_this_entrez)) {
    return(unique(BioMart_Entrez_type_only_df[["Gene_type"]][are_this_entrez]))
  } else {
    return(NA_integer_)
  }
})

BioMart_entrez_expanded_df <- data.frame(
  "Entrez_ID"      = rep(collected_entrezs_vec, lengths(BioMart_types_list)),
  "Original_type"  = unlist(BioMart_types_list, use.names = FALSE),
  stringsAsFactors = FALSE
)

new_order <- order(
  BioMart_entrez_expanded_df[["Entrez_ID"]],
  match(BioMart_entrez_expanded_df[["Original_type"]],
        names(BioMart_gene_type_map)
        )
)
BioMart_entrez_expanded_df <- BioMart_entrez_expanded_df[new_order, ]
rownames(BioMart_entrez_expanded_df) <- NULL


BioMart_entrez_matches <- match(BioMart_entrez_expanded_df[["Original_type"]],
                                names(BioMart_gene_type_map)
                                )
BioMart_entrez_expanded_df[["Recoded_type"]] <- BioMart_gene_type_map[BioMart_entrez_matches]





# Summarize BioMart's gene type information for Entrez IDs ----------------

entrezs_fac <- factor(BioMart_entrez_expanded_df[["Entrez_ID"]],
                      levels = unique(BioMart_entrez_expanded_df[["Entrez_ID"]])
                      )
stopifnot(identical(collected_entrezs_vec, as.integer(levels(entrezs_fac))))

original_types_collapsed <- tapply(BioMart_entrez_expanded_df[["Original_type"]],
                                   entrezs_fac,
                                   CollapseFunction
                                   )
recoded_types_collapsed <- tapply(BioMart_entrez_expanded_df[["Recoded_type"]],
                                  entrezs_fac,
                                  CollapseFunction
                                  )
entrez_single_type_vec <- tapply(BioMart_entrez_expanded_df[["Recoded_type"]],
                                 entrezs_fac,
                                 function(x) x[[1]]
                                 )




# Check for discrepancies between NCBI and BioMart ------------------------

are_identical <- mapply(identical, NCBI_types_recoded_vec, entrez_single_type_vec)
have_NA <- is.na(NCBI_types_recoded_vec) | is.na(entrez_single_type_vec)
are_identical <- are_identical & !(have_NA)
are_conflicting <- !(are_identical) & !(have_NA)

table(are_identical)
table(are_conflicting)




# Collate the data for Entrez IDs -----------------------------------------

entrez_to_gene_type_df <- data.frame(
  "Entrez_ID"              = collected_entrezs_vec,
  "Gene_symbol"            = MapToEntrezs(as.character(collected_entrezs_vec))[["Gene_symbol"]],
  "Consensus_type"         = ifelse(is.na(NCBI_types_recoded_vec),
                                    entrez_single_type_vec,
                                    NCBI_types_recoded_vec
                                    ),
  "Recoded_type_NCBI"      = NCBI_types_recoded_vec,
  "Recoded_type_BioMart"   = entrez_single_type_vec,
  "Is_conflicting"         = ifelse(have_NA, NA, are_conflicting),
  "Original_type_NCBI"     = NCBI_types_vec,
  "Original_types_BioMart" = original_types_collapsed,
  "Recoded_types_BioMart"  = recoded_types_collapsed,
  stringsAsFactors         = FALSE,
  row.names                = NULL
)




# Perform some checks -----------------------------------------------------

table(entrez_to_gene_type_df[["Is_consistent"]], useNA = "ifany")
entrez_to_gene_type_df[entrez_to_gene_type_df[["Is_conflicting"]] %in% TRUE, ]






# Save data ---------------------------------------------------------------

save(list = c("entrez_to_gene_type_df", "BioMart_ENSG_gene_type_df", "gencode_ENSG_gene_type_df"),
     file = file.path(general_RData_directory, "19) Compile the information on gene type.RData")
     )



