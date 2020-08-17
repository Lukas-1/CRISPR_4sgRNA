### 4th August 2020 ###



# Import packages and source code -----------------------------------------

library("TxDb.Hsapiens.UCSC.hg38.knownGene")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(general_RData_directory, "18) Process the annotations from GENCODE.RData"))
load(file.path(general_RData_directory, "19) Compile the information on gene type.RData"))






# Define functions --------------------------------------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End")

ResolveDuplicateFeatures <- function(locations_df) {

  if ("Ensembl_transcript_ID" %in% colnames(locations_df)) {
    transcript_vec <- locations_df[["Ensembl_transcript_ID"]]
  } else {
    transcript_vec <- rep(NA, nrow(locations_df))
  }

  new_order <- order(match(locations_df[["Chromosome"]],
                           c(paste0("chr", 1:23), c("Y", "Y", "M"))
                           ),
                     locations_df[["Chromosome"]],
                     match(locations_df[["Strand"]], c("+", "-")),
                     locations_df[["Start"]],
                     locations_df[["End"]],
                     GetMinEntrez(locations_df[["Entrez_ID"]]),
                     locations_df[["Ensembl_gene_ID"]],
                     transcript_vec
                     )

  locations_df <- locations_df[new_order, ]

  location_strings <- do.call(paste, c(locations_df[location_columns], sep = "_"))
  entrez_location_IDs <- ifelse(is.na(locations_df[["Entrez_ID"]]),
                                NA_character_,
                                paste0(locations_df[["Entrez_ID"]], "_", location_strings)
                                )
  ENSG_location_IDs <- ifelse(is.na(locations_df[["Ensembl_gene_ID"]]),
                              NA_character_,
                              paste0(locations_df[["Ensembl_gene_ID"]], "_", location_strings)
                              )

  location_fac <- factor(location_strings, levels = unique(location_strings))
  CheckThatFactorIsInOrder(location_fac)

  group_indices_list <- tapply(seq_along(location_fac),
                               location_fac,
                               function(x) {
                                 group_IDs_vec <- rep(NA_integer_, length(x))
                                 current_index <- 1L
                                 for (i in seq_along(x)) {
                                   if (!(is.na(group_IDs_vec[[i]]))) {
                                     next
                                   } else {
                                     this_entrez <- entrez_location_IDs[x[[i]]]
                                     this_ENSG <- ENSG_location_IDs[x[[i]]]
                                     are_this_entrez <- (entrez_location_IDs[x] == this_entrez) %in% TRUE
                                     are_this_ENSG <- (ENSG_location_IDs[x] == this_ENSG) %in% TRUE
                                     all_entrezs <- unique(c(entrez_location_IDs[x][are_this_ENSG], this_entrez))
                                     all_entrezs <- all_entrezs[!(is.na(all_entrezs))]
                                     all_ENSGs <- unique(c(ENSG_location_IDs[x][are_this_entrez], this_ENSG))
                                     all_ENSGs <- all_ENSGs[!(is.na(all_ENSGs))]
                                     are_shared_entrez <- entrez_location_IDs[x] %in% all_entrezs
                                     are_shared_ENSG <- ENSG_location_IDs[x] %in% all_ENSGs
                                     are_shared <- are_shared_entrez | are_shared_ENSG
                                     group_IDs_vec[are_shared] <- current_index
                                     current_index <- current_index + 1L
                                   }
                                 }
                                 return(group_IDs_vec)
                               }, simplify = FALSE)

  group_IDs_long <- paste0(location_strings, "|",
                           unlist(group_indices_list, use.names = FALSE)
                           )
  integer_group_IDs <- as.integer(factor(group_IDs_long, levels = unique(group_IDs_long)))

  assign("delete_integer_group_IDs", integer_group_IDs, envir = globalenv())

  sources_order <- c("TxDb", "BioMart", "Gencode")
  sources_vec <- tapply(locations_df[["Source"]],
                        integer_group_IDs,
                        function(x) {
                          unique_vec <- unique(x)
                          unique_vec <- unique_vec[order(match(unique_vec, sources_order))]
                          paste0(unique_vec, collapse = ", ")
                        })
  entrezs_vec <- tapply(locations_df[["Entrez_ID"]],
                        integer_group_IDs,
                        function(x) {
                          if (all(is.na(x))) {
                            return(NA_character_)
                          } else {
                            paste0(unique(x[!(is.na(x))]), collapse = ", ")
                          }
                        })
  ENSG_vec <- tapply(locations_df[["Ensembl_gene_ID"]],
                     integer_group_IDs,
                     function(x) {
                       if (all(is.na(x))) {
                         return(NA_character_)
                       } else {
                         paste0(unique(x[!(is.na(x))]), collapse = ", ")
                       }
                     })

  symbols_vec <- tapply(locations_df[["Original_symbol"]],
                        integer_group_IDs,
                        function(x) {
                          if (all(is.na(x))) {
                            return(NA_character_)
                          } else {
                            paste0(unique(x[!(is.na(x))]), collapse = ", ")
                          }
                        })

  full_locations_df <- locations_df

  are_duplicated <- duplicated(integer_group_IDs)
  locations_df <- locations_df[!(are_duplicated), ]

  are_problematic <- (ENSG_vec != locations_df[["Ensembl_gene_ID"]]) %in% TRUE

  if (any(are_problematic)) {
    message(paste0("Some locations (", sum(are_problematic), " in total) ",
                   " mapped to more than one Ensembl gene ID."
                   )
            )
    data.frame(locations_df[are_problematic, colnames(locations_df) != "Ensembl_gene_ID"],
               "Old_ENSG" = locations_df[["Ensembl_gene_ID"]][are_problematic],
               "New_ENSG" = ENSG_vec[are_problematic]
               )
    message("All transcript entries for problematic entries are listed below:")
    problematic_locations <- location_strings[!(are_duplicated)][are_problematic]
    print(full_locations_df[location_strings %in% problematic_locations, ])
  } else {
    message("No ambiguous Ensembl IDs (for the same locations) were found.")
  }

  are_different <- ((entrezs_vec != locations_df[["Entrez_ID"]]) %in% TRUE)
  data.frame("Old_entrez" = locations_df[["Entrez_ID"]][are_different],
             "New_entrez" = entrezs_vec[are_different]
             )

  locations_df[["Entrez_ID"]] <- entrezs_vec
  locations_df[["Ensembl_gene_ID"]] <- ENSG_vec
  locations_df[["Source"]] <- sources_vec
  locations_df[["Original_symbol"]] <- symbols_vec


  new_order <- order(GetMinEntrez(locations_df[["Entrez_ID"]]),
                     locations_df[["Ensembl_gene_ID"]],
                     match(locations_df[["Chromosome"]],
                           c(paste0("chr", 1:23), c("Y", "Y", "M"))
                           ),
                     locations_df[["Chromosome"]],
                     match(locations_df[["Strand"]], c("+", "-")),
                     locations_df[["Start"]],
                     locations_df[["End"]]
                     )
  locations_df <- locations_df[new_order, ]
  row.names(locations_df) <- NULL
  return(locations_df)
}






# Create data frames using TxDb.Hsapiens.UCSC.hg38.knownGene --------------

TxDb_genes_df <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                     single.strand.genes.only = FALSE
                                     ))
TxDb_transcripts_df <- as.data.frame(transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene))
TxDb_exons_df <- as.data.frame(exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))
TxDb_CDS_df <- as.data.frame(cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))

stopifnot(identical(sort(unique(TxDb_genes_df[["group_name"]])),
                    sort(unique(TxDb_transcripts_df[["group_name"]]))
                    ))
stopifnot(identical(sort(unique(TxDb_genes_df[["group_name"]])),
                    sort(unique(TxDb_transcripts_df[["group_name"]]))
                    ))
stopifnot(identical(sort(unique(TxDb_genes_df[["group_name"]])),
                    sort(unique(TxDb_exons_df[["group_name"]]))
                    ))








# Define column names -----------------------------------------------------

BioMart_columns <- c(
  "Entrez_ID", "ENSG", "ENST", "Gene_symbol",
  "Chromosome", "Strand",  "Gene_start", "Gene_end"
)
TxDb_columns <- c(
  "group_name", "seqnames", "strand", "start", "end"
)
location_columns <- c("Chromosome", "Strand", "Start", "End")
gencode_columns <- c(
  "Entrez_ID", "Ensembl_gene_ID", "Ensembl_transcript_ID", "Exon_ID",
  "Gene_symbol", location_columns
)

gene_columns <- c(
  "Entrez_ID", "Ensembl_gene_ID", "Original_symbol",
  location_columns
)








# Create a merged data frame of gene coordinates --------------------------

BioMart_genes_df <- unique(BioMart_tidied_df[, setdiff(BioMart_columns, "ENST")])
TxDb_genes_sel_df <- TxDb_genes_df[, TxDb_columns]
gencode_genes_df <- gencode_df[gencode_df[["Entry_type"]] == "gene",
                               setdiff(gencode_columns, c("Ensembl_transcript_ID", "Exon_ID"))
                               ]
row.names(BioMart_genes_df) <- NULL
row.names(gencode_genes_df) <- NULL

colnames(BioMart_genes_df) <- gene_columns
colnames(gencode_genes_df) <- gene_columns


colnames(TxDb_genes_sel_df) <- c("Entrez_ID", location_columns)

TxDb_genes_sel_df[["Original_symbol"]] <- NA_character_
TxDb_genes_sel_df[["Ensembl_gene_ID"]] <- NA_character_
TxDb_genes_sel_df[["Source"]] <- "TxDb"
BioMart_genes_df[["Source"]] <- "BioMart"
gencode_genes_df[["Source"]] <- "GENCODE"

gencode_genes_df <- gencode_genes_df[, c(gene_columns, "Source")]

gene_locations_df <- rbind.data.frame(
  TxDb_genes_sel_df, BioMart_genes_df, gencode_genes_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)

gene_locations_df <- ResolveDuplicateFeatures(gene_locations_df)





# Create a merged data frame of transcript coordinates --------------------

BioMart_transcripts_df <- BioMart_tidied_df[, BioMart_columns]
TxDb_transcripts_sel_df <- TxDb_transcripts_df[, TxDb_columns]
gencode_transcripts_df <- gencode_df[gencode_df[["Entry_type"]] == "transcript",
                                     setdiff(gencode_columns, "Exon_ID"),
                                     ]

row.names(gencode_transcripts_df) <- NULL

transcript_columns <- c(
  "Entrez_ID", "Ensembl_gene_ID", "Ensembl_transcript_ID", "Original_symbol",
  location_columns
)

colnames(BioMart_transcripts_df) <- transcript_columns
colnames(gencode_transcripts_df) <- transcript_columns

colnames(TxDb_transcripts_sel_df) <- c("Entrez_ID", location_columns)

TxDb_transcripts_sel_df[["Original_symbol"]] <- NA_character_
TxDb_transcripts_sel_df[["Ensembl_gene_ID"]] <- NA_character_
TxDb_transcripts_sel_df[["Ensembl_transcript_ID"]] <- NA_character_
TxDb_transcripts_sel_df[["Source"]] <- "TxDb"
BioMart_transcripts_df[["Source"]] <- "BioMart"
gencode_transcripts_df[["Source"]] <- "GENCODE"
gencode_transcripts_df <- gencode_transcripts_df[, c(transcript_columns, "Source")]

transcript_locations_df <- rbind.data.frame(
  TxDb_transcripts_sel_df, BioMart_transcripts_df, gencode_transcripts_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)

transcript_locations_df <- ResolveDuplicateFeatures(transcript_locations_df)






# Create a merged data frame of exon coordinates --------------------------

TxDb_exons_sel_df <- TxDb_exons_df[, TxDb_columns]
gencode_exons_df <- gencode_df[gencode_df[["Entry_type"]] == "exon",
                               gencode_columns
                               ]

row.names(gencode_exons_df) <- NULL

exon_columns <- c(
  "Entrez_ID", "Ensembl_gene_ID", "Ensembl_transcript_ID", "Exon_ID",
  "Original_symbol",
  location_columns
)

colnames(gencode_exons_df) <- exon_columns

colnames(TxDb_exons_sel_df) <- c("Entrez_ID", location_columns)

TxDb_exons_sel_df[["Original_symbol"]] <- NA_character_
TxDb_exons_sel_df[["Ensembl_gene_ID"]] <- NA_character_
TxDb_exons_sel_df[["Ensembl_transcript_ID"]] <- NA_character_
TxDb_exons_sel_df[["Exon_ID"]] <- NA_character_
TxDb_exons_sel_df[["Source"]] <- "TxDb"
gencode_exons_df[["Source"]] <- "GENCODE"
gencode_exons_df <- gencode_exons_df[, c(exon_columns, "Source")]

exon_locations_df <- rbind.data.frame(
  TxDb_exons_sel_df, gencode_exons_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)


exon_locations_df <- ResolveDuplicateFeatures(exon_locations_df)





# Create a merged data frame of CDS coordinates ---------------------------

TxDb_CDS_sel_df <- TxDb_exons_df[, TxDb_columns]
gencode_CDS_df <- gencode_df[gencode_df[["Entry_type"]] == "CDS",
                             gencode_columns
                             ]

row.names(gencode_CDS_df) <- NULL

CDS_columns <- c(
  "Entrez_ID", "Ensembl_gene_ID", "Ensembl_transcript_ID", "Exon_ID",
  "Original_symbol",
  location_columns
)

colnames(gencode_CDS_df) <- CDS_columns

colnames(TxDb_CDS_sel_df) <- c("Entrez_ID", location_columns)

TxDb_CDS_sel_df[["Original_symbol"]] <- NA_character_
TxDb_CDS_sel_df[["Ensembl_gene_ID"]] <- NA_character_
TxDb_CDS_sel_df[["Ensembl_transcript_ID"]] <- NA_character_
TxDb_CDS_sel_df[["Exon_ID"]] <- NA_character_
TxDb_CDS_sel_df[["Source"]] <- "TxDb"
gencode_CDS_df[["Source"]] <- "GENCODE"
gencode_CDS_df <- gencode_CDS_df[, c(CDS_columns, "Source")]

CDS_locations_df <- rbind.data.frame(
  TxDb_CDS_sel_df, gencode_CDS_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)


CDS_locations_df <- ResolveDuplicateFeatures(CDS_locations_df)






# Save data ---------------------------------------------------------------

save(list = c("gene_locations_df", "transcript_locations_df",
              "exon_locations_df", "CDS_locations_df"
              ),
     file = file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData")
     )






