### 4th August 2020 ###



# Import packages and source code -----------------------------------------

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("biomaRt")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder
source(file.path(general_functions_directory, "29) Determining gene types.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(general_RData_directory, "18) Process the annotations from GENCODE.RData"))
load(file.path(general_RData_directory, "19) Compile the information on gene type.RData"))






# Define functions --------------------------------------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End")
all_chromosomes <- paste0("chr", c(1:23, "Y", "Y", "M"))


UniqueNonNA <- function(input_vec) {
  unique(input_vec[!(is.na(input_vec))])
}



ResolveDuplicateFeaturesTwice <- function(locations_df) {
  locations_df <- ResolveDuplicateFeatures(locations_df)
  num_rows_first <- nrow(locations_df)
  locations_df <- ResolveDuplicateFeatures(locations_df)
  num_2nd_pass <- num_rows_first - nrow(locations_df)
  message(paste0(num_2nd_pass, " additional duplicates were resolved ",
                 "in a second pass."
                 ))
  locations_df[["Duplicate_category"]] <- AnnotateRemainingDuplicates(locations_df)
  return(locations_df)
}


ResolveDuplicateFeatures <- function(locations_df) {

  if ("Gene_symbols" %in% names(locations_df)) {
    symbol_column <- "Gene_symbols"
    entrez_column <- "Entrez_IDs"
    ENSG_column <- "Ensembl_gene_IDs"
  } else {
    symbol_column <- "Original_symbol"
    entrez_column <- "Entrez_ID"
    ENSG_column <- "Ensembl_gene_ID"
  }

  # Requires the CollapseNonNA function

  if ("Ensembl_transcript_ID" %in% names(locations_df)) {
    transcript_vec <- locations_df[["Ensembl_transcript_ID"]]
  } else {
    transcript_vec <- rep(NA, nrow(locations_df))
  }

  assign("delete_locations_df", locations_df, envir = globalenv())
  assign("delete_entrez_column", entrez_column, envir = globalenv())

  new_order <- order(match(locations_df[["Chromosome"]], all_chromosomes),
                     locations_df[["Chromosome"]],
                     match(locations_df[["Strand"]], c("+", "-")),
                     locations_df[["Start"]],
                     locations_df[["End"]],
                     GetMinEntrez(locations_df[[entrez_column]]),
                     locations_df[[ENSG_column]],
                     transcript_vec
                     )

  locations_df <- locations_df[new_order, ]

  location_strings <- do.call(paste, c(locations_df[location_columns], sep = "_"))
  entrez_location_IDs <- ifelse(is.na(locations_df[[entrez_column]]),
                                NA_character_,
                                paste0(locations_df[[entrez_column]], "_", location_strings)
                                )
  ENSG_location_IDs <- ifelse(is.na(locations_df[[ENSG_column]]),
                              NA_character_,
                              paste0(locations_df[[ENSG_column]], "_", location_strings)
                              )

  symbol_location_IDs <- ifelse(is.na(locations_df[[symbol_column]]),
                                NA_character_,
                                paste0(locations_df[[symbol_column]], "_", location_strings)
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
                                     this_ENSG   <- ENSG_location_IDs[x[[i]]]
                                     this_symbol <- symbol_location_IDs[x[[i]]]

                                     are_this_entrez <- (entrez_location_IDs[x] == this_entrez) %in% TRUE
                                     are_this_ENSG   <- (ENSG_location_IDs[x] == this_ENSG) %in% TRUE
                                     are_this_symbol <- (symbol_location_IDs[x] == this_symbol) %in% TRUE

                                      all_entrezs <- UniqueNonNA(c(entrez_location_IDs[x][are_this_ENSG], this_entrez))
                                     all_ENSGs   <- UniqueNonNA(c(ENSG_location_IDs[x][are_this_entrez], this_ENSG))

                                     are_shared_entrez <- entrez_location_IDs[x] %in% all_entrezs
                                     are_shared_ENSG   <- ENSG_location_IDs[x] %in% all_ENSGs

                                     if ((length(UniqueNonNA(entrez_location_IDs[x[are_this_symbol]])) == 1) &&
                                         (length(UniqueNonNA(ENSG_location_IDs[x[are_this_symbol]])) == 1)
                                         ) {
                                       are_shared_symbol <- are_this_symbol
                                     } else {
                                       are_shared_symbol <- rep(FALSE, length(x))
                                     }

                                     are_shared <- are_shared_entrez | are_shared_ENSG | are_shared_symbol

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

  sources_order <- c("TxDb", "BioMart", "Gencode")
  sources_vec <- tapply(locations_df[["Source"]],
                        integer_group_IDs,
                        function(x) CollapseNonNA(x, sources_order)
                        )

  entrezs_vec <- tapply(locations_df[[entrez_column]],     integer_group_IDs, CollapseNonNA)
  ENSG_vec    <- tapply(locations_df[[ENSG_column]],       integer_group_IDs, CollapseNonNA)
  symbols_vec <- tapply(locations_df[["Original_symbol"]], integer_group_IDs, CollapseNonNA)

  full_locations_df <- locations_df

  are_duplicated <- duplicated(integer_group_IDs)
  locations_df <- locations_df[!(are_duplicated), ]

  are_problematic <- (ENSG_vec != locations_df[[ENSG_column]]) %in% TRUE

  if (any(are_problematic)) {
    message(paste0("Some locations (", sum(are_problematic), " in total)",
                   " mapped to more than one Ensembl gene ID."
                   )
            )
    data.frame(locations_df[are_problematic, names(locations_df) != ENSG_column],
               "Old_ENSG" = locations_df[[ENSG_column]][are_problematic],
               "New_ENSG" = ENSG_vec[are_problematic],
               stringsAsFactors = FALSE
               )
    message("All transcript entries for problematic entries are listed below:")
    problematic_locations <- location_strings[!(are_duplicated)][are_problematic]
    print(full_locations_df[location_strings %in% problematic_locations, ])
  } else {
    message("No ambiguous Ensembl IDs (for the same locations) were found.")
  }

  are_different <- ((entrezs_vec != locations_df[[entrez_column]]) %in% TRUE)
  data.frame("Old_entrez" = locations_df[["Entrez_ID"]][are_different],
             "New_entrez" = entrezs_vec[are_different]
             )

  locations_df[["Entrez_IDs"]] <- entrezs_vec
  locations_df[["Ensembl_gene_IDs"]] <- ENSG_vec
  locations_df[["Source"]] <- sources_vec
  locations_df[["Original_symbol"]] <- symbols_vec

  new_order <- order(GetMinEntrez(locations_df[["Entrez_IDs"]]),
                     locations_df[["Ensembl_gene_IDs"]],
                     match(locations_df[["Chromosome"]],
                           c(paste0("chr", 1:23), c("Y", "Y", "M"))
                           ),
                     locations_df[["Chromosome"]],
                     match(locations_df[["Strand"]], c("+", "-")),
                     locations_df[["Start"]],
                     locations_df[["End"]]
                     )
  locations_df <- locations_df[new_order, ]

  locations_df[["Gene_types"]] <- GetGeneTypes(locations_df[[entrez_column]],
                                               locations_df[[ENSG_column]]
                                               )

  ENSG_to_symbol_vec <- TranslateENSGtoSymbol(locations_df[["Ensembl_gene_IDs"]])
  entrez_to_symbol_vec <- MapToEntrezs(entrez_IDs_vec = locations_df[["Entrez_IDs"]])[["Gene_symbol"]]
  locations_df[["Gene_symbols"]] <- ifelse(is.na(entrez_to_symbol_vec),
                                           ENSG_to_symbol_vec,
                                           entrez_to_symbol_vec
                                           )


  ncbi_vec <- ifelse(is.na(locations_df[["Entrez_IDs"]]),
                     NA_character_,
                     paste0("ncbi:", locations_df[["Entrez_IDs"]])
                     )
  ncbi_vec <- gsub(", ", ", ncbi:", ncbi_vec, fixed = TRUE)
  locations_df[["Gene_IDs"]] <- ifelse(is.na(ncbi_vec),
                                       locations_df[["Ensembl_gene_IDs"]],
                                       ncbi_vec
                                       )

  all_columns <- c("Gene_IDs", "Entrez_IDs", "Ensembl_gene_IDs",
                   "Gene_symbols", "Original_symbol",
                   "Gene_types", "Source",
                   "Ensembl_transcript_ID", "Exon_ID",
                   "Chromosome", "Strand", "Start", "End"
                   )
  columns_reordered <- intersect(all_columns, names(locations_df))
  locations_df <- locations_df[, columns_reordered]
  row.names(locations_df) <- NULL
  return(locations_df)
}


AnnotateRemainingDuplicates <- function(locations_df) {

  location_strings <- do.call(paste, c(locations_df[location_columns], sep = "_"))
  duplicated_locations <- unique(location_strings[duplicated(location_strings)])

  are_duplicated <- location_strings %in% duplicated_locations

  results_vec <- rep(NA, nrow(locations_df))

  for (location in duplicated_locations) {
    are_this_location <- location_strings[are_duplicated] %in% location
    sources <- locations_df[["Source"]][are_duplicated][are_this_location]
    are_TxDb <- sources == "TxDb"
    if (all(are_TxDb)) {
      replacement_vec <- "Other"
    } else {
      replacement_vec <- ifelse(are_TxDb, "Only TxDb", "Other")
    }
    results_vec[are_duplicated][are_this_location] <- replacement_vec
  }
  return(results_vec)
}





FixTxDb_columns <- function(TxDb_df, example_df) {
  names(TxDb_df) <- c("Entrez_ID", location_columns)
  TxDb_df[["Original_symbol"]] <- NA_character_
  TxDb_df[["Ensembl_gene_ID"]] <- NA_character_
  TxDb_df[["Source"]] <- "TxDb"

  if ("Ensembl_transcript_ID" %in% names(example_df)) {
    TxDb_df[["Ensembl_transcript_ID"]] <- NA_character_
  }
  if ("Exon_ID" %in% names(example_df)) {
    TxDb_df[["Exon_ID"]] <- NA_character_
  }
  TxDb_df <- TxDb_df[, union(gene_columns, names(TxDb_df))]
  return(TxDb_df)
}


GetUniqueIDs <- function(char_vec) {
  unique(unlist(strsplit(char_vec, ", ", fixed = TRUE)))
}


AreMissingInVec2 <- function(char_vec_1, char_vec_2) {
  unique_vec_1 <- GetUniqueIDs(char_vec_1)
  unique_vec_2 <- GetUniqueIDs(char_vec_2)
  missing_IDs <- setdiff(unique_vec_1, unique_vec_2)
  splits_list <- strsplit(char_vec_1, ", ", fixed = TRUE)
  have_missing_IDs <- vapply(splits_list,
                             function(x) any(x %in% missing_IDs),
                             logical(1)
                             )
  return(have_missing_IDs)
}



ENSGVectoSymbol <- function(ENSG_vec) {

  # use_host <- "www.ensembl.org"
  use_host <- "http://apr2020.archive.ensembl.org"
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = use_host))
  biomaRt_lookup_df <- getBM(
    filters    = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "external_gene_name"),
    values     = ENSG_vec,
    mart       = mart
  )
  biomaRt_package_matches <- match(ENSG_vec, biomaRt_lookup_df[["ensembl_gene_id"]])
  BioMart_export_matches  <- match(ENSG_vec, BioMart_ENSG_to_symbol_df[["ENSG"]])
  gencode_matches         <- match(ENSG_vec, gencode_ENSG_to_symbol_df[["Ensembl_gene_ID"]])

  results_df <- data.frame(
    "Ensembl_gene_ID"        = ENSG_vec,
    "biomaRt_package_symbol" = biomaRt_lookup_df[["external_gene_name"]][biomaRt_package_matches],
    "BioMart_export_symbol"  = BioMart_ENSG_to_symbol_df[["Gene_symbol"]][BioMart_export_matches],
    "GENCODE_symbol"         = gencode_ENSG_to_symbol_df[["Gene_symbol"]][gencode_matches],
    stringsAsFactors         = FALSE
  )
  return(results_df)
}



TranslateENSGtoSymbol <- function(comma_separated_ENSG_vec) {

  ENSG_splits <- strsplit(comma_separated_ENSG_vec, ", ", fixed = TRUE)
  expanded_df <- ExpandList(ENSG_splits)
  are_NA <- is.na(expanded_df[["Value"]])
  are_duplicated <- duplicated(expanded_df[["Value"]][!(are_NA)])

  unique_vec <- expanded_df[["Value"]][!(are_NA)][!(are_duplicated)]
  translated_df <- ENSGVectoSymbol(unique_vec)
  unique_ENSG_vec <- ifelse(is.na(translated_df[["biomaRt_package_symbol"]]),
                            ifelse(is.na(translated_df[["BioMart_export_symbol"]]),
                                   translated_df[["GENCODE_symbol"]],
                                   translated_df[["BioMart_export_symbol"]]
                                   ),
                            translated_df[["biomaRt_package_symbol"]]
                            )

  matches_vec <- match(expanded_df[["Value"]],
                       unique_vec
                       )

  expanded_df[["Symbol"]] <- unique_ENSG_vec[matches_vec]

  collapsed_symbols_vec <- tapply(expanded_df[["Symbol"]],
                                  expanded_df[["List_index"]],
                                  CollapseNonNA
                                  )
  return(collapsed_symbols_vec)
}





# Collect Ensembl ID-gene symbol mappings ---------------------------------

BioMart_ENSG_to_symbol_df <- unique(BioMart_tidied_df[, c("ENSG", "Gene_symbol")])
gencode_ENSG_to_symbol_df <- unique(gencode_df[, c("Ensembl_gene_ID", "Gene_symbol")])

stopifnot(!(anyDuplicated(BioMart_ENSG_to_symbol_df[["ENSG"]])))
stopifnot(!(anyDuplicated(gencode_ENSG_to_symbol_df[["Ensembl_gene_ID"]])))




# Create data frames using TxDb.Hsapiens.UCSC.hg38.knownGene --------------

MakeTxDbDf <- function(GRangesList_object) {
  results_df <- as.data.frame(GRangesList_object, stringsAsFactors = FALSE)
  results_df[["strand"]] <- as.character(results_df[["strand"]])
  results_df[["seqnames"]] <- as.character(results_df[["seqnames"]])
  return(results_df)
}

TxDb_genes_df <- MakeTxDbDf(genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                  single.strand.genes.only = FALSE
                                  ))
TxDb_transcripts_df <- MakeTxDbDf(transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene))
TxDb_exons_df <- MakeTxDbDf(exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))
TxDb_CDS_df <- MakeTxDbDf(cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))

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

exon_columns <- gencode_columns
exon_columns[exon_columns == "Gene_symbol"] <- "Original_symbol"

transcript_columns <- setdiff(exon_columns, "Exon_ID")
gene_columns <- setdiff(transcript_columns, "Ensembl_transcript_ID")






# Create a merged data frame of gene coordinates --------------------------

BioMart_genes_df <- unique(BioMart_tidied_df[, setdiff(BioMart_columns, "ENST")])
TxDb_genes_sel_df <- TxDb_genes_df[, TxDb_columns]
gencode_genes_df <- gencode_df[gencode_df[["Entry_type"]] == "gene",
                               setdiff(gencode_columns, c("Ensembl_transcript_ID", "Exon_ID"))
                               ]
row.names(BioMart_genes_df) <- NULL
row.names(gencode_genes_df) <- NULL

names(BioMart_genes_df) <- gene_columns
names(gencode_genes_df) <- gene_columns

TxDb_genes_sel_df <- FixTxDb_columns(TxDb_genes_sel_df, gencode_genes_df)

BioMart_genes_df[["Source"]] <- "BioMart"
gencode_genes_df[["Source"]] <- "GENCODE"

gencode_genes_df <- gencode_genes_df[, c(gene_columns, "Source")]

gene_locations_df <- rbind.data.frame(
  TxDb_genes_sel_df, BioMart_genes_df, gencode_genes_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)

gene_locations_df <- ResolveDuplicateFeaturesTwice(gene_locations_df)






# Create a merged data frame of transcript coordinates --------------------

BioMart_transcripts_df <- BioMart_tidied_df[, BioMart_columns]
TxDb_transcripts_sel_df <- TxDb_transcripts_df[, TxDb_columns]
gencode_transcripts_df <- gencode_df[gencode_df[["Entry_type"]] == "transcript",
                                     setdiff(gencode_columns, "Exon_ID"),
                                     ]

row.names(gencode_transcripts_df) <- NULL

names(BioMart_transcripts_df) <- transcript_columns
names(gencode_transcripts_df) <- transcript_columns

TxDb_transcripts_sel_df <- FixTxDb_columns(TxDb_transcripts_sel_df, gencode_transcripts_df)

BioMart_transcripts_df[["Source"]] <- "BioMart"
gencode_transcripts_df[["Source"]] <- "GENCODE"
gencode_transcripts_df <- gencode_transcripts_df[, c(transcript_columns, "Source")]

transcript_locations_df <- rbind.data.frame(
  TxDb_transcripts_sel_df, BioMart_transcripts_df, gencode_transcripts_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)

transcript_locations_df <- ResolveDuplicateFeaturesTwice(transcript_locations_df)





# Create a merged data frame of exon coordinates --------------------------

TxDb_exons_sel_df <- TxDb_exons_df[, TxDb_columns]
gencode_exons_df <- gencode_df[gencode_df[["Entry_type"]] == "exon",
                               gencode_columns
                               ]

row.names(gencode_exons_df) <- NULL

names(gencode_exons_df) <- exon_columns

TxDb_exons_sel_df <- FixTxDb_columns(TxDb_exons_sel_df, gencode_exons_df)

gencode_exons_df[["Source"]] <- "GENCODE"
gencode_exons_df <- gencode_exons_df[, c(exon_columns, "Source")]

exon_locations_df <- rbind.data.frame(
  TxDb_exons_sel_df, gencode_exons_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)

exon_locations_df <- ResolveDuplicateFeaturesTwice(exon_locations_df)





# Create a merged data frame of CDS coordinates ---------------------------

TxDb_CDS_sel_df <- TxDb_CDS_df[, TxDb_columns]
gencode_CDS_df <- gencode_df[gencode_df[["Entry_type"]] == "CDS",
                             gencode_columns
                             ]
row.names(gencode_CDS_df) <- NULL

names(gencode_CDS_df) <- exon_columns

TxDb_CDS_sel_df <- FixTxDb_columns(TxDb_CDS_sel_df, gencode_CDS_df)

gencode_CDS_df[["Source"]] <- "GENCODE"
gencode_CDS_df <- gencode_CDS_df[, c(exon_columns, "Source")]

CDS_locations_df <- rbind.data.frame(
  TxDb_CDS_sel_df, gencode_CDS_df,
  make.row.names = FALSE, stringsAsFactors = FALSE
)

CDS_locations_df <- ResolveDuplicateFeaturesTwice(CDS_locations_df)





# Explore the results -----------------------------------------------------

gene_unique_entrezs       <- GetUniqueIDs(gene_locations_df[["Entrez_IDs"]])
transcript_unique_entrezs <- GetUniqueIDs(transcript_locations_df[["Entrez_IDs"]])
exon_unique_entrezs       <- GetUniqueIDs(exon_locations_df[["Entrez_IDs"]])
CDS_unique_entrezs        <- GetUniqueIDs(CDS_locations_df[["Entrez_IDs"]])

length(gene_unique_entrezs)
length(transcript_unique_entrezs)
length(exon_unique_entrezs)
length(CDS_unique_entrezs)

have_missing_entrezs <- AreMissingInVec2(transcript_locations_df[["Entrez_IDs"]],
                                         exon_locations_df[["Entrez_IDs"]]
                                         )
have_multiple_entrezs <- grepl(", ", transcript_locations_df[["Entrez_IDs"]], fixed = TRUE)
are_on_standard_chromosome <- transcript_locations_df[["Chromosome"]] %in% all_chromosomes
stopifnot(!(any(have_missing_entrezs & are_on_standard_chromosome & !(have_multiple_entrezs))))
# ==> None of the Entrez IDs that are present in the transcript data,
#     but not in the exon data, seem particularly relevant.




# Collect CDSs (protein-coding genes) or exons (non-coding genes) ---------

have_missing_entrezs <- AreMissingInVec2(exon_locations_df[["Entrez_IDs"]],
                                         CDS_locations_df[["Entrez_IDs"]]
                                         )

have_missing_ENSG <- AreMissingInVec2(exon_locations_df[["Ensembl_gene_IDs"]],
                                      CDS_locations_df[["Ensembl_gene_IDs"]]
                                      )

are_to_keep <- have_missing_entrezs |
               (have_missing_ENSG & is.na(exon_locations_df[["Entrez_IDs"]]))

CDS_or_exon_locations_df <- rbind.data.frame(
  data.frame(
    CDS_locations_df,
    "Entry" = "CDS",
    stringsAsFactors = FALSE
  ),
  data.frame(
    exon_locations_df[are_to_keep, ],
    "Entry" = "exon",
    stringsAsFactors = FALSE
  ),
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

# Annotate duplicates that are present as both exons of ncRNAs and CDSs of protein-coding genes,
# while preserving the duplicate status based on exon annotations
CDS_or_exon_locations_df[["Duplicate_category"]] <- ifelse(is.na(CDS_or_exon_locations_df[["Duplicate_category"]]),
                                                           AnnotateRemainingDuplicates(CDS_or_exon_locations_df),
                                                           CDS_or_exon_locations_df[["Duplicate_category"]]
                                                           )




# Save data ---------------------------------------------------------------

save(list = c("gene_locations_df", "transcript_locations_df",
              "exon_locations_df", "CDS_locations_df",
              "CDS_or_exon_locations_df"
              ),
     file = file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData")
     )










