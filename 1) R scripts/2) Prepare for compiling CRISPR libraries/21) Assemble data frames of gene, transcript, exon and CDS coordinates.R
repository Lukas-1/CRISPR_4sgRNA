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
                                     are_identical_entrez <- (entrez_location_IDs[x] == entrez_location_IDs[x[[i]]]) %in% TRUE
                                     are_identical_ENSG <- (ENSG_location_IDs[x] == ENSG_location_IDs[x[[i]]]) %in% TRUE
                                     are_identical <- are_identical_entrez | are_identical_ENSG
                                     group_IDs_vec[are_identical] <- current_index
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











# Create a merged data frame of gene coordinates --------------------------

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



new_order <- order(match(gene_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   gene_locations_df[["Chromosome"]],
                   match(gene_locations_df[["Strand"]], c("+", "-")),
                   gene_locations_df[["Start"]],
                   gene_locations_df[["End"]],
                   as.integer(gene_locations_df[["Entrez_ID"]]),
                   gene_locations_df[["Ensembl_gene_ID"]]
                   )

gene_locations_df <- gene_locations_df[new_order, ]

location_strings <- do.call(paste, c(gene_locations_df[location_columns], sep = "_"))
entrez_location_IDs <- ifelse(is.na(gene_locations_df[["Entrez_ID"]]),
                              NA_character_,
                              paste0(gene_locations_df[["Entrez_ID"]], "_", location_strings)
                              )
ENSG_location_IDs <- ifelse(is.na(gene_locations_df[["Ensembl_gene_ID"]]),
                            NA_character_,
                            paste0(gene_locations_df[["Ensembl_gene_ID"]], "_", location_strings)
                            )

indices_list <- split(seq_len(nrow(gene_locations_df)),
                      location_strings
                      )
are_duplicated <- lengths(indices_list) > 1
all_have_entrezs <- vapply(indices_list,
                           function(x) !(anyNA(gene_locations_df[["Entrez_ID"]][x])),
                           logical(1)
                           )
all_have_ENSGs <- vapply(indices_list,
                         function(x) !(anyNA(gene_locations_df[["Ensembl_gene_ID"]][x])),
                         logical(1)
                         )

num_entrezs <- vapply(indices_list,
                      function(x) sum(!(is.na(unique(gene_locations_df[["Entrez_ID"]][x])))),
                      integer(1)
                      )
num_ENSGs <- vapply(indices_list,
                    function(x) sum(!(is.na(unique(gene_locations_df[["Ensembl_gene_ID"]][x])))),
                    integer(1)
                    )

group_IDs_vec <- rep(NA_integer_, nrow(gene_locations_df))
current_index <- 1L
for (i in seq_along(group_IDs_vec)) {
  if (!(is.na(group_IDs_vec[[i]]))) {
    next
  } else {
    are_identical_entrez <- (entrez_location_IDs == entrez_location_IDs[[i]]) %in% TRUE
    are_identical_ENSG <- (ENSG_location_IDs == ENSG_location_IDs[[i]]) %in% TRUE
    are_identical <- are_identical_entrez | are_identical_ENSG
    group_IDs_vec[are_identical] <- current_index
    current_index <- current_index + 1L
  }
}

sources_order <- c("TxDb", "BioMart", "Gencode")
sources_vec <- tapply(gene_locations_df[["Source"]],
                      group_IDs_vec,
                      function(x) {
                        unique_vec <- unique(x)
                        unique_vec <- unique_vec[order(match(unique_vec, sources_order))]
                        paste0(unique_vec, collapse = ", ")
                      })
entrezs_vec <- tapply(gene_locations_df[["Entrez_ID"]],
                      group_IDs_vec,
                      function(x) {
                        if (all(is.na(x))) {
                          return(NA_character_)
                        } else {
                          paste0(unique(x[!(is.na(x))]), collapse = ", ")
                        }
                      })
ENSG_vec <- tapply(gene_locations_df[["Ensembl_gene_ID"]],
                   group_IDs_vec,
                   function(x) {
                     if (all(is.na(x))) {
                       return(NA_character_)
                     } else {
                       unique(x[!(is.na(x))])
                     }
                   })


gene_locations_df <- gene_locations_df[!(duplicated(group_IDs_vec)), ]


stopifnot(!(any(ENSG_vec != gene_locations_df[["Ensembl_gene_ID"]]) %in% TRUE))

are_different <- ((entrezs_vec != gene_locations_df[["Entrez_ID"]]) %in% TRUE)
data.frame("Old_entrez" = gene_locations_df[["Entrez_ID"]][are_different],
           "New_entrez" = entrezs_vec[are_different]
           )

gene_locations_df[["Entrez_ID"]] <- entrezs_vec
gene_locations_df[["Ensembl_gene_ID"]] <- ENSG_vec
gene_locations_df[["Source"]] <- sources_vec

new_order <- order(GetMinEntrez(gene_locations_df[["Entrez_ID"]]),
                   gene_locations_df[["Ensembl_gene_ID"]],
                   match(gene_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   gene_locations_df[["Chromosome"]],
                   match(gene_locations_df[["Strand"]], c("+", "-")),
                   gene_locations_df[["Start"]],
                   gene_locations_df[["End"]]
                   )
gene_locations_df <- gene_locations_df[new_order, ]
row.names(gene_locations_df) <- NULL





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


new_order <- order(match(transcript_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   transcript_locations_df[["Chromosome"]],
                   match(transcript_locations_df[["Strand"]], c("+", "-")),
                   transcript_locations_df[["Start"]],
                   transcript_locations_df[["End"]],
                   GetMinEntrez(transcript_locations_df[["Entrez_ID"]]),
                   transcript_locations_df[["Ensembl_gene_ID"]],
                   transcript_locations_df[["Ensembl_transcript_ID"]]
                   )

transcript_locations_df <- transcript_locations_df[new_order, ]

location_strings <- do.call(paste, c(transcript_locations_df[location_columns], sep = "_"))
entrez_location_IDs <- ifelse(is.na(transcript_locations_df[["Entrez_ID"]]),
                              NA_character_,
                              paste0(transcript_locations_df[["Entrez_ID"]], "_", location_strings)
                              )
ENSG_location_IDs <- ifelse(is.na(transcript_locations_df[["Ensembl_gene_ID"]]),
                            NA_character_,
                            paste0(transcript_locations_df[["Ensembl_gene_ID"]], "_", location_strings)
                            )

indices_list <- split(seq_len(nrow(transcript_locations_df)),
                      location_strings
                      )
are_duplicated <- lengths(indices_list) > 1
all_have_entrezs <- vapply(indices_list,
                           function(x) !(anyNA(transcript_locations_df[["Entrez_ID"]][x])),
                           logical(1)
                           )
all_have_ENSGs <- vapply(indices_list,
                         function(x) !(anyNA(transcript_locations_df[["Ensembl_gene_ID"]][x])),
                         logical(1)
                         )

num_entrezs <- vapply(indices_list,
                      function(x) sum(!(is.na(unique(transcript_locations_df[["Entrez_ID"]][x])))),
                      integer(1)
                      )
num_ENSGs <- vapply(indices_list,
                    function(x) sum(!(is.na(unique(transcript_locations_df[["Ensembl_gene_ID"]][x])))),
                    integer(1)
                    )

group_IDs_vec <- rep(NA_integer_, nrow(transcript_locations_df))
current_index <- 1L
for (i in seq_along(group_IDs_vec)) {
  if (!(is.na(group_IDs_vec[[i]]))) {
    next
  } else {
    are_identical_entrez <- (entrez_location_IDs == entrez_location_IDs[[i]]) %in% TRUE
    are_identical_ENSG <- (ENSG_location_IDs == ENSG_location_IDs[[i]]) %in% TRUE
    are_identical <- are_identical_entrez | are_identical_ENSG
    group_IDs_vec[are_identical] <- current_index
    current_index <- current_index + 1L
  }
}

sources_order <- c("TxDb", "BioMart", "Gencode")
sources_vec <- tapply(transcript_locations_df[["Source"]],
                      group_IDs_vec,
                      function(x) {
                        unique_vec <- unique(x)
                        unique_vec <- unique_vec[order(match(unique_vec, sources_order))]
                        paste0(unique_vec, collapse = ", ")
                      })
entrezs_vec <- tapply(transcript_locations_df[["Entrez_ID"]],
                      group_IDs_vec,
                      function(x) {
                        if (all(is.na(x))) {
                          return(NA_character_)
                        } else {
                          paste0(unique(x[!(is.na(x))]), collapse = ", ")
                        }
                      })
ENSG_vec <- tapply(transcript_locations_df[["Ensembl_gene_ID"]],
                   group_IDs_vec,
                   function(x) {
                     if (all(is.na(x))) {
                       return(NA_character_)
                     } else {
                       paste0(unique(x[!(is.na(x))]), collapse = ", ")
                     }
                   })





full_transcript_locations_df <- transcript_locations_df


full_transcript_locations_df[full_transcript_locations_df[["Entrez_ID"]] %in% "221468", ]


are_duplicated <- duplicated(group_IDs_vec)
transcript_locations_df <- transcript_locations_df[!(are_duplicated), ]

are_problematic <- (ENSG_vec != transcript_locations_df[["Ensembl_gene_ID"]]) %in% TRUE

message(paste0("Some locations (", sum(are_problematic), " in total) ",
               " mapped to more than one Ensembl gene ID."
               )
        )
data.frame(transcript_locations_df[are_problematic, colnames(transcript_locations_df) != "Ensembl_gene_ID"],
           "Old_ENSG" = transcript_locations_df[["Ensembl_gene_ID"]][are_problematic],
           "New_ENSG" = ENSG_vec[are_problematic]
           )
message("All transcript entries for problematic entries are listed below:")
problematic_locations <- location_strings[!(are_duplicated)][are_problematic]
print(full_transcript_locations_df[location_strings %in% problematic_locations, ])


are_different <- ((entrezs_vec != transcript_locations_df[["Entrez_ID"]]) %in% TRUE)
data.frame("Old_entrez" = transcript_locations_df[["Entrez_ID"]][are_different],
           "New_entrez" = entrezs_vec[are_different]
           )


transcript_locations_df[["Entrez_ID"]] <- entrezs_vec
transcript_locations_df[["Ensembl_gene_ID"]] <- ENSG_vec
transcript_locations_df[["Source"]] <- sources_vec

new_order <- order(GetMinEntrez(transcript_locations_df[["Entrez_ID"]]),
                   transcript_locations_df[["Ensembl_gene_ID"]],
                   match(transcript_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   transcript_locations_df[["Chromosome"]],
                   match(transcript_locations_df[["Strand"]], c("+", "-")),
                   transcript_locations_df[["Start"]],
                   transcript_locations_df[["End"]]
                   )
transcript_locations_df <- transcript_locations_df[new_order, ]
row.names(transcript_locations_df) <- NULL




gooo


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


new_order <- order(match(exon_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   exon_locations_df[["Chromosome"]],
                   match(exon_locations_df[["Strand"]], c("+", "-")),
                   exon_locations_df[["Start"]],
                   exon_locations_df[["End"]],
                   GetMinEntrez(exon_locations_df[["Entrez_ID"]]),
                   exon_locations_df[["Ensembl_gene_ID"]],
                   exon_locations_df[["Ensembl_transcript_ID"]]
                   )

exon_locations_df <- exon_locations_df[new_order, ]

location_strings <- do.call(paste, c(exon_locations_df[location_columns], sep = "_"))
entrez_location_IDs <- ifelse(is.na(exon_locations_df[["Entrez_ID"]]),
                              NA_character_,
                              paste0(exon_locations_df[["Entrez_ID"]], "_", location_strings)
                              )
ENSG_location_IDs <- ifelse(is.na(exon_locations_df[["Ensembl_gene_ID"]]),
                            NA_character_,
                            paste0(exon_locations_df[["Ensembl_gene_ID"]], "_", location_strings)
                            )

indices_list <- split(seq_len(nrow(exon_locations_df)),
                      location_strings
                      )
are_duplicated <- lengths(indices_list) > 1
all_have_entrezs <- vapply(indices_list,
                           function(x) !(anyNA(exon_locations_df[["Entrez_ID"]][x])),
                           logical(1)
                           )
all_have_ENSGs <- vapply(indices_list,
                         function(x) !(anyNA(exon_locations_df[["Ensembl_gene_ID"]][x])),
                         logical(1)
                         )

num_entrezs <- vapply(indices_list,
                      function(x) sum(!(is.na(unique(exon_locations_df[["Entrez_ID"]][x])))),
                      integer(1)
                      )
num_ENSGs <- vapply(indices_list,
                    function(x) sum(!(is.na(unique(exon_locations_df[["Ensembl_gene_ID"]][x])))),
                    integer(1)
                    )

group_IDs_vec <- rep(NA_integer_, nrow(exon_locations_df))
current_index <- 1L

num_occurrences_vec <- table(location_strings)[location_strings]


location_fac <- factor(location_strings, levels = unique(location_strings))

group_indices_list <- tapply(seq_along(location_fac),
                             location_fac,
                             function(x) {
                               group_IDs_vec <- rep(NA_integer_, length(x))
                               current_index <- 1L
                               for (i in seq_along(x)) {
                                 if (!(is.na(group_IDs_vec[[i]]))) {
                                   next
                                 } else {
                                   are_identical_entrez <- (entrez_location_IDs[x] == entrez_location_IDs[x[[i]]]) %in% TRUE
                                   are_identical_ENSG <- (ENSG_location_IDs[x] == ENSG_location_IDs[x[[i]]]) %in% TRUE
                                   are_identical <- are_identical_entrez | are_identical_ENSG
                                   group_IDs_vec[are_identical] <- current_index
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
sources_vec <- tapply(exon_locations_df[["Source"]],
                      integer_group_IDs,
                      function(x) {
                        unique_vec <- unique(x)
                        unique_vec <- unique_vec[order(match(unique_vec, sources_order))]
                        paste0(unique_vec, collapse = ", ")
                      })
entrezs_vec <- tapply(exon_locations_df[["Entrez_ID"]],
                      integer_group_IDs,
                      function(x) {
                        if (all(is.na(x))) {
                          return(NA_character_)
                        } else {
                          paste0(unique(x[!(is.na(x))]), collapse = ", ")
                        }
                      })
ENSG_vec <- tapply(exon_locations_df[["Ensembl_gene_ID"]],
                   integer_group_IDs,
                   function(x) {
                     if (all(is.na(x))) {
                       return(NA_character_)
                     } else {
                       paste0(unique(x[!(is.na(x))]), collapse = ", ")
                     }
                   })





full_exon_locations_df <- exon_locations_df

are_duplicated <- duplicated(integer_group_IDs)
exon_locations_df <- exon_locations_df[!(are_duplicated), ]

are_problematic <- (ENSG_vec != exon_locations_df[["Ensembl_gene_ID"]]) %in% TRUE

message(paste0("Some locations (", sum(are_problematic), " in total) ",
               " mapped to more than one Ensembl gene ID."
               )
        )
data.frame(exon_locations_df[are_problematic, colnames(exon_locations_df) != "Ensembl_gene_ID"],
           "Old_ENSG" = exon_locations_df[["Ensembl_gene_ID"]][are_problematic],
           "New_ENSG" = ENSG_vec[are_problematic]
           )
message("All transcript entries for problematic entries are listed below:")
problematic_locations <- location_strings[!(are_duplicated)][are_problematic]
print(full_exon_locations_df[location_strings %in% problematic_locations, ])


are_different <- ((entrezs_vec != exon_locations_df[["Entrez_ID"]]) %in% TRUE)
data.frame("Old_entrez" = exon_locations_df[["Entrez_ID"]][are_different],
           "New_entrez" = entrezs_vec[are_different]
           )


exon_locations_df[["Entrez_ID"]] <- entrezs_vec
exon_locations_df[["Ensembl_gene_ID"]] <- ENSG_vec
exon_locations_df[["Source"]] <- sources_vec

new_order <- order(GetMinEntrez(exon_locations_df[["Entrez_ID"]]),
                   exon_locations_df[["Ensembl_gene_ID"]],
                   match(exon_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   exon_locations_df[["Chromosome"]],
                   match(exon_locations_df[["Strand"]], c("+", "-")),
                   exon_locations_df[["Start"]],
                   exon_locations_df[["End"]]
                   )
exon_locations_df <- exon_locations_df[new_order, ]
row.names(exon_locations_df) <- NULL





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


new_order <- order(match(CDS_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   CDS_locations_df[["Chromosome"]],
                   match(CDS_locations_df[["Strand"]], c("+", "-")),
                   CDS_locations_df[["Start"]],
                   CDS_locations_df[["End"]],
                   GetMinEntrez(CDS_locations_df[["Entrez_ID"]]),
                   CDS_locations_df[["Ensembl_gene_ID"]],
                   CDS_locations_df[["Ensembl_transcript_ID"]]
                   )

CDS_locations_df <- CDS_locations_df[new_order, ]

location_strings <- do.call(paste, c(CDS_locations_df[location_columns], sep = "_"))
entrez_location_IDs <- ifelse(is.na(CDS_locations_df[["Entrez_ID"]]),
                              NA_character_,
                              paste0(CDS_locations_df[["Entrez_ID"]], "_", location_strings)
                              )
ENSG_location_IDs <- ifelse(is.na(CDS_locations_df[["Ensembl_gene_ID"]]),
                            NA_character_,
                            paste0(CDS_locations_df[["Ensembl_gene_ID"]], "_", location_strings)
                            )

indices_list <- split(seq_len(nrow(CDS_locations_df)),
                      location_strings
                      )
# are_duplicated <- lengths(indices_list) > 1
# all_have_entrezs <- vapply(indices_list,
#                            function(x) !(anyNA(CDS_locations_df[["Entrez_ID"]][x])),
#                            logical(1)
#                            )
# all_have_ENSGs <- vapply(indices_list,
#                          function(x) !(anyNA(CDS_locations_df[["Ensembl_gene_ID"]][x])),
#                          logical(1)
#                          )
#
# num_entrezs <- vapply(indices_list,
#                       function(x) sum(!(is.na(unique(CDS_locations_df[["Entrez_ID"]][x])))),
#                       integer(1)
#                       )
# num_ENSGs <- vapply(indices_list,
#                     function(x) sum(!(is.na(unique(CDS_locations_df[["Ensembl_gene_ID"]][x])))),
#                     integer(1)
#                     )

group_IDs_vec <- rep(NA_integer_, nrow(CDS_locations_df))
current_index <- 1L
for (i in seq_along(group_IDs_vec)) {
  if (!(is.na(group_IDs_vec[[i]]))) {
    next
  } else {
    are_identical_entrez <- (entrez_location_IDs == entrez_location_IDs[[i]]) %in% TRUE
    are_identical_ENSG <- (ENSG_location_IDs == ENSG_location_IDs[[i]]) %in% TRUE
    are_identical <- are_identical_entrez | are_identical_ENSG
    group_IDs_vec[are_identical] <- current_index
    current_index <- current_index + 1L
  }
}

sources_order <- c("TxDb", "BioMart", "Gencode")
sources_vec <- tapply(CDS_locations_df[["Source"]],
                      group_IDs_vec,
                      function(x) {
                        unique_vec <- unique(x)
                        unique_vec <- unique_vec[order(match(unique_vec, sources_order))]
                        paste0(unique_vec, collapse = ", ")
                      })
entrezs_vec <- tapply(CDS_locations_df[["Entrez_ID"]],
                      group_IDs_vec,
                      function(x) {
                        if (all(is.na(x))) {
                          return(NA_character_)
                        } else {
                          paste0(unique(x[!(is.na(x))]), collapse = ", ")
                        }
                      })
ENSG_vec <- tapply(CDS_locations_df[["Ensembl_gene_ID"]],
                   group_IDs_vec,
                   function(x) {
                     if (all(is.na(x))) {
                       return(NA_character_)
                     } else {
                       paste0(unique(x[!(is.na(x))]), collapse = ", ")
                     }
                   })


full_CDS_locations_df <- CDS_locations_df

are_duplicated <- duplicated(group_IDs_vec)
CDS_locations_df <- CDS_locations_df[!(are_duplicated), ]

are_problematic <- (ENSG_vec != CDS_locations_df[["Ensembl_gene_ID"]]) %in% TRUE

message(paste0("Some locations (", sum(are_problematic), " in total) ",
               " mapped to more than one Ensembl gene ID."
               )
        )
data.frame(CDS_locations_df[are_problematic, colnames(CDS_locations_df) != "Ensembl_gene_ID"],
           "Old_ENSG" = CDS_locations_df[["Ensembl_gene_ID"]][are_problematic],
           "New_ENSG" = ENSG_vec[are_problematic]
           )
message("All transcript entries for problematic entries are listed below:")
problematic_locations <- location_strings[!(are_duplicated)][are_problematic]
print(full_CDS_locations_df[location_strings %in% problematic_locations, ])


are_different <- ((entrezs_vec != CDS_locations_df[["Entrez_ID"]]) %in% TRUE)
data.frame("Old_entrez" = CDS_locations_df[["Entrez_ID"]][are_different],
           "New_entrez" = entrezs_vec[are_different]
           )


CDS_locations_df[["Entrez_ID"]] <- entrezs_vec
CDS_locations_df[["Ensembl_gene_ID"]] <- ENSG_vec
CDS_locations_df[["Source"]] <- sources_vec

new_order <- order(GetMinEntrez(CDS_locations_df[["Entrez_ID"]]),
                   CDS_locations_df[["Ensembl_gene_ID"]],
                   match(CDS_locations_df[["Chromosome"]],
                         c(paste0("chr", 1:23), c("Y", "Y", "M"))
                         ),
                   CDS_locations_df[["Chromosome"]],
                   match(CDS_locations_df[["Strand"]], c("+", "-")),
                   CDS_locations_df[["Start"]],
                   CDS_locations_df[["End"]]
                   )
CDS_locations_df <- CDS_locations_df[new_order, ]
row.names(CDS_locations_df) <- NULL


