### 17th April 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))






# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")

HPA_directory            <- file.path(CRISPR_input_directory, "Human genome", "Human Protein Atlas")
surfaceome_directory     <- file.path(CRISPR_input_directory, "Human genome", "Surface proteome")

intermediate_directory   <- file.path(CRISPR_root_directory, "4) Intermediate files", "Annotation", "Surfaceome")
UniProt_input_directory  <- file.path(intermediate_directory, "1) Input to UniProt")
UniProt_output_directory <- file.path(intermediate_directory, "2) Output from UniProt")



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData"))




# Read in data ------------------------------------------------------------

full_HPA_df <- read.delim(file.path(HPA_directory, "proteinatlas_2021_03_09.tsv"),
                          stringsAsFactors = FALSE, check.names = FALSE
                          )

CD_df <- read.delim(file.path(surfaceome_directory, "HCDM_CD_list_2021-04-26.txt"),
                    strip.white = TRUE, header = FALSE,
                    stringsAsFactors = FALSE
                    )

ReadSurfaceomeExcel <- function(file_name, ...) {
  results_df <- data.frame(read_excel(file.path(surfaceome_directory, file_name),
                                      ...
                                      ),
                           stringsAsFactors = FALSE, check.names = FALSE
                           )
  names(results_df) <- sub("\\r\\n$", "", names(results_df))
  return(results_df)
}


mass_spec_A_df       <- ReadSurfaceomeExcel("2015 - A Mass Spectrometric-Derived Cell Surface Protein Atlas - File S1.xlsx",
                                            sheet = 1
                                            )

mass_spec_B_df       <- ReadSurfaceomeExcel("2015 - A Mass Spectrometric-Derived Cell Surface Protein Atlas - File S1.xlsx",
                                            sheet = 2
                                            )

silico_SURFY_df      <- ReadSurfaceomeExcel("2018 - The in silico human surfaceome - Dataset S1.xls",
                                            sheet = 3, skip = 1, guess_max = 10000
                                            )
silico_surfaceome_df <- ReadSurfaceomeExcel("2018 - The in silico human surfaceome - Dataset S1.xls",
                                            sheet = 7, skip = 1, guess_max = 10000
                                            )
silico_HeLa_df       <- ReadSurfaceomeExcel("2018 - The in silico human surfaceome - Dataset S1.xls",
                                            sheet = 8, skip = 1, guess_max = 10000
                                            )
silico_groups_df     <- ReadSurfaceomeExcel("2018 - The in silico human surfaceome - Dataset S1.xls",
                                            sheet = 10, skip = 1, guess_max = 10000
                                            )




# Export intermediate data ------------------------------------------------

stopifnot(all(silico_surfaceome_df[["UniProt name"]] %in% silico_SURFY_df[["UniProt name"]]))


WriteUniProtIDs <- function(IDs_vec, file_name) {
  write.table(IDs_vec,
              file = file.path(UniProt_input_directory, paste0(file_name, ".txt")),
              row.names = FALSE, quote = FALSE, col.names = FALSE
              )
}

WriteUniProtIDs(mass_spec_B_df[["ID_link"]],
                "2015 - A Mass Spectrometric-Derived Cell Surface Protein Atlas - UniProt accessions"
                )

WriteUniProtIDs(silico_SURFY_df[["UniProt name"]],
                "2018 - The in silico human surfaceome - SURFY - UniProt names"
                )

WriteUniProtIDs(silico_SURFY_df[["UniProt accession"]],
                "2018 - The in silico human surfaceome - SURFY - UniProt accessions"
                )

WriteUniProtIDs(silico_HeLa_df[["UniProt ID"]],
                "2018 - The in silico human surfaceome - HeLa - UniProt accessions"
                )




# Read in intermediate data -----------------------------------------------

mass_spec_entrezs_df  <- read.delim(file.path(UniProt_output_directory, "2015 - A Mass Spectrometric-Derived Cell Surface Protein Atlas - UniProt accessions to Entrez IDs.tab"),
                                    stringsAsFactors = FALSE, check.names = FALSE
                                    )
SURFY_name_entrezs_df <- read.delim(file.path(UniProt_output_directory, "2018 - The in silico human surfaceome - SURFY - UniProt names to Entrez IDs.tab"),
                                    stringsAsFactors = FALSE, check.names = FALSE
                                    )
SURFY_acc_entrezs_df  <- read.delim(file.path(UniProt_output_directory, "2018 - The in silico human surfaceome - SURFY - UniProt accessions to Entrez IDs.tab"),
                                    stringsAsFactors = FALSE, check.names = FALSE
                                    )
HeLa_entrezs_df       <- read.delim(file.path(UniProt_output_directory, "2018 - The in silico human surfaceome - HeLa - UniProt accessions to Entrez IDs.tab"),
                                    stringsAsFactors = FALSE, check.names = FALSE
                                    )



# Define functions --------------------------------------------------------

CategorizeHPAReliability <- function(char_vec) {
  ifelse(char_vec == "",
         NA,
         char_vec
         )
  results_fac <- factor(char_vec,
                        levels = c("Uncertain", "Approved", "Supported", "Enhanced"),
                        ordered = TRUE
                        )
  return(results_fac)
}


MapUniProtToEntrezs <- function(mappings_df, UniProt_IDs) {
  vapply(UniProt_IDs, function(x) {
    are_this_ID <- mappings_df[["From"]] %in% x
    if (!(any(are_this_ID))) {
      return(NA_character_)
    } else {
      paste0(sort(unique(mappings_df[["To"]][are_this_ID])), collapse = ", ")
    }
  }, "")
}


SplitEntrezs <- function(entrezs_vec) {
  results_vec <- unlist(strsplit(entrezs_vec, ", ", fixed = TRUE), use.names = FALSE)
  results_vec <- setdiff(unique(as.integer(results_vec)), NA_integer_)
  return(results_vec)
}



GetMSEntrezs <- function(UniProt_IDs) {
  matches_vec <- match(UniProt_IDs, mass_spec_AB_df[["ID_link"]])
  mapped_entrezs <- SplitEntrezs(mass_spec_AB_df[["Consensus_Entrez_ID"]][matches_vec])
  all_MS_entrezs <- ifelse(mass_spec_AB_df[["ENTREZ ac"]] == 0,
                           mass_spec_AB_df[["ENTREZ geneID"]],
                           mass_spec_AB_df[["ENTREZ ac"]]
                           )
  original_entrezs <- setdiff(all_MS_entrezs[matches_vec], 0L)
  results_list <- list(
    "only_original" = setdiff(original_entrezs, mapped_entrezs),
    "only_mapped"   = setdiff(mapped_entrezs, original_entrezs),
    "all_mapped"    = mapped_entrezs
  )
  return(results_list)
}



AddTranslatedSymbols <- function(input_df, symbol_column) {
  mapped_df <- MapToEntrezs(symbols_vec = input_df[[symbol_column]])
  mapped_df <- mapped_df[, c("Entrez_ID", "Gene_symbol", "Original_symbol")]
  names(mapped_df) <- c("Symbol_to_Entrez_IDs",
                        "Backtranslated_symbol",
                        "Original_symbol"
                        )
  results_df <- data.frame(
    input_df,
    mapped_df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  return(results_df)
}



AreConflictingEntrezs <- function(input_df) {
  are_identical <- mapply(identical,
                          input_df[["UniProt_acc_to_Entrez_IDs"]],
                          input_df[["Symbol_to_Entrez_IDs"]]
                          )
  are_conflicting <- !(are_identical) &
                     !(is.na(input_df[["UniProt_acc_to_Entrez_IDs"]])) &
                     !(is.na(input_df[["Symbol_to_Entrez_IDs"]]))
  return(are_conflicting)
}




AddConsensusEntrezs <- function(input_df) {
  are_conflicting <- AreConflictingEntrezs(input_df)
  uni_splits <- strsplit(input_df[["UniProt_acc_to_Entrez_IDs"]], ", ", fixed = TRUE)
  sym_splits <- strsplit(input_df[["Symbol_to_Entrez_IDs"]], ", ", fixed = TRUE)

  consensus_list <- lapply(which(are_conflicting),
                           function(x) {
                             if (any(uni_splits[[x]] %in% sym_splits[[x]])) {
                               intersect(uni_splits[[x]], sym_splits[[x]])
                             } else {
                               uni_splits[[x]]
                             }
                           })

  consensus_vec <- vapply(consensus_list, function(x) paste0(x, collapse = ", "), "")
  results_vec <- input_df[["UniProt_acc_to_Entrez_IDs"]]
  results_vec[are_conflicting] <- consensus_vec
  are_NA <- is.na(results_vec)
  results_vec[are_NA] <- input_df[["Symbol_to_Entrez_IDs"]][are_NA]

  input_df[["Consensus_Entrez_ID"]] <- results_vec
  return(input_df)
}





# Process the CD (Cluster of Differentiation) gene list -------------------

names(CD_df) <- c("CD_NAME", "NCBI_NAME", "GENE_NAME", "NCBI_OTHER_NAME")

CD_mapped_df <- MapToEntrezs(symbols_vec = CD_df[["NCBI_NAME"]])
combined_CD_df <- data.frame(CD_df, CD_mapped_df)

are_to_ignore <- (CD_df[["NCBI_NAME"]] == "carbohydrate") |
                 grepl("^see ", CD_df[["NCBI_NAME"]])

are_ambiguous <- !(CD_mapped_df[["Entrez_source"]] %in% c(1, 3))

combined_CD_df[are_ambiguous & !(are_to_ignore), ]

manual_mappings_vec <- c(
  "CD20"   = "MS4A1",
  "MCP"    = "CD46",
  "CD97"   = "ADGRE5",
  "KIR3DP1; KIR2DS6; KIRX" = "KIR3DP1",
  "SN"     = "SIGLEC1",
  "IL8RA"  = "CXCR1"
)

manual_entrezs_vec <- MapToEntrezs(symbols_vec = manual_mappings_vec)[["Entrez_ID"]]
names(manual_entrezs_vec) <- names(manual_mappings_vec)

manually_mapped_vec <- manual_mappings_vec[CD_df[["NCBI_NAME"]]]
stopifnot(identical(sum(!(is.na(manually_mapped_vec))),
                    length(manual_mappings_vec)
                    ))

CD_df[["Entrez_ID"]] <- ifelse(is.na(manually_mapped_vec),
                               CD_mapped_df[["Entrez_ID"]],
                               manually_mapped_vec
                               )



# Create a mapping from Entrez IDs to CD IDs ------------------------------

num_occurrences <- table(CD_df[["Entrez_ID"]])[CD_df[["Entrez_ID"]]]
CD_df[(num_occurrences > 1) %in% TRUE, ]

are_included <- !(is.na(CD_df[["Entrez_ID"]])) &
                (!(CD_df[["Entrez_ID"]] %in% "5788") |
                  (CD_df[["CD_NAME"]] %in% "CD45"))
CD_to_entrezs_list <- strsplit(CD_df[["Entrez_ID"]][are_included], ", ", fixed = TRUE)
names(CD_to_entrezs_list) <- CD_df[["CD_NAME"]][are_included]

entrez_to_CD_map <- rep(names(CD_to_entrezs_list),
                        lengths(CD_to_entrezs_list)
                        )
names(entrez_to_CD_map) <- unlist(CD_to_entrezs_list, use.names = FALSE)




# Explore the mass spectrometry-based surfaceome --------------------------

table(mass_spec_B_df[["CSPA category"]], useNA = "ifany")

stopifnot(all(mass_spec_B_df[["Organisme"]] == "Human"))





# Explore the predicted human surfaceome ----------------------------------

stopifnot(all(is.na(silico_SURFY_df[["Comment"]])))

table(silico_SURFY_df[["Surfaceome Label"]], useNA = "ifany")
table(silico_SURFY_df[["Surfaceome Label Source"]], useNA = "ifany")
table(silico_SURFY_df[["signalpeptide"]], useNA = "ifany")
table(silico_SURFY_df[["CSPA category"]], useNA = "ifany")
table(silico_SURFY_df[["MachineLearning FPR class (1=1%, 2=5%, 3=15%)"]], useNA = "ifany")

boxplot(silico_SURFY_df[["MachineLearning score"]] ~ silico_SURFY_df[["Surfaceome Label"]])

table(silico_surfaceome_df[["Surfaceome Label Source"]])
table(silico_surfaceome_df[["MachineLearning FPR class (1=1%, 2=5%, 3=15%)"]])

table(silico_groups_df[["Almen category"]])
table(silico_groups_df[["Almen subclass"]])

hist(silico_HeLa_df[["group probability"]])





# Combine data frames for the MS-based surfaceome -------------------------

mass_spec_B_mat <- as.matrix(mass_spec_B_df[, 7:ncol(mass_spec_B_df)])
stopifnot(all(mass_spec_B_mat %in% c(1, NA)))
mode(mass_spec_B_mat) <- "integer"

stopifnot(all(mass_spec_B_df[["ID_link"]] %in% mass_spec_A_df[["ID_link"]]))

mass_spec_AB_matches <- match(mass_spec_B_df[["ID_link"]],
                              mass_spec_A_df[["ID_link"]]
                              )
stopifnot(!(anyNA(mass_spec_AB_matches)))
mass_spec_AB_df <- mass_spec_A_df[mass_spec_AB_matches, ]
row.names(mass_spec_AB_df) <- NULL
mass_spec_AB_df[["CD"]] <- ifelse(mass_spec_AB_df[["CD"]] == "no",
                                  NA,
                                  mass_spec_AB_df[["CD"]]
                                  )

mass_spec_AB_df <- data.frame(
  mass_spec_AB_df[, !(names(mass_spec_AB_df) %in% c("Experiment tag", "protein probability", "num unique peps"))],
  mass_spec_B_df[, !(names(mass_spec_B_df) %in% colnames(mass_spec_B_mat))],
  mass_spec_B_mat,
  check.names = FALSE
)

are_duplicated <- duplicated(as.list(mass_spec_AB_df))
mass_spec_AB_df <- mass_spec_AB_df[, !(are_duplicated)]

stopifnot(identical(1L, unique(table(mass_spec_AB_df[["ID_link"]]))))

mass_spec_AB_df[["UniProt_acc_to_Entrez_IDs"]] <- MapUniProtToEntrezs(mass_spec_entrezs_df, mass_spec_AB_df[["ID_link"]])

mass_spec_AB_df <- AddTranslatedSymbols(mass_spec_AB_df, "ENTREZ gene symbol")
mass_spec_AB_df <- AddConsensusEntrezs(mass_spec_AB_df)




# Explore the combined data from the MS-based surfaceome ------------------

are_identical_entrezs <- mapply(identical,
                                mass_spec_AB_df[["ENTREZ ac"]],
                                mass_spec_AB_df[["ENTREZ geneID"]]
                                )
mass_spec_AB_df[!(are_identical_entrezs), ]


are_identical_CD <- mapply(identical,
                           mass_spec_AB_df[["CD"]],
                           mass_spec_AB_df[["CD.1"]]
                           )
mass_spec_AB_df[!(are_identical_CD), ]




# Define Entrez IDs for the in silico surfaceome --------------------------

silico_HeLa_df[["UniProt_acc_to_Entrez_IDs"]]        <- MapUniProtToEntrezs(HeLa_entrezs_df, silico_HeLa_df[["UniProt ID"]])
silico_surfaceome_df[["UniProt_acc_to_Entrez_IDs"]]  <- MapUniProtToEntrezs(SURFY_acc_entrezs_df, silico_surfaceome_df[["UniProt accession"]])
silico_surfaceome_df[["UniProt_name_to_Entrez_IDs"]] <- MapUniProtToEntrezs(SURFY_name_entrezs_df, silico_surfaceome_df[["UniProt name"]])

silico_surfaceome_df <- AddTranslatedSymbols(silico_surfaceome_df, "UniProt gene")
silico_surfaceome_df <- AddConsensusEntrezs(silico_surfaceome_df)




# Explore the mappings to Entrez IDs --------------------------------------

## For the MS-based data
table(table(mass_spec_entrezs_df[["From"]]))
table(table(mass_spec_entrezs_df[["To"]]))

table(mass_spec_B_df[["ENTREZ geneID"]] == 0)
stopifnot(!(anyNA(mass_spec_B_df[["ENTREZ geneID"]])))
table(table(mass_spec_B_df[["ENTREZ geneID"]]))

no_mapped_entrez <- is.na(mass_spec_AB_df[["UniProt_acc_to_Entrez_IDs"]])
table(!(no_mapped_entrez))
mass_spec_AB_df[no_mapped_entrez, ]
mass_spec_AB_df[AreConflictingEntrezs(mass_spec_AB_df), ]



## For the in silico data
table(table(SURFY_name_entrezs_df[["From"]]))
table(table(SURFY_name_entrezs_df[["To"]]))
table(table(silico_surfaceome_df[["GeneID"]]))

# The following Entrez IDs seem to have been rounded by Excel...
silico_surfaceome_df[silico_surfaceome_df[["GeneID"]] %in% "100000000", ]
silico_surfaceome_df[silico_surfaceome_df[["GeneID"]] %in% "101000000", ]
silico_surfaceome_df[silico_surfaceome_df[["GeneID"]] %in% "103000000", ]


no_mapped_entrez <- is.na(silico_surfaceome_df[["UniProt_acc_to_Entrez_IDs"]])
table(!(no_mapped_entrez))
head(silico_surfaceome_df[no_mapped_entrez, ])
table(!(is.na(silico_surfaceome_df[["UniProt_name_to_Entrez_IDs"]])))

silico_surfaceome_df[AreConflictingEntrezs(silico_surfaceome_df), ]




# Examine Entrez IDs for the MS-based surfaceome --------------------------

are_certain <- mass_spec_AB_df[["CSPA category"]] == "1 - high confidence"

mass_spec_original_entrezs <- setdiff(mass_spec_AB_df[["ENTREZ geneID"]][are_certain], 0)
mass_spec_mapped_entrezs <- SplitEntrezs(mass_spec_AB_df[["Consensus_Entrez_ID"]][are_certain])

only_original_entrezs <- setdiff(mass_spec_original_entrezs, mass_spec_mapped_entrezs)
only_mapped_entrezs   <- setdiff(mass_spec_mapped_entrezs, mass_spec_original_entrezs)

are_only_original <- mass_spec_AB_df[["ENTREZ geneID"]] %in% only_original_entrezs
mass_spec_AB_df[are_only_original, ]


are_only_mapped <- vapply(strsplit(mass_spec_AB_df[["Consensus_Entrez_ID"]], ", ", fixed = TRUE),
                          function(x) any(as.integer(x) %in% only_mapped_entrezs),
                          logical(1)
                          )
mass_spec_AB_df[are_only_mapped, ]




# Define Entrez IDs for the MS-based surfaceome ---------------------------

MS_surfaceome_entrezs <- union(mass_spec_original_entrezs, mass_spec_mapped_entrezs)

are_HeLa  <- mass_spec_AB_df[["HeLa"]] %in% 1
are_HEK   <- mass_spec_AB_df[["HEK"]] %in% 1
are_LN229 <- mass_spec_AB_df[["LN229"]] %in% 1

are_putative <- mass_spec_AB_df[["CSPA category"]] %in% c("1 - high confidence", "2 - putative")

fisher.test(are_certain, are_HeLa)
fisher.test(are_certain, are_HEK)
fisher.test(are_certain, are_LN229)

MS_HeLa_list  <- GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_putative & are_HeLa])
MS_HEK_list   <- GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_putative & are_HEK])
MS_LN229_list <- GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_putative & are_LN229])




# Define sets of genes for the in silico surfaceome -----------------------

silico_HeLa_entrezs       <- SplitEntrezs(silico_HeLa_df[["UniProt_acc_to_Entrez_IDs"]])
silico_surfaceome_entrezs <- SplitEntrezs(silico_surfaceome_df[["Consensus_Entrez_ID"]])




# Define the final list of surfaceome gene sets ---------------------------

surfaceome_list <- list(
  "SURFY_surfaceome_2018" = silico_surfaceome_entrezs,
  "CSC_surfaceome_2015"   = as.integer(MS_surfaceome_entrezs),
  "CSC_HeLa_2018"         = silico_HeLa_entrezs,
  "CSC_HeLa_2015"         = MS_HeLa_list[["all_mapped"]],
  "CSC_HEK_2015"          = MS_HEK_list[["all_mapped"]],
  "CSC_LN229_2015"        = MS_LN229_list[["all_mapped"]]
)

all_entrezs <- sort(unique(unlist(surfaceome_list)))

are_present_list <- lapply(surfaceome_list,
                           function(x) ifelse(all_entrezs %in% x,
                                              "Yes",
                                              ""
                                              )
                           )
surfaceome_mat <- do.call(cbind, are_present_list)

surfaceome_df <- data.frame(
  "Entrez_ID"   = all_entrezs,
  "Gene_symbol" = MapToEntrezs(all_entrezs)[["Gene_symbol"]],
  surfaceome_mat,
  stringsAsFactors = FALSE
)




# Tidy data from the Human Protein Atlas ----------------------------------

HPA_df <- data.frame(
  "Ensembl_gene_ID"      = full_HPA_df[["Ensembl"]],
  "Entrez_ID"            = NA,
  "Gene_symbol"          = full_HPA_df[["Gene"]],
  "Antibody_IDs"         = full_HPA_df[["Antibody"]],
  "Antibody_RRIDs"       = full_HPA_df[["Antibody RRID"]],
  "Reliability_IH"       = CategorizeHPAReliability(full_HPA_df[["Reliability (IH)"]]),
  "Reliability_IF"       = CategorizeHPAReliability(full_HPA_df[["Reliability (IF)"]]),
  "Subcellular_location" = gsub(",", ", ", full_HPA_df[["Subcellular location"]], fixed = TRUE),
  stringsAsFactors       = FALSE
)

HPA_df[["Antibody_RRIDs"]] <- gsub(": (?=$|,)", "", HPA_df[["Antibody_RRIDs"]], perl = TRUE)
HPA_df[["Antibody_RRIDs"]] <- gsub(": ", "/", HPA_df[["Antibody_RRIDs"]], fixed = TRUE)

HPA_mappings_df <- MapEnsemblIDs(HPA_df)
names(HPA_df)[names(HPA_df) == "Gene_symbol"] <- "Original_symbol"

HPA_df[["Entrez_ID"]] <- HPA_mappings_df[["Consensus_entrez"]]




# Save data ---------------------------------------------------------------

save(list = c("HPA_df", "surfaceome_df",
              "CD_df", "entrez_to_CD_map"
              ),
     file = file.path(general_RData_directory, "23) Compile data on protein localization.RData")
     )


