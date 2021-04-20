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
SURFY_name_entrezs_df  <- read.delim(file.path(UniProt_output_directory, "2018 - The in silico human surfaceome - SURFY - UniProt names to Entrez IDs.tab"),
                                    stringsAsFactors = FALSE, check.names = FALSE
                                    )
SURFY_acc_entrezs_df <- read.delim(file.path(UniProt_output_directory, "2018 - The in silico human surfaceome - SURFY - UniProt accessions to Entrez IDs.tab"),
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

#
# ProcessSurfaceEntrezs <- function(mappings_df, UniProt_IDs) {
#   entrez_IDs <- unique(mappings_df[["To"]])
#   results_list <- lapply(entrez_IDs, function(x) {
#     are_this_ID <- mappings_df[["To"]] == x
#     if (!(any(are_this_ID))) {
#       results_vec <- c("Num_IDs" = 0L, "Num_present" = 0L, "Num_absent" = 0L)
#     } else {
#       these_IDs <- mappings_df[["From"]][are_this_ID]
#       are_present <- these_IDs %in% UniProt_IDs
#       results_list <- list(
#         "Entrez_ID"   = x,
#         "Num_IDs"     = sum(are_this_ID),
#         "Num_present" = sum(are_present),
#         "Num_absent"  = sum(!(are_present)),
#         "Present_IDs" = paste0(these_IDs[are_present], collapse = ", "),
#         "Absent_IDs"  = paste0(these_IDs[!(are_present)], collapse = ", ")
#       )
#     }
#   })
#   results_df <- do.call(rbind.data.frame,
#                         c(results_list, list(stringsAsFactors = FALSE,
#                                              make.row.names = FALSE
#                                              )
#                           )
#                         )
#   return(results_df)
# }
#
#
#
# GetSurfaceEntrezs <- function(mappings_df, UniProt_IDs) {
#   entrezs_df <- ProcessSurfaceEntrezs(mappings_df, UniProt_IDs)
#   entrezs_vec <- entrezs_df[["Entrez_ID"]][entrezs_df[["Num_present"]] >= 1]
#   return(entrezs_vec)
# }



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

mass_spec_AB_df[["Mapped_Entrez_IDs"]] <- MapUniProtToEntrezs(mass_spec_entrezs_df, mass_spec_AB_df[["ID_link"]])





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

silico_HeLa_df[["Mapped_Entrez_IDs"]] <- MapUniProtToEntrezs(HeLa_entrezs_df, silico_HeLa_df[["UniProt ID"]])
silico_surfaceome_df[["Mapped_Entrez_IDs"]] <- MapUniProtToEntrezs(SURFY_acc_entrezs_df, silico_surfaceome_df[["UniProt accession"]])
silico_surfaceome_df[["Mapped_UniProt_name_to_Entrez_IDs"]] <- MapUniProtToEntrezs(SURFY_name_entrezs_df, silico_surfaceome_df[["UniProt name"]])




# Explore the mappings to Entrez IDs --------------------------------------

## For the MS-based data
table(table(mass_spec_entrezs_df[["From"]]))
table(table(mass_spec_entrezs_df[["To"]]))

table(mass_spec_B_df[["ENTREZ geneID"]] == 0)
stopifnot(!(anyNA(mass_spec_B_df[["ENTREZ geneID"]])))
table(table(mass_spec_B_df[["ENTREZ geneID"]]))

table(!(is.na(mass_spec_AB_df[["Mapped_Entrez_IDs"]])))
were_found <- mass_spec_B_df[["ID_link"]] %in% mass_spec_entrezs_df[["From"]]
table(were_found)
identical(were_found, !(is.na(mass_spec_AB_df[["Mapped_Entrez_IDs"]])))



## For the in silico data
table(table(SURFY_name_entrezs_df[["From"]]))
table(table(SURFY_name_entrezs_df[["To"]]))
table(table(silico_surfaceome_df[["GeneID"]]))

# The following Entrez IDs seem to have been rounded by Excel...
silico_surfaceome_df[silico_surfaceome_df[["GeneID"]] %in% "100000000", ]
silico_surfaceome_df[silico_surfaceome_df[["GeneID"]] %in% "101000000", ]

name_found <- silico_surfaceome_df[["UniProt name"]] %in% SURFY_name_entrezs_df[["From"]]
accession_found <- silico_surfaceome_df[["UniProt accession"]] %in% SURFY_acc_entrezs_df[["From"]]
table(name_found)
table(accession_found)
head(silico_surfaceome_df[!(accession_found), ])

identical(name_found, !(is.na(silico_surfaceome_df[["Mapped_UniProt_name_to_Entrez_IDs"]])))
identical(accession_found, !(is.na(silico_surfaceome_df[["Mapped_Entrez_IDs"]])))

table(name_found, !(is.na(silico_surfaceome_df[["Mapped_Entrez_IDs"]])))
table(accession_found, !(is.na(silico_surfaceome_df[["Mapped_UniProt_name_to_Entrez_IDs"]])))



HeLa_found <- silico_HeLa_df[["UniProt ID"]] %in% HeLa_entrezs_df[["From"]]
table(HeLa_found)
silico_HeLa_df[!(HeLa_found), ]

identical(HeLa_found, !(is.na(silico_HeLa_df[["Mapped_Entrez_IDs"]])))



# Examine Entrez IDs for the MS-based surfaceome --------------------------

are_certain <- mass_spec_AB_df[["CSPA category"]] == "1 - high confidence"
mass_spec_original_entrezs <- setdiff(mass_spec_AB_df[["ENTREZ geneID"]][are_certain], 0)

mass_spec_mapped_df <- ProcessSurfaceEntrezs(mass_spec_entrezs_df,
                                             mass_spec_AB_df[["ID_link"]][are_certain]
                                             )
mass_spec_mapped_entrezs <- mass_spec_mapped_df[["Entrez_ID"]][mass_spec_mapped_df[["Num_present"]] >= 1]

only_original_entrezs <- setdiff(mass_spec_original_entrezs, mass_spec_mapped_entrezs)
are_only_original <- mass_spec_AB_df[["ENTREZ geneID"]] %in% only_original_entrezs
mass_spec_AB_df[mass_spec_AB_df[["ENTREZ geneID"]] %in% only_original_entrezs, ]

only_mapped_entrezs <- setdiff(mass_spec_mapped_entrezs, mass_spec_original_entrezs)
only_mapped_IDs <- mass_spec_mapped_df[["Present_IDs"]][mass_spec_mapped_df[["Entrez_ID"]] %in% only_mapped_entrezs]
are_only_mapped <- mass_spec_AB_df[["ID_link"]] %in% only_mapped_IDs
mass_spec_AB_df[are_only_mapped, ]





# Define Entrez IDs for the MS-based surfaceome ---------------------------

MS_surfaceome_entrezs <- union(mass_spec_original_entrezs, mass_spec_mapped_entrezs)

mass_spec_mapped_entrezs_2 <- setdiff(as.integer(unlist(strsplit(mass_spec_AB_df[["Mapped_Entrez_IDs"]][are_certain], ", ", fixed = TRUE))), NA_integer_)
identical(sort(mass_spec_mapped_entrezs), sort(mass_spec_mapped_entrezs_2))



GetMSEntrezs <- function(UniProt_IDs) {
  mapped_entrezs <- GetSurfaceEntrezs(mass_spec_entrezs_df, UniProt_IDs)
  all_MS_entrezs <- ifelse(mass_spec_AB_df[["ENTREZ ac"]] == 0,
                           mass_spec_AB_df[["ENTREZ geneID"]],
                           mass_spec_AB_df[["ENTREZ ac"]]
                           )
  matches_vec <- match(UniProt_IDs, mass_spec_AB_df[["ID_link"]])
  original_entrezs <- setdiff(all_MS_entrezs[matches_vec], 0L)
  results_list <- list(
    "union"         = union(mapped_entrezs, original_entrezs),
    "only_original" = setdiff(original_entrezs, mapped_entrezs),
    "only_mapped"   = setdiff(mapped_entrezs, original_entrezs)
  )
  return(results_list)
}

stopifnot(identical(sort(GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_certain])[["union"]]),
                    sort(MS_surfaceome_entrezs)
                    ))

are_HeLa  <- mass_spec_AB_df[["HeLa"]] %in% 1
are_HEK   <- mass_spec_AB_df[["HEK"]] %in% 1
are_LN229 <- mass_spec_AB_df[["LN229"]] %in% 1

fisher.test(are_certain, are_HeLa)
fisher.test(are_certain, are_HEK)
fisher.test(are_certain, are_LN229)

are_putative <- mass_spec_AB_df[["CSPA category"]] %in% c("1 - high confidence", "2 - putative")

MS_HeLa_list  <- GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_putative & are_HeLa])
MS_HEK_list   <- GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_putative & are_HEK])
MS_LN229_list <- GetMSEntrezs(mass_spec_AB_df[["ID_link"]][are_putative & are_LN229])


HeLa_entrezs_2 <- setdiff(as.integer(unlist(strsplit(mass_spec_AB_df[["Mapped_Entrez_IDs"]][are_putative & are_HeLa], ", ", fixed = TRUE))), NA_integer_)
identical(sort(HeLa_entrezs_2), as.integer(sort(setdiff(MS_HeLa_list[["union"]], MS_HeLa_list[["only_original"]]))))

HEK_entrezs_2 <- setdiff(as.integer(unlist(strsplit(mass_spec_AB_df[["Mapped_Entrez_IDs"]][are_putative & are_HEK], ", ", fixed = TRUE))), NA_integer_)
identical(sort(HEK_entrezs_2), as.integer(sort(setdiff(MS_HEK_list[["union"]], MS_HEK_list[["only_original"]]))))

LN229_entrezs_2 <- setdiff(as.integer(unlist(strsplit(mass_spec_AB_df[["Mapped_Entrez_IDs"]][are_putative & are_LN229], ", ", fixed = TRUE))), NA_integer_)
identical(sort(LN229_entrezs_2), as.integer(sort(setdiff(MS_LN229_list[["union"]], MS_LN229_list[["only_original"]]))))





# Define Entrez IDs for the in silico surfaceome --------------------------

silico_HeLa_df[["Mapped_Entrez_IDs"]] <- MapUniProtToEntrezs(HeLa_entrezs_df, silico_HeLa_df[["UniProt ID"]])
silico_surfaceome_df[["Mapped_Entrez_IDs"]] <- MapUniProtToEntrezs(SURFY_acc_entrezs_df, silico_surfaceome_df[["UniProt accession"]])

silico_HeLa_entrezs <- GetSurfaceEntrezs(HeLa_entrezs_df, silico_HeLa_df[["UniProt ID"]])
silico_surfaceome_UniProt_names_to_entrezs <- GetSurfaceEntrezs(SURFY_name_entrezs_df, silico_surfaceome_df[["UniProt name"]])
silico_surfaceome_entrezs <- GetSurfaceEntrezs(SURFY_acc_entrezs_df, silico_surfaceome_df[["UniProt accession"]])

silico_surfaceome_entrezs_2 <- setdiff(as.integer(unlist(strsplit(silico_surfaceome_df[["Mapped_Entrez_IDs"]], ", ", fixed = TRUE))), NA_integer_)
identical(sort(silico_surfaceome_entrezs_2), as.integer(sort(silico_surfaceome_entrezs)))





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

HPA_mappings_df <- MapEnsemblIDs(HPA_df)
names(HPA_df)[names(HPA_df) == "Gene_symbol"] <- "Original_symbol"

HPA_df[["Entrez_ID"]] <- HPA_mappings_df[["Consensus_entrez"]]






# Save data ---------------------------------------------------------------

save(list = "HPA_df",
     file = file.path(general_RData_directory, "23) Compile data on protein localization.RData")
     )


