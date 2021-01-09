### 21st July 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "28) Merging the FANTOM5 and BioMart TSS data.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
human_genome_directory  <- file.path(CRISPR_root_directory, "2) Input data", "Human genome")
FANTOM5_input_directory <- file.path(human_genome_directory, "FANTOM5_liftover")
Ensembl_input_directory <- file.path(human_genome_directory, "Ensembl")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Read in data ------------------------------------------------------------

# The HGNC file was downloaded from https://www.genenames.org/download/custom/
HGNC_df <- read.table(file.path(human_genome_directory, "HGNC", "Custom_HGNC_download_including_NCBI_IDs__2020_03_11.tsv"),
                      sep = "\t", quote = "", fill = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE, header = TRUE, row.names = NULL
                      )

# The two FANTOM5 files were downloaded from: https://figshare.com/articles/Re-processing_of_the_data_generated_by_the_FANTOM5_project_hg38_v3_CAGE_peaks/4880063/4
# on 21 July 2019
FANTOM5_ann_df <- read.table(file.path(FANTOM5_input_directory, "hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE,
                             header = TRUE, row.names = NULL, fill = TRUE,
                             check.names = FALSE, comment.char = ""
                             )

FANTOM5_bed_df <- read.table(file.path(FANTOM5_input_directory, "hg38_liftover_CAGE_peaks_phase1and2.bed"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                             )

# The BioMart file was downloaded from https://www.ensembl.org/biomart/martview
BioMart_df <- read.table(file.path(human_genome_directory, "Ensembl", "BioMart_human_2020-03-25_mart_export.txt"),
                         sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                         )

# Read in a later version of the BioMart export file that also contains the gene type
BioMart_gene_type_df <- read.table(file.path(human_genome_directory, "Ensembl", "BioMart_human_2020-06-03_mart_export_with_gene_type.txt"),
                                   sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                                   )




# Define functions --------------------------------------------------------

StandardizeFANTOM5IDs <- function(char_vec) {
  results_vec <- gsub(" ", ", ", char_vec, fixed = TRUE)
  results_vec <- ifelse(char_vec == "", NA_character_, results_vec)
  return(results_vec)
}




# Make small adjustments to the data frames -------------------------------

names(FANTOM5_ann_df)[[1]]  <- sub("#", "", names(FANTOM5_ann_df)[[1]], fixed = TRUE)

names(FANTOM5_bed_df) <- c(
  "Chromosome",
  "Peak_start",
  "Peak_end",
  "Peak_ID",
  "Score",
  "Strand",
  "TSS_start",
  "TSS_stop",
  "Color_code"
)

### From the readme file:
# Description of the columns in CAGE peaks coordinates files
# - chromosome
# - start of CAGE peak region
# - end of CAGE peak region
# - name (ID) of the CAGE peak
# - score
# - strand of the CAGE peak
# - start of the representative TSS position
# - end of the representative TSS position (Note: end is always start+1)
# - rgb string for color coding (plus or minus strand only)






# Integrate "gene type" into the BioMart data frame -----------------------

IDs_older <- paste0(BioMart_df[["Gene stable ID"]], "__",
                    BioMart_df[["Transcript stable ID"]]
                    )
IDs_newer <- paste0(BioMart_gene_type_df[["Gene stable ID"]], "__",
                    BioMart_gene_type_df[["Transcript stable ID"]]
                    )
ID_matches <- match(IDs_older, IDs_newer)

for (column_name in c("Gene type", "Transcript type")) {
  BioMart_df[[column_name]] <- BioMart_gene_type_df[[column_name]][ID_matches]
}




# Check for unannotated FANTOM5 CAGE peaks --------------------------------

annotation_columns <- c("Transcript_name", "GeneID", "HGNC/MGI_ID",
                        "UniProt_ID", "Gene_name", "Gene_symbol",
                        "Gene_synonyms", "Gene_source"
                        )
annotation_mat <- as.matrix(FANTOM5_ann_df[, annotation_columns])

are_empty_mat <- apply(annotation_mat, 2, function(x) x == "")
are_all_empty <- apply(are_empty_mat, 1, all)
have_NA_distance <- is.na(FANTOM5_ann_df[["Distance"]])

stopifnot(identical(are_all_empty, have_NA_distance))





# Check FANTOM5 CAGE peak IDs ---------------------------------------------

table(FANTOM5_ann_df[["CAGE_Peak_ID"]] %in% FANTOM5_bed_df[["Peak_ID"]])
table(FANTOM5_bed_df[["Peak_ID"]] %in% FANTOM5_ann_df[["CAGE_Peak_ID"]])

any(duplicated(FANTOM5_ann_df[["CAGE_Peak_ID"]]))
any(duplicated(FANTOM5_bed_df[["Peak_ID"]]))




# Remove unannotated FANTOM5 CAGE peaks -----------------------------------

table(FANTOM5_bed_df[["Peak_ID"]] %in% FANTOM5_ann_df[["CAGE_Peak_ID"]])

FANTOM5_ann_df <- FANTOM5_ann_df[!(are_all_empty), ]
are_present <- FANTOM5_bed_df[["Peak_ID"]] %in% FANTOM5_ann_df[["CAGE_Peak_ID"]]
FANTOM5_bed_df <- FANTOM5_bed_df[are_present, ]

rownames(FANTOM5_ann_df) <- NULL
rownames(FANTOM5_bed_df) <- NULL





# Modify the HGNC data frame ----------------------------------------------

HGNC_df[["Entrez_ID"]] <- ifelse(is.na(HGNC_df[["NCBI Gene ID"]]),
                                 HGNC_df[["NCBI Gene ID(supplied by NCBI)"]],
                                 HGNC_df[["NCBI Gene ID"]]
                                 )




# Modify the FANTOM5 annotation data frame --------------------------------

HGNC_splits <- strsplit(FANTOM5_ann_df[["HGNC/MGI_ID"]], " ", fixed = TRUE)
HGNC_vec <- rep(NA_character_, length(HGNC_splits))
num_IDs <- lengths(HGNC_splits)
HGNC_matches <- match(FANTOM5_ann_df[["HGNC/MGI_ID"]], HGNC_df[["HGNC ID"]])
HGNC_vec[num_IDs == 1] <- as.character(HGNC_df[["Entrez_ID"]][HGNC_matches[num_IDs == 1]])
HGNC_vec[num_IDs > 1] <- vapply(HGNC_splits[num_IDs > 1], function(x) {
  my_matches <- match(x, HGNC_df[["HGNC ID"]])
  matched_IDs <- as.character(HGNC_df[["Entrez_ID"]][my_matches])
  are_NA <- is.na(matched_IDs)
  if (all(are_NA)) {
    return(NA_character_)
  } else {
    matched_IDs <- matched_IDs[!(are_NA)]
    return(paste0(matched_IDs, collapse = ", "))
  }
}, "")

FANTOM5_ann_df[["HGNC_entrez"]] <- HGNC_vec

FANTOM5_ann_df[["FANTOM5_entrez"]] <- StandardizeFANTOM5IDs(FANTOM5_ann_df[["GeneID"]])
FANTOM5_ann_df <- FANTOM5_ann_df[, names(FANTOM5_ann_df) != "GeneID"]

have_comma_FANTOM5 <- grepl(",", FANTOM5_ann_df[["FANTOM5_entrez"]], fixed = TRUE)
have_comma_HGNC <- grepl(",", FANTOM5_ann_df[["HGNC_entrez"]], fixed = TRUE)
table(have_comma_FANTOM5, have_comma_HGNC)

FANTOM5_ann_df[["Consensus_entrez"]] <- ifelse(is.na(FANTOM5_ann_df[["FANTOM5_entrez"]]),
                                               FANTOM5_ann_df[["HGNC_entrez"]],
                                               ifelse(have_comma_FANTOM5 & !(have_comma_HGNC),
                                                      ifelse(is.na(FANTOM5_ann_df[["HGNC_entrez"]]),
                                                             FANTOM5_ann_df[["FANTOM5_entrez"]],
                                                             FANTOM5_ann_df[["HGNC_entrez"]]
                                                             ),
                                                      FANTOM5_ann_df[["FANTOM5_entrez"]]
                                                      )
                                               )




# Tidy some additional FANTOM5 annotation columns -------------------------

# FANTOM5_ann_df[["HGNC_ID"]] <- StandardizeFANTOM5IDs(FANTOM5_ann_df[["HGNC/MGI_ID"]])
FANTOM5_ann_df[["HGNC_ID"]] <- gsub("HGNC:", "", FANTOM5_ann_df[["HGNC/MGI_ID"]], fixed = TRUE)

FANTOM5_ann_df <- FANTOM5_ann_df[, names(FANTOM5_ann_df) != "HGNC/MGI_ID"]

tidy_annotation_columns <- c(
  "Transcript_name", "HGNC_ID",
  "UniProt_ID", "Gene_symbol", "Gene_synonyms", "Gene_source"
)

for (column_name in tidy_annotation_columns) {
  FANTOM5_ann_df[[column_name]] <- StandardizeFANTOM5IDs(FANTOM5_ann_df[[column_name]])
}




# Build a combined FANTOM5 data frame -------------------------------------

FANTOM_ann_matches <- match(FANTOM5_bed_df[["Peak_ID"]], FANTOM5_ann_df[["CAGE_Peak_ID"]])

FANTOM5_df <- data.frame(
  FANTOM5_ann_df[FANTOM_ann_matches, c("CAGE_Peak_ID", "Consensus_entrez", "Gene_symbol")],
  FANTOM5_bed_df[, c("Chromosome", "Strand", "Peak_start", "Peak_end", "TSS_start", "TSS_stop")],
  FANTOM5_ann_df[FANTOM_ann_matches, "Distance", drop = FALSE],
  FANTOM5_bed_df["Score"],
  FANTOM5_ann_df[FANTOM_ann_matches, c("FANTOM5_entrez", tidy_annotation_columns)],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(FANTOM5_df)[names(FANTOM5_df) == "Consensus_entrez"] <- "Entrez_ID"

stopifnot(all(grepl("chr", FANTOM5_df[["Chromosome"]], fixed = TRUE)))






# Perform checks on the combined FANTOM5 data frame -----------------------

unique_IDs_FANTOM5_df <- unique(FANTOM5_df[, c("Entrez_ID", "Gene_symbol", "Chromosome")])

num_occurrences_gene_symbol <- table(unique_IDs_FANTOM5_df[["Gene_symbol"]])[unique_IDs_FANTOM5_df[["Gene_symbol"]]]
num_occurrences_entrez_ID   <- table(unique_IDs_FANTOM5_df[["Entrez_ID"]])[unique_IDs_FANTOM5_df[["Entrez_ID"]]]

duplicated_symbols_df <- unique_IDs_FANTOM5_df[(num_occurrences_gene_symbol > 1) %in% TRUE, ]
duplicated_entrezs_df <- unique_IDs_FANTOM5_df[(num_occurrences_entrez_ID > 1) %in% TRUE, ]

duplicated_symbols_df <- duplicated_symbols_df[order(duplicated_symbols_df[["Gene_symbol"]]), ]
duplicated_entrezs_df <- duplicated_entrezs_df[order(GetMinEntrez(duplicated_entrezs_df[["Entrez_ID"]])), ]

pasted_IDs <- ifelse(is.na(unique_IDs_FANTOM5_df[["Entrez_ID"]]) & is.na(unique_IDs_FANTOM5_df[["Gene_symbol"]]),
                     NA_character_,
                     paste0(unique_IDs_FANTOM5_df[["Entrez_ID"]], "__", unique_IDs_FANTOM5_df[["Gene_symbol"]])
                     )
num_occurrences_pasted_IDs <- table(pasted_IDs)[pasted_IDs]

my_order <- order(GetMinEntrez(unique_IDs_FANTOM5_df[["Entrez_ID"]]), pasted_IDs)
unique_IDs_FANTOM5_df[my_order, ][(num_occurrences_pasted_IDs[my_order] > 1) %in% TRUE, ]





# Filter and standardize the combined FANTOM5 data frame ------------------

FANTOM5_group_vec <- paste0(ifelse(is.na(FANTOM5_df[["Entrez_ID"]]),
                                   "",
                                   paste0(FANTOM5_df[["Entrez_ID"]], " | ")
                                   ),
                            FANTOM5_df[["Gene_symbol"]], " | ",
                            FANTOM5_df[["Chromosome"]]
                            )

have_a_gene_ID <- (!(is.na(FANTOM5_df[["Entrez_ID"]]))) |
                  (!(is.na(FANTOM5_df[["Gene_symbol"]])))

FANTOM5_df <- data.frame(
  "Group" = ifelse(have_a_gene_ID, FANTOM5_group_vec, NA_character_),
  FANTOM5_df,
  stringsAsFactors = FALSE,
  row.names = NULL
)

FANTOM5_filtered_df <- FANTOM5_df[have_a_gene_ID, ]
rownames(FANTOM5_filtered_df) <- NULL






# Construct a simplified BioMart data frame -------------------------------

chromosomes_vec <- TidyBioMartChromosomes(BioMart_df[["Chromosome/scaffold name"]])

BioMart_tidied_df <- BioMart_df

BioMart_tidied_df <- data.frame(
  "Group"             = paste0(ifelse(is.na(BioMart_tidied_df[["NCBI gene ID"]]),
                                      "",
                                      paste0(BioMart_tidied_df[["NCBI gene ID"]], " | ")
                                      ),
                               BioMart_tidied_df[["Gene name"]], " | ",
                               chromosomes_vec
                               ),
  "Entrez_ID"         = BioMart_tidied_df[["NCBI gene ID"]],
  "Gene_symbol"       = BioMart_tidied_df[["Gene name"]],
  "ENSG"              = BioMart_tidied_df[["Gene stable ID"]],
  "ENST"              = BioMart_tidied_df[["Transcript stable ID"]],
  "Chromosome"        = chromosomes_vec,
  "Strand"            = ifelse(BioMart_tidied_df[["Strand"]] == -1, "-", "+"),
  "Gene_start"        = BioMart_tidied_df[["Gene start (bp)"]],
  "Gene_end"          = BioMart_tidied_df[["Gene end (bp)"]],
  "Transcript_start"  = BioMart_tidied_df[["Transcript start (bp)"]],
  "Transcript_end"    = BioMart_tidied_df[["Transcript end (bp)"]],
  "TSS"               = BioMart_tidied_df[["Transcription start site (TSS)"]],
  "Transcript_length" = BioMart_tidied_df[["Transcript length (including UTRs and CDS)"]],
  "Gene_type"         = BioMart_tidied_df[["Gene type"]],
  "Transcript_type"   = BioMart_tidied_df[["Transcript type"]],
  stringsAsFactors    = FALSE
)

are_on_chromosome <- BioMart_df[["Chromosome/scaffold name"]] %in% c(as.character(1:22), "X", "Y", "MT")

BioMart_filtered_df <- BioMart_tidied_df[are_on_chromosome, ]
rownames(BioMart_filtered_df) <- NULL





# Merge the FANTOM5 and BioMart data frames -------------------------------

length(intersect(unique(FANTOM5_filtered_df[["Group"]]), unique(BioMart_filtered_df[["Group"]])))
length(unique(FANTOM5_filtered_df[["Group"]]))
length(unique(BioMart_filtered_df[["Group"]]))

combined_TSS_df <- MergeTSSData(FANTOM5_filtered_df, BioMart_filtered_df)






# Construct a data frame of data combined from FANTOM5 and BioMart --------

combined_TSS_mat <- as.matrix(combined_TSS_df[, c("TSS", "Score")])

combined_TSS_summary_list <- sapply(unique(combined_TSS_df[["Group"]]), function(x) {
  print(x)
  are_this_group <- combined_TSS_df[["Group"]] == x
  sub_mat <- combined_TSS_mat[are_this_group, , drop = FALSE]
  best_index <- which.max(sub_mat[, "Score"])
  first_index <- which(are_this_group)[[1]]
  results_list <- list(
    "Group"       = x,
    "Entrez_ID"   = combined_TSS_df[["Entrez_ID"]][[first_index]],
    "Gene_symbol" = combined_TSS_df[["Gene_symbol"]][[first_index]],
    "Chromosome"  = combined_TSS_df[["Chromosome"]][[first_index]],
    "Strand"      = if (length(best_index) == 1) combined_TSS_df[["Strand"]][are_this_group][[best_index]] else combined_TSS_df[["Strand"]][[first_index]],
    "Best_TSS"    = if (length(best_index) == 1) unname(sub_mat[best_index, "TSS"]) else NA_integer_,
    "First_TSS"   = min(sub_mat[, "TSS"]),
    "Last_TSS"    = max(sub_mat[, "TSS"])
  )
  return(results_list)

}, simplify = FALSE)

combined_TSS_summary_df <- do.call(rbind.data.frame, c(combined_TSS_summary_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))




# # Construct a summarized FANTOM5 data frame -------------------------------
#
# FANTOM5_filtered_mat <- as.matrix(FANTOM5_filtered_df[, c("Peak_start", "Peak_end", "TSS_start", "Score")])
#
# FANTOM5_summary_list <- sapply(unique(FANTOM5_filtered_df[["Group"]]), function(x) {
#   print(x)
#   are_this_group <- FANTOM5_filtered_df[["Group"]] == x
#   sub_mat <- FANTOM5_filtered_mat[are_this_group, , drop = FALSE]
#   best_index <- which.max(sub_mat[, "Score"])
#   first_index <- which(are_this_group)[[1]]
#   results_list <- list(
#     "Group"       = x,
#     "Entrez_ID"   = FANTOM5_filtered_df[["Entrez_ID"]][[first_index]],
#     "Gene_symbol" = FANTOM5_filtered_df[["Gene_symbol"]][[first_index]],
#     "Chromosome"  = FANTOM5_filtered_df[["Chromosome"]][[first_index]],
#     "Strand"      = FANTOM5_filtered_df[["Strand"]][are_this_group][[best_index]],
#     "Best_TSS"    = unname(sub_mat[best_index, "TSS_start"]),
#     "First_TSS"   = min(sub_mat[, "TSS_start"]),
#     "Last_TSS"    = max(sub_mat[, "TSS_start"]),
#     "First_peak"  = min(sub_mat[, "Peak_start"]),
#     "Last_peak"   = min(sub_mat[, "Peak_end"])
#   )
#   return(results_list)
#
# }, simplify = FALSE)
#
# FANTOM5_summary_df <- do.call(rbind.data.frame, c(FANTOM5_summary_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
#
#
#
#
# # Construct a summarized BioMart data frame -------------------------------
#
# BioMart_filtered_mat <- as.matrix(BioMart_filtered_df[, c("Gene_start", "Gene_end", "Transcript_start", "Transcript_end", "TSS", "Transcript_length")])
#
# BioMart_summary_list <- sapply(unique(BioMart_filtered_df[["Group"]]), function(x) {
#   print(x)
#   are_this_group <- BioMart_filtered_df[["Group"]] == x
#   sub_mat <- BioMart_filtered_mat[are_this_group, , drop = FALSE]
#   first_index <- which(are_this_group)[[1]]
#   results_list <- list(
#     "Group"                  = x,
#     "Entrez_ID"              = as.character(BioMart_filtered_df[["Entrez_ID"]][[first_index]]),
#     "Gene_symbol"            = BioMart_filtered_df[["Gene_symbol"]][[first_index]],
#     "Chromosome"             = BioMart_filtered_df[["Chromosome"]][[first_index]],
#     "Strand"                 = BioMart_filtered_df[["Strand"]][[first_index]],
#     "First_TSS"              = min(sub_mat[, "TSS"]),
#     "Last_TSS"               = max(sub_mat[, "TSS"]),
#     "First_transcript_start" = min(sub_mat[, "Transcript_start"]),
#     "Last_transcript_start"  = max(sub_mat[, "Transcript_end"]),
#     "First_gene_start"       = min(sub_mat[, "Gene_start"]),
#     "Last_gene_end"          = max(sub_mat[, "Gene_end"])
#   )
#   return(results_list)
#
# }, simplify = FALSE)
#
# BioMart_summary_df <- do.call(rbind.data.frame, c(BioMart_summary_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))




# Save data ---------------------------------------------------------------

save(list = c("combined_TSS_df", "BioMart_filtered_df", "FANTOM5_filtered_df",
              "BioMart_tidied_df", "FANTOM5_df",
              "combined_TSS_summary_df" #, "BioMart_summary_df", "FANTOM5_summary_df"
              ),
     file = file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData")
     )





