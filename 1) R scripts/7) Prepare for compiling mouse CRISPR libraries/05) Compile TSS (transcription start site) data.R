### 21st July 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "28) Merging the FANTOM5 and BioMart TSS data.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
mouse_genome_directory  <- file.path(CRISPR_root_directory, "2) Input data", "Mouse genome")
FANTOM5_input_directory <- file.path(mouse_genome_directory, "FANTOM5_liftover")
Ensembl_input_directory <- file.path(mouse_genome_directory, "Ensembl")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "6) Mouse - General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Read in data ------------------------------------------------------------

# The two FANTOM5 files were downloaded from: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/
# on 6 October 2020
FANTOM5_ann_df <- read.table(file.path(FANTOM5_input_directory, "mm10_liftover+new_CAGE_peaks_phase1and2_annot.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, fill = TRUE, check.names = FALSE, comment.char = ""
                             )

FANTOM5_bed_df <- read.table(file.path(FANTOM5_input_directory, "mm10_liftover+new_CAGE_peaks_phase1and2.bed"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                             )

# The BioMart file was downloaded from https://www.ensembl.org/biomart/martview
BioMart_df <- read.table(file.path(mouse_genome_directory, "Ensembl", "BioMart_mouse_2020-10-07_mart_export.txt"),
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






# Modify the FANTOM5 annotation data frame --------------------------------

MGI_splits <- strsplit(FANTOM5_ann_df[["HGNC/MGI_ID"]], " ", fixed = TRUE)
MGI_vec <- rep(NA_character_, length(MGI_splits))
num_IDs <- lengths(MGI_splits)
MGI_matches <- match(FANTOM5_ann_df[["HGNC/MGI_ID"]], names(mgi_to_entrez_list))
MGI_vec[num_IDs == 1] <- unlist(mgi_to_entrez_list, use.names = FALSE)[MGI_matches[num_IDs == 1]]
MGI_vec[num_IDs > 1] <- vapply(MGI_splits[num_IDs > 1], function(x) {
  my_matches <- match(x, names(mgi_to_entrez_list))
  matched_IDs <- unlist(mgi_to_entrez_list, use.names = FALSE)[my_matches]
  are_NA <- is.na(matched_IDs)
  if (all(are_NA)) {
    return(NA_character_)
  } else {
    matched_IDs <- matched_IDs[!(are_NA)]
    return(paste0(matched_IDs, collapse = ", "))
  }
}, "")

FANTOM5_ann_df[["MGI_entrez"]] <- MGI_vec

FANTOM5_ann_df[["FANTOM5_entrez"]] <- StandardizeFANTOM5IDs(FANTOM5_ann_df[["GeneID"]])

have_comma_FANTOM5 <- grepl(",", FANTOM5_ann_df[["FANTOM5_entrez"]], fixed = TRUE)
have_comma_MGI <- grepl(",", FANTOM5_ann_df[["MGI_entrez"]], fixed = TRUE)
table(have_comma_FANTOM5, have_comma_MGI)

FANTOM5_ann_df[["Consensus_entrez"]] <- ifelse(is.na(FANTOM5_ann_df[["FANTOM5_entrez"]]),
                                               FANTOM5_ann_df[["MGI_entrez"]],
                                               ifelse(have_comma_FANTOM5 & !(have_comma_MGI),
                                                      ifelse(is.na(FANTOM5_ann_df[["MGI_entrez"]]),
                                                             FANTOM5_ann_df[["FANTOM5_entrez"]],
                                                             FANTOM5_ann_df[["MGI_entrez"]]
                                                             ),
                                                      FANTOM5_ann_df[["FANTOM5_entrez"]]
                                                      )
                                               )




# Tidy some additional FANTOM5 annotation columns -------------------------

FANTOM5_ann_df[["MGI_ID"]] <- StandardizeFANTOM5IDs(FANTOM5_ann_df[["HGNC/MGI_ID"]])
FANTOM5_ann_df[["MGI_ID"]] <- sub("MGI:", "", FANTOM5_ann_df[["MGI_ID"]])

FANTOM5_ann_df[["UniProt_ID"]] <- StandardizeFANTOM5IDs(FANTOM5_ann_df[["HGNC/MGI_ID"]])

tidy_annotation_columns <- c(
  "Transcript_name", "HGNC/MGI_ID",
  "UniProt_ID", "Gene_name", "Gene_synonyms", "Gene_source"
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
names(FANTOM5_df)[names(FANTOM5_df) == "HGNC/MGI_ID"] <- "MGI_ID"

FANTOM5_df[["Gene_symbol"]] <- StandardizeFANTOM5IDs(FANTOM5_df[["Gene_symbol"]])

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
  "Group"             = paste0(ifelse(is.na(BioMart_tidied_df[["NCBI gene (formerly Entrezgene) ID"]]),
                                      "",
                                      paste0(BioMart_tidied_df[["NCBI gene (formerly Entrezgene) ID"]], " | ")
                                      ),
                               BioMart_tidied_df[["Gene name"]], " | ",
                               chromosomes_vec
                               ),
  "Entrez_ID"         = BioMart_tidied_df[["NCBI gene (formerly Entrezgene) ID"]],
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




# Save data ---------------------------------------------------------------

save(list = c("combined_TSS_df", "BioMart_filtered_df", "FANTOM5_filtered_df",
              "BioMart_tidied_df", "FANTOM5_df",
              "combined_TSS_summary_df" #, "BioMart_summary_df", "FANTOM5_summary_df"
              ),
     file = file.path(general_RData_directory, "05) Compile TSS (transcription start site) data.RData")
     )





