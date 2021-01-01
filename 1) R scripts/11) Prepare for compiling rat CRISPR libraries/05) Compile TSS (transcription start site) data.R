### 18th November 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "28) Merging the FANTOM5 and BioMart TSS data.R"))
source(file.path(general_functions_directory, "23) Translating between Ensembl IDs, gene symbols and Entrez IDs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
rat_genome_directory    <- file.path(CRISPR_root_directory, "2) Input data", "Rat genome")
FANTOM5_input_directory <- file.path(rat_genome_directory, "FANTOM5_liftover")
Ensembl_input_directory <- file.path(rat_genome_directory, "Ensembl")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "10) Rat - General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Rn.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Read in data ------------------------------------------------------------

# The two FANTOM5 files were downloaded from: https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/
# on 17 November 2020
FANTOM5_ann_df <- read.table(file.path(FANTOM5_input_directory, "rn6.cage_peak_ann.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE,
                             header = TRUE, row.names = NULL, fill = TRUE,
                             check.names = FALSE, comment.char = ""
                             )

FANTOM5_bed_df <- read.table(file.path(FANTOM5_input_directory, "rn6.cage_peak_coord.bed"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE,
                             header = FALSE, row.names = NULL
                             )

# The BioMart file was downloaded from https://www.ensembl.org/biomart/martview
BioMart_df <- read.table(file.path(rat_genome_directory, "Ensembl", "BioMart_rat_2020-11-17_mart_export.txt"),
                         sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                         )





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






# Check FANTOM5 CAGE peak IDs ---------------------------------------------

table(FANTOM5_ann_df[["CAGE_peaks_ID"]] %in% FANTOM5_bed_df[["Peak_ID"]])
table(FANTOM5_bed_df[["Peak_ID"]] %in% FANTOM5_ann_df[["CAGE_peaks_ID"]])

any(duplicated(FANTOM5_ann_df[["CAGE_peaks_ID"]]))
any(duplicated(FANTOM5_bed_df[["Peak_ID"]]))





# Modify the FANTOM5 annotation data frame --------------------------------

names(FANTOM5_ann_df)[names(FANTOM5_ann_df) == "Ensembl_gene_id"] <- "Ensembl_gene_ID"
names(FANTOM5_ann_df)[names(FANTOM5_ann_df) == "Ensembl_gene_name"] <- "Gene_symbol"

FANTOM5_ann_df[["Ensembl_gene_ID"]] <- sapply(strsplit(FANTOM5_ann_df[["Ensembl_gene_ID"]], ".", fixed = TRUE), "[[", 1)

FANTOM5_ann_df[["Ensembl_gene_ID"]] <- ifelse(FANTOM5_ann_df[["Ensembl_gene_ID"]] == "",
                                              NA_character_,
                                              FANTOM5_ann_df[["Ensembl_gene_ID"]]
                                              )

FANTOM5_ann_df[["Gene_symbol"]] <- ifelse(FANTOM5_ann_df[["Gene_symbol"]] == ".",
                                          NA_character_,
                                          FANTOM5_ann_df[["Gene_symbol"]]
                                          )


FANTOM5_ensembl_df <- MapEnsemblIDs(FANTOM5_ann_df, use_dataset = "rnorvegicus_gene_ensembl")
names(FANTOM5_ensembl_df)[names(FANTOM5_ensembl_df) == "OrgHs_entrez"] <- "OrgRn_entrez"
names(FANTOM5_ensembl_df)[names(FANTOM5_ensembl_df) == "OrgHs_num_mappings"] <- "OrgRn_num_mappings"




# Build a combined FANTOM5 data frame -------------------------------------

FANTOM_ann_matches <- match(FANTOM5_bed_df[["Peak_ID"]], FANTOM5_ann_df[["CAGE_peaks_ID"]])

FANTOM5_df <- data.frame(
  FANTOM5_ann_df[FANTOM_ann_matches, "CAGE_peaks_ID", drop = FALSE],
  "Entrez_ID" = FANTOM5_ensembl_df[["Consensus_entrez"]][FANTOM_ann_matches],
  FANTOM5_ensembl_df[FANTOM_ann_matches, "Gene_symbol", drop = FALSE],
  FANTOM5_bed_df[, c("Chromosome", "Strand", "Peak_start", "Peak_end", "TSS_start", "TSS_stop")],
  FANTOM5_bed_df["Score"],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

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





