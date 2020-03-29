### 21st July 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
FANTOM5_input_directory <- file.path(CRISPR_input_directory, "Human genome", "FANTOM5_liftover")
Ensembl_input_directory <- file.path(CRISPR_input_directory, "Human genome", "Ensembl")
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))





# Read in data ------------------------------------------------------------

# The two FANTOM5 files were downloaded from: https://figshare.com/articles/Re-processing_of_the_data_generated_by_the_FANTOM5_project_hg38_v3_CAGE_peaks/4880063/4
# on 21 July 2019

FANTOM5_ann_df <- read.table(file.path(FANTOM5_input_directory, "hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, fill = TRUE, check.names = FALSE, comment.char = ""
                             )

FANTOM5_bed_df <- read.table(file.path(FANTOM5_input_directory, "hg38_liftover_CAGE_peaks_phase1and2.bed"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                             )

# The BioMart file was downloaded from https://www.ensembl.org/biomart/martview
BioMart_df     <- read.table(file.path(CRISPR_input_directory, "Human genome", "Ensembl", "BioMart_human_2020-03-25_mart_export.txt"),
                             sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                             )





# Define functions --------------------------------------------------------

TSSWindows <- function(sorted_TSS_vec, total_range = (2500L * 1000L - 2L), TSS_range = 1202L, overlap = 30L) {
  start_vec <- sorted_TSS_vec[[1]] - (TSS_range / 2L)
  num_TSSs <- length(sorted_TSS_vec)
  if (num_TSSs > 1) {
    current_start <- start_vec
    current_max <- start_vec + total_range
    end_vec <- c()
    for (i in 2:num_TSSs) {
      supposed_end <- sorted_TSS_vec[[i]] + (TSS_range / 2L)
      if (supposed_end > current_max) {
        new_end <- sorted_TSS_vec[[i - 1]] + (TSS_range / 2L)
        current_start <- sorted_TSS_vec[[i]] - (TSS_range / 2L)
        if (new_end > current_start) {
          new_end <- new_end - overlap
        }
        current_max <- current_start + total_range
        end_vec <- c(end_vec, new_end)
        start_vec <- c(start_vec, current_start)
      }
    }
    end_vec <- c(end_vec, sorted_TSS_vec[[i]] + (TSS_range / 2L))
  } else {
    end_vec <- sorted_TSS_vec[[1]] + (TSS_range / 2)
  }
  results_mat <- cbind("Start" = start_vec, "End" = end_vec)
  return(results_mat)
}


ListToDf <- function(named_list) {
  if (is.null(names(named_list))) {
    stop("named_list must be a named list!")
  }
  results_df <- data.frame(
    "Names" = rep.int(names(named_list), vapply(named_list, nrow, integer(1))),
    do.call(rbind, named_list),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
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







# Check cage peak IDs -----------------------------------------------------

table(FANTOM5_ann_df[["CAGE_Peak_ID"]] %in% FANTOM5_bed_df[["Peak_ID"]])
table(FANTOM5_bed_df[["Peak_ID"]] %in% FANTOM5_ann_df[["CAGE_Peak_ID"]])

any(duplicated(FANTOM5_ann_df[["CAGE_Peak_ID"]]))
any(duplicated(FANTOM5_bed_df[["Peak_ID"]]))







# Build a combined FANTOM data frame --------------------------------------

FANTOM_ann_matches <- match(FANTOM5_bed_df[["Peak_ID"]], FANTOM5_ann_df[["CAGE_Peak_ID"]])

FANTOM5_df <- data.frame(
  FANTOM5_ann_df[FANTOM_ann_matches, c("CAGE_Peak_ID", "GeneID", "Gene_symbol")],
  FANTOM5_bed_df[, c("Chromosome", "Strand", "Peak_start", "Peak_end", "TSS_start", "TSS_stop")],
  FANTOM5_ann_df[FANTOM_ann_matches, "Distance", drop = FALSE],
  FANTOM5_bed_df["Score"],
  stringsAsFactors = FALSE
)

names(FANTOM5_df)[names(FANTOM5_df) == "GeneID"] <- "Entrez_ID"

stopifnot(all(grepl("chr", FANTOM5_df[["Chromosome"]], fixed = TRUE)))





# Perform checks on the combined FANTOM data frame ------------------------

unique_IDs_FANTOM5_df <- unique(FANTOM5_df[, c("Entrez_ID", "Gene_symbol", "Chromosome")])

num_occurrences_gene_symbol <- table(unique_IDs_FANTOM5_df[["Gene_symbol"]])[unique_IDs_FANTOM5_df[["Gene_symbol"]]]
num_occurences_entrez_ID    <- table(unique_IDs_FANTOM5_df[["Entrez_ID"]])[unique_IDs_FANTOM5_df[["Entrez_ID"]]]

unique_IDs_FANTOM5_df[(num_occurrences_gene_symbol > 1) %in% TRUE, ]
unique_IDs_FANTOM5_df[(num_occurences_entrez_ID > 1) %in% TRUE, ]

pasted_IDs <- paste0(unique_IDs_FANTOM5_df[["Entrez_ID"]], "__", unique_IDs_FANTOM5_df[["Gene_symbol"]])
num_occurrences_pasted_IDs <- table(pasted_IDs)[pasted_IDs]

my_order <- order(pasted_IDs)
unique_IDs_FANTOM5_df[my_order, ][num_occurrences_pasted_IDs[my_order] > 1, ]





# Filter and standardize the combined FANTOM5 data frame ------------------

FANTOM5_filtered_df <- FANTOM5_df[!((FANTOM5_df[["Entrez_ID"]] == "") & (FANTOM5_df[["Gene_symbol"]] == "")), ]

StandardizeIDs <- function(char_vec) {
  results_vec <- gsub(" ", ", ", char_vec, fixed = TRUE)
  results_vec <- ifelse(char_vec == "", NA_character_, results_vec)
  return(results_vec)
}

FANTOM5_filtered_df[["Entrez_ID"]] <- StandardizeIDs(FANTOM5_filtered_df[["Entrez_ID"]])
FANTOM5_filtered_df[["Gene_symbol"]] <- StandardizeIDs(FANTOM5_filtered_df[["Gene_symbol"]])

FANTOM5_filtered_df <- data.frame(
  "Group" = paste0(ifelse(FANTOM5_filtered_df[["Entrez_ID"]] == "",
                          "",
                          paste0(FANTOM5_filtered_df[["Entrez_ID"]], " | ")
                          ),
                   FANTOM5_filtered_df[["Gene_symbol"]], " | ",
                   FANTOM5_filtered_df[["Chromosome"]]
                   ),
  FANTOM5_filtered_df,
  stringsAsFactors = FALSE,
  row.names = NULL
)





# Construct a simplified BioMart data frame -------------------------------

are_on_chromosome <- BioMart_df[["Chromosome/scaffold name"]] %in% c(as.character(1:22), "X", "Y", "MT")

BioMart_filtered_df <- BioMart_df[are_on_chromosome, ]

BioMart_filtered_df[["Chromosome/scaffold name"]] <- ifelse(BioMart_filtered_df[["Chromosome/scaffold name"]] == "MT",
                                                            "M",
                                                            BioMart_filtered_df[["Chromosome/scaffold name"]]
                                                            )

BioMart_filtered_df <- data.frame(
  "Group"             = paste0(ifelse(is.na(BioMart_filtered_df[["NCBI gene ID"]]),
                                      "",
                                      paste0(BioMart_filtered_df[["NCBI gene ID"]], " | ")
                                      ),
                               BioMart_filtered_df[["Gene name"]], " | ", "chr",
                               BioMart_filtered_df[["Chromosome/scaffold name"]]
                               ),
  "Entrez_ID"         = BioMart_filtered_df[["NCBI gene ID"]],
  "Gene_symbol"       = BioMart_filtered_df[["Gene name"]],
  "ENSG"              = BioMart_filtered_df[["Gene stable ID"]],
  "ENST"              = BioMart_filtered_df[["Transcript stable ID"]],
  "Chromosome"        = paste0("chr", BioMart_filtered_df[["Chromosome/scaffold name"]]),
  "Strand"            = ifelse(BioMart_filtered_df[["Strand"]] == -1, "-", "+"),
  "Gene_start"        = BioMart_filtered_df[["Gene start (bp)"]],
  "Gene_end"          = BioMart_filtered_df[["Gene end (bp)"]],
  "Transcript_start"  = BioMart_filtered_df[["Transcript start (bp)"]],
  "Transcript_end"    = BioMart_filtered_df[["Transcript end (bp)"]],
  "TSS"               = BioMart_filtered_df[["Transcription start site (TSS)"]],
  "Transcript_length" = BioMart_filtered_df[["Transcript length (including UTRs and CDS)"]],
  stringsAsFactors    = FALSE
)





# Merge the FANTOM5 and BioMart data frames -------------------------------

length(intersect(unique(FANTOM5_filtered_df[["Group"]]), unique(BioMart_filtered_df[["Group"]])))
length(unique(FANTOM5_filtered_df[["Group"]]))
length(unique(BioMart_filtered_df[["Group"]]))

common_columns <- c("Group", "Entrez_ID", "Gene_symbol", "Chromosome", "Strand")


are_duplicated_FANTOM5 <- duplicated(FANTOM5_filtered_df[, c(common_columns, "TSS_start")])
are_duplicated_BioMart <- duplicated(BioMart_filtered_df[, c(common_columns, "TSS")])


table(are_duplicated_FANTOM5)
table(are_duplicated_BioMart)


combined_TSS_df <- rbind.data.frame(
  data.frame(
    "Source" = "FANTOM5",
    FANTOM5_filtered_df[!(are_duplicated_FANTOM5), common_columns],
    "TSS" = FANTOM5_filtered_df[["TSS_start"]][!(are_duplicated_FANTOM5)],
    FANTOM5_filtered_df[!(are_duplicated_FANTOM5), "Score", drop = FALSE],
    stringsAsFactors = FALSE
  ),
  unique(data.frame(
    "Source" = "BioMart",
    BioMart_filtered_df[!(are_duplicated_BioMart), c(common_columns, "TSS")],
    "Score" = NA_integer_,
    stringsAsFactors = FALSE
  ), MARGIN = 1),
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

entrez_to_symbols_vec <- MapToEntrezs(entrez_IDs_vec = combined_TSS_df[["Entrez_ID"]])[["Gene_symbol"]]

combined_TSS_df <- combined_TSS_df[order(GetMinEntrez(combined_TSS_df[["Entrez_ID"]]),
                                         !(mapply(identical, entrez_to_symbols_vec, combined_TSS_df[["Gene_symbol"]])),
                                         combined_TSS_df[["Group"]],
                                         combined_TSS_df[["TSS"]]
                                         ), ]

common_columns <- c(common_columns, "TSS")

are_duplicates_FANTOM5 <- duplicated(combined_TSS_df[, common_columns], fromLast = TRUE)
are_duplicates_BioMart <- duplicated(combined_TSS_df[, common_columns], fromLast = FALSE)

stopifnot(all(combined_TSS_df[["Source"]][are_duplicates_FANTOM5] == "FANTOM5"))
stopifnot(all(combined_TSS_df[["Source"]][are_duplicates_BioMart] == "BioMart"))

combined_TSS_df[["Source"]][are_duplicates_FANTOM5] <- "FANTOM5, BioMart"

combined_TSS_df <- combined_TSS_df[!(are_duplicates_BioMart), ]

row.names(combined_TSS_df) <- NULL








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



# Construct a summarized FANTOM5 data frame -------------------------------

FANTOM5_filtered_mat <- as.matrix(FANTOM5_filtered_df[, c("Peak_start", "Peak_end", "TSS_start", "Score")])

FANTOM5_summary_list <- sapply(unique(FANTOM5_filtered_df[["Group"]]), function(x) {
  print(x)
  are_this_group <- FANTOM5_filtered_df[["Group"]] == x
  sub_mat <- FANTOM5_filtered_mat[are_this_group, , drop = FALSE]
  best_index <- which.max(sub_mat[, "Score"])
  first_index <- which(are_this_group)[[1]]
  results_list <- list(
    "Group"       = x,
    "Entrez_ID"   = FANTOM5_filtered_df[["Entrez_ID"]][[first_index]],
    "Gene_symbol" = FANTOM5_filtered_df[["Gene_symbol"]][[first_index]],
    "Chromosome"  = FANTOM5_filtered_df[["Chromosome"]][[first_index]],
    "Strand"      = FANTOM5_filtered_df[["Strand"]][are_this_group][[best_index]],
    "Best_TSS"    = unname(sub_mat[best_index, "TSS_start"]),
    "First_TSS"   = min(sub_mat[, "TSS_start"]),
    "Last_TSS"    = max(sub_mat[, "TSS_start"]),
    "First_peak"  = min(sub_mat[, "Peak_start"]),
    "Last_peak"   = min(sub_mat[, "Peak_end"])
  )
  return(results_list)

}, simplify = FALSE)

FANTOM5_summary_df <- do.call(rbind.data.frame, c(FANTOM5_summary_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))




# Construct a summarized BioMart data frame -------------------------------

BioMart_filtered_mat <- as.matrix(BioMart_filtered_df[, c("Gene_start", "Gene_end", "Transcript_start", "Transcript_end", "TSS", "Transcript_length")])

BioMart_summary_list <- sapply(unique(BioMart_filtered_df[["Group"]]), function(x) {
  print(x)
  are_this_group <- BioMart_filtered_df[["Group"]] == x
  sub_mat <- BioMart_filtered_mat[are_this_group, , drop = FALSE]
  first_index <- which(are_this_group)[[1]]
  results_list <- list(
    "Group"                  = x,
    "Entrez_ID"              = as.character(BioMart_filtered_df[["Entrez_ID"]][[first_index]]),
    "Gene_symbol"            = BioMart_filtered_df[["Gene_symbol"]][[first_index]],
    "Chromosome"             = BioMart_filtered_df[["Chromosome"]][[first_index]],
    "Strand"                 = BioMart_filtered_df[["Strand"]][[first_index]],
    "First_TSS"              = min(sub_mat[, "TSS"]),
    "Last_TSS"               = max(sub_mat[, "TSS"]),
    "First_transcript_start" = min(sub_mat[, "Transcript_start"]),
    "Last_transcript_start"  = max(sub_mat[, "Transcript_end"]),
    "First_gene_start"       = min(sub_mat[, "Gene_start"]),
    "Last_gene_end"          = max(sub_mat[, "Gene_end"])
  )
  return(results_list)

}, simplify = FALSE)

BioMart_summary_df <- do.call(rbind.data.frame, c(BioMart_summary_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))




# Save data ---------------------------------------------------------------

save(list = c("combined_TSS_df", "BioMart_filtered_df", "FANTOM5_filtered_df",
              "combined_TSS_summary_df", "BioMart_summary_df", "FANTOM5_summary_df"
              ),
     file = file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData")
     )





