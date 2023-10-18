### 13th December 2019 ###




# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("liftOver")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
liftOver_input_directory <- file.path(CRISPR_input_directory, "Human genome", "liftOver")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

hg19tohg38_chain <- import.chain(file.path(liftOver_input_directory, "hg19ToHg38.over.chain"))

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData"))





# Load global variables ---------------------------------------------------

human_genes_hg19_GRanges <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)






# Process the genomic locations of sgRNAs from hCRISPRa-v2 ----------------

have_locations <- !(is.na(CRISPRa_df[["hCRISPRa_v2_ID"]])) & (CRISPRa_df[["Is_control"]] == "No")

hCRISPRa_v2_IDs_vec <- CRISPRa_df[["hCRISPRa_v2_ID"]][have_locations]

have_two_entries <- grepl(";", hCRISPRa_v2_IDs_vec, fixed = TRUE)
two_entries_splits <- strsplit(hCRISPRa_v2_IDs_vec[have_two_entries], "; ", fixed = TRUE)
hCRISPRa_v2_IDs_vec[have_two_entries] <- sapply(two_entries_splits, "[[", 1)

hCRISPRa_v2_IDs_second_half <- sapply(strsplit(hCRISPRa_v2_IDs_vec, "_[-+]_"), "[[", 2)
hCRISPRa_v2_location_vec <- as.integer(sapply(strsplit(hCRISPRa_v2_IDs_second_half, "-", fixed = TRUE), "[[", 1))

hCRISPRa_v2_ID_splits <- strsplit(hCRISPRa_v2_IDs_vec, "_", fixed = TRUE)
hCRISPRa_v2_strand_vec <- vapply(hCRISPRa_v2_ID_splits, function(x) x[x %in% c("-", "+")], "")




# Map Entrez IDs to (one or more) chromosomes -----------------------------

unique_entrez_ID_strings <- unique(CRISPRa_df[["Entrez_ID"]][have_locations])
expanded_entrez_IDs_df <- ExpandList(strsplit(unique_entrez_ID_strings, ", ", fixed = TRUE))
entrez_matches <- match(expanded_entrez_IDs_df[["Value"]], entrez_to_symbol_df[["Entrez_ID"]])
expanded_entrez_IDs_df[["Chromosome"]] <- entrez_to_symbol_df[["Chromosome"]][entrez_matches]
chromosomes_list <- split(expanded_entrez_IDs_df[["Chromosome"]], expanded_entrez_IDs_df[["List_index"]])
chromosomes_list <- lapply(chromosomes_list, function(x) sort(unique(unlist(strsplit(x, ", ", fixed = TRUE), use.names = FALSE))))

chromosomes_df <- ExpandList(chromosomes_list)
chromosomes_df[["Entrez_ID"]] <- unique_entrez_ID_strings[as.integer(names(chromosomes_list)[chromosomes_df[["List_index"]]])]
names(chromosomes_df)[[1]] <- "Chromosome"




# Create a data frame of possible locations -------------------------------

locations_short_df <- data.frame(
  CRISPRa_df[have_locations, c("Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "hCRISPRa_v2_transcript")],
  "Original_sequence" = CRISPRa_df[["sgRNA_sequence"]][have_locations],
  "Location"          = hCRISPRa_v2_location_vec,
  "Strand"            = c("+" = "-", "-" = "+")[hCRISPRa_v2_strand_vec],
  stringsAsFactors    = FALSE,
  row.names           = NULL
)

locations_df_list <- lapply(seq_len(nrow(locations_short_df)), function(x) {
  this_entrez <- locations_short_df[["Entrez_ID"]][[x]]
  if (is.na(this_entrez)) {
    results_df <- data.frame(locations_short_df[x, ], "Chromosome" = NA_character_, stringsAsFactors = FALSE)
  } else {
    chromosome_vec <- chromosomes_df[["Chromosome"]][chromosomes_df[["Entrez_ID"]] == this_entrez]
    results_df <- data.frame(locations_short_df[rep(x, length(chromosome_vec)), ],
                             "Chromosome" = chromosome_vec,
                             stringsAsFactors = FALSE
                             )
  }
  return(results_df)
})

locations_df <- do.call(rbind.data.frame, c(locations_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

# Convert the single location to start and end coordinates
are_neg_strand <- locations_df[["Strand"]] == "-"
locations_df[["Start"]] <- ifelse(are_neg_strand, locations_df[["Location"]] + 4L,  locations_df[["Location"]] - 21L)
locations_df[["End"]]   <- ifelse(are_neg_strand, locations_df[["Location"]] + 23L, locations_df[["Location"]] - 2L)

# Filter out chromosomal locations that are impossible (would lie outside the chromosome)
chromosome_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
chromosome_lengths <- chromosome_lengths[names(chromosome_lengths) %in% locations_df[["Chromosome"]]]
chromosome_lengths_vec <- chromosome_lengths[locations_df[["Chromosome"]]]
exceed_length <- locations_df[["Start"]] > chromosome_lengths_vec

# Exclude impossible locations, and also guides whose associated Entrez ID is "NA"
found_locations_df <- locations_df[exceed_length %in% FALSE, ]





# Annotate these locations (using a liftOver from hg19 to hg38) -----------

location_annotations_df <- LiftOverAndAnnotate(found_locations_df)
gene_overlaps_df <- FindOverlappingGenes(found_locations_df, gene_models_GRanges = human_genes_hg19_GRanges)

lifted_df <- data.frame(location_annotations_df, gene_overlaps_df, stringsAsFactors = FALSE, row.names = NULL)

are_confirmed <- toupper(substr(lifted_df[["Original_sequence"]], 2, 20)) == substr(lifted_df[["Sequence_hg19"]], 2, 20)

confirmed_lifted_df <- lifted_df[are_confirmed, ]



## Examine some problematic cases

selected_columns <- c("hCRISPRa_v2_transcript", "Gene_symbol", "Overlapping_symbols", "Original_symbol",
                      "Entrez_ID", "Overlapping_Entrez_IDs",
                      "Chromosome", "Strand", "Location",
                      "Original_sequence", "Sequence_hg19", "Sequence_liftOver",
                      "Num_overlapping_genes"
                      )

any_confirmed <- sapply(unique(lifted_df[["Combined_ID"]]), function(x) any(are_confirmed[lifted_df[["Combined_ID"]] == x]))

lifted_df[lifted_df[["Combined_ID"]] %in% names(any_confirmed)[!(any_confirmed)], selected_columns]

unique(CRISPRa_df[is.na(CRISPRa_df[["Entrez_ID"]]) & have_locations, c("Original_symbol", "hCRISPRa_v2_transcript")])









# Disambiguate Entrez IDs (based on the location of sgRNAs) ---------------

confirmed_lifted_df <- ReassignEntrezsByLocations(confirmed_lifted_df)


# Check the results of disambiguation
per_gene_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "hCRISPRa_v2_transcript", "Chromosome")
unique(confirmed_lifted_df[confirmed_lifted_df[["Entrez_assignment"]] %in% "Ambiguous chromosome, and ambiguous overlaps", per_gene_columns])
unique(confirmed_lifted_df[confirmed_lifted_df[["Entrez_assignment"]] %in% "Overlaps with gene", per_gene_columns])






# Put the location determined by liftOver into a single column ------------

liftOver_columns <- c("Chromosome_liftOver", "Strand_liftOver", "Start_liftOver", "End_liftOver")

confirmed_lifted_df[["Location_liftOver"]] <- MakeLocationStrings(setNames(confirmed_lifted_df[, liftOver_columns],
                                                                           c("Chromosome", "Strand", "Start", "End")
                                                                           )
                                                                  )



# Modify CRISPRa_df for sgRNAs from hCRISPRa-v2 ---------------------------

CRISPRa_df_IDs <- paste0(CRISPRa_df[["Combined_ID"]][have_locations], "__", CRISPRa_df[["sgRNA_sequence"]][have_locations])
lifted_df_IDs <- paste0(confirmed_lifted_df[["Combined_ID"]], "__", confirmed_lifted_df[["Original_sequence"]])

ID_matches <- match(CRISPRa_df_IDs, lifted_df_IDs)

unique(CRISPRa_df[have_locations, ][is.na(ID_matches), c("Entrez_ID", "Gene_symbol", "Original_symbol", "hCRISPRa_v2_transcript")])


new_entrezs_vec <- ifelse(is.na(ID_matches),
                          CRISPRa_df[["Entrez_ID"]][have_locations],
                          confirmed_lifted_df[["Entrez_ID"]][ID_matches]
                          )

CRISPRa_df[["Entrez_ID"]][have_locations] <- new_entrezs_vec
CRISPRa_df[["Combined_ID"]][have_locations] <- ifelse(is.na(CRISPRa_df[["Entrez_ID"]][have_locations]),
                                                    toupper(CRISPRa_df[["Original_symbol"]][have_locations]),
                                                    CRISPRa_df[["Entrez_ID"]][have_locations]
                                                    )

are_to_replace_symbols <- !(is.na(ID_matches)) & (confirmed_lifted_df[["Were_replaced"]][ID_matches] %in% TRUE)

if (any(are_to_replace_symbols)) {
  new_symbols_vec <- CRISPRa_df[["Gene_symbol"]][have_locations]
  new_symbols_vec[are_to_replace_symbols] <- MapToEntrezs(entrez_IDs_vec = new_entrezs_vec[are_to_replace_symbols])[["Gene_symbol"]]
  CRISPRa_df[["Gene_symbol"]][have_locations] <- new_symbols_vec
  CRISPRa_df[["Original_symbol"]][have_locations][are_to_replace_symbols]
}

new_sequences_vec <- ifelse(is.na(ID_matches),
                            CRISPRa_df[["sgRNA_sequence"]][have_locations],
                            paste0(substr(confirmed_lifted_df[["Sequence_hg19"]][ID_matches], 1, 1),
                                   substr(CRISPRa_df[["sgRNA_sequence"]][have_locations], 2, 20)
                                   )
                            )
CRISPRa_df[["sgRNA_sequence"]][have_locations] <- new_sequences_vec

CRISPRa_df[["Exchanged_5pG"]] <- NA_character_
CRISPRa_df[["Exchanged_5pG"]][have_locations] <- ifelse(is.na(ID_matches), "No", "Yes")

CRISPRa_df[["Location_liftOver"]] <- NA_character_
CRISPRa_df[["Location_liftOver"]][have_locations] <- confirmed_lifted_df[["Location_liftOver"]][ID_matches]


CRISPRa_df[["Original_PAM"]][have_locations] <- ifelse(is.na(CRISPRa_df[["Original_PAM"]][have_locations]),
                                                       confirmed_lifted_df[["PAM_liftOver"]][ID_matches],
                                                       CRISPRa_df[["Original_PAM"]][have_locations]
                                                       )





# Remove new duplicates that resulted from replacing 5' nucleotides -------

CRISPRa_df <- ResolveDuplicates(CRISPRa_df, concatenate_columns = c("Sublibrary", "hCRISPRa_v2_ID", "hCRISPRa_TSS_source"))





# Save data ---------------------------------------------------------------

save(list = "CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRa-v2 - CRISPRa_df.RData")
     )




