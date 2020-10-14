### 13th December 2019 ###




# Import packages and source code -----------------------------------------

library("BSgenome.Mmusculus.UCSC.mm10")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory  <- file.path(RData_directory, "7) Mouse - CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRa_RData_directory, "01) Compile predefined CRISPRa libraries - CRISPRa_df.RData"))





# Load global variables ---------------------------------------------------

mouse_genes_mm10_GRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
rm(human_genes_GRanges)





# Process the genomic locations of sgRNAs from mCRISPRa-v2 ----------------

have_locations <- !(is.na(CRISPRa_df[["mCRISPRa_v2_ID"]])) & (CRISPRa_df[["Is_control"]] == "No")

mCRISPRa_v2_IDs_vec <- CRISPRa_df[["mCRISPRa_v2_ID"]][have_locations]

have_two_entries <- grepl(";", mCRISPRa_v2_IDs_vec, fixed = TRUE)
two_entries_splits <- strsplit(mCRISPRa_v2_IDs_vec[have_two_entries], "; ", fixed = TRUE)
mCRISPRa_v2_IDs_vec[have_two_entries] <- sapply(two_entries_splits, "[[", 1)

mCRISPRa_v2_IDs_second_half <- sapply(strsplit(mCRISPRa_v2_IDs_vec, "_[-+]_"), "[[", 2)
mCRISPRa_v2_location_vec <- as.integer(sapply(strsplit(mCRISPRa_v2_IDs_second_half, "-", fixed = TRUE), "[[", 1))

mCRISPRa_v2_ID_splits <- strsplit(mCRISPRa_v2_IDs_vec, "_", fixed = TRUE)
mCRISPRa_v2_strand_vec <- vapply(mCRISPRa_v2_ID_splits, function(x) x[x %in% c("-", "+")], "")




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
  CRISPRa_df[have_locations, c("Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "mCRISPRa_v2_transcript")],
  "Original_sequence" = CRISPRa_df[["sgRNA_sequence"]][have_locations],
  "Location"          = mCRISPRa_v2_location_vec,
  "Strand"            = c("+" = "-", "-" = "+")[mCRISPRa_v2_strand_vec],
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
chromosome_lengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
chromosome_lengths <- chromosome_lengths[names(chromosome_lengths) %in% locations_df[["Chromosome"]]]
chromosome_lengths_vec <- chromosome_lengths[locations_df[["Chromosome"]]]
exceed_length <- locations_df[["Start"]] > chromosome_lengths_vec

# Exclude impossible locations, and also guides whose associated Entrez ID is "NA"
found_locations_df <- locations_df[exceed_length %in% FALSE, ]





# Annotate these locations (using a liftOver from hg19 to hg38) -----------

retrieved_df <- FindOverlappingGenes(found_locations_df, gene_models_GRanges = mouse_genes_mm10_GRanges, is_mouse = TRUE)
retrieved_df[["Retrieved_sequence"]] <- RetrieveSequences(retrieved_df, BSgenome.Mmusculus.UCSC.mm10)
retrieved_df[["Retrieved_PAM"]] <- GetNGGPAM(retrieved_df, BSgenome.Mmusculus.UCSC.mm10)

are_confirmed <- toupper(substr(retrieved_df[["Original_sequence"]], 2, 20)) == substr(retrieved_df[["Retrieved_sequence"]], 2, 20)

confirmed_df <- retrieved_df[are_confirmed, ]



## Examine some problematic cases

selected_columns <- c("mCRISPRa_v2_transcript", "Gene_symbol",
                      "Overlapping_symbols", "Original_symbol",
                      "Entrez_ID", "Overlapping_Entrez_IDs",
                      "Chromosome", "Strand", "Location",
                      "Original_sequence", "Retrieved_sequence",
                      "Num_overlapping_genes"
                      )

any_confirmed <- sapply(unique(retrieved_df[["Combined_ID"]]),
                        function(x) any(are_confirmed[retrieved_df[["Combined_ID"]] == x])
                        )

retrieved_df[retrieved_df[["Combined_ID"]] %in% names(any_confirmed)[!(any_confirmed)], selected_columns]

unique(CRISPRa_df[is.na(CRISPRa_df[["Entrez_ID"]]) & have_locations, c("Original_symbol", "mCRISPRa_v2_transcript")])









# Disambiguate Entrez IDs (based on the location of sgRNAs) ---------------

confirmed_df <- ReassignEntrezsByLocations(confirmed_df)


# Check the results of disambiguation
per_gene_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "mCRISPRa_v2_transcript", "Chromosome")
unique(confirmed_df[confirmed_df[["Entrez_assignment"]] %in% "Ambiguous chromosome, and ambiguous overlaps", per_gene_columns])
unique(confirmed_df[confirmed_df[["Entrez_assignment"]] %in% "Overlaps with gene", per_gene_columns])






# Modify CRISPRa_df for sgRNAs from mCRISPRa-v2 ---------------------------

CRISPRa_df_IDs <- paste0(CRISPRa_df[["Combined_ID"]][have_locations], "__", CRISPRa_df[["sgRNA_sequence"]][have_locations])
confirmed_df_IDs <- paste0(confirmed_df[["Combined_ID"]], "__", confirmed_df[["Original_sequence"]])

ID_matches <- match(CRISPRa_df_IDs, confirmed_df_IDs)

unique(CRISPRa_df[have_locations, ][is.na(ID_matches), c("Entrez_ID", "Gene_symbol", "Original_symbol", "mCRISPRa_v2_transcript")])


new_entrezs_vec <- ifelse(is.na(ID_matches),
                          CRISPRa_df[["Entrez_ID"]][have_locations],
                          confirmed_df[["Entrez_ID"]][ID_matches]
                          )

CRISPRa_df[["Entrez_ID"]][have_locations] <- new_entrezs_vec
CRISPRa_df[["Combined_ID"]][have_locations] <- ifelse(is.na(CRISPRa_df[["Entrez_ID"]][have_locations]),
                                                    toupper(CRISPRa_df[["Original_symbol"]][have_locations]),
                                                    CRISPRa_df[["Entrez_ID"]][have_locations]
                                                    )

are_to_replace_symbols <- !(is.na(ID_matches)) & (confirmed_df[["Were_replaced"]][ID_matches] %in% TRUE)

if (any(are_to_replace_symbols)) {
  new_symbols_vec <- CRISPRa_df[["Gene_symbol"]][have_locations]
  new_symbols_vec[are_to_replace_symbols] <- MapToEntrezs(entrez_IDs_vec = new_entrezs_vec[are_to_replace_symbols], is_mouse = TRUE)[["Gene_symbol"]]
  CRISPRa_df[["Gene_symbol"]][have_locations] <- new_symbols_vec
  CRISPRa_df[["Original_symbol"]][have_locations][are_to_replace_symbols]
}

new_sequences_vec <- ifelse(is.na(ID_matches),
                            CRISPRa_df[["sgRNA_sequence"]][have_locations],
                            paste0(substr(confirmed_df[["Retrieved_sequence"]][ID_matches], 1, 1),
                                   substr(CRISPRa_df[["sgRNA_sequence"]][have_locations], 2, 20)
                                   )
                            )
CRISPRa_df[["sgRNA_sequence"]][have_locations] <- new_sequences_vec

CRISPRa_df[["Exchanged_5pG"]] <- NA_character_
CRISPRa_df[["Exchanged_5pG"]][have_locations] <- ifelse(is.na(ID_matches), "No", "Yes")



CRISPRa_df[["Original_PAM"]][have_locations] <- ifelse(is.na(CRISPRa_df[["Original_PAM"]][have_locations]),
                                                       confirmed_df[["Retrieved_PAM"]][ID_matches],
                                                       CRISPRa_df[["Original_PAM"]][have_locations]
                                                       )





# Remove new duplicates that resulted from replacing 5' nucleotides -------

CRISPRa_df <- ResolveDuplicates(CRISPRa_df, concatenate_columns = c("Sublibrary", "mCRISPRa_v2_ID", "mCRISPRa_TSS_source"))





# Save data ---------------------------------------------------------------

save(list = "CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from mCRISPRa-v2 - CRISPRa_df.RData")
     )




