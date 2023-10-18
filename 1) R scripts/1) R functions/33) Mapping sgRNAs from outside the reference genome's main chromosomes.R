### 1st March 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))



# Load global variables ---------------------------------------------------

chromosome_names <- names(BSgenome.Hsapiens.UCSC.hg38)
message("Loading the human genome into RAM...")
chromosome_sequences_list <- sapply(chromosome_names, function(x) BSgenome.Hsapiens.UCSC.hg38[[x]], simplify = FALSE)



# Define functions --------------------------------------------------------

MapUnlocatedSgRNAs <- function(CRISPR_df) {

  ### Search alternative and unmapped chromosomal sequences for sgRNA matches
  have_no_location <- is.na(CRISPR_df[, "Chromosome"]) &
                      (CRISPR_df[, "Is_control"] %in% "No")
  unlocated_sgRNAs <- toupper(CRISPR_df[, "sgRNA_sequence"][have_no_location])
  # stopifnot(!(any(duplicated(unlocated_sgRNAs))))
  original_unlocated_df <- FindSequences(unlocated_sgRNAs)

  ## Summarize sgRNA matches
  unlocated_df <- original_unlocated_df[original_unlocated_df[["Num_MM"]] == 0, ]
  unlocated_df[["PAM"]] <- GetNGGPAM(unlocated_df)

  unlocated_search_df <- SummarizeFoundSequencesDf(unlocated_df, all_sequences = unlocated_sgRNAs)

  assign("delete_unlocated_search_df", unlocated_search_df, envir = globalenv())

  unlocated_CRISPR_df <- CRISPR_df[have_no_location, ]
  row.names(unlocated_CRISPR_df) <- NULL
  unlocated_CRISPR_df <- ExtendWithGenomeSearch(unlocated_CRISPR_df[, !(names(unlocated_CRISPR_df) %in% names(unlocated_search_df))],
                                                unlocated_search_df,
                                                allow_5pG = FALSE
                                                )

  ## Integrate data
  select_CRISPR_columns <- c(
    "Source", "Entrez_ID", "Gene_symbol", "Original_symbol",
    "sgRNA_sequence",
    "Entrez_chromosome", "Discordant_locations",
    "Chromosome", "Strand", "Start", "End", "PAM", "Num_0MM", "Is_control"
  )
  select_CRISPR_columns <- intersect(select_CRISPR_columns, names(CRISPR_df))
  select_search_columns <- c(
    "PAM", "Hits_chromosome", "Hits_strand",
    "Hits_start", "Hits_end", "Num_0MM", "Locations_0MM"
  )

  unlocated_CRISPR_df <- unlocated_CRISPR_df[, select_search_columns]
  for (search_column in c("PAM", "Num_0MM", "Locations_0MM")) {
    names(unlocated_CRISPR_df)[names(unlocated_CRISPR_df) == search_column] <- paste0("Hits_", search_column)
  }
  show_df <- data.frame(
    CRISPR_df[have_no_location, select_CRISPR_columns],
    unlocated_CRISPR_df, stringsAsFactors = FALSE
  )
  show_df <- show_df[show_df[["Is_control"]] %in% "No", ]
  show_df <- show_df[, names(show_df) != "Is_control"]

  ## Compile the locations
  hits_splits <- strsplit(show_df[, "Hits_Locations_0MM"], "; ", fixed = TRUE)
  first_hits <- sapply(hits_splits, "[[", 1)
  have_hits <- !(is.na(first_hits))
  first_hits_df <- LocationStringToDf(first_hits[have_hits])
  indices_vec <- rep(NA, length(first_hits))
  indices_vec[have_hits] <- seq_len(sum(have_hits))
  first_hits_df <- first_hits_df[indices_vec, ]

  location_columns <- c("Chromosome", "Strand", "Start", "End")
  lax_df <- CRISPR_df[, location_columns]
  for (location_column in location_columns) {
    lax_df[, location_column][have_no_location] <- first_hits_df[, location_column]
  }
  names(lax_df) <- paste0(names(lax_df), "_all")

  results_df <- data.frame(CRISPR_df, lax_df, stringsAsFactors = FALSE)
  return(results_df)
}


ReplaceUnlocatedgRNAs <- function(CRISPR_df) {
  for (location_column in c("Chromosome", "Strand", "Start", "End")) {
    all_column <- paste0(location_column, "_all")
    CRISPR_df[[location_column]] <- CRISPR_df[[all_column]]
    CRISPR_df[[all_column]] <- NULL
  }
  return(CRISPR_df)
}



