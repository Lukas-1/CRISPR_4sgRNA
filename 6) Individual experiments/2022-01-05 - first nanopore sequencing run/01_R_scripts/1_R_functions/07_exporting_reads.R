## 2022-04-14



# Load packages and source code -------------------------------------------

library("Biostrings")
CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
pacbio_plate1_directory <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
source(file.path(pacbio_plate1_directory, "1) R functions", "20) Summarizing data across wells.R"))
source(file.path(pacbio_plate1_directory, "1) R functions", "03) Creating reference plasmid sequences.R"))



# Define paths ------------------------------------------------------------

integrated_pacbio_RData_dir <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data", "3) R objects")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")



# Load data ---------------------------------------------------------------

load(file.path(integrated_pacbio_RData_dir, "04) Create reference sequences for each well - sg_sequences_df.RData"))



# Read in data ------------------------------------------------------------

amplicon_ref <- read.table(file.path(first_nanopore_dir, "02_input_data", "amplicon_4sg.txt"),
                           quote = "", stringsAsFactors = FALSE
                           )[, 1]

plasmid_gbk <- ReadInPlasmid_gbk(file.path(first_nanopore_dir, "02_input_data", "pYJA5_with_4sg.gbk"))




# Define functions --------------------------------------------------------

MakeQualityCategories <- function(reads_df, max_CCS = 11L) {
  num_passes <- seq_len(max_CCS)
  min_qualities <- vapply(num_passes, GetMinQuality, numeric(1))
  pass_filters_list <- lapply(seq_along(num_passes), function(x) {
    pass_CCS <- reads_df[, "Num_full_passes"] >= num_passes[[x]]
    pass_quality <- reads_df[, "Read_quality"] >= min_qualities[[x]]
    return(pass_CCS & pass_quality)
  })
  pass_mat <- do.call(cbind, rev(pass_filters_list))
  results_vec <- apply(pass_mat, 1, function(x) if (any(x)) abs((max_CCS + 1L) - which(x)[[1]]) else 0L)
  return(results_vec)
}


AbbreviateSymbols <- function(symbols_vec) {
  have_comma <- grepl(",", symbols_vec, fixed = TRUE)
  splits_list <- strsplit(symbols_vec[have_comma], ", ", fixed = TRUE)
  num_splits <- lengths(splits_list)
  final_vec <- vapply(seq_along(splits_list),
                      function(x) paste0(splits_list[[x]][[1]], "+", num_splits[[x]] - 1),
                      ""
                      )
  results_vec <- symbols_vec
  results_vec[have_comma] <- final_vec
  results_vec <- sub("Control_", "Control-", results_vec, fixed = TRUE)
  return(results_vec)
}


CreateReadDescriptors <- function(reads_df, mapped_df) {

  is_pacbio <- "Num_full_passes" %in% names(reads_df)
  if (is_pacbio) {
    qualities_vec <- MakeQualityCategories(reads_df)
  }

  symbol_mat <- as.matrix(mapped_df[, paste0("Symbol_sg", 1:4)])
  for (column_index in seq_len(ncol(symbol_mat))) {
    symbol_mat[, column_index] <- AbbreviateSymbols(symbol_mat[, column_index])
  }
  symbol_list <- as.list(data.frame(t(symbol_mat), stringsAsFactors = FALSE))
  symbol_vec <- vapply(symbol_list, function(x) {
    if (length(unique(x)) == 1) {
      paste0(x[[1]], "_all4")
    } else {
      paste0(x, collapse = "_")
    }
  }, "")

  read_numbers <- formatC(seq_len(nrow(reads_df)),
                          width = nchar(as.character(nrow(reads_df))),
                          format = "d", flag = "0"
                          )
  matches_vec <- match(seq_len(nrow(reads_df)), mapped_df[, "Read_number"])
  stopifnot(nrow(reads_df) >= max(mapped_df[, "Read_number"]))
  plasmid_vec <- ifelse(is.na(matches_vec), "_no_match", symbol_vec[matches_vec])
  results_vec <- paste0("R", read_numbers, "_",
                        if (is_pacbio) paste0("CCS", qualities_vec, "_") else "",
                        plasmid_vec
                        )
  return(results_vec)
}


SelectReads <- function(reads_df, reads_indices) {
  sequence_vec <- DNAStringSet(reads_df[, "Sequence"][reads_indices])
  quality_vec <- PhredQuality(reads_df[, "Quality"][reads_indices])
  results_vec <- QualityScaledBStringSet(sequence_vec, quality_vec)
  names(results_vec) <- reads_df[, "Read_string"][reads_indices]
  return(results_vec)
}



ExportReads <- function(reads_object, output_dir, file_name, max_length = 20000L) {
  are_too_long <- lengths(reads_object) > max_length
  stopifnot(!(all(are_too_long)))
  if (any(are_too_long)) {
    long_IDs <- names(reads_object)[are_too_long]
    write_message <- paste0("The following reads were too long (>",
                            max_length, " bp) and could not",
                            " be exported: ", paste0(long_IDs, collapse = ", ")
                            )
    write.table(write_message,
                file = file.path(output_dir, paste0(file_name, " - log.txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )
  }
  writeQualityScaledXStringSet(reads_object[!(are_too_long)],
                               filepath = file.path(output_dir, paste0(file_name, ".fastq"))
                               )
  return(invisible(NULL))
}



IndicesForGene <- function(entrez_or_symbol, mapped_df) {
  if (is.numeric(entrez_or_symbol)) {
    genes_df <- mapped_df[, paste0("Entrez_sg", 1:4)]
  } else {
    genes_df <- mapped_df[, paste0("Symbol_sg", 1:4)]
  }
  use_regex <- paste0("(^|, )", entrez_or_symbol, "($|, )")
  are_match_mat <- do.call(cbind, lapply(1:4, function(x) {
    grepl(use_regex, genes_df[, x])
  }))
  have_match <- rowSums(are_match_mat) >= 1
  which(have_match)
}


ExportReadsForGene <- function(entrez_or_symbol, reads_df, mapped_df, output_dir = fastq_dir) {
  gene_indices <- IndicesForGene(entrez_or_symbol, mapped_df)
  read_numbers <- mapped_df[, "Read_number"][gene_indices]
  gene_reads <- SelectReads(reads_df, read_numbers)
  gene_file_name <- paste0("Reads_for_", entrez_or_symbol)
  ExportReads(gene_reads, output_dir, gene_file_name)
}



TidyCRISPRaDf <- function(sg_df) {

  CRISPRa_df <- sg_df[sg_df[, "Modality"] %in% "CRISPRa", ]
  row.names(CRISPRa_df) <- NULL

  are_5plus <- (CRISPRa_df[, "Plate_4sg"] == 1) &
               (CRISPRa_df[, "Plate_name"] == "HA_5plus")
  plate_names <- ifelse(are_5plus, "HA5+", paste0("HA", CRISPRa_df[, "Plate_4sg"]))

  num_occurrences <- table(CRISPRa_df[, "Gene_symbol"])[CRISPRa_df[, "Gene_symbol"]]

  symbols_vec <- paste0(CRISPRa_df[, "Gene_symbol"], "_TSS", CRISPRa_df[, "TSS_number"])
  symbols_vec <- ifelse(num_occurrences == 1, CRISPRa_df[, "Gene_symbol"], symbols_vec)
  symbols_vec <- ifelse(is.na(symbols_vec),
                        sub("Control_", "Control-", CRISPRa_df[, "Target_ID"], fixed = TRUE),
                        symbols_vec
                        )

  CRISPRa_df[, "Plasmid_ID"] <- paste0(symbols_vec, "_",
                                       sub("_", "", plate_names, fixed = TRUE),
                                       "_Well", formatC(CRISPRa_df[, "Well_4sg"], width = 3, flag = "0")
                                       )
  CRISPRa_df[, "Num_plasmids_for_gene"] <- num_occurrences
  return(CRISPRa_df)
}




ExportReferences <- function(gene_symbol) {
  stopifnot(all(c("CRISPRa_df", "plasmid_lines_list", "plasmid_gbk") %in% ls(envir = globalenv())))
  stopifnot(identical(nrow(CRISPRa_df), length(plasmid_lines_list)))
  symbols_vec <- ifelse(is.na(CRISPRa_df[, "Gene_symbol"]),
                        CRISPRa_df[, "Plasmid_ID"],
                        CRISPRa_df[, "Gene_symbol"]
                        )
  are_this_symbol <- symbols_vec %in% gene_symbol
  plasmid_names <- CRISPRa_df[, "Plasmid_ID"][are_this_symbol]
  message("Exporting reference sequences for the following plasmids: ",
          paste0(plasmid_names, collapse = ", "), "!"
          )
  for (i in seq_along(plasmid_names)) {
    file_name <- paste0(plasmid_names[[i]], "_4sg.gbk")
    use_index <- which(are_this_symbol)[[i]]
    export_lines <- c(plasmid_gbk[["preceding"]],
                      plasmid_lines_list[[use_index]],
                      plasmid_gbk[["following"]]
                      )
    write.table(x         = export_lines,
                file      = file.path(reference_dir, file_name),
                col.names = FALSE,
                row.names = FALSE,
                quote     = FALSE
                )
  }
}






