## 2022-05-03




# Functions for annotating the library ------------------------------------

FindAffectedGenes <- function(CRISPR_df, sg_number) {
  use_columns <- c("Entrez_IDs", "Gene_symbols", paste0(c("Locations_0MM", "Num_0MM"), "_sg", sg_number))
  input_df <- CRISPR_df[, use_columns]
  names(input_df) <- c("Entrez_ID", "Gene_symbol", "Locations_0MM", "Num_0MM")
  for (location_column in c("Chromosome", "Strand", "Start", "End")) {
    input_df[, location_column] <- NA
  }
  input_df[, "Entrez_ID"] <- gsub(", ", "/", input_df[, "Entrez_ID"], fixed = TRUE)
  nearby_list <- AlignSummaryDf(FindNearbyTSSs, input_df, all_TSS_df)
  include_columns <- c("Affected_Entrez_IDs", "Affected_gene_symbols", "Affects_intended_main_TSS", "Intended_loci")
  affected_df <- nearby_list[["summary_df"]][, include_columns]
  names(affected_df) <- paste0(names(affected_df), paste0("_sg", sg_number))
  return(affected_df)
}



Disambiguate_IDs <- function(input_df, from_symbols_column, sg1_column, sg2_column) {
  from_symbols_splits <- strsplit(input_df[, from_symbols_column], ", ", fixed = TRUE)
  affected_sg1_splits <- strsplit(input_df[, sg1_column], "[,;] ")
  affected_sg2_splits <- strsplit(input_df[, sg2_column], "[,;] ")
  results_list <- mapply(function(x, y, z) {
    results_vec <- intersect(x, c(y, z))
    if (length(results_vec) == 0) {
      results_vec <- x
    }
    return(results_vec)
  }, from_symbols_splits, affected_sg1_splits, affected_sg2_splits, SIMPLIFY = FALSE)
  results_vec <- sapply(results_list, "[[", 1)
  return(results_vec)
}



GetGCcontent <- function(char_vec) {
  char_vec <- toupper(char_vec)
  if (all(substr(char_vec, 1, 1) == "G")) {
    char_vec <- substr(char_vec, 2, nchar(char_vec))
  }
  char_mat <- do.call(rbind, strsplit(char_vec, "", fixed = TRUE))
  are_GC_mat <- (char_mat == "G") | (char_mat == "C")
  num_GC_vec <- as.integer(rowSums(are_GC_mat))
  return(num_GC_vec)
}



ChooseBetweenMultiplePlasmids <- function(input_df) {
  num_plasmids <- table(input_df[, "Entrez_ID"])[as.character(input_df[, "Entrez_ID"])]
  are_first_vec <- ifelse(num_plasmids %in% 1, TRUE, NA)
  are_preferred_vec <- are_first_vec
  TSS_columns <- paste0("Affects_intended_main_TSS_sg", 1:2)
  target_main_TSS <- rowSums(as.matrix(input_df[, TSS_columns]))
  duplicated_entrezs <- unique(input_df[, "Entrez_ID"][(num_plasmids > 1) %in% TRUE])
  for (this_entrez in duplicated_entrezs) {
    are_this_entrez <- input_df[, "Entrez_ID"] %in% this_entrez
    this_seq <- seq_len(sum(are_this_entrez))
    are_first_vec[are_this_entrez] <- this_seq == 1
    are_preferred_vec[are_this_entrez] <- this_seq == which.max(target_main_TSS[are_this_entrez])
  }
  input_df[, "Num_plasmids_for_Entrez"] <- num_plasmids
  input_df[, "Is_first_plasmid"] <- are_first_vec
  input_df[, "Is_preferred_plasmid"] <- are_preferred_vec
  return(input_df)
}




# Functions for preparing input to GuideScan2 -----------------------------

MakeRangesDf <- function(input_vec, add5prime = 0L) {
  are_NA <- is.na(input_vec)
  splits_list <- strsplit(input_vec[!(are_NA)], ":", fixed = TRUE)
  assign("delete_splits_list", splits_list, envir = globalenv())
  start_end_vec <- sapply(splits_list, "[[", 2)
  start_end_splits <- strsplit(start_end_vec, "-", fixed = TRUE)
  results_df <- data.frame(
    "Chromosome" = sapply(splits_list, "[[", 1),
    "Start"      = as.integer(sapply(start_end_splits, "[[", 1)),
    "End"        = as.integer(sapply(start_end_splits, "[[", 2)),
    "Strand"     = sapply(splits_list, "[[", 3),
    stringsAsFactors = FALSE
  )
  indices_vec <- rep(NA, length(input_vec))
  indices_vec[!(are_NA)] <- seq_len(nrow(results_df))
  results_df <- results_df[indices_vec, ]
  row.names(results_df) <- NULL
  return(results_df)
}



GetFirstLocation <- function(input_df, sg_number) {
  intended_vec <- input_df[, paste0("Intended_loci_sg", sg_number)]
  intended_splits <- strsplit(intended_vec, ", ", fixed = TRUE)
  location_splits <- strsplit(input_df[, paste0("Locations_0MM_sg", sg_number)], "; ", fixed = TRUE)
  are_unavailable <- is.na(intended_vec)
  combined_splits <- ifelse(are_unavailable, location_splits, intended_splits)
  first_vec <- sapply(combined_splits, "[[", 1)
  are_to_replace <- are_unavailable & !(is.na(first_vec))
  reformat_splits <- strsplit(first_vec[are_to_replace], ":", fixed = TRUE)
  chr_splits <- strsplit(sapply(reformat_splits, "[[", 1), "(", fixed = TRUE)
  chr_vec <- sapply(chr_splits, "[[", 1)
  strand_vec <- sub(")", "", sapply(chr_splits, "[[", 2), fixed = TRUE)
  first_vec[are_to_replace] <- paste0(chr_vec, ":", sapply(reformat_splits, "[[", 2), ":", strand_vec)
  return(first_vec)
}



ExtractStrand <- function(char_vec) {
  splits_list <- strsplit(char_vec, ".23-", fixed = TRUE)
  first_half_vec <- sapply(splits_list, "[[", 1)
  splits_list <- strsplit(first_half_vec, "_", fixed = TRUE)
  use_indices <- lengths(splits_list) - 1L
  strand_vec <- mapply(function(x, y) x[[y]], splits_list, use_indices)
  strand_vec <- ifelse(strand_vec == "non-targeting", NA_character_, strand_vec)
  return(strand_vec)
}



ExpandPlasmidIDs <- function(input_df) {
  results_df <- data.frame(
    "Short_ID"       = NA,
    "Plasmid_number" = seq_len(nrow(input_df)),
    "Is_NT"          = input_df[, "gene"] == "negative_control",
    input_df,
    stringsAsFactors = FALSE
  )
  plasmid_numbers_vec <- gsub(" ", "0", format(results_df[, "Plasmid_number"]))
  symbols_vec <- ifelse(is.na(input_df[, "Gene_symbol"]),
                        gsub("_", "-", input_df[, "gene"], fixed = TRUE),
                        input_df[, "Gene_symbol"]
                        )
  short_IDs <- paste0(symbols_vec, "_", plasmid_numbers_vec)
  results_df[, "Short_ID"] <- short_IDs
  return(results_df)
}



DataForGuideScan2 <- function(input_df) {
  shared_df <- ExpandPlasmidIDs(input_df)
  sg1_df <- data.frame(
    shared_df,
    "Sg_number" = 1L,
    "Sequence"  = substr(input_df[, "protospacer_A"], 2, 20),
    "Strand"    = ExtractStrand(input_df[, "sgID_A"]),
    "Location"  = GetFirstLocation(input_df, 1),
    stringsAsFactors = FALSE
  )
  sg2_df <- data.frame(
    shared_df,
    "Sg_number" = 2L,
    "Sequence"  = substr(input_df[, "protospacer_B"], 2, 20),
    "Strand"    = ExtractStrand(input_df[, "sgID_B"]),
    "Location"  = GetFirstLocation(input_df, 2),
    stringsAsFactors = FALSE
  )
  combined_df <- rbind.data.frame(sg1_df, sg2_df, stringsAsFactors = FALSE,
                                  make.row.names = FALSE
                                  )
  combined_df <- combined_df[order(combined_df[, "Plasmid_number"]), ]
  row.names(combined_df) <- NULL
  combined_df[, "Sequence"] <- toupper(combined_df[, "Sequence"])
  return(combined_df)
}



Replace5PrimeG <- function(ranges_df) {
  are_pos <- ranges_df[, "Strand"] %in% "+"
  are_neg <- ranges_df[, "Strand"] %in% "-"
  ranges_df[, "Start"][are_pos] <- ranges_df[, "Start"][are_pos] - 1L
  ranges_df[, "End"][are_neg] <- ranges_df[, "End"][are_neg] + 1L
  return(ranges_df)
}



ReformatDf <- function(input_df, add5prime = FALSE) {
  are_valid <- !(is.na(input_df[, "Location"]))
  input_df <- input_df[are_valid, ]
  row.names(input_df) <- NULL
  ranges_df <- MakeRangesDf(input_df[, "Location"])
  if (add5prime) {
    ranges_df <- Replace5PrimeG(ranges_df)
    input_df[, "Sequence"] <- RetrieveSequences(ranges_df)
  }
  combined_IDs <- paste0(input_df[, "Short_ID"],
                         "_sg", input_df[, "Sg_number"]
                         )
  results_df <- data.frame(
    input_df[, names(input_df) != "Strand"],
    "Combined_ID" = combined_IDs,
    ranges_df,
    stringsAsFactors = FALSE
  )
  return(results_df)
}



FormatForGuideScan2 <- function(input_df, add5prime = FALSE) {
  if (!("Chromosome" %in% names(input_df))) {
    input_df <- ReformatDf(input_df, add5prime)
  }
  results_df <- data.frame(
    "id"         = input_df[, "Combined_ID"],
    "sequence"   = input_df[, "Sequence"],
    "pam"        = "NGG",
    "chromosome" = input_df[, "Chromosome"],
    "position"   = ifelse(input_df[, "Strand"] == "-",
                          input_df[, "Start"] - 3L,
                          input_df[, "Start"]
                          ),
    "sense"      = input_df[, "Strand"],
    stringsAsFactors = FALSE
  )
  return(results_df)
}



AddGuideScan2Output <- function(CRISPR_df, output_df, split_string = "_sg", ID_column = "Short_ID") {
  use_columns <- names(CRISPR_df)
  if (!(ID_column %in% names(CRISPR_df))) {
    CRISPR_df <- ExpandPlasmidIDs(CRISPR_df)
  }
  are_unique <- !(duplicated(output_df[, "id"]))
  ID_splits <- strsplit(output_df[, "id"][are_unique], split_string, fixed = TRUE)
  sg_numbers <- as.integer(sapply(ID_splits, "[[", 2))
  IDs_vec <- sapply(ID_splits, "[[", 1)
  spec_vec <- output_df[, "specificity"][are_unique]
  unique_sg_numbers <- sort(unique(sg_numbers))
  spec_mat <- do.call(cbind, lapply(unique_sg_numbers, function(x) {
    are_this_sg <- sg_numbers == x
    matches_vec <- match(CRISPR_df[, ID_column], IDs_vec[are_this_sg])
    spec_vec[are_this_sg][matches_vec]
  }))
  colnames(spec_mat) <- paste0("Specificity_sg", unique_sg_numbers)
  results_df <- data.frame(CRISPR_df[, use_columns],
                           spec_mat,
                           "Min_specificity" = apply(spec_mat, 1, min),
                           "Combined_specificity" = 1 / rowSums(1 / spec_mat)
                           )
  return(results_df)
}

