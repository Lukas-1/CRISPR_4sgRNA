### 7th August 2019 ####





# Import packages and source code -----------------------------------------

library("SNPlocs.Hsapiens.dbSNP151.GRCh38")
library("MafDb.1Kgenomes.phase1.GRCh38")
library("MafDb.1Kgenomes.phase3.GRCh38")
library("MafDb.gnomAD.r2.1.GRCh38")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R"))







# Define functions --------------------------------------------------------

cp <- function(p) {
  ### Adapted from https://stackoverflow.com/a/35119842
  p <- p[!(is.na(p))]
  if (length(p) == 0) {
    return(NA_real_)
  } else if (length(p) == 1) {
    return(p)
  } else {
    ev <- do.call(expand.grid, replicate(length(p), 0:1, simplify = FALSE))
    pe <- apply(ev, 1, function(x) prod(p * (x == 1) + (1 - p) * (x == 0)))
    cp_vec <- tapply(pe, rowSums(ev), sum)
    return(sum(cp_vec[-1]))
  }
}


NoNAmax <- function(numeric_vec) {
  if (all(is.na(numeric_vec))) {
    NA_real_
  } else {
    max(numeric_vec, na.rm = TRUE)
  }
}



FrequenciesFromMafDb <- function(ranges_df, SNP_package_name = "MafDb.1Kgenomes.phase1.GRCh38", frequency_cutoff = 0.001, round_frequencies = TRUE, columns_postfix = NULL) {

  CheckRangesDf(ranges_df)

  ### Construct a GRanges object from the positions_df data frame
  GRanges_object_sgRNAs <- GRanges(
    seqnames = sub("chr", "", ranges_df[, "Chromosome"], fixed = TRUE),
    ranges   = IRanges(start = ranges_df[, "Start"], end = ranges_df[, "End"]),
    strand   = ranges_df[, "Strand"]
  )

  ### Find overlaps with SNPs in the human genome ###
  message("Searching for SNPs within the indicated genomic ranges...")
  gc()
  GRanges_object_overlapping_SNPs <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP151.GRCh38, GRanges_object_sgRNAs)


  ### Extract data on SNP allele frequencies from the thousand genomes project ###
  message("Filtering for high-frequency SNPs...")
  GPos_object_frequencies_added <- suppressWarnings(gscores(get(SNP_package_name), GRanges_object_overlapping_SNPs))
  frequencies_added_df <- as.data.frame(mcols(GPos_object_frequencies_added), stringsAsFactors = FALSE)


  ### Filter out low-frequency minor alleles ###
  is_high_frequency <- ifelse(is.na(frequencies_added_df[, "AF"]), FALSE, frequencies_added_df[, "AF"] >= frequency_cutoff)
  high_frequencies_df <- frequencies_added_df[is_high_frequency, ]
  GPos_object_polymorphisms <- GPos_object_frequencies_added[is_high_frequency, ]


  ### Compile data on SNPs that overlap with sgRNAs ###
  message("Counting high-frequency SNPs within the regions...")
  SNP_overlap_matches_df <- as.data.frame(findOverlaps(GRanges_object_sgRNAs, GPos_object_polymorphisms))
  SNP_matches_list <- split(SNP_overlap_matches_df[, 2], factor(SNP_overlap_matches_df[, 1], levels = unique(SNP_overlap_matches_df[, 1])))

  match_matches_vec <- match(as.character(seq_len(nrow(ranges_df))), names(SNP_matches_list))
  num_SNPs_vec  <- lengths(SNP_matches_list)[match_matches_vec]
  num_SNPs_vec  <- ifelse(is.na(num_SNPs_vec), 0L, num_SNPs_vec)

  SNP_rsIDs_vec <- vapply(SNP_matches_list, function(x) paste0(high_frequencies_df[x, "RefSNP_id"], collapse = ", "), "")[match_matches_vec]
  SNP_AFs_vec <- vapply(SNP_matches_list, function(x) paste0(format(high_frequencies_df[x, "AF"], scientific = FALSE), collapse = ", "), "")[match_matches_vec]

  message("Calculating the probabilities of at least one high-frequency SNP being present within the regions...")
  max_AF_vec <- vapply(SNP_matches_list, function(x) NoNAmax(high_frequencies_df[x, "AF"]), numeric(1))[match_matches_vec]
  sum_AFs_vec <- vapply(SNP_matches_list, function(x) cp(high_frequencies_df[x, "AF"]), numeric(1))[match_matches_vec]

  if (round_frequencies) {
    sum_AFs_vec <- ifelse(sum_AFs_vec >= 0.1, signif(sum_AFs_vec, digits = 2), signif(sum_AFs_vec, digits = 1))
  }

  ### Compile results ###
  results_df <- data.frame(
    "num_SNPs"       = num_SNPs_vec,
    "SNP_IDs"        = SNP_rsIDs_vec,
    "SNP_AFs"        = SNP_AFs_vec,
    "SNP_AF_max"     = max_AF_vec,
    "SNP_AF_sum"     = sum_AFs_vec,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  if (!(is.null(columns_postfix))) {
    names(results_df) <- paste0(names(results_df), "_", columns_postfix)
  }
  return(results_df)
}




FrequenciesFromVCFFiles <- function(ranges_df, frequency_cutoff = 0.001, round_frequencies = TRUE, only_Kaviar = FALSE) {

  stopifnot("common_polymorphisms_df" %in% ls(envir = globalenv())) # requires common_polymorphisms_df in the global environment

  CheckRangesDf(ranges_df)

  GRanges_object_sgRNAs <- GRanges(
    seqnames = sub("chr", "", ranges_df[, "Chromosome"], fixed = TRUE),
    ranges   = IRanges(start = ranges_df[, "Start"], end = ranges_df[, "End"]),
    strand   = ranges_df[, "Strand"]
  )

  have_rsID_df <- common_polymorphisms_df[!(is.na(common_polymorphisms_df[, "rsID"])), ]

  exceed_cutoff <- ((1 - have_rsID_df[, "AF_Kaviar"]) >= frequency_cutoff) %in% TRUE
  if (!(only_Kaviar)) {
    exceed_cutoff <- exceed_cutoff |
                     (((1 - have_rsID_df[, "AF_1kGenomes"]) >= frequency_cutoff) %in% TRUE) |
                     (((1 - have_rsID_df[, "AF_TOPMED"]) >= frequency_cutoff) %in% TRUE)
  }

  GRanges_object_polymorphisms <- GRanges(
    seqnames = have_rsID_df[exceed_cutoff, "Chromosome"],
    ranges   = IRanges(start = have_rsID_df[exceed_cutoff, "Position"],
                       end = have_rsID_df[exceed_cutoff, "Position"] + nchar(have_rsID_df[exceed_cutoff, "Reference"]) - 1L
                       ),
    strand   = "*"
  )

  rsIDs_source_vec <- have_rsID_df[exceed_cutoff, "rsID"]
  AF_Kaviar_source_vec <- 1 - have_rsID_df[exceed_cutoff, "AF_Kaviar"]

  if (!(only_Kaviar)) {
    AF_TOPMED_source_vec <- 1 - have_rsID_df[exceed_cutoff, "AF_TOPMED"]
    AF_1k_source_vec     <- 1 - have_rsID_df[exceed_cutoff, "AF_1kGenomes"]
  }

  ### Compile data on SNPs that overlap with sgRNAs ###
  message("Counting high-frequency SNPs within the regions...")
  SNP_overlap_matches_df <- as.data.frame(findOverlaps(GRanges_object_sgRNAs, GRanges_object_polymorphisms))
  SNP_matches_list <- split(SNP_overlap_matches_df[, 2], factor(SNP_overlap_matches_df[, 1], levels = unique(SNP_overlap_matches_df[, 1])))

  match_matches_vec <- match(as.character(seq_len(nrow(ranges_df))), names(SNP_matches_list))
  num_SNPs_vec <- lengths(SNP_matches_list)[match_matches_vec]
  num_SNPs_vec <- ifelse(is.na(num_SNPs_vec), 0L, num_SNPs_vec)
  SNP_rsIDs_vec <- vapply(SNP_matches_list, function(x) paste0(rsIDs_source_vec[x], collapse = ", "), "")[match_matches_vec]

  message("Finding the highest-frequency SNP within the regions...")
  max_AF_Kaviar_vec <- vapply(SNP_matches_list, function(x) NoNAmax(AF_Kaviar_source_vec[x]), numeric(1))[match_matches_vec]
  if (!(only_Kaviar)) {
     max_AF_1kGenomes_vec <- vapply(SNP_matches_list, function(x) NoNAmax(AF_1k_source_vec[x]),     numeric(1))[match_matches_vec]
     max_AF_TOPMED_vec    <- vapply(SNP_matches_list, function(x) NoNAmax(AF_TOPMED_source_vec[x]), numeric(1))[match_matches_vec]
  }

  message("Calculating the probabilities of at least one high-frequency SNP being present within the regions...")
  sum_AFs_Kaviar_vec <- vapply(SNP_matches_list, function(x) cp(AF_Kaviar_source_vec[x]), numeric(1))[match_matches_vec]
  if (!(only_Kaviar)) {
     sum_AFs_1kGenomes_vec <- vapply(SNP_matches_list, function(x) cp(AF_1k_source_vec[x]),     numeric(1))[match_matches_vec]
     sum_AFs_TOPMED_vec    <- vapply(SNP_matches_list, function(x) cp(AF_TOPMED_source_vec[x]), numeric(1))[match_matches_vec]
  }

  if (round_frequencies) {
    round_columns <- c("AF_Kaviar_source_vec", "sum_AFs_Kaviar_vec")
    if (!(only_Kaviar)) {
      round_columns <- c(round_columns, c("AF_1k_source_vec", "AF_TOPMED_source_vec", "sum_AFs_1kGenomes_vec", "sum_AFs_TOPMED_vec"))
    }
    for (my_vec_name in round_columns) {
      new_vec <- ifelse(get(my_vec_name) >= 0.1, signif(get(my_vec_name), digits = 2), signif(get(my_vec_name), digits = 1))
      assign(my_vec_name, new_vec)
    }
  }

  message("Compiling results...")

  CompileAFs <- function(AF_source_vec) {
    SNP_AFs_vec <- vapply(SNP_matches_list, function(x) {
      if (all(is.na(AF_source_vec[x]))) {
        NA_character_
      } else {
        paste0(format(AF_source_vec[x][!is.na(AF_source_vec[x])], scientific = FALSE), collapse = ", ")
      }
      }, "")[match_matches_vec]
    return(SNP_AFs_vec)
  }

  SNP_AFs_Kaviar_vec <- CompileAFs(AF_Kaviar_source_vec)

  if (only_Kaviar) {
    results_df <- data.frame(
      "num_SNPs_vcf"         = num_SNPs_vec,
      "SNP_IDs_vcf"          = SNP_rsIDs_vec,
      "SNP_AFs_Kaviar"       = SNP_AFs_Kaviar_vec,
      "SNP_AF_max_Kaviar"    = max_AF_Kaviar_vec,
      "SNP_AF_sum_Kaviar"    = sum_AFs_Kaviar_vec,
      stringsAsFactors       = FALSE,
      row.names              = NULL
    )
  } else {
    SNP_AFs_1k_vec     <- CompileAFs(AF_1k_source_vec)
    SNP_AFs_TOPMED_vec <- CompileAFs(AF_TOPMED_source_vec)
    results_df <- data.frame(
      "num_SNPs_vcf"         = num_SNPs_vec,
      "SNP_IDs_vcf"          = SNP_rsIDs_vec,
      "SNP_AFs_1kGenomes"    = SNP_AFs_1k_vec,
      "SNP_AF_max_1kGenomes" = max_AF_1kGenomes_vec,
      "SNP_AF_sum_1kGenomes" = sum_AFs_1kGenomes_vec,
      "SNP_AFs_TOPMED"       = SNP_AFs_TOPMED_vec,
      "SNP_AF_max_TOPMED"    = max_AF_TOPMED_vec,
      "SNP_AF_sum_TOPMED"    = sum_AFs_TOPMED_vec,
      "SNP_AFs_Kaviar"       = SNP_AFs_Kaviar_vec,
      "SNP_AF_max_Kaviar"    = max_AF_Kaviar_vec,
      "SNP_AF_sum_Kaviar"    = sum_AFs_Kaviar_vec,
      stringsAsFactors       = FALSE,
      row.names              = NULL
    )
  }

  return(results_df)
}



AllPolymorphisms <- function(ranges_df, only_23bp_only_Kaviar = FALSE) {

  sgRNA_only_ranges_df <- ranges_df[, c("Chromosome", "Strand", "Start", "End")]
  sgRNA_plus_PAM_ranges_df <- sgRNA_only_ranges_df

  if (!(only_23bp_only_Kaviar)) {
    PAM_ranges_df <- sgRNA_only_ranges_df
    PAM_ranges_df[, "Start"] <- ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] + 2L, ranges_df[, "Start"] - 3L)
    PAM_ranges_df[, "End"]   <- ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] + 3L, ranges_df[, "Start"] - 2L)
  }

  sgRNA_plus_PAM_ranges_df[, "Start"] <- ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "Start"], ranges_df[, "Start"] - 3L)
  sgRNA_plus_PAM_ranges_df[, "End"]   <- ifelse(ranges_df[, "Strand"] == "+", ranges_df[, "End"] + 3L, ranges_df[, "End"])

  message("Finding overlaps between the full 23 nucleotide sequence (sgRNA + PAM) and polymorphisms in the human genome...")
  sgRNA_plus_PAM_polymorphisms_df <- AllPolymorphismsForRange(sgRNA_plus_PAM_ranges_df, columns_prefix = "all23", only_Kaviar = only_23bp_only_Kaviar)
  message("")

  if (!(only_23bp_only_Kaviar)) {

    message("Finding overlaps between the 20-nucleotide sgRNA sequences and polymorphisms in the human genome...")
    sgRNA_only_polymorphisms_df <- AllPolymorphismsForRange(sgRNA_only_ranges_df, columns_prefix = "sgRNA")
    message("")
    message("Finding overlaps between the 2-nucleotide (N)GG PAM sequences and polymorphisms in the human genome...")
    PAM_polymorphisms_df <- AllPolymorphismsForRange(PAM_ranges_df, columns_prefix = "PAM")
    message("")

    results_df <- data.frame(
      sgRNA_only_polymorphisms_df,
      PAM_polymorphisms_df,
      sgRNA_plus_PAM_polymorphisms_df,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    SNP_column_names <- grep("_SNP_", names(results_df), fixed = TRUE, value = TRUE)
    SNP_column_roots <- unique(sub("^(PAM|sgRNA|all23)_", "", SNP_column_names))

    SNP_AF_max_roots <- grep("_AF_max_", SNP_column_roots, fixed = TRUE, value = TRUE)
    NGG_AF_max_list <- sapply(SNP_AF_max_roots, function(root) {
      PAM_vec <- results_df[, paste0("PAM_", root)]
      sg_vec <- results_df[, paste0("sgRNA_", root)]
      results_vec <- vapply(seq_len(nrow(results_df)), function(x) {
        NoNAmax(c(PAM_vec[[x]], sg_vec[[x]]))
      }, numeric(1))
      return(results_vec)
    }, simplify = FALSE)
    NGG_AF_max_mat <- do.call(cbind, NGG_AF_max_list)

    SNP_rsID_roots <- grep("SNP_IDs_", SNP_column_roots, fixed = TRUE, value = TRUE)
    NGG_rsID_list <- sapply(SNP_rsID_roots, function(root) {
      PAM_list   <- strsplit(results_df[, paste0("PAM_", root)], ", ", fixed = TRUE)
      sg_list    <- strsplit(results_df[, paste0("sgRNA_", root)], ", ", fixed = TRUE)
      all23_list <- strsplit(results_df[, paste0("all23_", root)], ", ", fixed = TRUE)
      results_vec <- vapply(seq_len(nrow(results_df)), function(x) {
        rsIDs_vec <- intersect(all23_list[[x]], c(sg_list[[x]], PAM_list[[x]]))
        if (all(is.na(rsIDs_vec))) {
          return(NA_character_)
        } else {
          return(paste0(rsIDs_vec, collapse = ", "))
        }
      }, "")
      return(results_vec)
    }, simplify = FALSE)
    NGG_rsID_mat <- do.call(cbind, NGG_rsID_list)

    NGG_df <- data.frame(NGG_rsID_mat, NGG_AF_max_mat, stringsAsFactors = FALSE, row.names = NULL)
    NGG_df <- NGG_df[, order(match(names(NGG_df), SNP_column_roots))]
    names(NGG_df) <- paste0("all22_", names(NGG_df))
    results_df <- data.frame(results_df, NGG_df, stringsAsFactors = FALSE, row.names = NULL)

  } else {
    results_df <- sgRNA_plus_PAM_polymorphisms_df
  }

  return(results_df)
}



AllPolymorphismsForRange <- function(ranges_df, columns_prefix = NULL, only_Kaviar = FALSE) {

  freq_vcf <- FrequenciesFromVCFFiles(ranges_df, only_Kaviar = only_Kaviar)

  if (!(only_Kaviar)) {
    freq_1kGp1  <- FrequenciesFromMafDb(ranges_df, SNP_package_name = "MafDb.1Kgenomes.phase1.GRCh38", columns_postfix = "1kG_ph1")
    freq_1kGp3  <- FrequenciesFromMafDb(ranges_df, SNP_package_name = "MafDb.1Kgenomes.phase3.GRCh38", columns_postfix = "1kG_ph3")
    freq_gnomAD <- FrequenciesFromMafDb(ranges_df, SNP_package_name = "MafDb.gnomAD.r2.1.GRCh38",      columns_postfix = "gnomAD")

    sgRNA_polymorphisms_df <- data.frame(
      freq_vcf,
      freq_1kGp1,
      freq_1kGp3,
      freq_gnomAD,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    sgRNA_polymorphisms_df <- freq_vcf
  }

  if (!(is.null(columns_prefix))) {
    names(sgRNA_polymorphisms_df) <- paste0(columns_prefix, "_", names(sgRNA_polymorphisms_df))
  }
  return(sgRNA_polymorphisms_df)
}





