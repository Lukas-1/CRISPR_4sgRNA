### 27th September 2020 ###



# Import packages and source code -----------------------------------------

library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("data.table") # for data.table::frank




# Define folder paths -----------------------------------------------------

CRISPR_root_directory         <- "~/CRISPR"
general_RData_directory       <- file.path(CRISPR_root_directory, "3) RData files", "1) General")
human_genome_directory        <- file.path(CRISPR_root_directory, "2) Input data", "Human genome")
FANTOM5_input_directory       <- file.path(human_genome_directory, "FANTOM5_liftover")
Ensembl_input_directory       <- file.path(human_genome_directory, "Ensembl")

sub_project_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-27 - target HNRNPK")
intermediate_files_directory  <- file.path(sub_project_directory, "2) Intermediate files")
file_output_directory         <- file.path(sub_project_directory, "3) Output")






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
BioMart_df <- read.table(file.path(human_genome_directory, "Ensembl", "BioMart_human_2020-03-25_mart_export.txt"),
                         sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL, check.names = FALSE
                         )





# Define maps -------------------------------------------------------------

entrez_to_symbol_vec <- as.character(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])




# Define functions --------------------------------------------------------

MakeTxDbDf <- function(GRangesList_object, filter_entrezs = NULL) {

  TxDb_columns <- c(
    # "group"      = "Group",
    "group_name" = "Entrez_ID",
    "seqnames"   = "Chromosome",
    "strand"     = "Strand",
    "start"      = "Start",
    "end"        = "End",
    "width"      = "Width",
    "tx_name"    = "ENST_ID",
    "tx_id"      = "Tx_ID",
    "exon_id"    = "Exon_ID",
    "cds_id"     = "CDS_ID"
  )

  results_df <- as.data.frame(GRangesList_object, stringsAsFactors = FALSE)
  chosen_names <- intersect(names(TxDb_columns), names(results_df))
  results_df <- results_df[, chosen_names]
  names(results_df) <- TxDb_columns[names(results_df)]
  results_df[["Strand"]] <- as.character(results_df[["Strand"]])
  results_df[["Chromosome"]] <- as.character(results_df[["Chromosome"]])
  results_df <- results_df[order(as.integer(results_df[["Entrez_ID"]])), ]
  if (!(is.null(filter_entrezs))) {
    results_df <- results_df[results_df[["Entrez_ID"]] %in% filter_entrezs, ]
  }
  row.names(results_df) <- NULL
  return(results_df)
}


FilterHNRNPK <- function(input_df) {
  input_df[input_df[["Entrez_ID"]] %in% "3190", ]
}


SplitOffTargetsSummary <- function(off_targets_summary_vec) {
  offtargets_summary_splits <- strsplit(off_targets_summary_vec, "|", fixed = TRUE)
  results_df <- data.frame(
    "GuideScan_Num_2MM" = vapply(offtargets_summary_splits, function(x) if (all(is.na(x))) NA_integer_ else as.integer(sub("2:", "", x[[1]], fixed = TRUE)), integer(1)),
    "GuideScan_Num_3MM" = vapply(offtargets_summary_splits, function(x) if (all(is.na(x))) NA_integer_ else as.integer(sub("3:", "", x[[2]], fixed = TRUE)), integer(1)),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  return(results_df)
}




TidyGuideScanColumns <- function(guidescan_df) {
  results_df <- data.frame(
    "sgRNA_sequence"        = guidescan_df[["gRNA"]],
    "Chromosome"            = guidescan_df[["chromosome"]],
    "Strand"                = guidescan_df[["strand"]],
    "Start"                 = as.integer(guidescan_df[["target site start coordinate"]]),
    "End"                   = as.integer(guidescan_df[["target site end coordinate"]]),
    "GuideScan_efficiency"  = as.numeric(ifelse(guidescan_df[["cutting efficiency score"]] == "*",
                                                NA_character_,
                                                guidescan_df[["cutting efficiency score"]]
                                                )
                                         ),
    "GuideScan_specificity" = as.numeric(guidescan_df[["cutting specificity score"]]),
    SplitOffTargetsSummary(guidescan_df[["offtargets summary"]]),
    "GuideScan_Num_2or3MM"  = as.integer(guidescan_df[["offtargets sum"]]),
    "gRNA_label"            = guidescan_df[["gRNA label"]],
    stringsAsFactors        = FALSE,
    row.names               = NULL
  )
  # Make the GuideScan locations consistent with the locations returned by Biostrings::matchPattern
  results_df[["Start"]] <- ifelse(results_df[["Strand"]] %in% "-",
                                  results_df[["Start"]] + 3L,
                                  results_df[["Start"]]
                                  )
  results_df[["End"]] <- ifelse(results_df[["Strand"]] %in% "-",
                                results_df[["End"]] + 1L,
                                results_df[["End"]] - 2L
                                )

  results_df[["Cut_location"]] <- GetCutLocations(results_df)
  results_df[["sgRNA_string"]] <- MakesgRNAString(results_df)
  return(results_df)
}



GetCutLocations <- function(ranges_df) {
  ifelse(ranges_df[["Strand"]] == "+", ranges_df[["End"]] - 2L, ranges_df[["Start"]] + 3L)
}




MakesgRNAString <- function(use_df, chromosome = NULL) {
  if ("Chromosome" %in% names(use_df)) {
    chromosome <- use_df[["Chromosome"]]
  } else if (is.null(chromosome)) {
    stop("No chromosome provided!")
  }
  paste0(chromosome, ":",
         use_df[["Strand"]], ":",
         use_df[["Cut_location"]], "_",
         use_df[["sgRNA_sequence"]]
         )
}






ResolveMissingOffTargets <- function(CRISPR_df, use_for_zero = 0) {

  CFD_columns <- c("CRISPOR_4MM_specificity", "CRISPOR_3MM_specificity")

  lack_detailed_offtargets <- is.na(CRISPR_df[["CRISPOR_4MM_specificity"]]) &
                              !(is.na(CRISPR_df[["CRISPOR_CFD_specificity"]]))

  have_no_offtargets <- lack_detailed_offtargets & (CRISPR_df[["CRISPOR_CFD_specificity"]] %in% 100)
  if (any(have_no_offtargets)) {
    for (Num_MM_column in c(paste0("CRISPOR_Num_", 0:4, "MM"), "CRISPOR_Num_2or3MM")) {
      CRISPR_df[[Num_MM_column]][have_no_offtargets] <- 0L
    }
  }
  for (CFD_column in CFD_columns) {
    CRISPR_df[[CFD_column]][have_no_offtargets] <- 1
  }
  too_many_offtargets <- lack_detailed_offtargets & (CRISPR_df[["CRISPOR_CFD_specificity"]] %in% 0)
  for (CFD_column in CFD_columns) {
    CRISPR_df[[CFD_column]][too_many_offtargets] <- use_for_zero
  }
  return(CRISPR_df)
}



TidyCRISPORDf <- function(CRISPOR_df) {

  CRISPOR_columns <- c(
    "guideId"               = "gRNA_ID",
    "targetSeq"             = "Full_sequence",
    "mitSpecScore"          = "CRISPOR_MIT_specificity",
    "cfdSpecScore"          = "CRISPOR_CFD_specificity",
    "offtargetCount"        = "CRISPOR_off_target_count",
    "Doench '16-Score"      = "CRISPOR_Doench_efficacy",
    "Moreno-Mateos-Score"   = "CRISPOR_Moreno_Mateos",
    "Out-of-Frame-Score"    = "CRISPOR_out_of_frame",
    "Lindel-Score"          = "CRISPOR_lindel_score",
    "GrafEtAlStatus"        = "CRISPOR_Graf_status",
    "targetGenomeGeneLocus" = "CRISPOR_note"
  )

  chosen_names <- intersect(names(CRISPOR_columns), names(CRISPOR_df))
  results_df <- CRISPOR_df[, chosen_names]
  names(results_df) <- CRISPOR_columns[names(results_df)]

  stopifnot(typeof(results_df[["CRISPOR_MIT_specificity"]]) == "integer")
  stopifnot(typeof(results_df[["CRISPOR_CFD_specificity"]]) == "integer")

  results_df[["sgRNA_sequence"]] <- substr(results_df[["Full_sequence"]], 1, 20)
  results_df[["PAM"]] <- substr(results_df[["Full_sequence"]], 21, 23)

  # The following code may have to be changed if a sequence is submitted to
  # CRISPOR that does not lie on the minus strand.
  seq_ID <- unique(CRISPOR_df[["#seqId"]])
  stopifnot(length(seq_ID) == 1)
  seq_ID_splits <- strsplit(seq_ID, ":", fixed = TRUE)[[1]]
  seq_ID_range <- as.integer(strsplit(seq_ID_splits[[2]], "-", fixed = TRUE)[[1]])
  stopifnot(seq_ID_splits[[3]] == "-")
  results_df[["Strand"]] <- ifelse(grepl("forw", results_df[["gRNA_ID"]], fixed = TRUE), "-", "+")
  locations_vec <- as.integer(sub("(forw|rev)$", "", results_df[["gRNA_ID"]]))
  results_df[["Cut_location"]] <- ifelse(results_df[["Strand"]] == "-",
                                         seq_ID_range[[2]] - locations_vec + 5L,
                                         seq_ID_range[[2]] - locations_vec - 4L
                                         )
  results_df[["sgRNA_string"]] <- MakesgRNAString(results_df, seq_ID_splits[[1]])
  results_df <- results_df[, !(names(results_df) == "Full_sequence")]
  return(results_df)
}




SummarizeOfftargets <- function(offtargets_df) {
  offtargets_df[["mismatchCount"]] <- factor(offtargets_df[["mismatchCount"]],
                                             levels = 0:4, ordered = TRUE
                                             )
  offtargets_df[["cfdOfftargetScore"]] <- as.numeric(ifelse(offtargets_df[["cfdOfftargetScore"]] == "None",
                                                            NA,
                                                            offtargets_df[["cfdOfftargetScore"]]
                                                            )
                                                     )
  results_list <- tapply(seq_len(nrow(offtargets_df)),
                         factor(offtargets_df[["guideId"]],
                                levels = unique(offtargets_df[["guideId"]])
                                ),
                         function(x) {
                           mismatch_table <- as.integer(table(offtargets_df[["mismatchCount"]][x]))
                           names(mismatch_table) <- paste0("CRISPOR_Num_", 0:4, "MM")
                           specificity_unrounded <- 1 / (1 + sum(offtargets_df[["cfdOfftargetScore"]][x], na.rm = TRUE))
                           specificity_upto3MM <- 1 / (1 + sum(offtargets_df[["cfdOfftargetScore"]][x][offtargets_df[["mismatchCount"]][x] %in% 0:3], na.rm = TRUE))
                           return(c(
                             list("gRNA_ID" = offtargets_df[["guideId"]][[x[[1]]]]),
                             as.list(mismatch_table),
                             list(
                               "CRISPOR_3MM_specificity" = specificity_upto3MM,
                               "CRISPOR_4MM_specificity" = specificity_unrounded
                             )
                           ))
                         })

  results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  lack_detailed_offtargets <- is.na(results_df[["CRISPOR_4MM_specificity"]]) |
                              !(is.na(results_df[["CRISPOR_CFD_specificity"]]))
  stopifnot(!(any(lack_detailed_offtargets)))

  results_df[["CRISPOR_Num_2or3MM"]] <- as.integer(rowSums(as.matrix(results_df[, c("CRISPOR_Num_2MM", "CRISPOR_Num_3MM")])))
  return(results_df)
}



TidyGPPDf <- function(GPP_df) {
  GPP_columns <- c(
    "Pick Order"                   = "GPP_rank",
    "sgRNA Cut Position (1-based)" = "Cut_location",
    "Strand of sgRNA"              = "Strand",
    "sgRNA Sequence"               = "sgRNA_sequence",
    "PAM Sequence"                 = "PAM",
    "Other Target Matches"         = "GPP_note",
    "On-Target Efficacy Score"     = "Azimuth_2_efficiency_score"
  )
  chosen_names <- intersect(names(GPP_columns), names(GPP_df))
  results_df <- GPP_df[, chosen_names]
  names(results_df) <- GPP_columns[names(results_df)]
  results_df[["sgRNA_string"]] <- MakesgRNAString(results_df, "chr9")
  return(results_df)
}




CutLocationsToRanges <- function(sgRNA_df) {
  data.frame(
    sgRNA_df[, c("Chromosome", "Strand")],
    "Start" = sgRNA_df[["Cut_location"]] - 1L,
    "End"   = sgRNA_df[["Cut_location"]],
    stringsAsFactors = FALSE
  )
}


RangesDfToGRangesObject <- function(ranges_df) {
  GRanges_object <- GRanges(
    seqnames = ranges_df[["Chromosome"]],
    ranges   = IRanges(start = ranges_df[["Start"]], end = ranges_df[["End"]]),
    strand   = ranges_df[["Strand"]]
  )
  return(GRanges_object)
}





FilterOverlapping <- function(query_ranges_df, target_ranges_df) {

  query_GRanges <- RangesDfToGRangesObject(query_ranges_df)
  target_GRanges <- RangesDfToGRangesObject(target_ranges_df)
  hits_object <- findOverlaps(query_GRanges,
                              target_GRanges,
                              ignore.strand = TRUE,
                              select = "all"
                              )
  found_indices <- subjectHits(hits_object)
  return(target_ranges_df[found_indices, ])
}



ReturnOverlappingIDs <- function(query_ranges_df,
                                 target_ranges_df,
                                 entrez_ID_column = "Entrez_ID",
                                 return_symbols = TRUE
                                 ) {

  query_GRanges <- RangesDfToGRangesObject(query_ranges_df)
  target_GRanges <- RangesDfToGRangesObject(target_ranges_df)
  hits_object <- findOverlaps(query_GRanges,
                              target_GRanges,
                              ignore.strand = TRUE,
                              select = "all"
                              )

  query_fac <- factor(queryHits(hits_object), levels = seq_len(nrow(query_ranges_df)))

  target_IDs_vec <- target_ranges_df[[entrez_ID_column]]
  if (return_symbols) {
    target_IDs_vec <- entrez_to_symbol_vec[target_IDs_vec]
  }

  results_vec <- tapply(subjectHits(hits_object),
                        query_fac,
                        function(x) {
                          results_vec <- target_IDs_vec[x]
                          if (all(is.na(results_vec))) {
                            NA_character_
                          } else {
                            paste0(unique(results_vec[!(is.na(results_vec))]), collapse = ", ")
                          }
                        })
  return(results_vec)
}





AddOverlaps <- function(sgRNA_df) {
  guide_locations_df <- CutLocationsToRanges(sgRNA_df)
  results_df <- data.frame(
    sgRNA_df,
    "Gene_overlaps"       = ReturnOverlappingIDs(guide_locations_df, TxDb_genes_df),
    "Transcript_overlaps" = ReturnOverlappingIDs(guide_locations_df, TxDb_transcripts_df),
    "Exon_overlaps"       = ReturnOverlappingIDs(guide_locations_df, TxDb_exons_df),
    "CDS_overlaps"        = ReturnOverlappingIDs(guide_locations_df, TxDb_CDS_df),
    stringsAsFactors = FALSE
  )
  return(results_df)
}



SignedDistance <- function(query_df, target_df, return_ID_column = NULL) {

  query_GRanges_object <- RangesDfToGRangesObject(query_df)
  target_GRanges_object <- RangesDfToGRangesObject(target_df)

  hits_object <- distanceToNearest(query_GRanges_object, target_GRanges_object, ignore.strand = TRUE)

  stopifnot(identical(queryHits(hits_object), seq_len(nrow(query_df))))

  precede_vec <- precede(query_GRanges_object, target_GRanges_object, ignore.strand = TRUE)
  follow_vec <- follow(query_GRanges_object, target_GRanges_object, ignore.strand = TRUE)

  are_preceding <- mapply(identical, precede_vec, subjectHits(hits_object))

  target_strand <- unique(target_df[["Strand"]])
  stopifnot(identical(target_strand, "+") || identical(target_strand, "-"))
  if (target_strand == "-") {
    factor <- 1L
  } else {
    factor <- -1L
  }
  sign_vec <- ifelse(are_preceding, 1L, -1L) * factor

  distance_vec <- mcols(hits_object)[["distance"]] * sign_vec

  results_df <- data.frame("Distance" = distance_vec)

  if (!(is.null(return_ID_column))) {
    results_df <- data.frame(
      results_df,
      "Preceding_ID" = target_df[[return_ID_column]][precede_vec],
      "Following_ID" = target_df[[return_ID_column]][follow_vec],
      "Nearest_ID"   = target_df[[return_ID_column]][subjectHits(hits_object)],
      "Distance"     = distance_vec,
      stringsAsFactors = FALSE
    )
  }

  return(results_df)
}




GetDistances <- function(sgRNA_df, features_df, TSS_df) {

  sgRNA_locations_df <- CutLocationsToRanges(sgRNA_df)

  features_distance_df <- SignedDistance(sgRNA_locations_df, features_df, return_ID_column = "Feature_ID")

  if (identical(unique(TSS_df[["Strand"]]), "-")) {
    stopifnot(min(TSS_df[["Start"]]) > max(features_df[["End"]]))
  } else {
    stopifnot(max(TSS_df[["End"]]) < min(features_df[["Start"]]))
  }

  TSS_dist_mat <- vapply(seq_len(nrow(TSS_df)),
                         function(x) SignedDistance(sgRNA_locations_df, TSS_df[x, ])[["Distance"]],
                         integer(nrow(sgRNA_df))
                         )

  colnames(TSS_dist_mat) <- paste0("Distance_to_TSS", seq_len(ncol(TSS_dist_mat)))

  results_df <- data.frame(
    features_distance_df[, c("Nearest_ID", "Distance", "Preceding_ID", "Following_ID")],
    TSS_dist_mat,
    stringsAsFactors = FALSE
  )
  return(results_df)
}



TidyCDSDf <- function(CDS_df) {
  CDS_df <- CDS_df[order(CDS_df[["End"]], decreasing = TRUE), ]
  CDS_df[["Feature_ID"]] <- paste0("CDS", seq_len(nrow(CDS_df)))
  CDS_df <- CDS_df[, c("Feature_ID", "Chromosome", "Strand", "Start", "End")]
  return(CDS_df)
}



RankGuidesDf <- function(sgRNA_df, are_eligible = rep(TRUE, nrow(sgRNA_df))) {

  stopifnot(length(are_eligible) == nrow(sgRNA_df))
  are_specific <- sgRNA_df[["CRISPOR_3MM_specificity"]] > 0.2
  are_GrafOK <- sgRNA_df[["CRISPOR_Graf_status"]] %in% "GrafOK"
  have_polyT <- grepl("TTTT", sgRNA_df[["sgRNA_sequence"]], ignore.case = TRUE)

  are_eligible <- are_eligible & are_specific & are_GrafOK & !(have_polyT)

  results_columns <- c("Specificity_rank", "Efficacy_rank", "Combined_rank")
  ranks_mat <- matrix(nrow = nrow(sgRNA_df),
                      ncol = 3,
                      dimnames = list(NULL, results_columns)
                      )
  if (any(are_eligible)) {
    ranks_mat[are_eligible, "Specificity_rank"] <- rank(-(sgRNA_df[["CRISPOR_4MM_specificity"]][are_eligible]), ties.method = "first")
    ranks_mat[are_eligible, "Efficacy_rank"] <- rank(-(sgRNA_df[["Azimuth_2_efficiency_score"]][are_eligible]), ties.method = "first")
    combined_rank <- data.table::frank(
      list(rowMeans(ranks_mat[are_eligible, c("Specificity_rank", "Efficacy_rank")]),
           -(sgRNA_df[["CRISPOR_4MM_specificity"]][are_eligible]),
           -(sgRNA_df[["Azimuth_2_efficiency_score"]][are_eligible])
           ),
      ties.method = "first"
    )
    ranks_mat[are_eligible, "Combined_rank"] <- combined_rank
  }
  return(ranks_mat)
}



RetrieveSequence <- function(chromosome, strand, start, end) {
  library("BSgenome.Hsapiens.UCSC.hg38")
  result_string <- as.character(Views(BSgenome.Hsapiens.UCSC.hg38[[chromosome]], start, end))
  if (strand == "-") {
    result_string <- as.character(reverseComplement(DNAString(result_string)))
  } else if (strand != "+") {
    stop("Invalid value for the strand parameter!")
  }
  return(result_string)
}




# Create data frames using TxDb.Hsapiens.UCSC.hg38.knownGene --------------

TxDb_genes_df       <- MakeTxDbDf(genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                        single.strand.genes.only = FALSE
                                        ))
TxDb_transcripts_df <- MakeTxDbDf(transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene))
TxDb_exons_df       <- MakeTxDbDf(exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))
TxDb_CDS_df         <- MakeTxDbDf(cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))




# Build a combined FANTOM5 data frame -------------------------------------

names(FANTOM5_ann_df)[[1]] <- sub("#", "", names(FANTOM5_ann_df)[[1]], fixed = TRUE)

names(FANTOM5_bed_df) <- c(
  "Chromosome", "Peak_start", "Peak_end", "Peak_ID", "Score", "Strand",
  "TSS_start", "TSS_stop", "Color_code"
)

tidy_annotation_columns <- c(
  "Transcript_name", "HGNC/MGI_ID", "UniProt_ID", "Gene_name",
  "Gene_synonyms", "Gene_source"
)

FANTOM_ann_matches <- match(FANTOM5_bed_df[["Peak_ID"]], FANTOM5_ann_df[["CAGE_Peak_ID"]])

FANTOM5_df <- data.frame(
  FANTOM5_ann_df[FANTOM_ann_matches, c("CAGE_Peak_ID", "Gene_symbol")],
  FANTOM5_bed_df[, c("Chromosome", "Strand", "Peak_start", "Peak_end", "TSS_start", "TSS_stop")],
  FANTOM5_ann_df[FANTOM_ann_matches, "Distance", drop = FALSE],
  FANTOM5_bed_df["Score"],
  FANTOM5_ann_df[FANTOM_ann_matches, tidy_annotation_columns],
  check.names = FALSE,
  stringsAsFactors = FALSE
)




# Look up transcription start sites (TSSs) for HNRNPK ---------------------

HNRPNK_FANTOM5_df <- FANTOM5_df[FANTOM5_df[["Gene_symbol"]] %in% "HNRNPK", ]
HNRPNK_BioMart_df <- BioMart_df[BioMart_df[["NCBI gene ID"]] %in% "3190", ]
HNRPNK_BioMart_df[order(HNRPNK_BioMart_df[["Transcript start (bp)"]]), ]





# Look up the TSSs for the neighbouring RMI1 gene -------------------------

RMI1_FANTOM5_df <- FANTOM5_df[FANTOM5_df[["Gene_symbol"]] %in% "RMI1", ]
RMI1_BioMart_df <- BioMart_df[BioMart_df[["NCBI gene ID"]] %in% "80010", ]






# Define the search limits for gRNAs for HNRNPK ---------------------------

HNRNPK_TSS_indices <- order(HNRPNK_FANTOM5_df[["Score"]], decreasing = TRUE)[1:2]

HNRNPK_coordinates <- sort(
  c("gene_start"     = FilterHNRNPK(TxDb_genes_df)[["Start"]],
    "gene_end"       = FilterHNRNPK(TxDb_genes_df)[["End"]],
    "primary_TSS"    = HNRPNK_FANTOM5_df[["TSS_start"]][HNRNPK_TSS_indices[[1]]],
    "secondary_TSS"  = HNRPNK_FANTOM5_df[["TSS_start"]][HNRNPK_TSS_indices[[2]]]
    )
)


promoter_coordinates <- sort(
  c("ENSR00000236870_start" = 83974802,
    "ENSR00000236870_stop"  = 83987799,
    "ENSR00000884780_start" = 83964693,
    "ENSR00000884780_stop"  = 83966320
    )
)

sort(c(HNRNPK_coordinates, promoter_coordinates))


submit_df <- data.frame(
  "Chromosome" = FilterHNRNPK(TxDb_genes_df)[["Chromosome"]],
  "Start"      = FilterHNRNPK(TxDb_genes_df)[["Start"]] - 1000L,
  "End"        = FilterHNRNPK(TxDb_genes_df)[["End"]] + 1000L,
  "Names"      = ".",
  "Scores"     = ".",
  "Strand"     =  FilterHNRNPK(TxDb_genes_df)[["Strand"]],
  stringsAsFactors = FALSE
)




# Define a data frame of HNRNPK landmarks ---------------------------------

HNRNPK_TSS_df <- data.frame(
  "Feature_ID" = c("Primary TSS", "Secondary TSS"),
  HNRPNK_FANTOM5_df[HNRNPK_TSS_indices, c("Chromosome", "Strand")],
  "Start" = HNRPNK_FANTOM5_df[["TSS_start"]][HNRNPK_TSS_indices],
  "End" = HNRPNK_FANTOM5_df[["TSS_stop"]][HNRNPK_TSS_indices],
  stringsAsFactors = FALSE
)

HNRNPK_TSS_df <- HNRNPK_TSS_df[order(HNRNPK_TSS_df[["End"]], decreasing = TRUE), ]

HNRNPK_CDS_df <- TidyCDSDf(FilterHNRNPK(TxDb_CDS_df))





# Define a data frame of RMI1 hallmarks -----------------------------------

RMI1_TSS_df <- RMI1_FANTOM5_df[, c("Chromosome", "Strand", "TSS_start", "TSS_stop")]
names(RMI1_TSS_df)[names(RMI1_TSS_df) == "TSS_start"] <- "Start"
names(RMI1_TSS_df)[names(RMI1_TSS_df) == "TSS_stop"] <- "End"


RMI1_CDS_df <- TidyCDSDf(TxDb_CDS_df[TxDb_CDS_df[["Entrez_ID"]] %in% "80010", ])






# Filter for miRNA-7 exons ------------------------------------------------

miR7_df <- TxDb_exons_df[TxDb_exons_df[["Entrez_ID"]] %in% "407043", ]

args_list <- as.list(miR7_df[, 2:5])
names(args_list) <- tolower(names(args_list))
do.call(RetrieveSequence, args_list)







# Calculate the distance between HNRNPK and RM1 promoters -----------------

HNRNPK_TSS_df[["Start"]][HNRNPK_TSS_df[["Feature_ID"]] == "Primary TSS"] -
RMI1_TSS_df[["Start"]]

HNRNPK_TSS_df[["Start"]][HNRNPK_TSS_df[["Feature_ID"]] == "Secondary TSS"] -
RMI1_TSS_df[["Start"]]






# Export gene coordinates for external tools ------------------------------

write.table(paste0(submit_df[["Chromosome"]], ":", submit_df[["Start"]], "-", submit_df[["End"]]),
            file = file.path(intermediate_files_directory,
                             "GuideScan",
                             "HNRNPK_input_for_GuideScan.txt"
                             ),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )

write.table(submit_df,
            file = file.path(intermediate_files_directory,
                             "CRISPOR",
                             "HNRNPK_input_for_CRISPOR.bed"
                             ),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )

write.table(paste0("NC_000009.12:", submit_df[["Strand"]], ":", submit_df[["Start"]], "-", submit_df[["End"]]),
            file = file.path(intermediate_files_directory,
                             "GPP sgRNA designer",
                             "HNRNPK_input_for_GPP_sgRNA_designer.txt"
                             ),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )




# Read in data from external tools ----------------------------------------

ReadOutput <- function(sub_path, fill = FALSE) {
  read.table(file = file.path(intermediate_files_directory,
                              sub_path
                              ),
             sep = "\t", header = TRUE,
             row.names = NULL, quote = "", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE
             )
}

CRISPOR_guides_df     <- ReadOutput(file.path("CRISPOR", "CRISPOR_output__HRNPK__CRISPRko.tsv"))
CRISPOR_offtargets_df <- ReadOutput(file.path("CRISPOR", "CRISPOR_output__HRNPK__CRISPRko_offs.tsv"))

GPP_guides_df         <- ReadOutput(file.path("GPP sgRNA designer", "HNRNPK__GPP_output.txt"), fill = TRUE)

GuideScan_guides_df   <- read.csv(file = file.path(intermediate_files_directory,
                                                   "GuideScan",
                                                   "HNRNPK__GuideScan_output.csv"
                                                   ),
                                  skip = 1, header = TRUE, row.names = NULL,
                                  check.names = FALSE,
                                  quote = "\"", stringsAsFactors = FALSE,
                                  comment.char = ""
                                  )




# Tidy output from external tools -----------------------------------------

GuideScan_guides_df <- TidyGuideScanColumns(GuideScan_guides_df)

CRISPOR_offtargets_df <- SummarizeOfftargets(CRISPOR_offtargets_df)

CRISPOR_guides_df <- TidyCRISPORDf(CRISPOR_guides_df)

matches_vec <- match(CRISPOR_guides_df[["gRNA_ID"]],
                     CRISPOR_offtargets_df[["gRNA_ID"]]
                     )

CRISPOR_df <- data.frame(CRISPOR_guides_df,
                         CRISPOR_offtargets_df[matches_vec, names(CRISPOR_offtargets_df) != "gRNA_ID"],
                         stringsAsFactors = FALSE
                         )
CRISPOR_df <- ResolveMissingOffTargets(CRISPOR_df)

GPP_guides_df <- TidyGPPDf(GPP_guides_df)




# Combine data frames -----------------------------------------------------

CRISPOR_df <- CRISPOR_df[order(CRISPOR_df[["Cut_location"]], decreasing = TRUE), ]
row.names(CRISPOR_df) <- NULL

GPP_matches       <- match(CRISPOR_df[["sgRNA_string"]], GPP_guides_df[["sgRNA_string"]])
GuideScan_matches <- match(CRISPOR_df[["sgRNA_string"]], GuideScan_guides_df[["sgRNA_string"]])


GuideScan_columns <- c("GuideScan_efficiency", "GuideScan_specificity",
                       "GuideScan_Num_2MM", "GuideScan_Num_3MM"
                       )

guides_columns <- c(
  "Chromosome", "Strand", "Cut_location", "sgRNA_sequence", "PAM",
  "Azimuth_2_efficiency_score",
  "GuideScan_specificity", "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity",
  "CRISPOR_CFD_specificity",
  "CRISPOR_Num_0MM", "CRISPOR_Num_1MM",
  "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
  "GuideScan_Num_3MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
  "CRISPOR_Graf_status",
  "CRISPOR_note"
)

guides_df <- data.frame("Chromosome" = "chr9",
                        CRISPOR_df[, !(names(CRISPOR_df) %in% c("gRNA_ID", "sgRNA_string "))],
                        GPP_guides_df[GPP_matches, "Azimuth_2_efficiency_score", drop = FALSE],
                        GuideScan_guides_df[GuideScan_matches, GuideScan_columns],
                        stringsAsFactors = FALSE
                        )[, guides_columns]

guides_df <- AddOverlaps(guides_df)

HNRNPK_distances_df <- GetDistances(guides_df, HNRNPK_CDS_df, HNRNPK_TSS_df)
RMI1_distances_df <- GetDistances(guides_df, RMI1_CDS_df, RMI1_TSS_df)

HNRNPK_categories_vec <- vapply(seq_len(nrow(HNRNPK_distances_df)), function(x) {
  if (HNRNPK_distances_df[["Distance_to_TSS1"]][[x]] < 0) {
    "Upstream of the first (minor) TSS"
  } else if (HNRNPK_distances_df[["Distance_to_TSS2"]][[x]] < 0) {
    "Upstream of the main TSS"
  } else if (HNRNPK_distances_df[["Distance"]][[x]] == 0) {
    HNRNPK_distances_df[["Nearest_ID"]][[x]]
  } else if (is.na(HNRNPK_distances_df[["Preceding_ID"]][[x]])) {
    paste0("Upstream of ", HNRNPK_distances_df[["Nearest_ID"]][[x]])
  } else if (is.na(HNRNPK_distances_df[["Following_ID"]][[x]])) {
    paste0("Downstream of ", HNRNPK_distances_df[["Nearest_ID"]][[x]])
  } else {
    paste0("Between ", HNRNPK_distances_df[["Preceding_ID"]][[x]],
           " and ", HNRNPK_distances_df[["Following_ID"]][[x]]
           )
  }
}, "")



HNRNPK_supercategories_vec <- vapply(HNRNPK_categories_vec, function(x) {
  if (x  %in% c("Upstream of the first (minor) TSS", "Upstream of the main TSS")) {
    "Upstream of the main TSS"
  } else if (x == "Upstream of CDS1") {
    "5' UTR"
  } else if (x == "Downstream of CDS20") {
    "3' UTR"
  } else if (!(grepl("^Between CDS", x))) {
    "Coding sequence"
  } else {
    "Intron or UTR"
  }
}, "")


guides_df <- data.frame(
  guides_df,
  "Target_region"                   = HNRNPK_categories_vec,
  "Target_category"                 = HNRNPK_supercategories_vec,
  "Nearest_HNRNPK_CDS"              = HNRNPK_distances_df[["Nearest_ID"]],
  "Nearest_HNRNPK_CDS_distance"     = HNRNPK_distances_df[["Distance"]],
  "Distance_to_upstream_HNRNPK_TSS" = HNRNPK_distances_df[["Distance_to_TSS1"]],
  "Distance_to_main_HNRNPK_TSS"     = HNRNPK_distances_df[["Distance_to_TSS2"]],
  "Distance_to_RM1_TSS"             = RMI1_distances_df[["Distance_to_TSS1"]],
  stringsAsFactors = FALSE
)




# Identify prioritized guides ---------------------------------------------

target_CDS <- !(is.na(guides_df[["CDS_overlaps"]]))
target_other_genes <- (!(is.na(guides_df[["Transcript_overlaps"]]))) &
                      (guides_df[["Transcript_overlaps"]] != "HNRNPK")

are_distant <- (guides_df[["Target_region"]] == "Downstream of CDS20") &
               (is.na(guides_df[["Exon_overlaps"]]))

are_eligible <- !(target_CDS | target_other_genes | are_distant)

rank_mat <- matrix(nrow = nrow(guides_df), ncol = 3)
rank_mat[, ] <- -1L

for (sub_category in unique(guides_df[["Target_category"]])) {
  are_this_category <- guides_df[["Target_category"]] == sub_category
  this_mat <- RankGuidesDf(guides_df[are_this_category, ],
                           are_eligible = are_eligible[are_this_category]
                           )
  rank_mat[are_this_category, ] <- this_mat
  colnames(rank_mat) <- colnames(this_mat)
}

guides_df <- data.frame(guides_df, rank_mat)

guides_columns <- c(
  "sgRNA_sequence", "PAM", "Chromosome", "Strand", "Cut_location",
  "Target_region", "CRISPOR_note",
  "Target_category",
  "Combined_rank", "Specificity_rank", "Efficacy_rank",
  "Azimuth_2_efficiency_score", "GuideScan_specificity", "CRISPOR_3MM_specificity",
  "CRISPOR_4MM_specificity", "CRISPOR_CFD_specificity", "CRISPOR_Num_0MM",
  "CRISPOR_Num_1MM", "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
  "GuideScan_Num_3MM", "CRISPOR_Num_3MM", "CRISPOR_Num_4MM",
  "CRISPOR_Graf_status",
  "Nearest_HNRNPK_CDS", "Nearest_HNRNPK_CDS_distance",
  "Distance_to_upstream_HNRNPK_TSS", "Distance_to_main_HNRNPK_TSS", "Distance_to_RM1_TSS",
  "Gene_overlaps", "Transcript_overlaps", "Exon_overlaps", "CDS_overlaps"
)

guides_df <- guides_df[, guides_columns]





# Pick the top 4 guides ---------------------------------------------------

are_top25 <- ifelse(guides_df[["Target_category"]] == "Intron or UTR",
                    guides_df[["Combined_rank"]] %in% 1:10,
                    guides_df[["Combined_rank"]] %in% 1:5
                    )

selected_guides_df <- guides_df[are_top25, ][guides_df[["CRISPOR_3MM_specificity"]][are_top25] > 0.5, ]

new_order <- order(match(selected_guides_df[["Target_region"]], selected_guides_df[["Target_region"]]),
                   -(selected_guides_df[["CRISPOR_4MM_specificity"]])
                   )
selected_guides_df <- selected_guides_df[new_order, ]
# selected_guides_df <- selected_guides_df[!(duplicated(selected_guides_df[["Target_region"]])), ]

exclude_regions <- c("Upstream of the first (minor) TSS", "Downstream of CDS20")
selected_guides_df <- selected_guides_df[!(selected_guides_df[["Target_region"]] %in% exclude_regions), ]
row.names(selected_guides_df) <- NULL

select_indices <- c(2, 4, 8, 12)





# Check for additional highly specific guides -----------------------------

are_selected <- are_eligible &
                (guides_df[["Target_category"]] %in% "Intron or UTR") &
                (guides_df[["CRISPOR_3MM_specificity"]] > 0.55) &
                (guides_df[["Azimuth_2_efficiency_score"]] > 0.5) &
                (!(is.na(guides_df[["Combined_rank"]])))

new_order <- order(guides_df[["CRISPOR_4MM_specificity"]][are_selected],
                   decreasing = TRUE
                   )
downstream_guides_df <- guides_df[which(are_selected)[new_order], ]
row.names(downstream_guides_df) <- NULL

are_far_downstream <- !(downstream_guides_df[["Target_region"]] %in% selected_guides_df[["Target_region"]])





# Add two more guides to four already selected ----------------------------

selected_guides_df <- rbind.data.frame(selected_guides_df[select_indices, ],
                                       downstream_guides_df[which(are_far_downstream)[c(3, 1)], ],
                                       stringsAsFactors = FALSE
                                       )




# Export data -------------------------------------------------------------

ExportDf <- function(export_df, file_name) {
  have_NA_columns <- vapply(export_df, anyNA, logical(1))
  for (i in which(have_NA_columns)) {
    export_df[[i]] <- ifelse(is.na(export_df[[i]]),
                             "",
                             export_df[[i]]
                             )
  }
  write.table(export_df,
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )
}

ExportDf(guides_df[are_eligible, ], "HNRNPK_guides")
ExportDf(selected_guides_df, "HNRNPK_top6_nonCDS")



