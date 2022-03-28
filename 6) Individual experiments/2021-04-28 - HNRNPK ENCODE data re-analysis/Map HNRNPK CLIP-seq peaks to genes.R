### 1st May 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
# source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R")) # For CheckThatFactorIsInOrder
source(file.path(general_functions_directory, "29) Determining gene types.R"))
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")
file_directory          <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-28 - HNRNPK ENCODE data re-analysis")
file_input_directory    <- file.path(file_directory, "1) Input")
file_output_directory   <- file.path(file_directory, "2) Output")
RData_directory         <- file.path(file_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))





# Read in data ------------------------------------------------------------

ReadNarrowPeaks <- function(file_name) {
  peaks_df <- read.delim(file.path(file_input_directory, file_name),
                         header = FALSE, stringsAsFactors = FALSE
                         )
  # See https://genome.ucsc.edu/FAQ/FAQformat.html#format12
  narrowPeak_columns <- c(
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "signalValue", "pValue", "qValue", "peak"
  )
  names(peaks_df) <- narrowPeak_columns

  stopifnot(identical(unique(peaks_df[["qValue"]]), -1L))
  stopifnot(identical(unique(peaks_df[["peak"]]), -1L))

  peaks_df <- peaks_df[, !(names(peaks_df) %in% c("qValue", "peak"))]
  return(peaks_df)
}

HepG2_IDR_df <- ReadNarrowPeaks("206.01v02.IDR.out.0102merged.bed.blacklist_removed.bed.narrowPeak.encode.bed")
K562_IDR_df  <- ReadNarrowPeaks("326.01v02.IDR.out.0102merged.bed.blacklist_removed.bed.narrowPeak.encode.bed")





# Define functions --------------------------------------------------------

ReformatNarrowPeaks <- function(input_df) {
  results_df <- data.frame(
    "Index"        = seq_len(nrow(input_df)),
    "Chromosome"   = input_df[, "chrom"],
    "Strand"       = input_df[, "strand"],
    "Start"        = input_df[, "chromStart"],
    "End"          = input_df[, "chromEnd"],
    "Enrichment"   = input_df[, "signalValue"],
    "P_value"      = 10^(-(input_df[, "pValue"])),
    stringsAsFactors = FALSE
  )
  return(results_df)
}


AnnotatePeaks <- function(peaks_df, genes_df) {

  ## Prepare genes_df for use in the following functions
  rename_gene_columns <- c(
    "Gene_IDs" = "Affected_gene_IDs",
    "Source"   = "Gene_source"
  )
  for (column_name in names(rename_gene_columns)) {
    names(genes_df)[names(genes_df) == column_name] <- rename_gene_columns[[column_name]]
  }


  ## Prepare peaks_df for use in the following functions

  peaks_df[["Locus"]] <- MakeLocationStrings(peaks_df)
  peaks_df[["Guide_locus"]] <- peaks_df[["Locus"]]


  ## Calculate overlaps

  combined_df <- FindOverlappingHits(peaks_df,
                                     genes_df,
                                     retain_columns = c("Locus"),
                                     ignore_strand = FALSE
                                     )

  ## Expand the data frame
  location_columns <- c("Chromosome", "Strand", "Start", "End")
  full_combined_df <- AlignHits(peaks_df[, !(names(peaks_df) %in% location_columns)],
                                combined_df
                                )


  ## Prepare the data frame for summarization

  full_combined_df[["Intended_Entrez_ID"]] <- ""
  full_combined_df[["Intended_gene_symbol"]] <- ""
  full_combined_df[["Num_loci"]] <- 1L
  full_combined_df <- AddEntrezIDAvailable(full_combined_df, genes_df)


  ## Summarize all hits per peak

  hits_df <- SummarizeFullDf(full_combined_df, tolerate_num_affected = TRUE)


  ## Separate Ensembl gene IDs from Entrez gene IDs

  split_gene_IDs <- strsplit(hits_df[["Affected_gene_IDs"]], ", ", fixed = TRUE)

  entrez_IDs_vec <- vapply(split_gene_IDs, function(x) {
    are_ENSG <- grepl("^ENSG", x)
    results_vec <- x[!(are_ENSG)]
    stopifnot(all((results_vec == "") | !(is.na(as.integer(results_vec)))))
    paste0(results_vec, collapse = ", ")
  }, "")

  ENSG_IDs_vec <- vapply(split_gene_IDs, function(x) {
    are_ENSG <- grepl("^ENSG", x)
    if (!(any(are_ENSG))) {
      return(NA_character_)
    } else {
      paste0(x[are_ENSG], collapse = ", ")
    }
  }, "")



  ## Replace empty strings with NAs

  for (column_name in c("Affected_gene_IDs", "Affected_gene_symbols")) {
    hits_df
  }
  entrez_IDs_vec <- ifelse(entrez_IDs_vec == "", NA, entrez_IDs_vec)


  ## Assign gene types

  gene_types_vec <- GetGeneTypes(entrezs_vec = entrez_IDs_vec,
                                 ENSG_vec = ENSG_IDs_vec
                                 )


  ## Assemble the final data frame

  results_df <- data.frame(
    "Peak_number" = peaks_df[, "Index"],
    peaks_df[, c("Chromosome", "Strand", "Start", "End", "Enrichment", "P_value")],
    "Gene_IDs" = hits_df[, "Affected_gene_IDs"],
    "Entrez_gene_IDs" = entrez_IDs_vec,
    "Gene_symbols" = hits_df[, "Affected_gene_symbols"],
    "Number_of_genes" = hits_df[, "Num_affected_genes"],
    "Gene_types" = gene_types_vec,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



WriteTable <- function(export_df, file_path) {
  are_character <- vapply(export_df, is.character, logical(1))
  for (i in which(are_character)) {
    export_df[[i]] <- ifelse(is.na(export_df[[i]]), "", export_df[[i]])
  }
  write.table(export_df,
              file      = file_path,
              sep       = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote     = FALSE
              )
}




# Explore data ------------------------------------------------------------

table(HepG2_IDR_df[["name"]])
table(K562_IDR_df[["name"]])




# Prepare CLIP-seq data frames --------------------------------------------

HepG2_peaks_df <- ReformatNarrowPeaks(HepG2_IDR_df)
K562_peaks_df <- ReformatNarrowPeaks(K562_IDR_df)




# Prepare the transcript locations data frame -----------------------------

transcripts_df <- transcript_locations_df
transcripts_df[["Gene_types"]][transcripts_df[["Gene_types"]] %in% "protein-coding, ncRNA"] <- "protein-coding"
transcripts_df <- PrepareGenesDf(transcripts_df)


### Correct missing Entrez IDs

have_entrez <- !(is.na(transcripts_df[["Entrez_IDs"]]))

have_entrez_ENSG <- setdiff(transcripts_df[["Ensembl_gene_IDs"]][have_entrez], NA)

are_problematic <- !(have_entrez) & (transcripts_df[["Ensembl_gene_IDs"]] %in% have_entrez_ENSG)

table(transcripts_df[["Source"]][are_problematic])

problematic_gene_IDs <- transcripts_df[["Gene_IDs"]][are_problematic]

for (problematic_ID in problematic_gene_IDs) {
  are_this_ID <- transcripts_df[["Gene_IDs"]] %in% problematic_ID
  stopifnot(all(is.na(transcripts_df[["Entrez_IDs"]][are_this_ID])))
  this_ID_ENSG <- unique(transcripts_df[["Ensembl_gene_IDs"]][are_this_ID])
  stopifnot(length(this_ID_ENSG) == 1)
  are_this_ENSG <- transcripts_df[["Ensembl_gene_IDs"]] %in% this_ID_ENSG
  this_ENSG_entrez <- setdiff(transcripts_df[["Entrez_IDs"]][are_this_ENSG], NA)
  stopifnot(length(this_ENSG_entrez) == 1)
  stopifnot(!(grepl(",", this_ENSG_entrez, fixed = TRUE)))
  this_ENSG_gene_ID <- unique(transcripts_df[["Gene_IDs"]][are_this_ENSG & have_entrez])
  stopifnot(length(this_ENSG_gene_ID) == 1)
  stopifnot(!(grepl(",", this_ENSG_entrez, fixed = TRUE)))
  transcripts_df[["Gene_IDs"]][are_this_ENSG] <- this_ENSG_gene_ID
  transcripts_df[["Entrez_IDs"]][are_this_ENSG] <- this_ENSG_entrez
}


### Manually correct some genes that still have conflicting gene symbols

change_symbols <- c(
  "AC139887.2" = "LOC100129917",
  "CASC19"     = "PCAT2",
  "AC104211.1" = "FLJ46284",
  "AL117329.1" = "LOC284930"
)

use_transcripts_df <- transcripts_df

for (i in seq_along(change_symbols)) {
  is_to_change <- use_transcripts_df[["Gene_symbols"]] %in% names(change_symbols)[[i]]
  stopifnot(sum(is_to_change) == 1)
  use_transcripts_df[["Gene_symbols"]][is_to_change] <- change_symbols[[i]]
}




# Find overlaps with peaks ------------------------------------------------

HepG2_CLIPseq_df <- AnnotatePeaks(HepG2_peaks_df, use_transcripts_df)
K562_CLIPseq_df  <- AnnotatePeaks(K562_peaks_df, use_transcripts_df)




# Export data -------------------------------------------------------------

WriteTable(HepG2_CLIPseq_df,
           file.path(file_output_directory, "HepG2_CLIPseq_peaks.tsv")
           )

WriteTable(K562_CLIPseq_df,
           file.path(file_output_directory, "K562_CLIPseq_peaks.tsv")
           )




# Save data ---------------------------------------------------------------

save(list = c("HepG2_CLIPseq_df", "K562_CLIPseq_df", "transcripts_df"),
     file = file.path(RData_directory, "Map HNRNPK CLIP-seq peaks to genes.RData")
     )


