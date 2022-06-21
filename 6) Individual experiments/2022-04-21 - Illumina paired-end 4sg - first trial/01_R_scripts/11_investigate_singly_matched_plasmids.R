### 2022-06-20



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts", "1) R functions")
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define paths ------------------------------------------------------------

experiments_directory   <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir             <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 4sg - first trial")
rdata_dir               <- file.path(project_dir, "03_R_objects")
tables_dir              <- file.path(project_dir, "04_output_data", "Tables")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__lumi_df.RData"))

load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))



# Define functions --------------------------------------------------------

ExpandUnmatchedReads <- function(input_df, unmatched_sg) {
  stopifnot(unmatched_sg %in% 1:2)
  if (unmatched_sg == 2) {
    matched_sg <- 1L
  } else {
    matched_sg <- 2L
  }
  are_eligible <- (!(is.na(input_df[, paste0("Plasmid_sg", matched_sg)]))) &
                  is.na(input_df[, paste0("Plasmid_sg", unmatched_sg)])
  plasmids_vec <- input_df[, paste0("Plasmid_sg", matched_sg)][are_eligible]
  reads_vec <- input_df[, paste0("Aligned_read_sg", unmatched_sg)][are_eligible]
  plasmid_ID_splits <- strsplit(plasmids_vec, ", ", fixed = TRUE)
  expanded_df <- data.frame(
    "Plasmid_ID"            = unlist(plasmid_ID_splits),
    "Num_possible_plasmids" = rep(lengths(plasmid_ID_splits), lengths(plasmid_ID_splits)),
    "Sequence"              = rep(reads_vec, lengths(plasmid_ID_splits)),
    stringsAsFactors        = FALSE
  )
  return(expanded_df)
}


TabulateSequences <- function(input_df, plasmid_ID, fraction_cutoff = NULL) {
  are_this_plasmid <- input_df[, "Plasmid_ID"] == plasmid_ID
  if (!(any(are_this_plasmid))) {
    stop("No reads were found for this plasmid!")
  }
  sg_table <- table(input_df[, "Sequence"][are_this_plasmid])
  num_occurrences <- as.integer(sg_table)
  fractions_vec <- num_occurrences / sum(are_this_plasmid)
  new_order <- order(num_occurrences, decreasing = TRUE)
  if (!(is.null(fraction_cutoff))) {
    exceed_cutoff <- fractions_vec[new_order] >= fraction_cutoff
    if (!(any(exceed_cutoff))) {
      new_order <- new_order[[1]]
    } else {
      new_order <- new_order[exceed_cutoff]
    }
  }
  results_df <- data.frame(
    "Plasmid_ID"      = plasmid_ID,
    "Sequence"        = names(sg_table)[new_order],
    "Num_occurrences" = num_occurrences[new_order],
    "Fraction"        = fractions_vec[new_order],
    stringsAsFactors  = FALSE
  )
  return(results_df)
}


TabulateAllPlasmids <- function(plasmids_vec, expanded_df) {
  num_total <- length(plasmids_vec)
  digits_width <- nchar(as.character(num_total))
  only_sg1_df_list <- vector(mode = "list", length = num_total)
  results_df_list <- lapply(seq_along(plasmids_vec), function(i) {
    if ((i %in% 1) || ((i %% 50) == 0)) {
      message(paste0("Processing plasmid #", format(i, width = digits_width),
                     " out of ", num_total, " (", plasmids_vec[[i]], ")..."
                     ))
    }
    current_df <- TabulateSequences(expanded_df, plasmids_vec[[i]], fraction_cutoff = 0.2)
    current_df[, "Num_sequences"] <- nrow(current_df)
    return(current_df)
  })
  results_df <- do.call(rbind.data.frame, c(results_df_list, stringsAsFactors = FALSE))
  return(results_df)
}



# Identify plasmids for which only one sgRNA was matched ------------------

only_has_sg1 <- (counts_df[, "Data_contains_sg1"] & (!(counts_df[, "Data_contains_sg2"])))
only_has_sg2 <- (counts_df[, "Data_contains_sg2"] & (!(counts_df[, "Data_contains_sg1"])))

counts_switch_vec <- counts_df[, "Sum_MaySwitch_xMM"] - counts_df[, "Sum_NoEvidentSwitch_xMM"]
fraction_switch_vec <- counts_switch_vec / counts_df[, "Sum_MaySwitch_xMM"]
are_all_switched <- fraction_switch_vec == 1

only_sg1_plasmids <- counts_df[, "Plasmid_ID"][only_has_sg1 & (!(are_all_switched))]
only_sg2_plasmids <- counts_df[, "Plasmid_ID"][only_has_sg2 & (!(are_all_switched))]



# Tabulate the most common sequences for the second sgRNA -----------------

only_sg1_expanded_df <- ExpandUnmatchedReads(lumi_df, 2)
only_sg2_expanded_df <- ExpandUnmatchedReads(lumi_df, 1)

only_sg1_df <- TabulateAllPlasmids(only_sg1_plasmids, only_sg1_expanded_df)
only_sg2_df <- TabulateAllPlasmids(only_sg2_plasmids, only_sg2_expanded_df)
single_sg_df <- rbind.data.frame(
  data.frame(only_sg1_df[1], "Found_sg" = 1L, only_sg1_df[, 2:ncol(only_sg1_df)]),
  data.frame(only_sg2_df[1], "Found_sg" = 2L, only_sg2_df[, 2:ncol(only_sg2_df)]),
  stringsAsFactors = FALSE
)



# Search the human genome for matches to sgRNA sequences ------------------

unique_sequences <- unique(single_sg_df[, "Sequence"])
stopifnot(!(any(grepl("N", unique_sequences, fixed = TRUE))))

sequences_df <- FindSequences(unique_sequences, max.mismatch = 0)
sequences_df[["PAM"]] <- GetNGGPAM(sequences_df)



# Summarize sgRNA matches -------------------------------------------------

genome_search_df <- SummarizeFoundSequencesDf(sequences_df, all_sequences = unique_sequences)

matches_vec <- match(single_sg_df[, "Sequence"], genome_search_df[, "Sequence"])

use_columns <- c("Num_0MM", "Locations_0MM")
loci_df <- genome_search_df[matches_vec, use_columns]
row.names(loci_df) <- NULL



# Find all TSSs targeted by each sgRNA ------------------------------------

matches_vec <- match(single_sg_df[, "Plasmid_ID"], CRISPRoff_df[, "Plasmid_ID"])
single_df <- data.frame(CRISPRoff_df[matches_vec, c("Plasmid_ID", "Gene_symbol", "Entrez_ID", "sgID_AB")],
                        single_sg_df[, names(single_sg_df) != "Plasmid_ID"],
                        loci_df,
                        stringsAsFactors = FALSE,
                        row.names = NULL
                        )
location_columns <- c("Chromosome", "Strand", "Start", "End")
for (location_column in location_columns) {
  single_df[, location_column] <- NA
}
nearby_list <- AlignSummaryDf(FindNearbyTSSs, single_df, all_TSS_df)
use_columns <- c("Affects_intended_gene", "Unintended_Entrez_IDs", "Unintended_gene_symbols")
single_df <- data.frame(single_df[, !(names(single_df) %in% location_columns)],
                          nearby_list[["summary_df"]][,  use_columns],
                          stringsAsFactors = FALSE
                          )



# Add information on the reference sequences ------------------------------

sg1_vec <- substr(CRISPRoff_df[, "protospacer_A"][matches_vec], 2, 20)
sg2_vec <- substr(CRISPRoff_df[, "protospacer_B"][matches_vec], 2, 20)
single_df <- data.frame(
  single_df[, 1:4],
  "Found_sg"             = ifelse(single_df[, "Found_sg"] == 1, "A", "B"),
  "Found_sg_reference"   = ifelse(single_df[, "Found_sg"] == 1, sg1_vec, sg2_vec),
  "Missing_sg_reference" = ifelse(single_df[, "Found_sg"] == 1, sg2_vec, sg1_vec),
  "Observed_sequence"    = single_df[, "Sequence"],
  single_df["Affects_intended_gene"],
  single_df[, setdiff(7:ncol(single_df), which(names(single_df) == "Affects_intended_gene"))],
  stringsAsFactors = FALSE
)



# Remove alternative sequences that are likely sequencing errors ----------

plasmids_vec <- single_df[, "Plasmid_ID"]
plasmids_fac <- factor(plasmids_vec, levels = unique(plasmids_vec))
single_df_list <- split(single_df, plasmids_fac)
single_df_list <- lapply(single_df_list, function(x) {
  do_target <- x[, "Affects_intended_gene"] %in% TRUE
  have_matches <- x[, "Num_0MM"] > 0
  if (!(any(have_matches))) {
    max_fraction <- max(x[, "Fraction"])
    x <- x[x[, "Fraction"] == max_fraction, ]
  } else if (any(do_target) && (!(all(do_target)))) {
    use_cutoff <- min(x[, "Fraction"][do_target])
    x <- x[x[, "Fraction"] >= use_cutoff, ]
  } else {
    if (any(have_matches) && (!(all(have_matches)))) {
      x <- x[have_matches, ]
    }
  }
  x[, "Num_sequences"] <- nrow(x)
  return(x)
})

single_df <- do.call(rbind.data.frame,
                     c(single_df_list,
                       stringsAsFactors = FALSE,
                       make.row.names = FALSE
                     ))



# Export data -------------------------------------------------------------

export_df <- single_df
export_df[, "Affects_intended_gene"] <- ifelse(export_df[, "Affects_intended_gene"],
                                               "Yes", "No"
                                               )
write.table(export_df, file = file.path(tables_dir, "Only_one_sgRNA_found.tsv"),
            sep = "\t", na = "", row.names = FALSE, quote = FALSE
            )



# Save data ---------------------------------------------------------------

save(list = "single_df",
     file = file.path(rdata_dir, "11_investigate_singly_matched_plasmids.RData")
     )



