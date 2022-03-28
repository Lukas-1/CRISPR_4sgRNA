## 2022-02-17


# Load packages and source code -------------------------------------------

root_directory        <- "~/CRISPR/6) Individual experiments"
R_functions_directory <- file.path(root_directory, "2020-08-29 - PacBio - first 384-well plate/1) R functions")
source(file.path(R_functions_directory, "02) Analyzing reads.R"))



# Define paths ------------------------------------------------------------

root_directory <- file.path(root_directory, "2022-01-05 - first nanopore sequencing run")
rdata_dir <- file.path(root_directory, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "06_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "07_look_up_aligned_sgRNAs.RData"))





# Prepare the CRISPR library for identifying non-unique plasmids ----------

upper_library_mat <- toupper(as.matrix(sg_sequences_df[, paste0("Sequence_sg", 1:4)]))
plasmids_list_list <- lapply(1:4, function(x) {
  split(sg_sequences_df[, "Plasmid_ID"],
        factor(upper_library_mat[, x], levels = unique(upper_library_mat[, x]))
        )
})
are_unique_library_mat <- do.call(cbind, lapply(1:4, function(x) {
  seq_vec <- upper_library_mat[, x]
  !(duplicated(seq_vec) | duplicated(seq_vec, fromLast = TRUE))
}))

sg_combos_vec <- do.call(paste, c(as.list(data.frame(upper_library_mat)), sep = "_"))
combo_plasmids_list <- split(sg_sequences_df[, "Plasmid_ID"], sg_combos_vec)
combo_plasmids_vec <- vapply(combo_plasmids_list, paste0, collapse = ", ", "")




# Identify sequences that map to more than one plasmid --------------------

sequences_mat <- do.call(cbind, lapply(1:4, function(x) {
  ifelse(matched_df[, paste0("Num_MM_sg", x)] %in% 0,
         matched_df[, paste0("Aligned_read_sg", x)],
         matched_df[, paste0("Correct_sgRNA_sg", x)]
         )
}))
colnames(sequences_mat) <- paste0("Sequence_sg", 1:4)

have_guide <- rowSums(is.na(sequences_mat)) != 4
have_guide_mat <- sequences_mat[have_guide, ]

all_unique_mat <- do.call(cbind, lapply(1:4, function(x) {
  matches_vec <- match(have_guide_mat[, x], upper_library_mat[, x])
  are_unique_library_mat[matches_vec, x]
}))
are_unique <- rowSums(is.na(all_unique_mat) | all_unique_mat) == 4




# Map reads to plasmids ---------------------------------------------------

possible_plasmids_list_list <- lapply(1:4, function(x) {
  lapply(have_guide_mat[!(are_unique), x], function(y) {
    if (is.na(y)) {
      NA
    } else {
      plasmids_list_list[[x]][[y]]
    }
  })
})

short_NA_mat <- do.call(cbind, lapply(1:4, function(x) {
  is.na(possible_plasmids_list_list[[x]])
}))

IntersectLists <- function(index_1, index_2) {
  both_non_NA <- !(short_NA_mat[, index_1] | short_NA_mat[, index_2])
  intersect_list <- mapply(intersect,
                           possible_plasmids_list_list[[index_1]][both_non_NA],
                           possible_plasmids_list_list[[index_2]][both_non_NA]
                           )
  have_intersect <- lengths(intersect_list) > 0
  possible_plasmids_list_list[[index_1]][both_non_NA][have_intersect] <- intersect_list[have_intersect]
  possible_plasmids_list_list[[index_2]][both_non_NA][have_intersect] <- intersect_list[have_intersect]
  return(possible_plasmids_list_list)
}

old_possible_plasmids_list_list <- NULL

for (h in 1:10) {
  if (identical(old_possible_plasmids_list_list, possible_plasmids_list_list)) {
    break
  }
  message("Iteration #", h, "...")
  old_possible_plasmids_list_list <- possible_plasmids_list_list
  for (i in 1:3) {
    possible_plasmids_list_list <- IntersectLists(i, i + 1)
  }
  for (i in 1:2) {
    possible_plasmids_list_list <- IntersectLists(i, i + 2)
  }
  possible_plasmids_list_list <- IntersectLists(1, 4)
}




# Collapse lists of plasmids into comma-separated vectors -----------------

plasmids_mat <- do.call(cbind, lapply(1:4, function(x) {
  plasmids_vec <- rep(NA, nrow(have_guide_mat))
  matches_vec <- match(have_guide_mat[are_unique, x], upper_library_mat[, x])
  plasmids_vec[are_unique] <- sg_sequences_df[, "Plasmid_ID"][matches_vec]
  collapsed_vec <- vapply(possible_plasmids_list_list[[x]],
                          function(y) if (all(is.na(y))) NA_character_ else paste0(y, collapse = ", "),
                          ""
                          )
  plasmids_vec[!(are_unique)] <- collapsed_vec
  return(plasmids_vec)
}))
colnames(plasmids_mat) <- paste0("Plasmid_sg", 1:4)





# Map plasmids to genes ---------------------------------------------------

MakeGenesMat <- function(genes_column) {
  results_mat <- do.call(cbind, lapply(1:4, function(x) {
    plasmid_splits <- strsplit(plasmids_mat[, x], ", ", fixed = TRUE)
    long_vec <- unlist(plasmid_splits, use.names = FALSE)
    indices_vec <- rep(seq_along(plasmid_splits), lengths(plasmid_splits))
    matches_vec <- match(long_vec, sg_sequences_df[, "Plasmid_ID"])
    lookup_vec <- sg_sequences_df[, genes_column]
    if (genes_column == "Entrez_ID") {
      lookup_vec <- ifelse(is.na(lookup_vec),
                           sg_sequences_df[, "Gene_symbol"],
                           lookup_vec
                           )
    }
    genes_vec <- lookup_vec[matches_vec]
    results_vec <- rep(NA, length(plasmid_splits))
    genes_list <- split(genes_vec, indices_vec)
    are_NA <- is.na(genes_list)
    results_vec[!(are_NA)] <- vapply(genes_list[!(are_NA)],
                                     function(x) paste0(unique(x), collapse = ", "),
                                     ""
                                     )
    return(results_vec)
  }))
  return(results_mat)
}

symbols_mat <- MakeGenesMat("Gene_symbol")
colnames(symbols_mat) <- paste0("Symbol_sg", 1:4)
entrezs_mat <- MakeGenesMat("Entrez_ID")
colnames(entrezs_mat) <- paste0("Entrez_sg", 1:4)



# Identify template switches ----------------------------------------------

multiple_list <- as.list(data.frame(t(plasmids_mat), stringsAsFactors = FALSE))
unique_list <- lapply(multiple_list, function(x) unique(x[!(is.na(x))]))

rle_list <- unique_list
have_switch <- lengths(unique_list) > 1
rle_list[have_switch] <- lapply(multiple_list[have_switch],
                                function(x) rle(x[!(is.na(x))])[["values"]]
                                )
switch_back <- lengths(rle_list) != lengths(unique_list)
switch_twice <- lengths(rle_list) > (lengths(unique_list) + 1)



# Prepare additional relevant data ----------------------------------------

num_MM_mat <- as.matrix(matched_df[, paste0("Num_MM_sg", 1:4)])[have_guide, ]
rownames(num_MM_mat) <- NULL

qual_string_mat <- as.matrix(matched_df[, paste0("Quality_sg", 1:4)])[have_guide, ]
qual_string_mat <- gsub(" ", "", qual_string_mat, fixed = TRUE)
NA_mat <- is.na(have_guide_mat)
qualities_vec <- GetMeanQuality(qual_string_mat[!(NA_mat)])
qualities_mat <- matrix(nrow = nrow(NA_mat), ncol = 4)
qualities_mat[!(NA_mat)] <- qualities_vec
colnames(qualities_mat) <- paste0("Quality_sg", 1:4)



# Create a data frame combining all relevant data -------------------------

nano_df <- data.frame(
  "Read_number" = which(have_guide),
  matched_df[have_guide, 1:4],
  have_guide_mat, qualities_mat, num_MM_mat,
  plasmids_mat,
  symbols_mat, entrezs_mat,
  "Num_matched_sgRNAs"    = rowSums(!(is.na(have_guide_mat))),
  "Num_targeted_plasmids" = lengths(unique_list),
  "Num_template_switches" = lengths(rle_list) - 1L,
  "Switch_sg1_to_sg2"     = plasmids_mat[, 1] != plasmids_mat[, 2],
  "Switch_sg2_to_sg3"     = plasmids_mat[, 2] != plasmids_mat[, 3],
  "Switch_sg3_to_sg4"     = plasmids_mat[, 3] != plasmids_mat[, 4],
  "Num_switch_backs"      = ifelse(switch_twice, 2L,
                                   ifelse(switch_back == 1, 1L, 0L)
                                   ),
  stringsAsFactors = FALSE
)



# Obtain read counts per plasmid ------------------------------------------

GetCounts <- function(library_plasmids,
                      reads_df,
                      exclude_template_switch = FALSE,
                      min_num_guides = 1L
                      ) {

  plasmids_mat <- as.matrix(reads_df[, paste0("Plasmid_sg", 1:4)])
  if (exclude_template_switch || (min_num_guides >= 2)) {
    are_selected <- rep(TRUE, nrow(reads_df))
    if (exclude_template_switch) {
      are_selected <- reads_df[, "Num_template_switches"] == 0
    }
    are_selected[are_selected] <- reads_df[, "Num_matched_sgRNAs"][are_selected] >= min_num_guides
    plasmids_mat <- plasmids_mat[are_selected, ]
  }

  multiple_list <- as.list(data.frame(t(plasmids_mat), stringsAsFactors = FALSE))
  unique_list <- lapply(multiple_list, function(x) unique(x[!(is.na(x))]))

  non_unique_vec <- table(unlist(unique_list, use.names = FALSE))
  unique_plasmids_list <- strsplit(names(non_unique_vec), ", ", fixed = TRUE)
  plasmids_vec <- unlist(unique_plasmids_list)
  lengths_vec <- lengths(unique_plasmids_list)
  counts_vec <- rep(non_unique_vec, lengths_vec)
  weights_vec <- 1 / rep(lengths_vec, lengths_vec)
  counts_vec <- weights_vec * counts_vec

  matches_vec <- match(library_plasmids, plasmids_vec)
  matched_counts_vec <- counts_vec[matches_vec]
  matched_counts_vec[is.na(matched_counts_vec)] <- 0

  return(matched_counts_vec)
}


counts_all4_vec       <- GetCounts(sg_sequences_df[, "Plasmid_ID"], nano_df,
                                   exclude_template_switch = TRUE,
                                   min_num_guides = 4
                                   )
counts_2ormore_vec    <- GetCounts(sg_sequences_df[, "Plasmid_ID"], nano_df,
                                   exclude_template_switch = TRUE,
                                   min_num_guides = 2
                                   )
counts_noswitch_vec   <- GetCounts(sg_sequences_df[, "Plasmid_ID"], nano_df,
                                   exclude_template_switch = TRUE
                                   )
counts_unselected_vec <- GetCounts(sg_sequences_df[, "Plasmid_ID"], nano_df,
                                   exclude_template_switch = FALSE
                                   )

use_columns <- c("Plasmid_ID", "Gene_symbol", "Entrez_ID",
                 "TSS_ID", "Is_main_TSS", "Plate_ID",
                 "Well_number", "Is_obsolete"
                 )
counts_df <- data.frame(
  "Count_perfect"    = counts_all4_vec,
  "Count_2_filtered" = counts_2ormore_vec,
  "Count_filtered"   = counts_noswitch_vec,
  "Count_unfiltered" = counts_unselected_vec,
  sg_sequences_df[, use_columns],
  stringsAsFactors = FALSE
)




# Explore summary statistics of the Nanopore sequencing results -----------

## Prepare for exploring summary statistics
nano4_df <- nano_df[nano_df[, "Num_matched_sgRNAs"] == 4, ]
row.names(nano4_df) <- NULL
perfect_df <- nano4_df[nano4_df[["Num_template_switches"]] == 0, ]
row.names(perfect_df) <- NULL

## 36.0% of aligned sgRNA sequences could be assigned to an sgRNA from the library.
num_not_NA_total <- sum(!(is.na(sequences_mat)))
num_not_NA_total / (nrow(sequences_mat) * 4)

## Of these, 80.8% were a perfect match, and 19.2% showed a (unique) match with one mismatched base.
## In absolute terms, this corresponds to 29.0% of aligned sequences with a
## perfect match, and 6.9% with one mismatch.
num_not_NA_filtered <- sum(!(is.na(num_MM_mat)))
sum(num_MM_mat, na.rm = TRUE) / num_not_NA_filtered
sum(num_MM_mat == 0, na.rm = TRUE) / (nrow(sequences_mat) * 4)
sum(num_MM_mat, na.rm = TRUE) / (nrow(sequences_mat) * 4)

## 71.9% of reads have at least one gRNA that could be assigned to an sgRNA from the library.
nrow(have_guide_mat) / nrow(sequences_mat)

## Of these, 37.1% had only one matching sgRNA, 33.3% had exactly two matching sgRNAs,
## 22.1% had exactly 3 matching sgRNAs, and 7.5% had 4 matching sgRNAs.
table(nano_df[, "Num_matched_sgRNAs"]) / nrow(nano_df)

## For 5.4% of reads, all 4 sgRNA sequences could be assigned to an sgRNA from the library.
nrow(nano4_df) / nrow(sequences_mat)

## 63.5% of reads with 4 sgRNAs showed a template switch.
## 21.7% showed two or more template switches!
## However, only 6.3% targeted 3 or more plasmids, indicating that switches to another plasmid,
## and then back to the first plasmid, were quite common. This affected 16.5% of reads.
table(nano4_df[, "Num_template_switches"]) / nrow(nano4_df)
sum(nano4_df[, "Num_template_switches"] >= 1) / nrow(nano4_df)
sum(nano4_df[, "Num_template_switches"] >= 2) / nrow(nano4_df)
sum(nano4_df[, "Num_targeted_plasmids"] >= 3) / nrow(nano4_df)
sum(nano4_df[, "Num_switch_backs"] > 0) / nrow(nano4_df)

## The likelihood of a template switch between two neighbouring sgRNAs is 39.7%
## between sg1 and sg2, 22.8% between sg2 and sg3, and 25.5% between sg3 and
## sg4. It makes sense that the switch rate between sg1 and sg2 would be higher,
## since the stretch of DNA between these two sgRNAs is longer (and contains the
## trimethoprim resistance element encoding for dihydrofolate reductase).
table(nano4_df[, "Switch_sg1_to_sg2"]) / nrow(nano4_df)
table(nano4_df[, "Switch_sg2_to_sg3"]) / nrow(nano4_df)
table(nano4_df[, "Switch_sg3_to_sg4"]) / nrow(nano4_df)

## Indeed, for 14.1% of plasmids where sg1 and sg4 matched one plasmid,
## sg2 and/or sg3 matched another plasmid. This has implications for strategies
## where only two sgRNAs are sequenced. Even when targeting sg1 and sg4,
## the prevalence of template switching would still have to be reduced.
sg1_matches_sg4 <- (nano4_df[, "Plasmid_sg1"] == nano4_df[, "Plasmid_sg4"])
double_switch <- sg1_matches_sg4 & (nano4_df[, "Num_template_switches"] > 0)
table(double_switch) / nrow(nano4_df)
sum(nano4_df[, "Plasmid_sg1"] == nano4_df[, "Plasmid_sg4"]) / nrow(nano4_df)
sum(double_switch) / sum(sg1_matches_sg4)


## When including only ideal reads (with four recognized sgRNAs that all map
## to the same plasmid(s)), 83.4% of plasmids in the library were represented
## by at least one read.

perfect_vec <- do.call(paste, c(as.list(data.frame(nano4_df[, paste0("Sequence_sg", 1:4)])), sep = "_"))
table(names(combo_plasmids_vec) %in% perfect_vec) / length(combo_plasmids_vec)

are_present <- counts_df[, "Count_perfect"] > 0
sum(are_present) / nrow(counts_df)
median(counts_df[, "Count_perfect"])

## When including only reads where a) at least two sequences match sgRNAs in
## library and b) no template switches are observed, 99.3% of plasmids are
## represented at least once, and 11.4% of plasmids are represented
## by at least 100 reads.
sum(counts_df[, "Count_2_filtered"] > 0) / nrow(counts_df)
sum(counts_df[, "Count_2_filtered"] >= 100) / nrow(counts_df)
median(counts_df[, "Count_2_filtered"])

## When including only reads where no template switches are observed, 99.7% of plasmids are
## represented at least once, and 44% of plasmids are represented
## by at least 100 reads.
sum(counts_df[, "Count_filtered"] > 0) / nrow(counts_df)
sum(counts_df[, "Count_filtered"] >= 100) / nrow(counts_df)
median(counts_df[, "Count_filtered"])

## When including all reads, 99.8% of plasmids are represented at least once,
## 70.5% of plasmids are represented by at least 100 reads, and
## 88.9% of plasmids are represented by at least 50 reads.
sum(counts_df[, "Count_unfiltered"] > 0) / nrow(counts_df)
sum(counts_df[, "Count_unfiltered"] >= 100) / nrow(counts_df)
sum(counts_df[, "Count_unfiltered"] >= 50) / nrow(counts_df)
median(counts_df[, "Count_unfiltered"])




sg2_matches_sg3 <- (nano4_df[, "Plasmid_sg2"] == nano4_df[, "Plasmid_sg3"])
switch_outside <- sg2_matches_sg3 & (nano4_df[, "Num_template_switches"] > 0)
table(switch_outside) / nrow(nano4_df)
table(nano4_df[, "Num_template_switches"][switch_outside])
sum(nano4_df[, "Plasmid_sg2"] == nano4_df[, "Plasmid_sg3"]) / nrow(nano4_df)
sum(switch_outside) / sum(sg2_matches_sg3)


sg3_matches_sg4 <- (nano4_df[, "Plasmid_sg3"] == nano4_df[, "Plasmid_sg4"])
switch_outside <- sg3_matches_sg4 & (nano4_df[, "Num_template_switches"] > 0)
table(switch_outside) / nrow(nano4_df)
sum(nano4_df[, "Plasmid_sg3"] == nano4_df[, "Plasmid_sg4"]) / nrow(nano4_df)
sum(switch_outside) / sum(sg3_matches_sg4)


sum(nano4_df[, "Plasmid_sg1"] != nano4_df[, "Plasmid_sg2"]) / nrow(nano4_df)



# Save data ---------------------------------------------------------------

save(nano_df,
     file = file.path(rdata_dir, "08_assign_sgRNAs_to_plasmids.RData")
     )


