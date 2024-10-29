## 2024-06-10


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "06_assigning_sgRNAs_to_plasmids.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2024-05-10 - prepooled vs postpooled - Nanopore")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "05_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "07_assign_sgRNAs_to_plasmids__num_reads_mat.RData"))
load(file.path(rdata_dir, "07_assign_sgRNAs_to_plasmids__nano_df.RData"))



# Define functions --------------------------------------------------------

TwoWithoutSwitchOr3OrMore <- function(char_vec) {
  if (sum(is.na(char_vec)) == 2) {
    AtLeastNumGuides(char_vec, num_guides = 2)
  } else {
    AtLeastNumGuides(char_vec, num_guides = 3)
  }
}


CustomMakeCountsDf <- function(sg_df, mapped_df) {

  plasmids_mat <- as.matrix(mapped_df[, paste0("Plasmid_sg", 1:4)])
  two_or_more_mat <- t(apply(plasmids_mat, 1, AtLeastNumGuides, num_guides = 2))
  two_filtered_mat <- t(apply(plasmids_mat, 1, TwoWithoutSwitchOr3OrMore))
  three_or_more_mat <- t(apply(plasmids_mat, 1, AtLeastNumGuides, num_guides = 3))

  have_no_switch  <- mapped_df[, "Num_template_switches"] == 0
  no_switch_all_4 <- have_no_switch & (mapped_df[, "Num_matched_sgRNAs"] == 4)

  counts_all4_vec           <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                         plasmids_mat[no_switch_all_4, ],
                                         )
  counts_3ormore_vec        <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                         three_or_more_mat,
                                         )
  counts_2_unknown_or_3_vec <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                         two_filtered_mat
                                         )
  counts_2ormore_vec        <- GetCounts(sg_sequences_df[, "Plasmid_ID"],
                                         two_or_more_mat
                                         )

  use_columns <- c("Plasmid_ID", "Gene_symbol", "Entrez_ID",
                   "TSS_ID", "Is_main_TSS", "Plate_ID",
                   "Well_number", "Is_obsolete"
                   )
  counts_df <- data.frame(
    "Count_perfect"        = counts_all4_vec,
    "Count_3_or_more"      = counts_3ormore_vec,
    "Count_2_unknown_or_3" = counts_2_unknown_or_3_vec,
    "Count_2_or_more"      = counts_2ormore_vec,
    sg_sequences_df[, use_columns],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(counts_df)
}



# Obtain read counts per plasmid -----------------------------------------

samples_fac <- interaction(nano_df[, "Condition"], nano_df[, "Replicate"], sep = "_", lex.order = TRUE)
levels(samples_fac) <- paste0(rep(c("prepool_T0", "prepool_T12", "postpool_T0", "postpool_T12"), each = 2),
                              rep(paste0("_rep", 1:2), times = 4)
                              )

counts_df_list <- sapply(levels(samples_fac), function(x) {
  are_this_level <- samples_fac == x
  message("Producing read counts from ", sum(are_this_level), " reads for the sample '", x, "'...")
  use_df <- nano_df[are_this_level, ]
  row.names(use_df) <- NULL
  CustomMakeCountsDf(sg_sequences_df, use_df)
}, simplify = FALSE)



# Combine count data ------------------------------------------------------

prefix_vec <- c(
  "all4"           = "perfect",
  "3plus"          = "3_or_more",
  "2strictOR3plus" = "2_unknown_or_3",
  "2plus"          = "2_or_more"
)

counts_mat <- do.call(cbind, lapply(names(prefix_vec), function(x) {
  select_column <- prefix_vec[[x]]
  counts_mat <- do.call(cbind, lapply(counts_df_list, function(x) x[, paste0("Count_", select_column)]))
  colnames(counts_mat) <- paste0("Count_", x, "_", names(counts_df_list))
  return(counts_mat)
}))

counts_df <- data.frame(counts_df_list[[1]][, !(startsWith(names(counts_df_list[[1]]), "Count_"))], counts_mat)



# Count the number of reads in each category per sample -------------------

num_mapped_mat <- do.call(cbind, lapply(1:4, function(x) tapply(!(is.na(nano_df[, paste0("Plasmid_sg", x)])), samples_fac, sum)))
colnames(num_mapped_mat) <- paste0("Num_mapped_sg", 1:4)

num_1MM_mat <- do.call(cbind, lapply(1:4, function(x) tapply(nano_df[, paste0("Num_MM_sg", x)], samples_fac, sum, na.rm = TRUE)))
colnames(num_1MM_mat) <- paste0("Num_1MM_sg", 1:4)

num_matched_mat <- do.call(cbind, lapply(1:4, function(x) tapply(nano_df[, "Num_matched_sgRNAs"], samples_fac, function(y) sum(y == x))))
colnames(num_matched_mat) <- c("Num_1_sgRNA", paste0("Num_", 2:4, "_sgRNAs"))

num_any_mapped_mat <- cbind(
  num_reads_mat[, c("Total_num_reads", "Num_mapped_reads")],
  num_mapped_mat,
  num_1MM_mat,
  num_matched_mat
)



# Produce per-sample metrics for fully mapped reads -----------------------

have_all4 <- nano_df[, "Num_matched_sgRNAs"] == 4


num_mismatched_all4_vec <- rowSums(nano_df[have_all4, paste0("Num_MM_sg", 1:4)])
num_mismatches_mat <- do.call(cbind, lapply(0:4, function(x) tapply(num_mismatched_all4_vec, samples_fac[have_all4], function(y) sum(y == x))))
colnames(num_mismatches_mat) <- c("Num_0_mismatched_sgRNAs", "Num_1_mismatched_sgRNA", paste0("Num_", 2:4, "_mismatched_sgRNAs"))


num_targeted_mat <- do.call(cbind, lapply(1:4, function(x) {
  tapply(nano_df[, "Num_targeted_plasmids"][have_all4], samples_fac[have_all4], function(y) sum(y == x))
}))
colnames(num_targeted_mat) <- c("Num_1_targeted_plasmid", paste0("Num_", 2:4, "_targeted_plasmids"))


num_switches_mat <- do.call(cbind, lapply(0:3, function(x) {
  tapply(nano_df[, "Num_template_switches"][have_all4], samples_fac[have_all4], function(y) sum(y == x))
}))
colnames(num_switches_mat) <- c("Num_0_template_switches", "Num_1_template_switch", paste0("Num_", 2:3, "_template_switches"))


num_switchbacks_mat <- do.call(cbind, lapply(0:2, function(x) {
  tapply(nano_df[, "Num_switch_backs"][have_all4], samples_fac[have_all4], function(y) sum(y == x))
}))
colnames(num_switchbacks_mat) <- c("Num_0_switch_backs", "Num_1_switch_back", "Num_2_switch_backs")



GetNumMatchingMat <- function(input_df, samples_fac, sg1, sg2) {
  stopifnot(all(c(sg1, sg2) %in% 1:4), nrow(input_df) == length(samples_fac))
  have_all4 <- input_df[, "Num_matched_sgRNAs"] == 4
  are_eligible <- have_all4
  are_eligible[have_all4] <- nano_df[, paste0("Plasmid_sg", sg1)][have_all4] == nano_df[, paste0("Plasmid_sg", sg2)][have_all4]
  correct_plasmids_vec <- nano_df[, paste0("Plasmid_sg", sg1)][are_eligible]
  other_sgs <- setdiff(1:4, c(sg1, sg2))
  first_other_matches_vec <- nano_df[, paste0("Plasmid_sg", other_sgs[[1]])][are_eligible] == correct_plasmids_vec
  second_other_matches_vec <- nano_df[, paste0("Plasmid_sg", other_sgs[[2]])][are_eligible] == correct_plasmids_vec
  num_matching_vec <- ifelse(first_other_matches_vec & second_other_matches_vec,
                             4L,
                             ifelse(first_other_matches_vec | second_other_matches_vec, 3L, 2L)
                             )
  num_matching_mat <- do.call(cbind, lapply(2:4, function(x) {
    tapply(num_matching_vec, samples_fac[are_eligible], function(y) sum(y == x))
  }))
  colnames(num_matching_mat) <- paste0("Num_", 2:4, "_matching_sg", sg1, "and", "sg", sg2)
  return(num_matching_mat)
}

num_all4_mapped_mat <- cbind(
  "Num_perfect_reads" = tapply(have_all4, samples_fac, sum),
  num_mismatches_mat,
  num_targeted_mat,
  num_switches_mat,
  num_switchbacks_mat,
  "Num_switch_sg1_to_sg2" = tapply(nano_df[, "Switch_sg1_to_sg2"][have_all4], samples_fac[have_all4], sum),
  "Num_switch_sg2_to_sg3" = tapply(nano_df[, "Switch_sg2_to_sg3"][have_all4], samples_fac[have_all4], sum),
  "Num_switch_sg3_to_sg4" = tapply(nano_df[, "Switch_sg3_to_sg4"][have_all4], samples_fac[have_all4], sum),
  GetNumMatchingMat(nano_df, samples_fac, 1, 2),
  GetNumMatchingMat(nano_df, samples_fac, 2, 3),
  GetNumMatchingMat(nano_df, samples_fac, 3, 4)
)



# Save data ---------------------------------------------------------------

save(list = c("counts_df", "num_any_mapped_mat", "num_all4_mapped_mat"),
     file = file.path(rdata_dir, "08_produce_read_counts_and_metrics.RData")
     )

