### 7th November 2022 ###


# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(general_RData_directory, "22) Compile data on essential genes - essential_df.RData"))



# Prepare TSS_df ----------------------------------------------------------

TSS_df <- PrepareTSSDf(all_TSS_df,
                       only_consistent_chromosomes = TRUE,
                       only_protein_coding         = TRUE,
                       only_best_TSS               = TRUE,
                       check_entrezs_and_symbols   = TRUE
                       )
include_chromosomes <- c(paste0("chr", 1:23), "chrX")
are_included <- (!(is.na(TSS_df[, "Score"]))) &
                (TSS_df[, "Chromosome"] %in% include_chromosomes)
TSS_df <- TSS_df[are_included, ]
row.names(TSS_df) <- NULL



# Tidy TSS_df -------------------------------------------------------------

are_all_identical <- sapply(names(TSS_df), function(x) length(unique(TSS_df[, x])) == 1)
omit_columns <- c(names(TSS_df)[are_all_identical],
                  c("Start", "End", "Original_symbols", "Source")
                  )
TSS_df <- TSS_df[, !(names(TSS_df) %in% omit_columns)]

names(TSS_df)[names(TSS_df) == "Entrez_IDs"] <- "Entrez_ID"
names(TSS_df)[names(TSS_df) == "Gene_symbols"] <- "Gene_symbol"

TSS_df[, "Entrez_ID"] <- as.integer(TSS_df[, "Entrez_ID"])



# Order TSS_df ------------------------------------------------------------

new_order <- order(match(TSS_df[, "Chromosome"], include_chromosomes),
                   TSS_df[, "TSS"]
                   )
TSS_df <- TSS_df[new_order, ]
row.names(TSS_df) <- NULL

strands_fac <- factor(TSS_df[, "Strand"])
TSS_df_list_list <- lapply(split(TSS_df, strands_fac), function(x) {
  split(x, factor(x[, "Chromosome"], levels = include_chromosomes))
})



# Extend TSS_df with gene essentiality data -------------------------------

essential_fraction <- essential_df[["CRISPR_num_essential"]] /
                      essential_df[["CRISPR_num_cell_lines"]]

are_essential <- essential_fraction > 0.95
are_non_essential <- (essential_fraction == 0) &
                     (essential_df[, "CRISPR_mean_probability"] <= 0.03)
three_categories <- ifelse(are_essential, "Essential",
                           ifelse(are_non_essential, "Non-essential", "Other")
                           )

matches_vec <- match(TSS_df[, "Entrez_ID"], essential_df[, "Entrez_ID"])
TSS_df[, "Essentiality"] <- three_categories[matches_vec]



# Identify all gene pairs -------------------------------------------------

max_distance <- 20000L
mat_list <- lapply(c("-", "+"), function(current_strand) {
  print(current_strand)
  if (current_strand == "+") {
    opposite_strand <- "-"
  } else {
    opposite_strand <- "+"
  }
  mat_list <- lapply(include_chromosomes, function(current_chromosome) {
    print(current_chromosome)
    current_df <- TSS_df_list_list[[current_strand]][[current_chromosome]]
    opposite_df <- TSS_df_list_list[[opposite_strand]][[current_chromosome]]
    mat_list <- lapply(seq_len(nrow(current_df)), function(x) {
      distances_vec <- (opposite_df[["TSS"]] - current_df[["TSS"]][[x]])
      are_within_range <- (distances_vec > 0) & (distances_vec <= max_distance)
      if (any(are_within_range)) {
        num_genes <- sum(are_within_range)
        results_mat <- cbind(
          "Entrez_ID_1" = rep(current_df[["Entrez_ID"]][[x]], num_genes),
          "Entrez_ID_2" = opposite_df[, "Entrez_ID"][are_within_range]
        )
      } else {
        results_mat <- NULL
      }
      return(results_mat)
    })
    results_mat <- do.call(rbind, mat_list)
    return(results_mat)
  })
  results_mat <- do.call(rbind, mat_list)
  return(results_mat)
})
all_mat <- do.call(rbind, mat_list)



# Annotate gene pairs -----------------------------------------------------

matches_1 <- match(all_mat[, "Entrez_ID_1"], TSS_df[, "Entrez_ID"])
matches_2 <- match(all_mat[, "Entrez_ID_2"], TSS_df[, "Entrez_ID"])
include_columns <- c("Gene_symbol", "Chromosome", "Strand", "TSS", "Essentiality")
df_1 <- TSS_df[matches_1, include_columns]
names(df_1) <- paste0(include_columns, "_1")
df_2 <- TSS_df[matches_2, include_columns]
names(df_2) <- paste0(include_columns, "_2")

pairs_df <- data.frame(
  "Entrez_ID_1"   = all_mat[, "Entrez_ID_1"],
  "Entrez_ID_2"   = all_mat[, "Entrez_ID_2"],
  df_1,
  df_2,
  row.names = NULL
)
stripped_columns <- sapply(strsplit(names(pairs_df), "_", fixed = TRUE), "[[", 1)
pairs_df <- pairs_df[, order(match(stripped_columns, stripped_columns))]

stopifnot(identical(pairs_df[, "Chromosome_1"], pairs_df[, "Chromosome_2"]))
stopifnot(!(any(pairs_df[, "Strand_1"] == pairs_df[, "Strand_2"])))

pairs_df <- pairs_df[, names(pairs_df) != "Chromosome_2"]
names(pairs_df)[names(pairs_df) == "Chromosome_1"] <- "Chromosome"

pairs_df[, "Distance"] <- pairs_df[, "TSS_2"] - pairs_df[, "TSS_1"]



# Re-order pairs_df -------------------------------------------------------

column_index <- which(names(pairs_df) == "TSS_2")
indices_vec <- unique(c(which(seq_len(ncol(pairs_df)) <= column_index),
                        which(names(pairs_df) == "Distance"),
                        which(seq_len(ncol(pairs_df)) > column_index)
                        ))
pairs_df <- pairs_df[, indices_vec]

new_order <- order(match(pairs_df[, "Chromosome"], include_chromosomes),
                   pairs_df[, "TSS_1"]
                   )
pairs_df <- pairs_df[new_order, ]

row.names(pairs_df) <- NULL



# Classify essentiality of pairs ------------------------------------------

three_classes <- ifelse(((pairs_df[, "Essentiality_1"] == "Essential") &
                         (pairs_df[, "Essentiality_2"] == "Non-essential")) |
                        ((pairs_df[, "Essentiality_1"] == "Non-essential") &
                         (pairs_df[, "Essentiality_2"] == "Essential")),
                        "Discordant",
                        ifelse(((pairs_df[, "Essentiality_1"] == "Non-essential") &
                                (pairs_df[, "Essentiality_2"] == "Non-essential")),
                               "Both non-essential", "Mixed"
                               )
                        )
pairs_df[, "Combination"] <- three_classes



# Save data ---------------------------------------------------------------

bidirectional_df <- pairs_df

save(bidirectional_df,
     file = file.path(general_RData_directory,
                      "24) Enumerate pairs of genes at bidirectional promoters.RData"
                      )
     )




