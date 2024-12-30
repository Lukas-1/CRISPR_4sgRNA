## 2024-11-18


# Preparatory functions ---------------------------------------------------

GetSgPairsMat <- function() {
  sg_pairs_mat <- combn(1:4, 2)
  colnames(sg_pairs_mat) <- paste0("sg", sg_pairs_mat[1, ], "_sg", sg_pairs_mat[2, ])
  pairs_in_order <- c("sg1_sg2", "sg2_sg3", "sg3_sg4",
                      "sg1_sg3", "sg2_sg4",
                      "sg1_sg4"
                      )
  sg_pairs_mat <- sg_pairs_mat[, pairs_in_order]
  return(sg_pairs_mat)
}



# Functions for computing error rates -------------------------------------

CategorDfToMat <- function(categor_df, use_column) {
  unique_features <- unique(categor_df[, "Feature"])
  num_features <- length(unique_features)
  return(matrix(categor_df[, use_column],
                ncol = num_features,
                nrow = nrow(categor_df) / num_features,
                byrow = TRUE,
                dimnames = list(NULL, unique_features)
                ))
}


CompareSwitchedNonswitched <- function(categor_df,
                                       mapped_df,
                                       sg_A,
                                       sg_B,
                                       include_reads = "all",
                                       only_deletions = FALSE,
                                       k = 10000
                                       ) {

  if (!(include_reads %in% c("all", "fully mapped", "minimally mapped"))) {
    stop("An invalid value was provided for the 'include_reads' parameter!")
  }
  CheckThatIntegerVectorIsInOrder(categor_df[, "Read_number"])
  stopifnot(all(c(sg_A, sg_B) %in% 1:4))

  read_numbers <- unique(categor_df[, "Read_number"])
  matches_vec <- match(read_numbers, mapped_df[, "Read_number"])
  stopifnot(!(anyNA(matches_vec)))

  ## Determine which reads to include (e.g. reads where both sg_A or sg_B are mapped)
  plasmid_columns <- paste0("Plasmid_sg", sg_A:sg_B)
  plasmid_mat <- as.matrix(mapped_df[matches_vec, plasmid_columns])
  are_included <- (!(is.na(plasmid_mat[, 1]))) & (!(is.na(plasmid_mat[, ncol(plasmid_mat)])))
  if (include_reads == "fully mapped") {
    are_included <- are_included & (mapped_df[matches_vec, "Num_matched_sgRNAs"] == 4)
  } else if (include_reads == "minimally mapped") {
    are_included <- are_included & (mapped_df[matches_vec, "Num_matched_sgRNAs"] == 2)
  }
  plasmid_mat <- plasmid_mat[are_included, ]

  ## Identify template switches
  are_switched <- apply(plasmid_mat, 1, function(x) length(unique(x[!(is.na(x))])) > 1)
  message("For this analysis, ",
          sum(are_included), " reads (out of a total of ", length(read_numbers),
          ") were included.\nOf these, ", sum(are_switched), " (",
          format(sum(are_switched) / length(are_switched) * 100, digits = 1, nsmall = 1),
          "%) featured a template switch\nbetween sg", sg_A, " and sg", sg_B,
          ", whereas ", sum(!(are_switched)), " (",
          format(sum(!(are_switched)) / length(are_switched) * 100, digits = 1, nsmall = 1),
          "%) did not.\n"
          )

  ## Produce a logical matrix (is correct? => TRUE/FALSE) for each read and feature
  if (only_deletions) {
    are_incorrect_mat <- CategorDfToMat(categor_df, "Mostly_deleted")
  } else {
    are_incorrect_mat <- !(CategorDfToMat(categor_df, "Is_correct"))
  }
  stopifnot(nrow(are_incorrect_mat) == length(read_numbers))
  are_incorrect_mat <- are_incorrect_mat[are_included, ]

  ## Calculate error rates for switched and non-switched reads
  errors_df <- data.frame(
    "Feature"                   = colnames(are_incorrect_mat),
    "Num_incorrect_nonswitched" = as.integer(colSums(are_incorrect_mat[!(are_switched), ])),
    "Num_incorrect_switched"    = as.integer(colSums(are_incorrect_mat[are_switched, ])),
    row.names = NULL
  )
  errors_df[, "Fraction_incorrect_nonswitched"] <- errors_df[, "Num_incorrect_nonswitched"] / sum(!(are_switched))
  errors_df[, "Fraction_incorrect_switched"] <- errors_df[, "Num_incorrect_switched"] / sum(are_switched)
  errors_df[, "Delta_fraction"] <- errors_df[, "Fraction_incorrect_switched"] - errors_df[, "Fraction_incorrect_nonswitched"]

  ## Bootstrap p values for error rates
  set.seed(1)
  resampled_delta_mat <- t(vapply(seq_len(k), function(x) {
    random_are_switched <- sample(are_switched)
    (colSums(are_incorrect_mat[!(random_are_switched), ]) / sum(!(random_are_switched))) -
    colSums(are_incorrect_mat[random_are_switched, ]) / sum(random_are_switched)
  }, numeric(nrow(errors_df))))

  p_values_one_sided <- vapply(seq_len(nrow(errors_df)), function(x) {
    sum(resampled_delta_mat[, x] >= errors_df[, "Delta_fraction"][[x]])
  }, integer(1))
  p_values_two_sided <- vapply(seq_len(nrow(errors_df)), function(x) {
    sum(abs(resampled_delta_mat[, x]) >= abs(errors_df[, "Delta_fraction"][[x]]))
  }, integer(1))

  errors_df[, "P_one_sided"] <- p_values_one_sided / k
  errors_df[, "P_two_sided"] <- p_values_two_sided / k
  return(errors_df)
}


