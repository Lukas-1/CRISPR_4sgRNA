### 13 July 2020 ###



# Define functions --------------------------------------------------------


TidyBioMartChromosomes <- function(chromosomes_vec) {
  results_vec <- ifelse(chromosomes_vec == "MT", "M", chromosomes_vec)
  are_on_chromosome <- results_vec %in% c(as.character(1:22), "X", "Y", "M")
  results_vec <- paste0(ifelse(are_on_chromosome, "chr",  ""),
                        results_vec
                        )
  return(results_vec)
}



MergeTSSData <- function(FANTOM5_df,
                         BioMart_df,
                         common_columns = c("Group",
                                            "Entrez_ID",
                                            "Gene_symbol",
                                            "Chromosome",
                                            "Strand"
                                            ),
                         use_legacy_TSS_order = TRUE
                         ) {

  are_duplicated_FANTOM5 <- duplicated(FANTOM5_df[, c(common_columns, "TSS_start")])
  are_duplicated_BioMart <- duplicated(BioMart_df[, c(common_columns, "TSS")])

  table(are_duplicated_FANTOM5)
  table(are_duplicated_BioMart)

  combined_TSS_df <- rbind.data.frame(
    data.frame(
      "Source" = "FANTOM5",
      FANTOM5_df[!(are_duplicated_FANTOM5), common_columns],
      "TSS" = FANTOM5_df[["TSS_start"]][!(are_duplicated_FANTOM5)],
      FANTOM5_df[!(are_duplicated_FANTOM5), "Score", drop = FALSE],
      stringsAsFactors = FALSE
    ),
    unique(data.frame(
      "Source" = "BioMart",
      BioMart_df[!(are_duplicated_BioMart), c(common_columns, "TSS")],
      "Score" = NA_integer_,
      stringsAsFactors = FALSE
    ), MARGIN = 1),
    stringsAsFactors = FALSE,
    make.row.names = FALSE
  )

  if (use_legacy_TSS_order) {
    TSS_order_vec <- combined_TSS_df[["TSS"]]
  } else {
    TSS_order_vec <- ifelse(combined_TSS_df[["Strand"]] == "+",
                            combined_TSS_df[["TSS"]],
                            -(combined_TSS_df[["TSS"]])
                            )
  }

  if (all(c("Gene_symbol", "Group") %in% common_columns)) {
    entrez_to_symbols_vec <- MapToEntrezs(entrez_IDs_vec = combined_TSS_df[["Entrez_ID"]])[["Gene_symbol"]]

    combined_TSS_df <- combined_TSS_df[order(GetMinEntrez(combined_TSS_df[["Entrez_ID"]]),
                                             !(mapply(identical, entrez_to_symbols_vec, combined_TSS_df[["Gene_symbol"]])),
                                             combined_TSS_df[["Group"]],
                                             TSS_order_vec
                                             ), ]
  } else {
    combined_TSS_df <- combined_TSS_df[order(GetMinEntrez(combined_TSS_df[["Entrez_ID"]]),
                                             combined_TSS_df[["Entrez_ID"]],
                                             TSS_order_vec
                                             ), ]
  }

  common_columns <- c(common_columns, "TSS")

  are_duplicates_FANTOM5 <- duplicated(combined_TSS_df[, common_columns], fromLast = TRUE)
  are_duplicates_BioMart <- duplicated(combined_TSS_df[, common_columns], fromLast = FALSE)

  stopifnot(all(combined_TSS_df[["Source"]][are_duplicates_FANTOM5] == "FANTOM5"))
  stopifnot(all(combined_TSS_df[["Source"]][are_duplicates_BioMart] == "BioMart"))

  combined_TSS_df[["Source"]][are_duplicates_FANTOM5] <- "FANTOM5, BioMart"

  combined_TSS_df <- combined_TSS_df[!(are_duplicates_BioMart), ]

  row.names(combined_TSS_df) <- NULL
  return(combined_TSS_df)
}


