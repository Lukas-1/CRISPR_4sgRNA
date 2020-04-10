### 9th April 2020 ###



# Import packages and source code -----------------------------------------

library("vioplot")
library("eulerr")
library("gridExtra")
library("ggplot2")
library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R")) # For GetMinEntrez
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))




# Define maps -------------------------------------------------------------

categorical_4sg_columns <- c("Are_overlapping", "Have_homologies")

CFD_specificity_scores <- c("GuideScan_specificity", "CRISPOR_3MM_specificity",
                            "CRISPOR_4MM_specificity"
                            )

range_0_to_1_columns <- CFD_specificity_scores

range_0_to_100_columns <- c("CRISPOR_Doench_efficacy", "GuideScan_efficiency",
                            "CRISPOR_CFD_specificity",  "CRISPOR_MIT_specificity",
                            grep("SNP_AF_(sum|max)_", SNP_column_names, value = TRUE)
                            )

core_numeric_column_labels <- c(
  "all22_SNP_AF_max_Kaviar"   = "Most frequent polymorphism within the sgRNA/PAM region (Kaviar database)",
  "CRISPOR_Doench_efficacy"   = "Efficacy score (Rule Set 2) from CRISPOR",
  "GuideScan_efficiency"      = "Efficacy score (Rule Set 2) from GuideScan",
  "GuideScan_specificity"     = "GuideScan specificity score",
  "CRISPOR_3MM_specificity"   = "CRISPOR specificity score (up to 3 mismatches)",
  "CRISPOR_4MM_specificity"   = "CRISPOR specificity score (up to 4 mismatches)",
  "CRISPOR_CFD_specificity"   = "CRISPOR CFD specificity score",
  "CRISPOR_MIT_specificity"   = "CRISPOR MIT specificity score",
  "Deviation_from_TSS_window" = "Deviation from the optimal window near the TSS (-75 to -150 bp)"
)

numeric_column_labels <- c(
  core_numeric_column_labels,
  "GuideScan_Num_2MM"    = "GuideScan \u2013 number of 2MM sites",
  "GuideScan_Num_3MM"    = "GuideScan \u2013 number of 3MM sites",
  "GuideScan_Num_2or3MM" = "GuideScan \u2013 number of 2MM or 3MM sites",
  "CRISPOR_Num_2MM"      = "CRISPOR \u2013 number of 2MM sites",
  "CRISPOR_Num_3MM"      = "CRISPOR \u2013 number of 3MM sites",
  "CRISPOR_Num_2or3MM"   = "CRISPOR \u2013 number of 2MM or 3MM sites",
  "CRISPOR_Num_4MM"      = "CRISPOR \u2013 number of 4MM sites"
)




specificity_score_cutoffs <- c(
  "GuideScan_specificity"   = 0.2,
  "CRISPOR_3MM_specificity" = 0.2,
  "CRISPOR_4MM_specificity" = 0.04,
  "CRISPOR_CFD_specificity" = 50,
  "CRISPOR_MIT_specificity" = 50
)

categorical_columns <- c(
  "Do_not_meet_criteria",
  "Not_mapped",
  "Poly_T",
  "Non_canonical_PAM",
  "CRISPOR_Graf_status",
  names(core_numeric_column_labels)
)


sublibraries_order <- c(
  "Calabrese_1to3",
  "Calabrese_4to6",
  "hCRISPRa-v2_1to5",
  "hCRISPRa-v2_6to10",
  "Dolcetto_1to3",
  "Dolcetto_4to6",
  "hCRISPRi-v2.1_1to5",
  "hCRISPRi-v2,1_6to10",
  "Brunello",
  "TKOv3",
  "GPP_1to10",
  "GPP_top10",
  "GPP_11to50",
  "GPP_rest",
  "None"
)



sublibraries_labels <- c(
  "Calabrese_1to3"    = "# 1-3",
  "Calabrese_4to6"    = "# 4-6",
  "hCRISPRa-v2_1to5"  = "Top 5",
  "hCRISPRa-v2_6to10" = "Supp 5",
  "Dolcetto_1to3"     = "# 1-3",
  "Dolcetto_4to6"     = "# 4-6",
  "hCRISPRi-v2_1to5"  = "Top 5",
  "hCRISPRi-v2_6to10" = "Supp 5",
  "GPP_1to10"         = "Top 10",
  "GPP_top10"         = "Top 10",
  "GPP_11to50"        = "# 11-50",
  "GPP_rest"          = "Other",
  "None"              = ""
)





args_list <- list(
  "B1) sources - filtered" = list(
    "show_sublibraries"      = FALSE,
    "filter_top4"            = TRUE
  ),
  "B2) sources - filtered (4-guide combination)" = list(
    "aggregate_scores"       = TRUE
  ),
  "B3) sources - unfiltered - only top 4" = list(
    "show_sublibraries"      = FALSE,
    "filter_top4"            = TRUE,
    "filter_complete_genes"  = FALSE,
    "filter_complete_scores" = FALSE
  ),
  "B4) sources - unfiltered (4-guide combination)" = list(
    "aggregate_scores"       = TRUE,
    "show_sublibraries"      = FALSE,
    "filter_complete_genes"  = FALSE,
    "filter_complete_scores" = TRUE
  ),
  "B5) sources - unfiltered - all guides" = list(
    "show_sublibraries"      = FALSE,
    "filter_top4"            = FALSE
  ),
  "C1) sources - rest vs. 4sg - filtered" = list(
    "show_rest_v_4sg"        = TRUE,
    "filter_top4"            = TRUE
  ),
  "C2) sources - rest vs. 4sg - unfiltered - only top 4" = list(
    "show_rest_v_4sg"        = TRUE,
    "filter_top4"            = TRUE,
    "filter_complete_genes"  = FALSE,
    "filter_complete_scores" = FALSE
  ),
  "C3) sources - rest vs. 4sg - unfiltered - all guides" = list(
    "show_rest_v_4sg"        = TRUE,
    "filter_top4"            = FALSE
  )
)


CRISPRa_args_list <- list(
  "D1) sources - sublibraries - filtered" = list(
    "show_sublibraries"      = TRUE,
    "filter_top4"            = TRUE
  ),
  "D2) sources - sublibraries - unfiltered - only top 4" = list(
    "show_sublibraries"      = TRUE,
    "filter_top4"            = TRUE,
    "filter_complete_genes"  = FALSE,
    "filter_complete_scores" = FALSE
  ),
  "D3) sources - sublibraries - unfiltered - all guides" = list(
    "show_sublibraries"      = TRUE,
    "filter_top4"            = FALSE
  ),
  "D4) sources - sublibraries (GPP collapsed) - unfiltered" = list(
    "show_sublibraries"      = TRUE,
    "filter_top4"            = FALSE,
    "collapse_GPP"           = TRUE
  )
)


CRISPRko_args_list <- list(
  "D) sources - sublibraries (unfiltered)" = list(
    "show_sublibraries" = TRUE,
    "filter_top4"       = FALSE
  )
)







# Settings for generating multi-plot layouts ------------------------------

use_layout_mat <- rbind(c(3, 3), c(1, 2), c(4, 4))
layout_widths <- c(2, 6)
layout_heights <- c(1, 19, 0.3)
pdf_width <- 8
pdf_height <- 4







# General helper functions ------------------------------------------------

FilterCRISPRDf <- function(CRISPR_df) {
  CRISPR_df[(CRISPR_df[["Is_control"]] == "No") & (CRISPR_df[["Source"]] != "Curated"), ]
}

palify_cache_101 <- list()
Palify <- function(myhex, fraction_pale = 0.5) {
  if (myhex %in% names(palify_cache_101)) {
    color_vec <- palify_cache_101[[myhex]]
  } else {
    color_vec <- colorRampPalette(c(myhex, "#FFFFFF"))(101)
    palify_cache_101[[myhex]] <- color_vec
    assign("palify_cache_101", palify_cache_101, envir = globalenv())
  }
  color_vec[[round(fraction_pale * 100) + 1]]
}






# Helper functions for filtering CRISPRa sgRNAs ---------------------------

# show_columns <- c(
#   "Entrez_ID", "Gene_symbol", "hCRISPRa_v2_transcript", "hCRISPRa_v2_rank", "Rank",
#   "Num_TSSs", "TSS_number", "Allocated_TSS", "TSS_ID", "AltTSS_ID",
#   "Entrez_chromosome", "Chromosome", "Start", "Best_TSS"
# )



GetMainTSS <- function(CRISPRa_df) {

  are_not_controls <- CRISPRa_df[["Is_control"]] == "No"
  altTSS_ID_fac <- factor(CRISPRa_df[["AltTSS_ID"]][are_not_controls],
                          levels = unique(CRISPRa_df[["AltTSS_ID"]][are_not_controls])
                          )
  CheckThatFactorIsInOrder(altTSS_ID_fac)

  combined_ID_vec <- tapply(CRISPRa_df[["Combined_ID"]][are_not_controls], altTSS_ID_fac, unique)

  min_distance_from_TSS_vec <- tapply(CRISPRa_df[["Distance_from_TSS"]][are_not_controls],
                                      altTSS_ID_fac,
                                      function(x) if (all(is.na(x))) NA_integer_ else min(abs(x), na.rm = TRUE)
                                      )

  combined_ID_fac <- factor(combined_ID_vec, levels = unique(combined_ID_vec))
  CheckThatFactorIsInOrder(combined_ID_fac)


  altTSS_ID_table <- table(altTSS_ID_fac)
  best_TSS_ID_vec <- tapply(seq_along(combined_ID_fac),
                            combined_ID_fac,
                            function(x) {
                              altTSS_IDs <- names(min_distance_from_TSS_vec)[x]
                              if (length(x) == 1) {
                                best_TSS <- names(min_distance_from_TSS_vec)[x]
                              } else {
                                distances <- min_distance_from_TSS_vec[x]
                                if (all(is.na(distances))) {
                                  common_ID_frequencies <- altTSS_ID_table[altTSS_IDs]
                                  best_TSS <- altTSS_IDs[common_ID_frequencies == max(common_ID_frequencies)][[1]]
                                } else {
                                  best_TSS <- altTSS_IDs[which.min(distances)]
                                }
                                return(best_TSS)
                              }
                            })

  best_TSS_expanded_vec <- best_TSS_ID_vec[match(CRISPRa_df[["Combined_ID"]][are_not_controls], names(best_TSS_ID_vec))]

  best_TSS_including_controls_vec <- rep(NA_character_, nrow(CRISPRa_df))
  best_TSS_including_controls_vec[are_not_controls] <- best_TSS_expanded_vec

  return(best_TSS_including_controls_vec)
}








# Helper functions for vertically aligning text on plots ------------------

StripExpression <- function(my_expression) {
  if (is.character(my_expression)) {
    literal_string <- paste0("\"", capture.output(cat(my_expression)), "\"")
  } else {
    literal_string <- capture.output(my_expression)
    if (substr(literal_string, 1, 11) == "expression(") {
      literal_string <- substr(literal_string, 12, nchar(literal_string) - 1)
    }
  }
  return(literal_string)
}

ConcatenateExpressions <- function(expression_list, my_sep = "  \u2013  ") {
  literal_strings <- sapply(expression_list, StripExpression)
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}

VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("g")), use_expression, expression(phantom("h")))
  return(ConcatenateExpressions(my_list, my_sep = ""))
}

Embolden <- function(my_string) {
  parse(text = paste0("bold(", StripExpression(my_string), ")"))
}

AdjustTextVec <- function(char_vec) {
  sapply(char_vec, function(x) Embolden(VerticalAdjust(x)))
}





# Helper functions for plotting -------------------------------------------

SubgroupColorsDf <- function() {
  dark_colors <- vapply(c("Blues", "Greens", "Reds", "Purples", "Greys"),
                        function(x) colorRampPalette(brewer.pal(9, x))(18)[[17]],
                        ""
                        )
  names(dark_colors) <- substr(names(dark_colors), 1, nchar(names(dark_colors)) - 1)
  dark_colors <- c(
    dark_colors[1:3],
    "BlueGreen" = "#006969",
    dark_colors["Purple"],
    "Brown" = "#581200",
    dark_colors["Grey"]
  )
  colors_df <- data.frame(
    "Color_name" = names(dark_colors),
    "Lightened"  = FALSE,
    "Pale"       = vapply(dark_colors, Palify, fraction_pale = 0.9, ""),
    "Medium"     = vapply(dark_colors, Palify, fraction_pale = 0.775, ""),
    "Dark"       = dark_colors,
    stringsAsFactors = FALSE
  )
  colors_df <- colors_df[rep(seq_len(nrow(colors_df)), each = 2), ]
  rownames(colors_df) <- NULL
  colors_df[["Lightened"]] <- rep_len(c(TRUE, FALSE), length.out = nrow(colors_df))
  for (column_name in c("Pale", "Medium", "Dark")) {
    colors_df[[column_name]][colors_df[["Lightened"]]] <- vapply(colors_df[[column_name]][colors_df[["Lightened"]]], Palify, fraction_pale = 0.5, "")
  }
  return(colors_df)
}



FilterOriginalColorsDf <- function(colors_df, plot_df, show_rest_v_4sg = FALSE, omit_labels_for_single_subgroup = TRUE) {

  colors_df <- colors_df[colors_df[["Color_name"]] %in% c("Blue", "Green", "Red", "Purple"), ]
  groups_subgroups_vec <- paste0(plot_df[["Group"]], "__", plot_df[["Subgroup"]])
  unique_df <- plot_df[!(duplicated(groups_subgroups_vec)), c("Group", "Subgroup")]

  if (show_rest_v_4sg) {
    colors_df <- colors_df[-(nrow(colors_df) - 1), ]
  } else {
    new_order <- order(match(colors_df[["Color_name"]], colors_df[["Color_name"]]),
                       colors_df[["Lightened"]]
                       )
    colors_df <- colors_df[new_order, ]
    unique_groups_table <- table(factor(unique_df[["Group"]], levels = unique(unique_df[["Group"]])))
    are_to_keep <- unlist(lapply(unique_groups_table, function(x) {
      if (x == 1) {
        c(TRUE, FALSE)
      } else if (x == 2) {
        c(TRUE, TRUE)
      } else {
        stop("Unexpected number of subgroups!")
      }
    }), use.names = FALSE)
    colors_df <- colors_df[are_to_keep, ]
  }

  colors_df[["Group_label"]] <- as.character(unique_df[["Group"]])
  if (show_rest_v_4sg) {
    colors_df[["Group_label"]] <- ifelse(colors_df[["Group_label"]] == "4sg", "All 4sg", colors_df[["Group_label"]])
    colors_df[["Subgroup_label"]] <- as.character(unique_df[["Subgroup"]])
  } else {
    colors_df[["Subgroup_label"]] <- sublibraries_labels[as.character(unique_df[["Subgroup"]])]
  }

  if (omit_labels_for_single_subgroup) {
    num_occurrences_vec <- table(colors_df[["Color_name"]])[colors_df[["Color_name"]]]
    colors_df[["Subgroup_label"]] <- ifelse(num_occurrences_vec > 1, colors_df[["Subgroup_label"]], "")
  }

  are_new_color <- rep.int(NA, nrow(colors_df))
  last_color <- NA
  for (i in seq_len(nrow(colors_df))) {
    this_color <- colors_df[["Color_name"]][[i]]
    are_new_color[[i]] <- !(identical(this_color, last_color))
    last_color <- this_color
  }
  colors_df[["Are_new_color"]] <- are_new_color

  rownames(colors_df) <- NULL
  return(colors_df)
}






ExpandedSubgroupsDf <- function(CRISPR_df, collapse_GPP = FALSE, show_rest_v_4sg = FALSE) {

  sources_vec <- as.character(ReformatSourceToFactor(CRISPR_df[["Source"]]))
  sources_splits <- strsplit(sources_vec, ", ", fixed = TRUE)

  CRISPR_df_expanded <- CRISPR_df[rep(seq_len(nrow(CRISPR_df)), lengths(sources_splits)), ]

  sources_vec <- unlist(sources_splits, use.names = FALSE)
  sources_fac <- factor(sources_vec, levels = intersect(libraries_order, sources_vec))

  if (show_rest_v_4sg) {
    subgroups_fac <- factor(ifelse(CRISPR_df_expanded[["Rank"]] %in% 1:4, "4sg", "Rest"), levels = c("Rest", "4sg"))
  } else {
    if ("hCRISPRa-v2" %in% sources_vec) {
      sources_sublibraries <- ifelse(sources_vec == "hCRISPRa-v2",
                                     paste0("hCRISPRa-v2_", ifelse(CRISPR_df_expanded[["hCRISPRa_v2_rank"]] %in% 1:5, "1to5", "6to10")),
                                     ifelse(sources_vec == "Calabrese",
                                            paste0("Calabrese_", ifelse(CRISPR_df_expanded[["Calabrese_rank"]] %in% "1/2/3", "1to3", "4to6")),
                                            ifelse(sources_vec == "GPP",
                                                   paste0("GPP_", ifelse(CRISPR_df_expanded[["GPP_rank"]] %in% 1:10, "top10", "rest")),
                                                   sources_vec
                                                   )
                                            )
                                     )
    } else {
      sources_sublibraries <- ifelse(sources_vec == "GPP",
                                     paste0("GPP_", ifelse(CRISPR_df_expanded[["GPP_rank"]] %in% 1:10, "1to10", "11to50")),
                                     "None"
                                     )

    }
    if (collapse_GPP) {
      sources_sublibraries <- ifelse(sources_vec == "GPP", "None", sources_sublibraries)
    }
    subgroups_fac <- factor(sources_sublibraries, levels = intersect(sublibraries_order, sources_sublibraries))
  }


  plot_df <- data.frame(
    CRISPR_df_expanded,
    "Group"    = sources_fac,
    "Subgroup" = subgroups_fac,
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[order(plot_df[["Group"]], plot_df[["Subgroup"]]), ]
  rownames(plot_df) <- NULL

  are_top4 <- CRISPR_df[["Rank"]] %in% 1:4
  top4_df <- data.frame(
    CRISPR_df[are_top4, ],
    "Group"    = "4sg",
    "Subgroup" = "None",
    stringsAsFactors = FALSE
  )
  results_df <- rbind.data.frame(plot_df, top4_df, make.row.names = FALSE, stringsAsFactors = FALSE)
  return(results_df)
}



FilterCompleteTop4 <- function(CRISPR_df) {
  num_occurrences_vec <- table(CRISPR_df[["Combined_ID"]])[CRISPR_df[["Combined_ID"]]]
  results_df <- CRISPR_df[num_occurrences_vec >= 4, ]
  return(results_df)
}



FilterTop4 <- function(expanded_CRISPR_df,
                       filter_complete_genes  = TRUE,
                       filter_complete_scores = FALSE,
                       data_column            = NULL,
                       show_sublibraries      = FALSE
                       ) {

  assign("delete_expanded_CRISPR_df",     expanded_CRISPR_df,     envir = globalenv())
  assign("delete_filter_complete_genes",  filter_complete_genes,  envir = globalenv())
  assign("delete_filter_complete_scores", filter_complete_scores, envir = globalenv())
  assign("delete_data_column",            data_column,            envir = globalenv())

  if (filter_complete_scores && is.null(data_column)) {
    stop("If 'filter_complete_scores' is TRUE, the 'data_column' argument must be specified!")
  }

  ## Distinguish between CRISPRa and CRISPRko
  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(expanded_CRISPR_df)
  if (is_CRISPRko && !("Entrez_source_Brunello" %in% colnames(expanded_CRISPR_df))) {
    stop("No identifying column names were found in expanded_CRISPR_df!")
  }

  ## Filter for valid Entrez IDs
  have_entrez <- (expanded_CRISPR_df[["Entrez_ID"]] == expanded_CRISPR_df[["Combined_ID"]]) %in% TRUE
  are_protein_coding <- expanded_CRISPR_df[["Entrez_ID"]] %in% collected_entrez_IDs
  are_eligible <- have_entrez & are_protein_coding

  if (!(is_CRISPRko) && show_sublibraries) {
    are_eligible <- are_eligible & !(expanded_CRISPR_df[["Subgroup"]] == "hCRISPRa-v2_6to10") # Otherwise, a tiny group of "supp 5" genes is displayed in the plot
  }
  eligible_entrezs_vec <- expanded_CRISPR_df[["Entrez_ID"]][are_eligible]


  ## Filter for Entrez IDs with at least 4 guides in each of the libraries
  if (filter_complete_genes) {
    library_counts_mat <- as.matrix(table(factor(eligible_entrezs_vec, levels = unique(eligible_entrezs_vec)),
                                          expanded_CRISPR_df[["Group"]][are_eligible]
                                          )
                                    )
    are_complete_entrezs <- rowSums(library_counts_mat >= 4) == 4
    are_complete <- expanded_CRISPR_df[["Entrez_ID"]] %in% rownames(library_counts_mat)[are_complete_entrezs]
    are_eligible <- are_eligible & are_complete
  }
  filtered_df <- expanded_CRISPR_df[are_eligible, ]


  ## Choose the top 4 guides from GPP
  GPP_df <- filtered_df[filtered_df[["Group"]] == "GPP", ]
  GPP_df <- GPP_df[order(as.integer(GPP_df[["Entrez_ID"]])), ]

  GPP_entrezs_fac <- factor(GPP_df[["Entrez_ID"]], levels = unique(GPP_df[["Entrez_ID"]]))

  CheckThatFactorIsInOrder(GPP_entrezs_fac)

  GPP_rank_list <- split(GPP_df[["GPP_rank"]], GPP_entrezs_fac)
  GPP_are_chosen_list <- lapply(GPP_rank_list, function(x) rank(x, ties.method = "first") %in% 1:4)
  GPP_are_chosen <- unlist(GPP_are_chosen_list, use.names = FALSE)
  GPP_chosen_df <- GPP_df[GPP_are_chosen, ]
  if (filter_complete_genes) {
    stopifnot(all(table(GPP_chosen_df[["Combined_ID"]]) == 4))
  }


  ## Select the 4sg guides
  FourSg_df <- filtered_df[filtered_df[["Group"]] == "4sg", ]
  assign("delete_FourSg_df", FourSg_df, envir = globalenv())
  if (!(filter_complete_genes)) {
    FourSg_df <- FilterCompleteTop4(FourSg_df)
  }
  stopifnot(all(table(FourSg_df[["Combined_ID"]]) == 4))


  if (!(is_CRISPRko)) {

    ## Choose the guides for Calabrese
    Doench_df <- filtered_df[filtered_df[["Group"]] %in% c("Dolcetto", "Calabrese"), ]
    Doench_df <- Doench_df[order(as.integer(Doench_df[["Entrez_ID"]])), ]

    if (!(filter_complete_genes)) {
      Doench_df <- FilterCompleteTop4(Doench_df)
    }

    Doench_entrezs_fac <- factor(Doench_df[["Entrez_ID"]], levels = unique(Doench_df[["Entrez_ID"]]))

    CheckThatFactorIsInOrder(Doench_entrezs_fac)

    rank_column <- grep("^(Calabrese|Dolcetto)_rank", colnames(Doench_df))
    Doench_rank_list <- split(Doench_df[[rank_column]], Doench_entrezs_fac)

    set.seed(1)
    Doench_are_chosen_list <- lapply(Doench_rank_list, function(x) {
      are_SetA <- x == "1/2/3"
      num_SetA <- sum(are_SetA)
      if (num_SetA > 4) {
        choose_indices <- sample(which(are_SetA), size = 4)
        are_chosen <- seq_along(x) %in% choose_indices
      } else {
        num_needed <- 4 - num_SetA
        choose_indices <- which(!(are_SetA))[sample.int(n = num_needed, size = num_needed)]
        are_chosen <- are_SetA | (seq_along(x) %in% choose_indices)
      }
      stopifnot(sum(are_chosen) == 4)
      return(are_chosen)
    })
    Doench_are_chosen <- unlist(Doench_are_chosen_list, use.names = FALSE)
    Doench_chosen_df <- Doench_chosen_df[Doench_are_chosen, ]


    ## Choose the guides for hCRISPRa-v2
    hCRISPR_v2_df <- filtered_df[filtered_df[["Group"]] %in% c("hCRISPRa-v2", "hCRISPRa-v2.1", "hCRISPRa-v2.0"), ]

    hCRISPR_v2_df <- hCRISPR_v2_df[order(match(hCRISPR_v2_df[["Combined_ID"]], hCRISPR_v2_df[["Combined_ID"]])), ]
    hCRISPR_v2_combined_IDs <- factor(hCRISPR_v2_df[["Combined_ID"]], levels = unique(hCRISPR_v2_df[["Combined_ID"]]))
    rank_column <- grep("_v2_rank", colnames(hCRISPR_v2_df), fixed = TRUE)
    hCRISPR_v2_rank_list <- split(as.integer(hCRISPR_v2_df[[rank_column]]), hCRISPR_v2_combined_IDs)
    are_top4 <- unlist(lapply(hCRISPR_v2_rank_list, function(x) {
      top_4_ranks <- sort(unique(x))[1:4]
      x %in% top_4_ranks
    }))
    hCRISPRa_v2_chosen_df <- hCRISPR_v2_df[are_top4, ]
    hCRISPRa_v2_chosen_df <- ChooseOriginalTop4(hCRISPRa_v2_chosen_df)


    ## Combine the final data frame for CRISPRa
    if (filter_complete_genes) {
      stopifnot(length(unique(c(nrow(Calabrese_chosen_df), nrow(hCRISPRa_v2_chosen_df), nrow(GPP_chosen_df), nrow(FourSg_df)))) == 1)
    }

    combined_chosen_df <- rbind.data.frame(
      Calabrese_chosen_df,
      hCRISPRa_v2_chosen_df,
      GPP_chosen_df,
      FourSg_df,
      stringsAsFactors = FALSE,
      make.row.names = FALSE
    )

  } else {
    Brunello_chosen_df <- ChooseOriginalTop4(filtered_df[filtered_df[["Group"]] == "Brunello", ])
    TKOv3_chosen_df <- ChooseOriginalTop4(filtered_df[filtered_df[["Group"]] == "TKOv3", ])
    if (filter_complete_genes) {
      stopifnot(length(unique(c(nrow(Brunello_chosen_df), nrow(TKOv3_chosen_df), nrow(GPP_chosen_df), nrow(FourSg_df)))) == 1)
    }
    combined_chosen_df <- rbind.data.frame(
      Brunello_chosen_df,
      TKOv3_chosen_df,
      GPP_chosen_df,
      FourSg_df,
      stringsAsFactors = FALSE,
      make.row.names = FALSE
    )
  }

  if (filter_complete_scores) {

    tabulate_fac <- factor(combined_chosen_df[["Entrez_ID"]], levels = unique(combined_chosen_df[["Entrez_ID"]]))
    if (!(filter_complete_genes)) {
      tabulate_fac <- droplevels(interaction(tabulate_fac, combined_chosen_df[["Group"]], lex.order = TRUE))
    }
    assign("delete_combined_chosen_df", combined_chosen_df, envir = globalenv())
    assign("delete_tabulate_fac",       tabulate_fac, envir = globalenv())
    assign("delete_data_column",       data_column, envir = globalenv())
    if (data_column == "Have_homologies") {
      num_present_table <- table(tabulate_fac)
    } else {
      if (data_column == "Are_overlapping") {
        check_column <- "Cut_location"
      } else {
        check_column <- data_column
      }
      num_present_table <- table(tabulate_fac[!(is.na(combined_chosen_df[[check_column]]))])
    }
    num_occurrences_vec <- num_present_table[as.character(tabulate_fac)]
    are_complete_entries <- num_occurrences_vec == max(num_present_table)
    combined_chosen_df <- combined_chosen_df[are_complete_entries, ]
  }

  ## Re-order the final data frame
  new_order <- order(combined_chosen_df[["Group"]], combined_chosen_df[["Subgroup"]])
  combined_chosen_df <- combined_chosen_df[new_order, ]
  rownames(combined_chosen_df) <- NULL
  return(combined_chosen_df)
}



ChooseOriginalTop4 <- function(CRISPR_df) {

  if ("Original_order" %in% colnames(CRISPR_df)) {
    new_order <- order(match(CRISPR_df[["Entrez_ID"]], CRISPR_df[["Entrez_ID"]]),
                       CRISPR_df[["Original_order"]]
                       )
    CRISPR_df <- CRISPR_df[new_order, ]
  }

  num_occurrences_vec <- table(CRISPR_df[["Entrez_ID"]])[CRISPR_df[["Entrez_ID"]]]
  fourguides_df <- CRISPR_df[num_occurrences_vec == 4, ]
  other_df <- CRISPR_df[num_occurrences_vec > 4, ]

  # print(head(other_df[, c("Entrez_ID", "Gene_symbol", "Original_symbol", "Source", "sgRNA_sequence")], 10))

  other_order <- order(as.integer(other_df[["Entrez_ID"]]),
                       other_df[["Original_symbol"]] != ""
                       )
  other_df <- other_df[other_order, ]
  other_df_splits <- split(other_df, factor(other_df[["Entrez_ID"]], levels = unique(other_df[["Entrez_ID"]])))
  other_df_list <- lapply(other_df_splits, function(x) x[1:4, ])
  other_df <- do.call(rbind.data.frame, c(other_df_list, stringsAsFactors = FALSE, make.row.names = FALSE))

  chosen_df <- rbind.data.frame(fourguides_df, other_df, stringsAsFactors = FALSE, make.row.names = FALSE)

  stopifnot(all(table(chosen_df[["Entrez_ID"]]) == 4))
  return(chosen_df)
}






MutuallyExclusiveSubgroupsDf <- function(CRISPR_df, use_column) {

  ## Rename and re-order the groups
  groups_fac <- ReformatSourceToFactor(CRISPR_df[["Source"]])

  ## Split into subgroups
  chosen_fac <- factor(ifelse(CRISPR_df[["Rank"]] %in% 1:4, "4sg", "Rest"), levels = c("Rest", "4sg"))
  subgroups_fac <- interaction(groups_fac, chosen_fac, sep = "_", lex.order = TRUE)

  ## Build a data frame for plotting
  plot_df <- data.frame(
    "Data"     = CRISPR_df[[use_column]],
    "Group"    = groups_fac,
    "Subgroup" = subgroups_fac,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[order(plot_df[["Subgroup"]]), ]
  rownames(plot_df) <- NULL
  return(plot_df)
}



SourcesSubgroupLabels <- function(plot_df, x_positions, colors_df, show_sublibraries = TRUE, show_rest_v_4sg = FALSE, extra_space = FALSE) {

  ## Draw the group labels
  if (extra_space) {
    first_line_factor  <- 0.105
  } else {
    first_line_factor  <- 0.045
  }

  if (show_sublibraries || show_rest_v_4sg) {
    if (extra_space) {
      second_line_factor <- 0.18
    } else {
      second_line_factor <- 0.105
    }
  } else {
    if (extra_space) {
      second_line_factor <- 0.12
    } else {
      second_line_factor <- 0.06
    }
  }

  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  if (show_sublibraries || show_rest_v_4sg) {
    ## Draw the subgroup labels
    if (show_sublibraries) {
      text_colors <- colors_df[["Dark"]]
    } else {
      text_colors <- ifelse(colors_df[["Subgroup_label"]] == "Rest",
                            Palify("#000000", fraction_pale = 0.5),
                            "#000000"
                            )
    }
    text(x      = x_positions,
         y      = par("usr")[[3]] - (plot_height * first_line_factor),
         labels = AdjustTextVec(colors_df[["Subgroup_label"]]),
         adj    = c(0.5, 0.5),
         cex    = 0.7,
         col    = text_colors,
         font   = 2,
         xpd    = NA
         )
  }

  ## Draw the group labels
  have_subgroups <- colors_df[["Subgroup_label"]] != ""
  group_colors_fac <- factor(colors_df[["Color_name"]], levels = unique(colors_df[["Color_name"]]))
  if (show_sublibraries) {
    line_factor_vec <- ifelse(tapply(have_subgroups, group_colors_fac, any),
                              second_line_factor,
                              (first_line_factor + second_line_factor) / 2
                              )
  } else {
    line_factor_vec <- second_line_factor
  }

  text(x      = tapply(x_positions, group_colors_fac, mean),
       y      = par("usr")[[3]] - (plot_height * line_factor_vec),
       labels = AdjustTextVec(tapply(as.character(colors_df[["Group_label"]]), group_colors_fac, unique)),
       adj    = c(0.5, 0.5),
       cex    = 0.8,
       col    = vapply(split(colors_df, group_colors_fac),
                       function(x) x[["Dark"]][order(x[["Lightened"]])][[1]], #colorRampPalette(x)(3)[[2]]),
                       ""
                       ),
       font   = 2,
       xpd    = NA
       )

  if (show_sublibraries && !(show_rest_v_4sg)) {
    ## Draw line segments
    intermediate_line_factor <- mean(c(rep.int(first_line_factor, 5), rep.int(second_line_factor, 4)))
    for (group in unique(colors_df[["Group_label"]][have_subgroups])) {
      are_this_group <- colors_df[["Group_label"]] == group
      segments(x0  = min(x_positions[are_this_group]),
               x1  = max(x_positions[are_this_group]),
               y0  = par("usr")[[3]] - (plot_height * intermediate_line_factor),
               col = colors_df[["Medium"]][are_this_group & !(colors_df[["Lightened"]])],
               xpd = NA
               )
    }
  }
  return(invisible(NULL))
}



UniqueSequencesSubgroupLabels <- function(plot_df, x_positions, colors_df, extra_space = FALSE, only_chosen = FALSE) {
  if (extra_space) {
    first_line_factor  <- 0.065
    second_line_factor <- 0.13
  } else {
    first_line_factor  <- 0.035
    second_line_factor <- 0.1
  }
  if (only_chosen) {
    second_line_factor <- second_line_factor <- 0.04
  }
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  num_subgroups <- nrow(colors_df)

  ## Draw the small subgroup labels (Rest/4sg)
  if (!(only_chosen)) {
    text(x      = x_positions,
         y      = par("usr")[[3]] - (plot_height * first_line_factor),
         labels = c("Rest", "4sg"),
         adj    = c(0.5, 1),
         cex    = 0.75,
         col    = rep_len(c(Palify("#000000", fraction_pale = 0.5), "#000000"),
                          length.out = num_subgroups
                          ),
         font   = 2,
         xpd    = NA
         )
  }

  ## Draw the group labels
  group_labels <- sub(", ", ",\n", levels(plot_df[["Group"]]))
  group_labels[[length(group_labels)]] <- "All 3\nsources"
  old_lheight <- par("lheight" = 1.15)
  colors_vec <- colors_df[["Dark"]]
  if (!(only_chosen)) {
    colors_vec <- colors_vec[!(colors_df[["Lightened"]])]
    x_positions <- tapply(x_positions, rep(seq_len(num_subgroups / 2), each = 2), mean)
  }
  text(x      = x_positions,
       y      = par("usr")[[3]] - (plot_height * second_line_factor),
       labels = group_labels,
       adj    = c(0.5, 1),
       cex    = 0.8,
       col    = colors_vec,
       font   = 2,
       xpd    = NA
       )
  par(old_lheight)
  return(invisible(NULL))
}





# Functions for plotting numerical data -----------------------------------

FixNumericColumns <- function(CRISPR_df, y_column) {
  if (grepl("SNP", y_column, fixed = TRUE)) {
    CRISPR_df[[y_column]] <- ifelse(is.na(CRISPR_df[[y_column]]) & !(is.na(CRISPR_df[["Start"]])),
                                    0,
                                    CRISPR_df[[y_column]]
                                    ) * 100
  } else if (y_column == "Deviation_from_TSS_window") {
    distance_vec <- CRISPR_df[["Distance_from_TSS"]]
    are_too_near <- distance_vec > -75
    are_too_far <- distance_vec < -150
    deviation_vec <- ifelse(!(are_too_near) & !(are_too_far),
                            0L,
                            ifelse(are_too_near,
                                   abs(-75L - distance_vec),
                                   abs(-150L - distance_vec)
                                   )
                            )
    deviation_vec <- ifelse(deviation_vec > 1000, 1000L, deviation_vec)
    CRISPR_df[["Deviation_from_TSS_window"]] <- deviation_vec
  }
  return(CRISPR_df)
}



GetAxisLimits <- function(numeric_vec, column_name = NULL, provide_other_limits = FALSE, extend_range_fraction = 0.02) {
  max_value <- max(numeric_vec, na.rm = TRUE)
  min_value <- min(numeric_vec, na.rm = TRUE)
  span <- max_value - min_value
  variable_limits <- c(min_value - (span * 0.02), max_value + (span * 0.02))
  if (!(is.null(column_name)) && (column_name %in% range_0_to_1_columns)) {
    limits <- c(-0.02, 1.02)
    if (min_value > 0.02) {
      limits[[1]] <- 0
    }
    if (max_value < 0.98) {
      if ((max_value < 0.3) && (column_name == "CRISPOR_4MM_specificity")) {
        limits <- variable_limits
      } else {
        limits[[2]] <- 1
      }
    }
  } else if (!(is.null(column_name)) && (column_name %in% range_0_to_100_columns)) {
    limits <- c(-2, 102)
    if (min_value > 2) {
      limits[[1]] <- 0
    }
    if (max_value < 98) {
      limits[[2]] <- 100
    }
  } else {
    if (provide_other_limits) {
      limits <- variable_limits
    } else {
      return(NULL)
    }
  }
  return(limits)
}


DrawViolinGridAndAxes <- function(y_column,
                                  show_title       = TRUE,
                                  aggregate_scores = aggregate_scores,
                                  title_cex        = 0.9,
                                  title_line       = 1.3
                                  ) {
  tick_locations <- axTicks(2)
  tick_distance <- tick_locations[[2]] - tick_locations[[1]]
  grid_locations <- seq(from = tick_locations[[1]] - (tick_distance / 2),
                        to   = tick_locations[[length(tick_locations)]] + (tick_distance / 2),
                        by   = tick_distance
                        )
  grid_locations <- grid_locations[(grid_locations > par("usr")[[3]]) & (grid_locations < par("usr")[[4]])]

  abline(h = grid_locations, col = "gray95", lwd = 0.5)

  ticks_for_grid <- tick_locations
  if (y_column %in% c("CRISPOR_3MM_specificity", "GuideScan_specificity")) {
    are_0point2 <- tick_locations == 0.2
    if (any(are_0point2)) {
      abline(h = 0.2, col = "gray50", lwd = 0.5)
      ticks_for_grid <- ticks_for_grid[!(are_0point2)]
    }
  }
  abline(h = ticks_for_grid, col = "gray90", lwd = 0.5)

  tick_labels <- format(tick_locations)
  if ((y_column == "Deviation_from_TSS_window") && (tick_labels[[length(tick_labels)]] == "1000")) {
    tick_labels <- sapply(tick_labels, as.expression)
    tick_labels[[length(tick_labels)]] <- bquote("" >= 1000)
  }
  if (grepl("SNP", y_column, fixed = TRUE)) {
    tick_labels <- paste0(tick_labels, "%")
  }
  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.3,
       cex.axis = 0.8,
       lwd      = 0.75
       )
  title_text <- numeric_column_labels[[y_column]]
  if (aggregate_scores) {
    title_text <- paste0("Aggregate ", title_text)
  }
  if (show_title) {
    title(title_text, cex.main = title_cex, line = 1.3)
  }
  return(invisible(title_text))
}




PlotViolin <- function(plot_df,
                       x_positions,
                       x_limits,
                       colors_df,
                       y_column_name,
                       show_title         = TRUE,
                       large_count_labels = FALSE,
                       aggregate_scores   = FALSE,
                       title_cex          = 0.9,
                       title_line         = 1.3
                       ) {

  assign("delete_plot_df",       plot_df,       envir = globalenv())
  assign("delete_x_positions",   x_positions,   envir = globalenv())
  assign("delete_colors_df",     colors_df,     envir = globalenv())
  assign("delete_y_column_name", y_column_name, envir = globalenv())

  stopifnot(all(c("Numeric_data", "Groups_factor", "Point_colors") %in% colnames(plot_df)))
  stopifnot(all(c("Pale", "Medium", "Dark") %in% colnames(colors_df)))

  y_limits <- GetAxisLimits(plot_df[["Numeric_data"]], y_column_name, provide_other_limits = TRUE)
  plot(1,
       xlim = x_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  title_text <- DrawViolinGridAndAxes(y_column_name,
                                      show_title       = show_title,
                                      aggregate_scores = aggregate_scores,
                                      title_cex        = title_cex,
                                      title_line       = title_line
                                      )

  jittered_vec <- x_positions[as.integer(plot_df[["Groups_factor"]])] +
                  rnorm(n = nrow(plot_df), mean = 0, sd = 0.03)
  points_alpha <- 0.1
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

  ## Draw the violin plots
  vioplot(plot_df[["Numeric_data"]] ~ plot_df[["Groups_factor"]],
          add      = TRUE,
          at       = x_positions,
          pchMed   = NA,
          drawRect = FALSE,
          col      = colors_df[["Medium"]],
          border   = NA,
          wex      = 0.75,
          axes     = FALSE
          )

  ## Draw the jittered points
  points(x   = jittered_vec,
         y   = plot_df[["Numeric_data"]],
         cex = 0.5,
         col = paste0(plot_df[["Point_colors"]], alpha_hex),
         pch = 16
         )

  ## Draw the superimposed boxplots
  boxplot(plot_df[["Numeric_data"]] ~ plot_df[["Groups_factor"]],
          add       = TRUE,
          at        = x_positions,
          cex       = 0.2,
          boxwex    = 0.275,
          outline   = FALSE,
          names     = rep.int("", length(x_positions)),
          whisklty  = "blank",
          staplewex = 0,
          whisklwd  = 0,
          staplelty = 0,
          col       = colors_df[["Pale"]],
          border    = colors_df[["Dark"]],
          axes      = FALSE,
          lwd       = 0.75
          )

  box(xpd = NA, lwd = 0.75)

  ## Draw the subgroup numbers
  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  if (large_count_labels) {
    count_plot_fraction <- 0.0275
    count_cex <- 0.7
  } else {
    count_plot_fraction <- 0.0125
    count_cex <- 0.5
  }

  text(x      = x_positions,
       y      = par("usr")[[3]] - (plot_height * count_plot_fraction),
       labels = table(plot_df[["Groups_factor"]][!(is.na(plot_df[["Numeric_data"]]))]),
       adj    = c(0.5, 1),
       cex    = count_cex,
       col    = "gray40",
       xpd    = NA
       )

  return(invisible(title_text))
}



ViolinBox_UniqueTwoGroups <- function(CRISPR_df, y_column, show_title = TRUE) {

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), y_column)

  ## Prepare the group colors
  colors_df <- SubgroupColorsDf()
  colors_df <- colors_df[colors_df[["Color_name"]] == "Grey", ]

  ## Build a data frame for plotting
  plot_df <- data.frame(
    "Numeric_data" = CRISPR_df[[y_column]],
    "Chosen"       = CRISPR_df[["Rank"]] %in% 1:4,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[order(plot_df[["Chosen"]]), ]
  rownames(plot_df) <- NULL
  plot_df[["Point_colors"]] <- rep.int(colors_df[["Dark"]], as.integer(table(plot_df[["Chosen"]])))

  plot_df[["Groups_factor"]] <- factor(plot_df[["Chosen"]])


  ## Daw the violins
  x_positions <- 1:2
  num_groups <- 2
  x_limits <- c(0.5 - (num_groups * 0.04), num_groups + 0.5 + (num_groups * 0.04))
  old_mar <- par(mar = c(5, 4, 3.3, 0.2) + 0.1)
  title_text <- PlotViolin(plot_df, x_positions, x_limits, colors_df, y_column, show_title = show_title)


  ## Draw the group labels
  old_lheight <- par("lheight" = 1.15)
  text(x      = x_positions,
       y      = par("usr")[[3]] - ((par("usr")[[4]] - par("usr")[[3]]) * 0.085),
       labels = c("Rest", "4sg"),
       adj    = c(0.5, 1),
       cex    = 0.9,
       col    =  "black",
       font   = 2,
       xpd    = NA
       )
  par(old_lheight)

  ## Final steps
  box(xpd = NA, lwd = 0.75)
  par(old_mar)
  results_list <- list("plot_df" = plot_df, "title_text" = title_text)
  return(invisible(results_list))
}



### Delete this!!

# CRISPR_df <- merged_CRISPRko_df
# y_column <- "GuideScan_specificity"
# show_rest_v_4sg        = FALSE
# aggregate_scores       = TRUE
# show_sublibraries      = !(show_rest_v_4sg || aggregate_scores)
# filter_top4            = aggregate_scores
# filter_complete_genes  = FALSE
# filter_complete_scores = FALSE
# collapse_GPP           = filter_top4
#
# CRISPR_df              <- merged_replaced_CRISPRa_df
# y_column               <- "Deviation_from_TSS_window"
# aggregate_scores       = FALSE
# filter_top4            = aggregate_scores
# show_sublibraries      = FALSE
# filter_complete_genes  = FALSE
# filter_complete_scores = FALSE
# show_rest_v_4sg        = FALSE
# collapse_GPP           = aggregate_scores





ViolinBox_Sources <- function(CRISPR_df,
                              y_column,
                              show_rest_v_4sg        = FALSE,
                              aggregate_scores       = FALSE,
                              show_sublibraries      = !(show_rest_v_4sg || aggregate_scores),
                              filter_top4            = aggregate_scores,
                              filter_complete_genes  = TRUE,
                              filter_complete_scores = TRUE,
                              collapse_GPP           = filter_top4
                              ) {

  if (show_sublibraries && show_rest_v_4sg) {
    stop("The 'show_sublibraries' and 'show_rest_v_4sg' arguments may not both be TRUE!")
  }
  if (aggregate_scores) {
    if (show_sublibraries) {
      stop("If 'aggregate_scores' is TRUE, the 'show_sublibraries' argument cannot also be TRUE!")
    }
    if (!(filter_top4 && filter_complete_scores)) {
      stop("If 'aggregate_scores' is TRUE, 'filter_top4' and 'filter_complete_scores' must also be TRUE!")
    }
    if (!(y_column %in% CFD_specificity_scores)) {
      stop("Only specificity scores can be aggregated!")
    }
  }

  ## Set up the data frames for plotting
  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), y_column)

  if (filter_top4) {
    is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
    if (!(is_CRISPRko)) {
      are_main_transcript <- (GetMainTSS(CRISPR_df) == CRISPR_df[["AltTSS_ID"]]) %in% TRUE
      CRISPR_df <- CRISPR_df[are_main_transcript, ]
    }
  }

  plot_df <- ExpandedSubgroupsDf(CRISPR_df, collapse_GPP = collapse_GPP, show_rest_v_4sg = show_rest_v_4sg)

  if (filter_top4) {
    plot_df <- FilterTop4(plot_df,
                          filter_complete_genes  = filter_complete_genes,
                          filter_complete_scores = filter_complete_scores,
                          data_column            = y_column,
                          show_sublibraries      = show_sublibraries
                          )
  }
  assign("delete_plot_df", plot_df, envir = globalenv())
  colors_df <- FilterOriginalColorsDf(SubgroupColorsDf(), plot_df, show_rest_v_4sg = show_rest_v_4sg)

  plot_df <- plot_df[, c("Combined_ID", "Entrez_ID", "Gene_symbol", "Group", "Subgroup", y_column)]


  if (show_sublibraries || show_rest_v_4sg) {
    plot_df[["Groups_factor"]] <- droplevels(interaction(plot_df[["Group"]], plot_df[["Subgroup"]], lex.order = TRUE))
  } else {
    colors_df <- colors_df[!(colors_df[["Lightened"]]), ]
    plot_df[["Groups_factor"]] <- plot_df[["Group"]]
  }

  if (aggregate_scores) {
    assign("delete_plot_df", plot_df, envir = globalenv())
    plot_df <- plot_df[order(plot_df[["Group"]], as.integer(plot_df[["Entrez_ID"]])), ]
    gene_IDs_fac <- factor(plot_df[["Entrez_ID"]], levels = unique(plot_df[["Entrez_ID"]]))
    interaction_fac <- interaction(plot_df[["Groups_factor"]], gene_IDs_fac, lex.order = TRUE)
    if (!(filter_complete_genes)) {
      interaction_fac <- droplevels(interaction_fac)
    }
    stopifnot(all(table(interaction_fac) == 4))
    CheckThatFactorIsInOrder(interaction_fac)
    aggregate_scores_vec <- tapply(plot_df[[y_column]], interaction_fac, AggregateSpecificityScores)
    are_to_keep <- rep_len(c(TRUE, FALSE, FALSE, FALSE), length.out = nrow(plot_df))
    plot_df <- plot_df[are_to_keep, ]
    plot_df[[y_column]] <- aggregate_scores_vec
  }

  plot_df[["Point_colors"]] <- rep(colors_df[["Dark"]], as.integer(table(plot_df[["Groups_factor"]])))
  plot_df[["Numeric_data"]] <- plot_df[[y_column]]


  ## Determine the spacing and the x positions
  spaces_vec <- ifelse(colors_df[["Are_new_color"]][-1], 1.15, 0.72)
  x_positions <- cumsum(c(1, spaces_vec))

  num_subgroups <- nlevels(plot_df[["Groups_factor"]])
  x_limits <- c((x_positions[[1]] - 0.5) - (num_subgroups * 0.04),
                x_positions[[length(x_positions)]] + 0.5 + (num_subgroups * 0.04)
                )

  ## Draw the violins
  old_mar <- par(mar = c(4.4, 4, 3.5, 3) + 0.1)
  PlotViolin(plot_df, x_positions, x_limits, colors_df, y_column,
             aggregate_scores = aggregate_scores,
             title_cex = 0.8, title_line = 1.5
             )

  SourcesSubgroupLabels(plot_df, x_positions, colors_df,
                        show_sublibraries = show_sublibraries,
                        show_rest_v_4sg = show_rest_v_4sg, extra_space = TRUE
                        )

  ## Final steps
  par(old_mar)
  return(invisible(NULL))
}






ViolinBox_UniqueLibraries <- function(CRISPR_df, y_column, show_title = TRUE) {

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), y_column)

  plot_df <- MutuallyExclusiveSubgroupsDf(CRISPR_df, y_column)

  num_subgroups <- nlevels(plot_df[["Subgroup"]])


  ## Draw the violins
  colors_df <- SubgroupColorsDf()

  plot_df[["Point_colors"]] <- rep(colors_df[["Dark"]], as.integer(table(plot_df[["Subgroup"]])))
  plot_df[["Groups_factor"]] <- plot_df[["Subgroup"]]
  plot_df[["Numeric_data"]] <- plot_df[["Data"]]

  x_positions <- seq_len(num_subgroups) + (rep_len(c(1, -1), length.out = num_subgroups) * 0.18)
  x_limits <- c(0.85 - (num_subgroups * 0.04), num_subgroups + 0.15 + (num_subgroups * 0.04))
  old_mar <- par(mar = c(5, 4, 3.3, 3) + 0.1)
  title_text <- PlotViolin(plot_df, x_positions, x_limits, colors_df, y_column, show_title = show_title)


  ## Draw the group labels
  UniqueSequencesSubgroupLabels(plot_df, x_positions, colors_df, extra_space = TRUE)


  ## Final steps
  par(old_mar)
  results_list <- list("plot_df" = plot_df, "title_text" = title_text)
  return(invisible(results_list))
}




UniquePointsBoxPlots <- function(CRISPR_df) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
  if (is_CRISPRko) {
    numeric_column_labels <- numeric_column_labels[names(numeric_column_labels) != "Deviation_from_TSS_window"]
  }

  for (make_PDF in c(FALSE, TRUE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "Box plots - A) unique points.pdf"),
          width = pdf_width, height = pdf_height
          )
    }
    for (y_column in names(numeric_column_labels)) {
      layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
      ViolinBox_UniqueTwoGroups(CRISPR_df, y_column, show_title = FALSE)
      ViolinBox_results <- ViolinBox_UniqueLibraries(CRISPR_df, y_column, show_title = FALSE)
      OuterTitleForLayout(ViolinBox_results[["title_text"]], extra_space = TRUE)
      layout(1)
    }
    if (make_PDF) {
      dev.off()
    }
  }
  numeric_seq <- seq_along(numeric_column_labels)
  file_numbers <- FormatFixedWidthInteger(numeric_seq)
  for (i in seq_along(numeric_column_labels)) {
    file_name <- paste0(file_numbers[[i]], ") ", names(numeric_column_labels)[[i]])
    png(file = file.path(output_plots_directory, "Box plots - A) unique points", paste0("Box plots - ", file_name, ".png")),
        width = pdf_width, height = pdf_height - 0.14, units = "in", res = 600
        )
    layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
    ViolinBox_UniqueTwoGroups(CRISPR_df, names(numeric_column_labels)[[i]], show_title = FALSE)
    ViolinBox_results <- ViolinBox_UniqueLibraries(CRISPR_df, names(numeric_column_labels)[[i]], show_title = FALSE)
    OuterTitleForLayout(ViolinBox_results[["title_text"]])
    dev.off()
  }
}



SourcesBoxPlots <- function(CRISPR_df) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% names(my_df)
  if (is_CRISPRko) {
    args_list <- c(args_list, CRISPRko_args_list)
    numeric_column_labels <- numeric_column_labels[names(numeric_column_labels) != "Deviation_from_TSS_window"]
  } else {
    args_list <- c(args_list, CRISPRa_args_list)
  }
  use_width <- pdf_width * 0.85
  use_height <- pdf_height * 1.2

  for (make_PNG in c(FALSE, TRUE)) {
    for (make_PDF in c(FALSE, TRUE)) {
      if (make_PNG && make_PDF) {
        next
      }
      for (file_name in names(args_list)) {
        if (make_PDF) {
          PDF_file_name <- paste0("Box plots - ", file_name, ".pdf")
          pdf(file = file.path(output_plots_directory, PDF_file_name),
              width = use_width, height = use_height
              )
        }
        numeric_seq <- seq_along(numeric_column_labels)
        file_numbers <- FormatFixedWidthInteger(numeric_seq)
        for (i in numeric_seq) {
          y_column <- names(numeric_column_labels)[[i]]
          if (identical(args_list[[file_name]][["aggregate_scores"]], TRUE) &&
              !(y_column %in% CFD_specificity_scores)
              ) {
            next
          }
          if (make_PNG) {
            folder_path <- file.path(output_plots_directory, paste0("Box plots - ", file_name))
            if (!(dir.exists(folder_path))) {
              dir.create(folder_path)
            }
            PNG_file_name <- paste0("Box plots - ", file_name, " - ",
                                    file_numbers[[i]], ") ", y_column, ".png"
                                    )
            png(file = file.path(folder_path, PNG_file_name),
                width = use_width, height = use_height, units = "in", res = 600
                )
          }
          do.call(ViolinBox_Sources, c(list(y_column = y_column, CRISPR_df = CRISPR_df), args_list[[file_name]]))
          if (make_PNG) {
            dev.off()
          }
        }
        if (make_PDF) {
          dev.off()
        }
      }
    }
  }
}






# Functions for plotting categorical data ---------------------------------

MakeCategoricalList <- function(CRISPR_df, use_column, use_cutoff) {
  assign("delete_use_column", use_column, envir = globalenv())
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  assign("delete_use_cutoff", use_cutoff, envir = globalenv())
  title_cex <- 1
  if (use_column == "Are_overlapping") {
    logical_vec <- CRISPR_df[["Num_overlaps"]] > 0
    use_title <- "Overlap with other guides (fewer than 50 base pairs apart)"
  } else if (use_column == "Not_mapped") {
    logical_vec <- is.na(CRISPR_df[["Start"]])
    use_title <- "Cannot be mapped to a specific location in the genome"
  } else if (use_column == "Do_not_meet_criteria") {
    logical_vec <- !(MeetCriteria(CRISPR_df))
    use_title <- "Do not meet all criteria"
  } else if (use_column == "Deviation_from_TSS_window") {
    logical_vec <- CRISPR_df[["Deviation_from_TSS_window"]] != 0
    use_title <- "Lie outside the optimal window around the TSS (-150 to -75 bp)"
    use_column <- "Outside_optimal_TSS_window"
  } else if (use_column == "CRISPOR_Graf_status") {
    logical_vec <- CRISPR_df[["CRISPOR_Graf_status"]] != "GrafOK"
    use_title <- "Violate criteria of Graf et al."
    use_column <- "Graf_not_OK"
  } else if (use_column == "Poly_T") {
    logical_vec <- grepl("TTTT", CRISPR_df[["sgRNA_sequence"]], ignore.case = TRUE)
    use_title <- "Contain TTTT sequence"
  } else if (use_column == "Non_canonical_PAM") {
    logical_vec <- substr(CRISPR_df[["PAM"]], 2, 3) != "GG"
    use_title <- "Non-canonical PAM (i.e. not NGG)"
  } else if (use_column %in% grep("_SNP_AF_(max|sum)_", colnames(CRISPR_df), value = TRUE)) {
    if (is.null(use_cutoff)) {
      use_cutoff <- SNP_frequency_cutoff * 100
    }
    logical_vec <- CRISPR_df[[use_column]] > use_cutoff
    use_title <- paste0("Overlap with a genetic polymorphism (frequency >",
                        use_cutoff, "%)"
                        )
    use_column <- "Overlap_with_SNP"
  } else if (use_column %in% names(specificity_score_cutoffs)) {
    if (is.null(use_cutoff)) {
      use_cutoff <- specificity_score_cutoffs[[use_column]]
    }
    logical_vec <- CRISPR_df[[use_column]] < use_cutoff
    use_title <- paste0("Unspecific (<", use_cutoff, ") according to the ",
                        numeric_column_labels[[use_column]]
                        )
    title_cex <- 0.9
    use_column <- "Are_unspecific"
  } else if (use_column %in% c("GuideScan_efficiency", "CRISPOR_Doench_efficacy")) {
    if (is.null(use_cutoff)) {
      use_cutoff <- 50
    }
    logical_vec <- CRISPR_df[[use_column]] < use_cutoff
    if (use_column == "GuideScan_efficiency") {
      via_string  <- "GuideScan"
    } else if (use_column == "CRISPOR_Doench_efficacy") {
      via_string <- "CRISPOR"
    }
    use_title <- paste0("Suboptimal efficacy (<", use_cutoff,
                        ") according to the Rule Set 2 score (via ", via_string, ")"
                        )
    title_cex <- 0.9
    use_column <- "Is_inefficient"
  } else {
    stop("Unknown column selected!")
  }
  results_list <- list(
    "logical_vec" = logical_vec,
    "use_column"  = use_column,
    "title_text"  = use_title,
    "title_cex"   = title_cex
  )
  return(results_list)
}



BarPlotXlimits <- function(spaces_vec, extend_factor = 0.04, side_space = NULL) {
  bar_width <- 1
  half_bar <- bar_width / 2
  bar_positions <- cumsum(spaces_vec + bar_width) - half_bar
  x_min <- min(bar_positions) - half_bar
  x_max <- max(bar_positions) + half_bar
  if (is.null(side_space)) {
    x_span <- x_max - x_min
    side_space <- x_span * extend_factor
  }
  x_limits <- c(x_min - side_space, x_max + side_space)
  return(x_limits)
}


DrawBarplotGridAndAxis <- function() {
  segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]], y0 = seq(0, 1, by = 0.2),   col = "gray95", lwd = 0.5)
  segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]], y0 = seq(0.1, 1, by = 0.2), col = "gray95", lwd = 0.5)
  axis(2, labels = paste0(seq(0, 100, by = 20), "%"), at = seq(0, 1, by = 0.2),
       mgp = c(3, 0.45, 0), las = 1, cex.axis = 0.8, lwd = 0.75, tcl = -0.3
       )
  return(invisible(NULL))
}



AnnotateBars <- function(bar_positions_vec, counts_mat, proportions_mat, text_cex, smaller_plot_factor = FALSE) {

  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  if (smaller_plot_factor) {
    percentage_plot_factor <- 0.03
    fractions_plot_factor <- 0.0675
  } else {
    percentage_plot_factor <- 0.0275
    fractions_plot_factor <- 0.075
  }

  ## Annotate the bars with percentages
  percent_vec <- proportions_mat[2, ] * 100
  old_scipen <- options(scipen = 1)
  proportions_vec <- ifelse(percent_vec < 1,
                            signif(percent_vec, digits = 1),
                            ifelse(percent_vec > 99.5,
                                   signif(percent_vec, digits = 3),
                                   signif(percent_vec, digits = 2)
                                   )
                            )

  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * percentage_plot_factor),
       labels = paste0(proportions_vec, "%"),
       adj    = c(0.5, 1),
       cex    = text_cex * 1.3,
       col    = "black",
       font   = 2,
       xpd    = NA
       )
  options(old_scipen)

  ## Annotate the bars with fractions (counts)
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * fractions_plot_factor),
       labels = colSums(counts_mat),
       adj    = c(0.5, 0.5),
       cex    = text_cex,
       col    = "black",
       xpd    = NA
       )
  line_widths <- strwidth(colSums(counts_mat), cex = text_cex)
  segments(x0  = bar_positions_vec - (line_widths / 2),
           x1  = bar_positions_vec + (line_widths / 2),
           y0  = par("usr")[[4]] + (plot_height * (fractions_plot_factor + 0.0125)),
           col = "black",
           xpd = NA,
           lwd = 0.3
           )
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * (fractions_plot_factor + 2 * 0.0125)),
       labels = counts_mat[2, ],
       adj    = c(0.5, 0.5),
       cex    = text_cex,
       col    = "black",
       xpd    = NA
       )
  return(invisible(NULL))
}




PlotBars <- function(proportions_mat, spaces_vec, x_limits, colors_df) {

  plot(1,
       xlim = x_limits,
       ylim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       ann  = FALSE,
       type = "n"
       )
  DrawBarplotGridAndAxis()

  ## Draw the bar plots()
  bar_positions <- barplot(height    = proportions_mat[1, ],
                           col       = colors_df[["Pale"]],
                           ylim      = c(0, 1),
                           las       = 1,
                           tcl       = -0.4,
                           border    = NA,
                           space     = spaces_vec,
                           names.arg = rep.int("", nrow(colors_df)),
                           axes      = FALSE,
                           add       = TRUE
                           )[, 1]

  barplot(height    = proportions_mat[2, ],
          offset    = proportions_mat[1, ],
          col       = colors_df[["Dark"]],
          border    = NA,
          space     = spaces_vec,
          names.arg = rep.int("", nrow(colors_df)),
          axes      = FALSE,
          add       = TRUE
          )
  # box(bty = "l", xpd = NA, lwd = 0.75)
  return(invisible(bar_positions))
}





BarPlot_UniqueTwoGroups <- function(CRISPR_df,
                                    use_column,
                                    use_cutoff     = NULL,
                                    show_title     = TRUE,
                                    only_chosen    = FALSE
                                    ) {

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), use_column)

  ## Prepare the data
  categorical_list <- MakeCategoricalList(CRISPR_df, use_column, use_cutoff)
  colors_df <- SubgroupColorsDf()
  colors_df <- colors_df[(colors_df[["Color_name"]] == "Grey") & !(colors_df[["Lightened"]]), ]

  if (only_chosen) {
    are_top_4 <- CRISPR_df[["Rank"]] %in% 1:4
    CRISPR_df <- CRISPR_df[are_top_4, ]
    categorical_list[["logical_vec"]] <- categorical_list[["logical_vec"]][are_top_4]
    num_subgroups <- 1L
    x_spaces <- NULL
  } else {
    num_subgroups <- 2L
    x_spaces <- 0.5
  }

  plot_df <- data.frame(
    "Data"   = categorical_list[["logical_vec"]],
    "Chosen" = CRISPR_df[["Rank"]] %in% 1:4,
    stringsAsFactors = FALSE
  )
  counts_mat <- as.matrix(table(factor(plot_df[["Data"]], levels = c(FALSE, TRUE)), plot_df[["Chosen"]]))
  proportions_mat <- prop.table(counts_mat, margin = 2)


  ## Prepare the plot region
  if (only_chosen) {
    x_limits <- c(-0.34, 1.74)
  } else {
    x_limits <- c(0.1, 3.4)
  }
  old_mar <- par(mar = c(4, 4, 5, 0.2) + 0.1)

  if (show_title) {
    title(categorical_list[["title_text"]], cex.main = categorical_list[["title_cex"]], line = 3.15)
  }

  bar_positions <- PlotBars(proportions_mat, x_spaces, x_limits, colors_df)

  ## Draw the group labels
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  old_lheight <- par("lheight" = 1.1)
  if (only_chosen) {
    group_labels <- "All chosen guides"
  } else {
    group_labels <- c("Rest", "4sg")
  }
  text(x      = bar_positions,
       y      = par("usr")[[3]] - (plot_height * 0.055),
       labels = group_labels,
       adj    = c(0.5, 1),
       cex    = 0.9,
       col    =  "black",
       font   = 2,
       xpd    = NA
       )
  par(old_lheight)
  AnnotateBars(bar_positions, counts_mat, proportions_mat, text_cex = 0.5)

  ## Final steps
  par(old_mar)
  results_list <- list("plot_df" = plot_df, categorical_list[c("title_text", "title_cex")])
  return(invisible(results_list))
}




# CRISPR_df <- merged_replaced_CRISPRa_df
# use_column <- "Have_homologies"
# use_cutoff             = NULL
# show_rest_v_4sg        = FALSE
# show_sublibraries      = FALSE
# filter_top4            = use_column %in% categorical_4sg_columns
# filter_complete_genes  = FALSE
# filter_complete_scores = TRUE
# collapse_GPP           = filter_top4




BarPlot_Sources <- function(CRISPR_df,
                            use_column,
                            use_cutoff             = NULL,
                            show_rest_v_4sg        = FALSE,
                            show_sublibraries      = !(show_rest_v_4sg || (use_column %in% categorical_4sg_columns)),
                            filter_top4            = use_column %in% categorical_4sg_columns,
                            filter_complete_genes  = TRUE,
                            filter_complete_scores = TRUE,
                            collapse_GPP           = filter_top4
                            ) {

  assign("delete_use_column",             use_column,             envir = globalenv())
  assign("delete_use_cutoff",             use_cutoff,             envir = globalenv())
  assign("delete_show_rest_v_4sg",        show_rest_v_4sg,        envir = globalenv())
  assign("delete_show_sublibraries",      show_sublibraries,      envir = globalenv())
  assign("delete_filter_top4",            filter_top4,            envir = globalenv())
  assign("delete_filter_complete_genes",  filter_complete_genes,  envir = globalenv())
  assign("delete_filter_complete_scores", filter_complete_scores, envir = globalenv())
  assign("delete_collapse_GPP",           collapse_GPP,           envir = globalenv())

  if (show_sublibraries && show_rest_v_4sg) {
    stop("The 'show_sublibraries' and 'show_rest_v_4sg' arguments may not both be TRUE!")
  }
  if (!(filter_top4) && (use_column == "Are_overlapping")) {
    warning("It makes little sense to graph the fraction of overlapping guides across libraries unless the top 4 guides are selected!")
  }
  if (!(filter_top4) && (use_column == "Have_homologies")) {
    stop("The 'filter_top4' argument must be TRUE if the 'Have_homologies' parameter is selected!")
  }

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), use_column)

  if (filter_top4) {
    is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
    if (!(is_CRISPRko)) {
      are_main_transcript <- (GetMainTSS(CRISPR_df) == CRISPR_df[["AltTSS_ID"]]) %in% TRUE
      CRISPR_df <- CRISPR_df[are_main_transcript, ]
    }
  }

  plot_df <- ExpandedSubgroupsDf(CRISPR_df, collapse_GPP = collapse_GPP, show_rest_v_4sg = show_rest_v_4sg)

  if (use_column == "Have_homologies") {
    categorical_list <- list(
      "logical_vec" = NULL,
      "use_column"  = use_column,
      "title_text"  = "Top 4 guides share identical sub-sequences (8 base pairs or longer)",
      "title_cex"   = 0.9
    )
  } else {
    categorical_list <- MakeCategoricalList(plot_df, use_column, use_cutoff)
  }

  new_column <- categorical_list[["use_column"]]

  plot_df[[new_column]] <- categorical_list[["logical_vec"]]

  if (filter_top4) {
    plot_df <- FilterTop4(plot_df,
                          filter_complete_genes  = filter_complete_genes,
                          filter_complete_scores = filter_complete_scores,
                          data_column            = new_column,
                          show_sublibraries      = show_sublibraries
                          )
  }
  colors_df <- FilterOriginalColorsDf(SubgroupColorsDf(), plot_df, show_rest_v_4sg = show_rest_v_4sg)

  if (show_sublibraries || show_rest_v_4sg) {
    plot_df[["Groups_factor"]] <- droplevels(interaction(plot_df[["Group"]], plot_df[["Subgroup"]], lex.order = TRUE))
  } else {
    colors_df <- colors_df[!(colors_df[["Lightened"]]), ]
    plot_df[["Groups_factor"]] <- plot_df[["Group"]]
  }

  if (new_column %in% categorical_4sg_columns) {
    if (use_column == "Are_overlapping") {
      new_order <- order(plot_df[["Group"]],
                         GetMinEntrez(plot_df[["Entrez_ID"]]),
                         plot_df[["Cut_location"]]
                         )
    } else {
      new_order <- order(plot_df[["Group"]],
                         GetMinEntrez(plot_df[["Entrez_ID"]])
                         )
    }
    plot_df <- plot_df[new_order, ]
    gene_IDs_fac <- factor(plot_df[["Entrez_ID"]], levels = unique(plot_df[["Entrez_ID"]]))
    interaction_fac <- interaction(plot_df[["Groups_factor"]], gene_IDs_fac, lex.order = TRUE)
    if (!(filter_complete_genes)) {
      interaction_fac <- droplevels(interaction_fac)
    }
    stopifnot(all(table(interaction_fac) == 4))
    CheckThatFactorIsInOrder(interaction_fac)

    # stopifnot(all(tapply(plot_df[["Chromosome"]], interaction_fac, function(x) length(unique(x[!(is.na(x))]))) <= 1))
    are_to_keep <- rep_len(c(TRUE, FALSE, FALSE, FALSE), length.out = nrow(plot_df))
    if (use_column == "Are_overlapping") {
      guide_list <- lapply(seq_len(4), function(x) x == seq_len(4))
      min_space <- 50L
      num_overlaps_vec <- tapply(plot_df[["Cut_location"]],
                                 interaction_fac,
                                 function(x) {
                                   overlap_numbers_vec <- vapply(guide_list,
                                                                 function(y) sum(abs(x[y] - x[!(y)]) < min_space),
                                                                 integer(1)
                                                                 )
                                   total_overlaps <- sum(overlap_numbers_vec) / 2
                                   return(total_overlaps)
                                 })
      plot_df <- plot_df[are_to_keep, ]
      plot_df[["Are_overlapping"]] <- num_overlaps_vec != 0
    } else if (use_column == "Have_homologies") {
      longest_shared_vec <- tapply(plot_df[["sgRNA_sequence"]], interaction_fac, LongestSharedSubsequence)
      plot_df <- plot_df[are_to_keep, ]
      plot_df[["Have_homologies"]] <- longest_shared_vec >= 8
    }
  }

  counts_mat <- as.matrix(table(factor(plot_df[[new_column]], levels = c(FALSE, TRUE)), plot_df[["Groups_factor"]]))
  proportions_mat <- prop.table(counts_mat, margin = 2)

  spaces_vec <- ifelse(colors_df[["Are_new_color"]], 1.3, 0.3)

  x_limits <- BarPlotXlimits(spaces_vec, side_space = max(spaces_vec))

  old_mar <- par(mar = c(4, 4, 6, 3) + 0.1)

  bar_positions <- PlotBars(proportions_mat, spaces_vec, x_limits, colors_df)

  title(categorical_list[["title_text"]], cex.main = categorical_list[["title_cex"]] * 0.9, line = 3.625)

  SourcesSubgroupLabels(plot_df, bar_positions, colors_df, show_sublibraries = show_sublibraries, show_rest_v_4sg = show_rest_v_4sg)

  AnnotateBars(bar_positions, counts_mat, proportions_mat, text_cex = 0.4, smaller_plot_factor  = TRUE)

  ## Final steps
  par(old_mar)
  return(invisible(NULL))
}




BarPlot_UniqueLibraries <- function(CRISPR_df,
                                    use_column,
                                    use_cutoff     = NULL,
                                    show_title     = TRUE,
                                    only_chosen    = FALSE
                                    ) {

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), use_column)

  ## Prepare the data
  categorical_list <- MakeCategoricalList(CRISPR_df, use_column, use_cutoff)
  CRISPR_df[[categorical_list[["use_column"]]]] <- categorical_list[["logical_vec"]]
  if (only_chosen) {
    CRISPR_df <- CRISPR_df[CRISPR_df[["Rank"]] %in% 1:4, ]
  }
  plot_df <- MutuallyExclusiveSubgroupsDf(CRISPR_df, categorical_list[["use_column"]])
  colors_df <- SubgroupColorsDf()
  if (only_chosen) {
    colors_df <- colors_df[!(colors_df[["Lightened"]]), -2]
    group_column <- "Group"
  } else {
    group_column <- "Subgroup"
  }
  counts_mat <- as.matrix(table(factor(plot_df[["Data"]], levels = c(FALSE, TRUE)), plot_df[[group_column]]))
  proportions_mat <- prop.table(counts_mat, margin = 2)
  num_subgroups <- nlevels(plot_df[[group_column]])
  if (only_chosen) {
    spaces_vec <- 0.8
  } else {
    spaces_vec <- rep.int(0.8, num_subgroups) + (rep_len(c(1, -1), length.out = num_subgroups) * 0.5)
  }
  ## Prepare the plot region
  if (only_chosen && (num_subgroups == 7)) {
    x_limits <- c(0.028, 13.372)
  } else {
    x_limits <- BarPlotXlimits(spaces_vec)
  }

  old_mar <- par(mar = c(4, 4, 5, 3) + 0.1)

  bar_positions <- PlotBars(proportions_mat, spaces_vec, x_limits, colors_df)

  if (show_title) {
    title(categorical_list[["title_text"]], cex.main = categorical_list[["title_cex"]], line = 3.15)
  }

  UniqueSequencesSubgroupLabels(plot_df, bar_positions, colors_df, only_chosen = only_chosen)
  AnnotateBars(bar_positions, counts_mat, proportions_mat, text_cex = 0.4)

  ## Final steps
  par(old_mar)
  results_list <- c(list("plot_df" = plot_df), categorical_list[c("title_text", "title_cex")])
  return(invisible(results_list))
}






UniqueSequencesBarPlots <- function(CRISPR_df) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
  if (is_CRISPRko) {
    categorical_columns <- setdiff(categorical_columns, "Deviation_from_TSS_window")
  }

  use_width <- pdf_width * 0.85
  use_height <- pdf_height

  for (only_chosen in c(TRUE, FALSE)) {
    if (only_chosen) {
      folder_name <- "Bar charts - A2) unique sequences - only 4sg"
      use_columns <- c("Are_overlapping", categorical_columns)
    } else {
      folder_name <- "Bar charts - A1) unique sequences - all guides"
      use_columns <- categorical_columns
    }
    for (make_PNG in c(FALSE, TRUE)) {
      for (make_PDF in c(FALSE, TRUE)) {
        if (make_PNG && make_PDF) {
          next
        }
        if (make_PDF) {
          pdf(file = file.path(output_plots_directory,
                               paste0(folder_name, ".pdf")
                               ),
              width = use_width, height = use_height
              )
        }
        categorical_seq <- seq_along(use_columns)
        file_numbers <- FormatFixedWidthInteger(categorical_seq)
        for (i in categorical_seq) {
          use_column <- use_columns[[i]]
          if (make_PNG) {
            folder_path <- file.path(output_plots_directory, folder_name)
            if (!(dir.exists(folder_path))) {
              dir.create(folder_path)
            }
            PNG_file_name <- paste0(folder_name, " - ",
                                    file_numbers[[i]], ") ", use_column, ".png"
                                    )
            png(file = file.path(folder_path, PNG_file_name),
                width = use_width, height = use_height, units = "in", res = 600
                )
          }

          layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
          BarPlot_UniqueTwoGroups(CRISPR_df, use_column, show_title = FALSE, only_chosen = only_chosen)
          StackedBarPlot_results <- BarPlot_UniqueLibraries(CRISPR_df, use_column, show_title = FALSE, only_chosen = only_chosen)
          OuterTitleForLayout(StackedBarPlot_results[["title_text"]], extra_space = FALSE)
          layout(1)

          if (make_PNG) {
            dev.off()
          }
        }
        if (make_PDF) {
          dev.off()
        }
      }
    }
  }
}




SourcesBarPlots <- function(CRISPR_df) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% colnames(CRISPR_df)
  if (is_CRISPRko) {
    categorical_columns <- setdiff(categorical_columns, "Deviation_from_TSS_window")
  }

  if (is_CRISPRko) {
    args_list <- c(args_list, CRISPRko_args_list)
  } else {
    args_list <- c(args_list, CRISPRa_args_list)
  }
  args_list <- lapply(args_list, function(x) x[names(x) != "aggregate_scores"])

  use_width <- pdf_width * 0.85
  use_height <- pdf_height * 1.425

  for (make_PNG in c(FALSE, TRUE)) {
    for (make_PDF in c(FALSE, TRUE)) {
      if (make_PNG && make_PDF) {
        next
      }
      for (file_name in names(args_list)) {
        if (grepl("4-guide combination", file_name, fixed = TRUE)) {
          use_columns <- categorical_4sg_columns
        } else {
          use_columns <- categorical_columns
        }
        if (make_PDF) {
          PDF_file_name <- paste0("Bar charts - ", file_name, ".pdf")
          pdf(file = file.path(output_plots_directory, PDF_file_name),
              width = use_width, height = use_height
          )
        }
        categorical_seq <- seq_along(use_columns)
        file_numbers <- FormatFixedWidthInteger(categorical_seq)

        for (i in categorical_seq) {
          use_column <- use_columns[[i]]
          if ((length(args_list[[file_name]]) == 0) &&
              !(use_column %in% categorical_4sg_columns)
              ) {
            next
          }
          if (make_PNG) {
            folder_path <- file.path(output_plots_directory, paste0("Bar charts - ", file_name))
            if (!(dir.exists(folder_path))) {
              dir.create(folder_path)
            }
            PNG_file_name <- paste0("Bar charts - ", file_name, " - ",
                                    file_numbers[[i]], ") ", use_column, ".png"
                                    )
            png(file = file.path(folder_path, PNG_file_name),
                width = use_width, height = use_height, units = "in", res = 600
                )
          }
          do.call(BarPlot_Sources, c(list(CRISPR_df = CRISPR_df, use_column = use_column), args_list[[file_name]]))
          if (make_PNG) {
            dev.off()
          }
        }
        if (make_PDF) {
          dev.off()
        }
      }
    }
  }
}












# Functions for generating histograms -------------------------------------

DrawHistogram <- function(overview_df, column_name) {
  hist_results <- hist(overview_df[[column_name]], breaks = 100, plot = FALSE)
  y_max <- max(hist_results[["density"]])
  hist(overview_df[[column_name]],
       las      = 1,
       mgp      = c(2.2, 0.45, 0),
       tcl      = -0.3,
       col      = brewer.pal(9, "Blues")[[8]],
       border   = NA,
       breaks   = 100,
       xlim     = c(-0.02, 1),
       ylim     = c(0 - (y_max * 0.02), y_max + (y_max * 0.02)),
       main     = paste0("Aggregate ", numeric_column_labels[[column_name]]),
       cex.main = 1,
       xlab     = "",
       xaxs     = "i",
       yaxs     = "i",
       freq     = FALSE
       )
  box(bty = "l")
  return(invisible(NULL))
}


SharedSubsequencesBarplot <- function(overview_df) {
  longest_subsequence_table <- table(factor(sgRNAs_overview_df[["Longest_subsequence"]], levels = 1:20))
  bar_positions <- barplot(longest_subsequence_table,
                           las       = 1,
                           mgp       = c(3, 0.45, 0),
                           tcl       = -0.3,
                           col       = brewer.pal(9, "Blues")[[8]],
                           main      = "Length of the longest shared subsequence",
                           cex.main  = 1,
                           names.arg = "",
                           border    = NA,
                           space     = 0.4,
                           ylim      = c(0, 10000),
                           xlim      = c(0, 28.4),
                           xaxs      = "i",
                           ylab      = "Count"
                           )
  box(xpd = NA, bty = "l")
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  text(x      = bar_positions[, 1],
       y      = par("usr")[[3]] - (plot_height * 0.1),
       labels = as.character(1:20),
       cex    = 0.8,
       font   = 2,
       xpd    = NA
       )
  text(x      = bar_positions[, 1],
       y      = par("usr")[[3]] - (plot_height * 0.03),
       labels = longest_subsequence_table,
       cex    = 0.5,
       font   = 2,
       xpd    = NA,
       col    = "gray50"
       )
  text(x      = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) / 2),
       y      = par("usr")[[3]] - (plot_height * 0.2),
       labels = "Number of base pairs",
       xpd    = NA
       )
  return(invisible(NULL))
}


Plot4sgData <- function(overview_df) {
  for (make_PDF in c(TRUE, FALSE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "Histograms - 4sg combination.pdf"),
          width = pdf_width, height = pdf_height * 1.3
          )
    }
    par("oma" = c(0, 1, 0, 1))
    SharedSubsequencesBarplot(overview_df)
    for (column_name in CFD_specificity_scores) {
      DrawHistogram(overview_df, column_name)
    }
    if (make_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}






# Functions for generating Venn diagrams ----------------------------------

PlotVennDiagrams <- function(CRISPR_df) {
  for (make_PDF in c(FALSE, TRUE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "Venn diagrams.pdf"),
          width = pdf_width * 1.1, height = pdf_height * 1.05
          )
    }
    CRISPR_df <- FilterCRISPRDf(CRISPR_df)
    are_chosen <- CRISPR_df[["Rank"]] %in% 1:4
    all_sources_fac <- ReformatSourceToFactor(CRISPR_df[["Source"]])

    euler_not_chosen <- PlotVennDiagram(all_sources_fac[!(are_chosen)])
    euler_4sg <- PlotVennDiagram(all_sources_fac[are_chosen])
    empty_plot <- ggplot() + theme_void()
    title_rest <- ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = 'bold("Rest")', parse = TRUE, size = 5)
    title_4sg  <- ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = 'bold("4sg")',  parse = TRUE, size = 5)
    grid.arrange(empty_plot, empty_plot, empty_plot, empty_plot, empty_plot,
                 empty_plot, title_rest, empty_plot, title_4sg, empty_plot,
                 empty_plot, euler_not_chosen, empty_plot, euler_4sg, empty_plot,
                 nrow = 3, widths = c(0.1, 1, 0.1, 1, 0.1), heights = c(0.075, 0.075, 1)
                 )
    if (make_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}


PlotVennDiagram <- function(sources_fac, draw_plot = FALSE) {
  sources_table <- table(sources_fac)
  sources_names <- gsub(", ", "&", names(sources_table), fixed = TRUE)
  sources_table <- as.integer(sources_table)
  names(sources_table) <- sources_names
  eulerr_options("padding" = grid::unit(0.25, "lines"))
  euler_plot <- plot(eulerr::euler(sources_table, shape = "circle"),
                     fills = c(brewer.pal(9, "Blues")[[3]], brewer.pal(9, "Greens")[[3]], brewer.pal(9, "Reds")[[3]]),
                     edges = FALSE,
                     labels = list(cex = 1),
                     quantities = list(font = 2, cex = 0.4)
                     )
  if (draw_plot) {
    print(euler_plot)
  }
  return(invisible(euler_plot))
}





# Functions for producing scatter plots -----------------------------------

ConvertCFDScores <- function(numeric_vec) {
  1 / (1 + (10000 / numeric_vec) - 100)
}


GaussianJitter <- function(numeric_vec) {
  numeric_vec +
  (rnorm(n = length(numeric_vec), mean = 0, sd = 0.05) *
  ((max(numeric_vec, na.rm = TRUE) - min(numeric_vec, na.rm = TRUE)) / 50))
}


ScatterPlot <- function(CRISPR_df,
                        x_column,
                        y_column,
                        identical_axes               = FALSE,
                        mark_diagonal                = identical_axes,
                        convert_CRISPOR_to_GuideScan = FALSE,
                        convert_GuideScan_to_CRISPOR = FALSE,
                        point_cex                    = 0.5,
                        point_alpha                  = 0.6,
                        custom_axis_limits           = NULL,
                        show_title                   = NULL,
                        add_jitter                   = FALSE,
                        make_PNG                     = FALSE,
                        only_top4                    = FALSE
                        ) {

  assign("delete_CRISPR_df",      CRISPR_df,      envir = globalenv())
  assign("delete_x_column",       x_column,       envir = globalenv())
  assign("delete_y_column",       y_column,       envir = globalenv())
  assign("delete_identical_axes", identical_axes, envir = globalenv())
  assign("delete_only_top4",      only_top4,      envir = globalenv())

  CRISPR_df <- FilterCRISPRDf(CRISPR_df)

  if (only_top4) {
    CRISPR_df <- CRISPR_df[CRISPR_df[["Rank"]] %in% 1:4, ]
  }

  if (convert_CRISPOR_to_GuideScan) {
    CRISPR_df[["CRISPOR_CFD_specificity"]] <- ConvertCFDScores(CRISPR_df[["CRISPOR_CFD_specificity"]])
  }
  if (convert_GuideScan_to_CRISPOR) {
    CRISPR_df[["GuideScan_specificity"]] <- (100 / (100 + ((1 / CRISPR_df[["GuideScan_specificity"]]) - 1))) * 100
  }

  if (make_PNG) {
    if (!("current_number" %in% ls(envir = globalenv()))) {
      current_number <- 1L
    }
    number_string <- formatC(current_number, width = 2, flag = "0")
    file_name <- paste0("Scatter plots - ",
                        number_string, ") ", gsub(":", " - ", show_title), ".png"
                        )
    if (only_top4) {
      sub_folder <- "Scatter plots - B) only 4sg"
    } else {
      sub_folder <- "Scatter plots - A) all guides"
    }
    png(file = file.path(output_plots_directory, sub_folder, file_name),
        width = 5.75, height = 5.75, units = "in", res = 600
        )
    current_number <- current_number + 1L
    assign("current_number", current_number, envir = globalenv())
  }

  old_par <- par(mar = rep.int(5, 4))
  x_vec <- CRISPR_df[[x_column]]
  y_vec <- CRISPR_df[[y_column]]

  if (add_jitter) {
    x_vec <- GaussianJitter(x_vec)
    y_vec <- GaussianJitter(y_vec)
  }

  if (identical_axes) {
    if (is.null(custom_axis_limits)) {
      x_axis_limits <- GetAxisLimits(x_vec, x_column, provide_other_limits = FALSE)
      stopifnot(identical(x_axis_limits, GetAxisLimits(y_vec, y_column, provide_other_limits = FALSE)))
      if (is.null(x_axis_limits)) {
        are_not_NA <- (!(is.na(x_vec))) & (!(is.na(y_vec)))
        x_axis_limits <- GetAxisLimits(c(x_vec[are_not_NA], y_vec[are_not_NA]), provide_other_limits = TRUE)
      }
    } else {
      x_axis_limits <- custom_axis_limits
    }
    y_axis_limits <- x_axis_limits
  } else {
    x_axis_limits <- GetAxisLimits(x_vec, x_column, provide_other_limits = TRUE)
    y_axis_limits <- GetAxisLimits(y_vec, y_column, provide_other_limits = TRUE)
  }

  plot(x_vec,
       y_vec,
       las  = 1,
       mgp  = c(2.7, 0.5, 0),
       tcl  = -0.4,
       type = "n",
       xlab = numeric_column_labels[[x_column]],
       ylab = numeric_column_labels[[y_column]],
       xlim = x_axis_limits,
       ylim = y_axis_limits,
       xaxs = "i",
       yaxs = "i"
       )

  assign("delete_mark_diagonal", mark_diagonal, envir = globalenv())
  if (mark_diagonal) {
    abline(a = 0, b = 1, col = "gray88", lwd = 0.5)
  }

  line_color <- "gray88"
  line_type <- "solid"
  line_width <- 0.75
  if (x_column %in% c("GuideScan_specificity", "CRISPOR_3MM_specificity")) {
    abline(v = 0.2, col = line_color, lty = line_type, lwd = line_width)
  }
  if (x_column %in% c("CRISPOR_CFD_specificity")) {
    abline(v = 80, col = line_color, lty = line_type, lwd = line_width)
  }
  if (y_column %in% c("GuideScan_specificity", "CRISPOR_3MM_specificity")) {
    abline(h = 0.2, col = line_color, lty = line_type, lwd = line_width)
  }
  if (y_column %in% c("CRISPOR_CFD_specificity")) {
    abline(h = 80, col = line_color, lty = line_type, lwd = line_width)
  }
  box()

  alpha_hex <- substr(rgb(1, 1, 1, point_alpha), 8, 9)
  points(x_vec,
         y_vec,
         col = paste0(brewer.pal(9, "Blues")[[7]], alpha_hex),
         pch = 16,
         cex = point_cex
         )
  box()
  if (!(is.null(show_title))) {
    title(show_title, cex.main = par("cex") * 0.9)
  }
  assign("delete_x_vec", x_vec, envir = globalenv())
  assign("delete_y_vec", y_vec, envir = globalenv())
  par(old_par)

  if (make_PNG) {
    dev.off()
  }
  return(invisible(NULL))
}





MakeScatterPlots <- function(CRISPR_df) {

  for (i in 1:3) {
    if (i == 1) {
      make_PDF <- FALSE
      make_PNG <- FALSE
    } else if (i == 2) {
      make_PDF <- FALSE
      make_PNG <- TRUE
    } else {
      make_PDF <- TRUE
      make_PNG <- FALSE
    }
    for (only_top4 in c(TRUE, FALSE)) {
      for (only_selected in c(FALSE)) {
        if ((only_selected) && !(make_PDF)) {
          next
        }
        assign("current_number", 1L, envir = globalenv())
        if (make_PDF) {
          if (only_top4) {
            append_to_filename <- " - B) only 4sg"
          } else {
            append_to_filename <- " - A) all guides"
          }
          if (only_selected) {
            append_to_filename <- paste0(append_to_filename, " - selected")
          }
          plot_dimensions <- 5.75
          pdf(file = file.path(output_plots_directory, paste0("Scatter plots", append_to_filename, ".pdf")),
              width = plot_dimensions, height = plot_dimensions
              )
        }
        ScatterPlot(CRISPR_df, "GuideScan_specificity", "GuideScan_efficiency",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "Efficacy vs. specificity (GuideScan)",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )
        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_3MM_specificity", "CRISPOR_Doench_efficacy",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title = "Efficacy vs. specificity (CRISPOR)",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_Doench_efficacy",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title = "Efficacy vs. specificity (CRISPOR, up to 4MM)",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_Doench_efficacy", "GuideScan_efficiency",
                      point_alpha = 0.2, point_cex = 0.4, add_jitter = TRUE,
                      show_title = "Original vs. updated Doench efficacy scores",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
        }

        ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                    show_title = "Specificity score vs. number of off-targets (GuideScan)",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )

        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_Num_2or3MM", "CRISPOR_3MM_specificity",
                      point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                      show_title = "Specificity score vs. number of off-targets (CRISPOR)",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
          if (FALSE) {
            ScatterPlot(CRISPR_df, "GuideScan_specificity", "CRISPOR_CFD_specificity",
                        point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                        make_PNG = make_PNG, only_top4 = only_top4
                        )
            ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_CFD_specificity", # This is just for confirmation
                        convert_GuideScan_to_CRISPOR = TRUE, point_alpha = 0.2, point_cex = 0.4,
                        make_PNG = make_PNG, only_top4 = only_top4
                        )
          }
        }

        ScatterPlot(CRISPR_df, "CRISPOR_Num_2or3MM", "GuideScan_Num_2or3MM",
                    point_alpha = 0.5, point_cex = 0.2, identical_axes = TRUE,
                    show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )
        ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "CRISPOR_Num_2or3MM",
                    point_alpha = 0.5, point_cex = 0.2,
                    identical_axes = TRUE, custom_axis_limits = c(0, 1000),
                    show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in)",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )
        ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "CRISPOR_Num_2or3MM",
                    point_alpha = 0.5, point_cex = 0.2,
                    identical_axes = TRUE, custom_axis_limits = c(0, 200),
                    show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in more)",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )

        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_Num_2MM", "GuideScan_Num_2MM",
                      point_alpha = 0.5, point_cex = 0.2, identical_axes = TRUE,
                      show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
          ScatterPlot(CRISPR_df, "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
                      point_alpha = 0.5, point_cex = 0.2,
                      identical_axes = TRUE,
                      custom_axis_limits = c(0, 300),
                      show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in)",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
          ScatterPlot(CRISPR_df, "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
                      point_alpha = 0.5, point_cex = 0.2,
                      identical_axes = TRUE, custom_axis_limits = c(0, 50),
                      show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in more)",
                      add_jitter = TRUE,
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_CFD_specificity", "GuideScan_specificity",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title  = "Specificity \u2013 GuideScan score vs. original CRISPOR CFD score",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
        }

        ScatterPlot(CRISPR_df, "CRISPOR_3MM_specificity", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4, identical_axes = TRUE,
                    show_title = "Specificity score \u2013 GuideScan vs. CRISPOR",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )

        ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "Specificity score \u2013 GuideScan vs. CRISPOR (4MM)",
                    make_PNG = make_PNG, only_top4 = only_top4
                    )

        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_3MM_specificity",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title = "CRISPOR specificity scores: 3MM vs. 4MM",
                      make_PNG = make_PNG, only_top4 = only_top4
                      )
        }
        if (make_PDF) {
          dev.off()
        }
      }
    }
  }
}






# Functions for generating multi-plot layouts -----------------------------

MakeEmptyPlot <- function() {
  plot(1, xlim = c(0, 1), ylim = c(0, 1), type = "n", ann = FALSE, axes = FALSE)
}


OuterTitleForLayout <- function(title_text, title_cex = 1, extra_space = FALSE) {
  old_mar <- par(mar = rep.int(0, 4))
  if (extra_space) {
    title_y <- -0.57
  } else {
    title_y <- -0.8
  }
  MakeEmptyPlot()
  text(x      = 0.5,
       y      = title_y,
       labels = title_text,
       cex    = title_cex,
       font   = 2,
       xpd    = NA
       )
  MakeEmptyPlot()
  par(old_mar)
}






