### 9th April 2020 ###



# Import packages and source code -----------------------------------------

library("rcartocolor")
library("png")
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
  "Deviation_from_TSS_window" = "Deviation from the optimal window near the TSS"
)

numeric_column_labels <- c(
  core_numeric_column_labels,
  "Restricted_distance_from_TSS" = "Position relative to the TSS",
  "GuideScan_Num_2MM"            = "GuideScan \u2013 number of 2MM sites",
  "GuideScan_Num_3MM"            = "GuideScan \u2013 number of 3MM sites",
  "GuideScan_Num_2or3MM"         = "GuideScan \u2013 number of 2MM or 3MM sites",
  "CRISPOR_Num_2MM"              = "CRISPOR \u2013 number of 2MM sites",
  "CRISPOR_Num_3MM"              = "CRISPOR \u2013 number of 3MM sites",
  "CRISPOR_Num_2or3MM"           = "CRISPOR \u2013 number of 2MM or 3MM sites",
  "CRISPOR_Num_4MM"              = "CRISPOR \u2013 number of 4MM sites"
)


aggregate_sum_column_labels <- c(
  "all22_SNP_AF_max_Kaviar" = "Expected to overlap with a genetic polymorphism"
)




specificity_score_cutoffs <- c(
  "GuideScan_specificity"   = 0.2,
  "CRISPOR_3MM_specificity" = 0.2,
  "CRISPOR_4MM_specificity" = 0.04,
  "CRISPOR_CFD_specificity" = 50,
  "CRISPOR_MIT_specificity" = 50
)


unintended_target_labels <- c(
  "Does_not_affect_intended_gene"          = "Do not seem to target the intended gene",

  "Affects_any_unintended_gene"            = "Affect a gene other than the intended gene",
  "Affects_unintended_protein_gene"        = "Affect an unintended protein-coding gene",
  "Affects_any_genes_at_other_loci"        = "Affect an unintended gene at another locus",
  "Affects_protein_genes_at_other_loci"    = "Affect an unintended protein-coding gene at another locus",

  "Affects_unintended_main_TSS"            = "Affect a gene other than the intended gene (only main TSS)",
  "Affects_unintended_protein_main_TSS"    = "Affect an unintended protein-coding gene (only main TSS)",
  "Affects_main_TSS_at_other_loci"         = "Affect an unintended gene at another locus (only main TSS)",
  "Affects_protein_main_TSS_at_other_loci" = "Affect an unintended protein-coding gene at another locus (only main TSS)"
)


categorical_columns <- c(
  names(unintended_target_labels),
  "Do_not_meet_criteria",
  "Not_mapped",
  "Poly_T",
  "Non_canonical_PAM",
  "CRISPOR_Graf_status",
  "Expected_all22_SNP_AF_max_Kaviar",
  names(core_numeric_column_labels)
)


TSS_columns <- c(
  "Deviation_from_TSS_window",
  "Restricted_distance_from_TSS",
  "Affects_unintended_main_TSS",
  "Affects_unintended_protein_main_TSS",
  "Affects_main_TSS_at_other_loci",
  "Affects_protein_main_TSS_at_other_loci"
)



sublibraries_order <- c(
  "Calabrese_1to3",
  "Calabrese_4to6",
  "hCRISPRa-v2_1to5",
  "hCRISPRa-v2_6to10",
  "Dolcetto_1to3",
  "Dolcetto_4to6",
  "hCRISPRi-v2.1_1to5",
  "hCRISPRi-v2.1_6to10",
  "Brunello",
  "TKOv3",
  "GPP_1to10",
  "GPP_top10",
  "GPP_11to50",
  "GPP_rest",
  "None"
)



sublibraries_labels <- c(
  "Calabrese_1to3"      = "# 1-3",
  "Calabrese_4to6"      = "# 4-6",
  "hCRISPRa-v2_1to5"    = "Top 5",
  "hCRISPRa-v2_6to10"   = "Supp 5",
  "Dolcetto_1to3"       = "# 1-3",
  "Dolcetto_4to6"       = "# 4-6",
  "hCRISPRi-v2.1_1to5"  = "Top 5",
  "hCRISPRi-v2.1_6to10" = "Supp 5",
  "GPP_1to10"           = "Top 10",
  "GPP_top10"           = "Top 10",
  "GPP_11to50"          = "# 11-50",
  "GPP_rest"            = "Other",
  "None"                = ""
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




# Helper functions for addint annotation for unintended targets  ----------

AddOtherTargetBooleans <- function(CRISPR_df,
                                   all_genes_df,
                                   protein_genes_df,
                                   all_genes_main_TSS_df = NULL,
                                   protein_genes_main_TSS_df = NULL
                                   ) {
  stopifnot(length(unique(nrow(CRISPR_df), nrow(all_genes_df), nrow(protein_genes_df))) == 1)

  CRISPR_df[["Does_not_affect_intended_gene"]]         <- !(all_genes_df[["Affects_intended_gene"]])

  CRISPR_df[["Affects_any_unintended_gene"]]           <- all_genes_df[["Affects_unintended_gene"]]
  CRISPR_df[["Affects_unintended_protein_gene"]]       <- protein_genes_df[["Affects_unintended_gene"]]

  CRISPR_df[["Affects_any_genes_at_other_loci"]]       <- all_genes_df[["Affects_genes_at_other_loci"]]
  CRISPR_df[["Affects_protein_genes_at_other_loci"]]   <- protein_genes_df[["Affects_genes_at_other_loci"]]

  if (!(is.null(all_genes_main_TSS_df))) {
    CRISPR_df[["Affects_unintended_main_TSS"]]            <- all_genes_main_TSS_df[["Affects_unintended_gene"]]
    CRISPR_df[["Affects_unintended_protein_main_TSS"]]    <- protein_genes_main_TSS_df[["Affects_unintended_gene"]]
    CRISPR_df[["Affects_main_TSS_at_other_loci"]]         <- all_genes_main_TSS_df[["Affects_genes_at_other_loci"]]
    CRISPR_df[["Affects_protein_main_TSS_at_other_loci"]] <- protein_genes_main_TSS_df[["Affects_genes_at_other_loci"]]
  }

  return(CRISPR_df)
}






# Helper functions specific to CRISPRa or CRISPRi -------------------------

IsCRISPRa <- function(CRISPR_df) {
  "hCRISPRa_v2_rank" %in% names(CRISPR_df)
}

GetTSSWindow <- function(is_CRISPRa) {
  if (is_CRISPRa) {
    lower_bound <- -150
    higher_bound <- -75
  } else {
    lower_bound <- +25
    higher_bound <- +75
  }
  return(c(lower_bound, higher_bound))
}


TSSWindowString <- function(is_CRISPRa) {
  TSS_window <- GetTSSWindow(is_CRISPRa)
  TSS_strings <- formatC(TSS_window, format = "d", flag = "+")
  result_string <- paste0(" (", TSS_strings[[1]], " to ", TSS_strings[[2]],
                          " bp)"
                          )
  return(result_string)
}




GetMainTSS <- function(CRISPR_df) {

  are_not_controls <- CRISPR_df[["Is_control"]] == "No"
  altTSS_ID_fac <- factor(CRISPR_df[["AltTSS_ID"]][are_not_controls],
                          levels = unique(CRISPR_df[["AltTSS_ID"]][are_not_controls])
                          )
  CheckThatFactorIsInOrder(altTSS_ID_fac)

  combined_ID_vec <- tapply(CRISPR_df[["Combined_ID"]][are_not_controls], altTSS_ID_fac, unique)

  min_distance_from_TSS_vec <- tapply(CRISPR_df[["Distance_from_TSS"]][are_not_controls],
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

  best_TSS_expanded_vec <- best_TSS_ID_vec[match(CRISPR_df[["Combined_ID"]][are_not_controls], names(best_TSS_ID_vec))]

  best_TSS_including_controls_vec <- rep(NA_character_, nrow(CRISPR_df))
  best_TSS_including_controls_vec[are_not_controls] <- best_TSS_expanded_vec

  return(best_TSS_including_controls_vec)
}






# Helper functions for vertically aligning text on plots ------------------

StripExpression <- function(my_expression) {
  assign("delete_my_expression", my_expression, envir = globalenv())
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
  assign("delete_expression_list", expression_list, envir = globalenv())
  literal_strings <- sapply(expression_list, StripExpression)
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}

VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("gh")), use_expression, expression(phantom("gh")))
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
    }), use.names = FALSE
    )
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

  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())

  sources_vec <- as.character(ReformatSourceToFactor(CRISPR_df[["Source"]]))
  sources_splits <- strsplit(sources_vec, ", ", fixed = TRUE)

  CRISPR_df_expanded <- CRISPR_df[rep(seq_len(nrow(CRISPR_df)), lengths(sources_splits)), ]

  sources_vec <- unlist(sources_splits, use.names = FALSE)
  sources_fac <- factor(sources_vec, levels = intersect(libraries_order, sources_vec))

  is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df_expanded)

  if (show_rest_v_4sg) {
    subgroups_fac <- factor(ifelse(CRISPR_df_expanded[["Rank"]] %in% 1:4, "4sg", "Rest"), levels = c("Rest", "4sg"))
  } else {

    if (is_CRISPRko) {
      sources_sublibraries <- ifelse(sources_vec == "GPP",
                                     paste0("GPP_", ifelse(CRISPR_df_expanded[["GPP_rank"]] %in% 1:10, "1to10", "11to50")),
                                     "None"
                                     )

    } else {
      hCRISPR_rank_column <- grep("v2_rank", names(CRISPR_df_expanded), fixed = TRUE, value = TRUE)
      Doench_rank_column <- grep("(Dolcetto|Calabrese)_rank", names(CRISPR_df_expanded), value = TRUE)
      are_hCRISPR <- sources_vec %in% c("hCRISPRa-v2", "hCRISPRi-v2.1")
      are_Doench <- sources_vec %in% c("Calabrese", "Dolcetto")
      sources_sublibraries <- ifelse(are_hCRISPR,
                                     paste0(sources_vec, "_", ifelse(CRISPR_df_expanded[[hCRISPR_rank_column]] %in% 1:5, "1to5", "6to10")),
                                     ifelse(are_Doench,
                                            paste0(sources_vec, "_", ifelse(CRISPR_df_expanded[[Doench_rank_column]] %in% "1/2/3", "1to3", "4to6")),
                                            ifelse(sources_vec == "GPP",
                                                   paste0("GPP_", ifelse(CRISPR_df_expanded[["GPP_rank"]] %in% 1:10, "top10", "rest")),
                                                   sources_vec
                                                   )
                                            )
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

  are_valid_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list, show_messages = FALSE)

  top4_df <- data.frame(
    CRISPR_df[are_valid_4sg, ],
    "Group"    = "4sg",
    "Subgroup" = "None",
    stringsAsFactors = FALSE
  )
  assign("delete_m_top4_df", top4_df, envir = globalenv())
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
  is_CRISPRko <- "Entrez_source_Brunello" %in% names(expanded_CRISPR_df)
  if (is_CRISPRko && !("Entrez_source_Brunello" %in% names(expanded_CRISPR_df))) {
    stop("No identifying column names were found in expanded_CRISPR_df!")
  }

  ## Filter for valid Entrez IDs
  have_entrez <- (expanded_CRISPR_df[["Entrez_ID"]] == expanded_CRISPR_df[["Combined_ID"]]) %in% TRUE
  are_protein_coding <- expanded_CRISPR_df[["Entrez_ID"]] %in% collected_entrez_IDs
  are_eligible <- have_entrez & are_protein_coding

  if (!(is_CRISPRko) && show_sublibraries) {
    are_eligible <- are_eligible & !(grepl("_6to10", expanded_CRISPR_df[["Subgroup"]], fixed = TRUE)) # Otherwise, a tiny group of "supp 5" genes is displayed in the plot
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

    ## Choose the guides for the John Doench libraries
    Doench_df <- filtered_df[filtered_df[["Group"]] %in% c("Dolcetto", "Calabrese"), ]
    Doench_df <- Doench_df[order(as.integer(Doench_df[["Entrez_ID"]])), ]

    if (!(filter_complete_genes)) {
      Doench_df <- FilterCompleteTop4(Doench_df)
    }

    Doench_entrezs_fac <- factor(Doench_df[["Entrez_ID"]], levels = unique(Doench_df[["Entrez_ID"]]))

    CheckThatFactorIsInOrder(Doench_entrezs_fac)

    rank_column <- grep("^(Calabrese|Dolcetto)_rank", names(Doench_df))
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
    Doench_chosen_df <- Doench_df[Doench_are_chosen, ]

    ## Choose the guides for the hCRISPR-v2 libraries
    hCRISPR_v2_df <- filtered_df[filtered_df[["Group"]] %in% c("hCRISPRa-v2", "hCRISPRi-v2.1", "hCRISPRi-v2.0"), ]

    hCRISPR_v2_df <- hCRISPR_v2_df[order(match(hCRISPR_v2_df[["Combined_ID"]], hCRISPR_v2_df[["Combined_ID"]])), ]
    hCRISPR_v2_combined_IDs <- factor(hCRISPR_v2_df[["Combined_ID"]], levels = unique(hCRISPR_v2_df[["Combined_ID"]]))
    rank_column <- grep("_v2_rank", names(hCRISPR_v2_df), fixed = TRUE)
    hCRISPR_v2_rank_list <- split(as.integer(hCRISPR_v2_df[[rank_column]]), hCRISPR_v2_combined_IDs)
    are_top4 <- unlist(lapply(hCRISPR_v2_rank_list, function(x) {
      top_4_ranks <- sort(unique(x))[1:4]
      x %in% top_4_ranks
    }))
    hCRISPR_v2_chosen_df <- hCRISPR_v2_df[are_top4, ]
    hCRISPR_v2_chosen_df <- ChooseOriginalTop4(hCRISPR_v2_chosen_df)


    ## Combine the final data frame for CRISPRa
    if (filter_complete_genes) {
      stopifnot(length(unique(c(nrow(Doench_chosen_df), nrow(hCRISPR_v2_chosen_df), nrow(GPP_chosen_df), nrow(FourSg_df)))) == 1)
    }

    combined_chosen_df <- rbind.data.frame(
      Doench_chosen_df,
      hCRISPR_v2_chosen_df,
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

  if ("Original_order" %in% names(CRISPR_df)) {
    new_order <- order(match(CRISPR_df[["Entrez_ID"]], CRISPR_df[["Entrez_ID"]]),
                       CRISPR_df[["Original_order"]]
                       )
    CRISPR_df <- CRISPR_df[new_order, ]
  }

  num_occurrences_vec <- table(CRISPR_df[["Entrez_ID"]])[CRISPR_df[["Entrez_ID"]]]
  fourguides_df <- CRISPR_df[num_occurrences_vec == 4, ]
  other_df <- CRISPR_df[num_occurrences_vec > 4, ]

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





SourcesSubgroupLabels <- function(plot_df,
                                  x_positions,
                                  colors_df,
                                  show_sublibraries = TRUE,
                                  show_rest_v_4sg   = FALSE,
                                  extra_space       = FALSE
                                  ) {

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
      second_line_factor <- 0.09
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
    assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
    assign("delete_y_column", y_column, envir = globalenv())
    CRISPR_df[[y_column]] <- ifelse(is.na(CRISPR_df[[y_column]]) & !(is.na(CRISPR_df[["Start"]])),
                                    0,
                                    CRISPR_df[[y_column]]
                                    ) * 100
  } else if (y_column == "Deviation_from_TSS_window") {
    distance_vec <- CRISPR_df[["Distance_from_TSS"]]
    is_CRISPRa <- IsCRISPRa(CRISPR_df)
    TSS_window <- GetTSSWindow(is_CRISPRa)
    are_below_lower_bound <- distance_vec < TSS_window[[1]]
    are_above_higher_bound <- distance_vec > TSS_window[[2]]
    deviation_vec <- ifelse(!(are_below_lower_bound) & !(are_above_higher_bound),
                            0L,
                            ifelse(are_above_higher_bound,
                                   abs(TSS_window[[2]] - distance_vec),
                                   abs(TSS_window[[1]] - distance_vec)
                                   )
                            )
    deviation_vec <- ifelse(deviation_vec > 1000, 1000L, deviation_vec)
    CRISPR_df[["Deviation_from_TSS_window"]] <- deviation_vec
  } else if (y_column == "Restricted_distance_from_TSS") {
    distance_vec <- CRISPR_df[["Distance_from_TSS"]]
    distance_vec <- ifelse(distance_vec < (-1000), -1000L, distance_vec)
    distance_vec <- ifelse(distance_vec > 1000, 1000L, distance_vec)
    CRISPR_df[["Restricted_distance_from_TSS"]] <- distance_vec
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
                                  is_CRISPRa,
                                  show_title            = TRUE,
                                  aggregate_scores      = aggregate_scores,
                                  title_cex             = 0.9,
                                  title_line            = 1.3,
                                  no_outside_annotation = FALSE,
                                  use_raster_array      = NULL
                                  ) {
  tick_locations <- axTicks(2)
  if (is.null(use_raster_array)) {
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
  }

  if (!(no_outside_annotation)) {
    tick_labels <- format(tick_locations)
    if ((y_column == "Deviation_from_TSS_window") &&
        (tick_labels[[length(tick_labels)]] == "1000")
        ) {
      tick_labels <- sapply(tick_labels, as.expression)
      tick_labels[[length(tick_labels)]] <- bquote("" >= 1000)
    } else if ((y_column == "Restricted_distance_from_TSS") &&
               (tick_labels[[length(tick_labels)]] == " 1000") &&
               (tick_labels[[1]] == "-1000")) {
      tick_labels <- sapply(tick_labels, as.expression)
      tick_labels[[length(tick_labels)]] <- bquote("" >= 1000)
      tick_labels[[1]] <- bquote("" <= "\u22121000")
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
  }

  title_text <- numeric_column_labels[[y_column]]
  if (aggregate_scores) {
    title_text <- paste0("Aggregate ", title_text)
  } else if (y_column == "Deviation_from_TSS_window") {
    title_text <- paste0(title_text, TSSWindowString(is_CRISPRa))
  }
  if (show_title && !(no_outside_annotation)) {
    title(title_text, cex.main = title_cex, line = title_line)
  }
  return(invisible(title_text))
}






PlotViolin <- function(plot_df,
                       x_positions,
                       x_limits,
                       colors_df,
                       y_column_name,
                       is_CRISPRa,
                       show_title            = TRUE,
                       large_count_labels    = FALSE,
                       aggregate_scores      = FALSE,
                       title_cex             = 0.9,
                       title_line            = 1.3,
                       no_outside_annotation = FALSE,
                       use_raster_array      = NULL,
                       point_cex             = 0.5
                       ) {

  stopifnot(all(c("Numeric_data", "Groups_factor", "Point_colors") %in% names(plot_df)))
  stopifnot(all(c("Pale", "Medium", "Dark") %in% names(colors_df)))

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
                                      is_CRISPRa,
                                      show_title            = show_title,
                                      aggregate_scores      = aggregate_scores,
                                      title_cex             = title_cex,
                                      title_line            = title_line,
                                      no_outside_annotation = no_outside_annotation,
                                      use_raster_array      = use_raster_array
                                      )

  if (is.null(use_raster_array)) {

    set.seed(1)
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
           cex = point_cex * par("cex"),
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
  } else {
    rasterImage(use_raster_array,
                xleft = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                )
  }

  if (!(no_outside_annotation)) {
    box(xpd = NA, lwd = 0.75)
  }

  ## Draw the subgroup numbers

  if (!(no_outside_annotation)) {
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
  }
  return(invisible(title_text))
}



ViolinBox_UniqueTwoGroups <- function(CRISPR_df,
                                      y_column,
                                      show_title            = TRUE,
                                      no_outside_annotation = FALSE,
                                      use_raster_array      = NULL,
                                      point_cex             = 0.5
                                      ) {

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
  if (no_outside_annotation) {
    old_mar <- par(mai = rep(0, 4))
  } else {
    old_mar <- par(mai = (c(1, 0.8, 0.66, 0.04) + 0.02) * 0.66)
  }

  title_text <- PlotViolin(plot_df, x_positions, x_limits, colors_df,
                           y_column, IsCRISPRa(CRISPR_df),
                           show_title = show_title,
                           no_outside_annotation = no_outside_annotation,
                           use_raster_array = use_raster_array,
                           point_cex = point_cex
                           )

  ## Draw the group labels
  if (!(no_outside_annotation)) {
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
  }

  ## Final steps
  par(old_mar)
  results_list <- list("plot_df" = plot_df, "title_text" = title_text)
  return(invisible(results_list))
}






ViolinBox_UniqueLibraries <- function(CRISPR_df,
                                      y_column,
                                      show_title            = TRUE,
                                      no_outside_annotation = FALSE,
                                      use_raster_array      = NULL,
                                      point_cex             = 0.5
                                      ) {

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
  if (no_outside_annotation) {
    old_mar <- par(mai = rep(0, 4))
  } else {
    old_mar <- par(mai = (c(1, 0.8, 0.66, 0.6) + 0.02) * 0.66)
  }

  title_text <- PlotViolin(plot_df, x_positions, x_limits, colors_df,
                           y_column, IsCRISPRa(CRISPR_df),
                           show_title = show_title,
                           no_outside_annotation = no_outside_annotation,
                           use_raster_array = use_raster_array,
                           point_cex = point_cex
                           )


  ## Draw the group labels
  if (!(no_outside_annotation)) {
    UniqueSequencesSubgroupLabels(plot_df, x_positions, colors_df, extra_space = TRUE)
  }

  ## Final steps
  par(old_mar)
  results_list <- list("plot_df" = plot_df, "title_text" = title_text)
  return(invisible(results_list))
}




UniquePointsBoxPlots <- function(CRISPR_df, embed_raster_within_PDFs = TRUE) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    numeric_column_labels <- numeric_column_labels[!(names(numeric_column_labels) %in% TSS_columns)]
  }

  main_folder_path <- file.path(output_plots_directory, "Box plots", "Box plots - A) unique points")

  for (make_PDF in c(FALSE, TRUE)) {
    if (make_PDF) {
      pdf(file = file.path(main_folder_path, "Box plots - A) unique points.pdf"),
          width = pdf_width, height = pdf_height
          )
    }
    for (y_column in names(numeric_column_labels)) {

      if (make_PDF && embed_raster_within_PDFs) {
        PDF_device <- dev.cur()
        temp_left_path <- file.path(main_folder_path, "temp_left.png")
        temp_right_path <- file.path(main_folder_path, "temp_right.png")

        left_width_fraction <- layout_widths[[1]] / sum(layout_widths)
        right_width_fraction <- layout_widths[[2]] / sum(layout_widths)
        height_fraction  <- layout_heights[[2]] / sum(layout_heights)

        temp_left_width <- pdf_width * left_width_fraction - (0.88 * 0.66)
        temp_right_width <- pdf_width * right_width_fraction - (1.44 * 0.66)
        temp_height <- pdf_height * height_fraction - (1.7 * 0.66)

        png(file = temp_left_path,
            width = temp_left_width, height = temp_height,
            units = "in", res = 600
            )
        ViolinBox_UniqueTwoGroups(CRISPR_df, y_column, show_title = FALSE,
                                  no_outside_annotation = TRUE,
                                  point_cex = 0.4
                                  )
        dev.off()

        png(file = temp_right_path,
            width = temp_right_width, height = temp_height,
            units = "in", res = 600
            )
        ViolinBox_UniqueLibraries(CRISPR_df, y_column, show_title = FALSE,
                                  no_outside_annotation = TRUE,
                                  point_cex = 0.4
                                  )
        dev.off()

        left_raster_array <- readPNG(temp_left_path)
        right_raster_array <- readPNG(temp_right_path)

        file.remove(temp_left_path)
        file.remove(temp_right_path)

        dev.set(PDF_device)

        layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
        ViolinBox_UniqueTwoGroups(CRISPR_df, y_column, show_title = FALSE,
                                  use_raster_array = left_raster_array
                                  )
        ViolinBox_results <- ViolinBox_UniqueLibraries(CRISPR_df, y_column, show_title = FALSE,
                                                       use_raster_array = right_raster_array
                                                       )
        OuterTitleForLayout(ViolinBox_results[["title_text"]], extra_space = TRUE)
        layout(1)

      } else {
        layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
        ViolinBox_UniqueTwoGroups(CRISPR_df, y_column, show_title = FALSE)
        ViolinBox_results <- ViolinBox_UniqueLibraries(CRISPR_df, y_column, show_title = FALSE)
        OuterTitleForLayout(ViolinBox_results[["title_text"]], extra_space = TRUE)
        layout(1)
      }
    }
    if (make_PDF) {
      dev.off()
    }
  }
  numeric_seq <- seq_along(numeric_column_labels)
  file_numbers <- FormatFixedWidthInteger(numeric_seq)
  for (i in seq_along(numeric_column_labels)) {
    file_name <- paste0(file_numbers[[i]], ") ", names(numeric_column_labels)[[i]])
    png(file = file.path(main_folder_path, "PNGs", paste0("Box plots - ", file_name, ".png")),
        width = pdf_width, height = pdf_height - 0.14, units = "in", res = 600
        )
    layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
    ViolinBox_UniqueTwoGroups(CRISPR_df, names(numeric_column_labels)[[i]], show_title = FALSE)
    ViolinBox_results <- ViolinBox_UniqueLibraries(CRISPR_df, names(numeric_column_labels)[[i]], show_title = FALSE)
    OuterTitleForLayout(ViolinBox_results[["title_text"]])
    dev.off()
  }
}






ViolinBox_Sources <- function(CRISPR_df,
                              y_column,
                              show_rest_v_4sg        = FALSE,
                              aggregate_scores       = FALSE,
                              show_sublibraries      = !(show_rest_v_4sg || aggregate_scores),
                              filter_top4            = aggregate_scores,
                              filter_complete_genes  = TRUE,
                              filter_complete_scores = TRUE,
                              collapse_GPP           = filter_top4,
                              no_outside_annotation  = FALSE,
                              use_raster_array       = NULL,
                              draw_plot              = FALSE
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
    is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df)
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
  colors_df <- FilterOriginalColorsDf(SubgroupColorsDf(), plot_df, show_rest_v_4sg = show_rest_v_4sg)

  plot_df_columns <- c("Combined_ID", "Entrez_ID", "Gene_symbol", "Group", "Subgroup")
  if (!(aggregate_scores)) {
    plot_df_columns <- c(plot_df_columns, "sgRNA_sequence")
    if (y_column %in% c("Deviation_from_TSS_window", "Restricted_distance_from_TSS")) {
      plot_df_columns <- c(plot_df_columns, "Distance_from_TSS")
    }
  }
  plot_df_columns <- c(plot_df_columns, y_column)
  plot_df <- plot_df[, plot_df_columns]

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

  if (draw_plot) {
    ## Determine the spacing and the x positions
    spaces_vec <- ifelse(colors_df[["Are_new_color"]][-1], 1.15, 0.72)
    x_positions <- cumsum(c(1, spaces_vec))

    num_subgroups <- nlevels(plot_df[["Groups_factor"]])
    x_limits <- c((x_positions[[1]] - 0.5) - (num_subgroups * 0.04),
                  x_positions[[length(x_positions)]] + 0.5 + (num_subgroups * 0.04)
                  )

    ## Draw the violins
    if (no_outside_annotation) {
      old_mar <- par(mar = rep(0, 4))
    } else {
      old_mar <- par(mar = c(4.4, 4, 3.5, 3) + 0.1)
    }
    PlotViolin(plot_df, x_positions, x_limits, colors_df, y_column,
               IsCRISPRa(CRISPR_df),
               aggregate_scores = aggregate_scores,
               title_cex = 0.8, title_line = 1.55,
               no_outside_annotation = no_outside_annotation,
               use_raster_array = use_raster_array
               )

    if (!(no_outside_annotation)) {
      SourcesSubgroupLabels(plot_df, x_positions, colors_df,
                            show_sublibraries = show_sublibraries,
                            show_rest_v_4sg = show_rest_v_4sg,
                            extra_space = TRUE
                            )
    }
    par(old_mar)
  }

  return(invisible(plot_df))
}







SourcesBoxPlots <- function(CRISPR_df, embed_raster_within_PDFs = TRUE) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    args_list <- c(args_list, CRISPRko_args_list)
    numeric_column_labels <- numeric_column_labels[!(names(numeric_column_labels) %in% TSS_columns)]
  } else {
    args_list <- c(args_list, CRISPRa_args_list)
  }
  use_width <- pdf_width * 0.85
  use_height <- pdf_height * 1.2

  export_raw_data_columns <- c(
    "CRISPOR_Doench_efficacy", "GuideScan_efficiency",
    "GuideScan_specificity", "CRISPOR_3MM_specificity",
    "CRISPOR_4MM_specificity"
  )
  raw_data_exported_vec_list <- sapply(names(args_list)[1:5],
                                       function(x) sapply(export_raw_data_columns, function(y) FALSE),
                                       simplify = FALSE
                                       )

  for (make_PNG in c(FALSE, TRUE)) {
    for (make_PDF in c(FALSE, TRUE)) {
      if (make_PNG && make_PDF) {
        next
      }
      for (file_name in names(args_list)) {
        only_four_groups <- file_name %in% names(args_list)[1:5]
        main_folder_path <- file.path(output_plots_directory,
                                      "Box plots",
                                      paste0("Box plots - ", file_name)
                                      )
        PNG_folder_path <- file.path(main_folder_path, "PNGs")
        all_folder_paths <- c(main_folder_path, PNG_folder_path)
        if (only_four_groups) {
          raw_data_path <- file.path(main_folder_path, "Raw data")
          all_folder_paths <- c(all_folder_paths, raw_data_path)
        }
        for (folder_path in all_folder_paths){
          if (!(dir.exists(folder_path))) {
            dir.create(folder_path)
          }
        }

        if (only_four_groups) { # Make the plot a bit less wide, in order to compensate for only 4 groups
          current_width <- use_width - 0.9
          current_height <- use_height + 0.55
        } else {
          current_width <- use_width
          current_height <- use_height
        }

        if (make_PDF) {
          PDF_file_name <- paste0("Box plots - ", file_name, ".pdf")

          pdf(file = file.path(main_folder_path, PDF_file_name),
              width = current_width, height = current_height
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

          if (make_PDF && embed_raster_within_PDFs) {
            PDF_device <- dev.cur()
            temp_path <- file.path(main_folder_path, "temp.png")
            temp_width <- current_width - 1.44 # subtract the margin
            temp_height <- current_height - 1.62 # subtract the margin

            png(file   = temp_path,
                width  = temp_width,
                height = temp_height,
                units  = "in",
                res    = 600
                )
            plot_df <- do.call(ViolinBox_Sources,
                               c(list(CRISPR_df             = CRISPR_df,
                                      y_column              = y_column,
                                      no_outside_annotation = TRUE
                                      ),
                                 args_list[[file_name]]
                                 )
                               )
            dev.off()
            raster_array <- readPNG(temp_path)
            file.remove(temp_path)
            dev.set(PDF_device)
            do.call(ViolinBox_Sources,
                    c(list(CRISPR_df        = CRISPR_df,
                           y_column         = y_column,
                           use_raster_array = raster_array
                           ),
                      args_list[[file_name]]
                      )
                    )
          } else {
            if (make_PNG) {
              folder_path <- file.path(output_plots_directory, paste0("Box plots - ", file_name))
              if (!(dir.exists(folder_path))) {
                dir.create(folder_path)
              }
              PNG_file_name <- paste0("Box plots - ", file_name, " - ",
                                      file_numbers[[i]], ") ", y_column, ".png"
                                      )
              png(file = file.path(PNG_folder_path, PNG_file_name),
                  width = current_width, height = current_height,
                  units = "in", res = 600
                  )
            }
            plot_df <- do.call(ViolinBox_Sources, c(list(y_column = y_column, CRISPR_df = CRISPR_df), args_list[[file_name]]))
            if (make_PNG) {
              dev.off()
            }
          }
          if (only_four_groups && (y_column %in% export_raw_data_columns) &&
              !(raw_data_exported_vec_list[[file_name]][[y_column]])
              ) {
            plot_df <- TidyRawData(plot_df)
            raw_data_file_name <- paste0("Box plots - raw data - ", file_name, " - ",
                                         file_numbers[[i]], ") ", y_column,
                                         ".csv"
                                         )
            write.csv(plot_df,
                      file = file.path(raw_data_path, raw_data_file_name),
                      quote = FALSE, row.names = FALSE
                      )
            raw_data_exported_vec_list[[file_name]][[y_column]] <- TRUE
          }
        }
        if (make_PDF) {
          dev.off()
        }
      }
    }
  }
}




TidyRawData <- function(plot_df) {
  if (identical(plot_df[["Combined_ID"]], plot_df[["Entrez_ID"]])) {
    plot_df <- plot_df[, names(plot_df) != "Combined_ID"]
  }
  results_df <- plot_df[, !(duplicated(as.list(plot_df)))]
  if (length(unique(results_df[["Subgroup"]])) == 1) {
    results_df[["Subgroup"]] <- NULL
  }
  results_df <- results_df[, !(names(results_df) %in% "Point_colors")]
  return(results_df)
}




# Functions for plotting categorical data ---------------------------------

MakeCategoricalList <- function(CRISPR_df, use_column, use_cutoff) {
  assign("delete_use_column", use_column, envir = globalenv())
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  assign("delete_use_cutoff", use_cutoff, envir = globalenv())
  title_cex <- 1
  if (use_column %in% names(unintended_target_labels)) {
    logical_vec <- CRISPR_df[[use_column]] > 0
    use_title <- unintended_target_labels[[use_column]]
  } else if (use_column == "Are_overlapping") {
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
    is_CRISPRa <- IsCRISPRa(CRISPR_df)
    use_title <- paste0("Lie outside the optimal window around the TSS",
                        TSSWindowString(is_CRISPRa)
                        )
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
  } else if (use_column %in% grep("_SNP_AF_(max|sum)_", names(CRISPR_df), value = TRUE)) {
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
                        ") according to the Rule Set 2 score (via ",
                        via_string, ")"
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


RoundSmallPercentages <- function(numeric_vec) {
  if (all(numeric_vec < 1)) {
    results_vec <- ifelse(numeric_vec < 0.1,
                          signif(numeric_vec, digits = 1),
                          signif(numeric_vec, digits = 2)
                          )
  } else if (all(numeric_vec < 10)) {
    results_vec <- ifelse(numeric_vec < 0.1,
                          signif(numeric_vec, digits = 1),
                          round(numeric_vec, digits = 2)
                          )
  } else {
    results_vec <- round(numeric_vec, digits = 1)
  }
  return(results_vec)
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

  old_scipen <- options(scipen = 1)

  ## Annotate the bars with percentages
  if (nrow(proportions_mat) == 1) {
    percent_vec <- proportions_mat[1, ] * 100
    proportions_vec <- format(RoundSmallPercentages(percent_vec))
  } else {
    percent_vec <- proportions_mat[2, ] * 100
    proportions_vec <- ifelse(percent_vec < 1,
                              signif(percent_vec, digits = 1),
                              ifelse(percent_vec > 99.5,
                                     signif(percent_vec, digits = 3),
                                     signif(percent_vec, digits = 2)
                                     )
                              )
  }

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

  if (nrow(proportions_mat) == 1) {
    dividends <- RoundSmallPercentages(counts_mat[1, ])
    divisors <- counts_mat[2, ]
  } else {
    dividends <- counts_mat[2, ]
    divisors <-  colSums(counts_mat)
  }

  ## Annotate the bars with fractions (counts)
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * fractions_plot_factor),
       labels = divisors,
       adj    = c(0.5, 0.5),
       cex    = text_cex,
       col    = "black",
       xpd    = NA
       )
  line_widths <- strwidth(divisors, cex = text_cex)
  segments(x0  = bar_positions_vec - (line_widths / 2),
           x1  = bar_positions_vec + (line_widths / 2),
           y0  = par("usr")[[4]] + (plot_height * (fractions_plot_factor + 0.0125)),
           col = "black",
           xpd = NA,
           lwd = 0.3
           )
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * (fractions_plot_factor + 2 * 0.0125)),
       labels = dividends,
       adj    = c(0.5, 0.5),
       cex    = text_cex,
       col    = "black",
       xpd    = NA
       )
  return(invisible(NULL))
}




PlotBars <- function(bars_mat,
                     spaces_vec,
                     x_limits,
                     colors_df,
                     create_plot = TRUE,
                     use_y_limits = NULL
                     ) {

  if (nrow(bars_mat) == 1) {
    if (is.null(use_y_limits)) {
      use_y_limits <- c(0, max(bars_mat[1, ] * 1.04))
    }
    y_axis_pos <- pretty(use_y_limits)
    y_limits <- c(y_axis_pos[[1]], y_axis_pos[[length(y_axis_pos)]])
    first_colors <- colors_df[["Dark"]]
  } else {
    y_limits <- c(0, 1)
    first_colors <- colors_df[["Pale"]]
  }

  if (create_plot) {
    plot(1,
         xlim = x_limits,
         ylim = y_limits,
         xaxs = "i",
         yaxs = "i",
         axes = FALSE,
         ann  = FALSE,
         type = "n"
         )

    if (nrow(bars_mat) == 1) {
      segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]], y0 = y_axis_pos,
               col = "gray95", lwd = 0.5
               )
      if (all(bars_mat[1, ] < 1)) {
        tick_labels <- paste0(y_axis_pos * 100, "%")
      } else {
        tick_labels <- y_axis_pos
      }
      axis(side      = 2,
           at        = y_axis_pos,
           labels    = tick_labels,
           lwd       = 0.75,
           las       = 1,
           mgp       = c(3, 0.45, 0),
           tcl       = -0.3,
           cex.axis  = 0.8
           )
    } else {
      DrawBarplotGridAndAxis()
    }
  }

  ## Draw the bar plots()
  bar_positions <- barplot(height    = bars_mat[1, ],
                           col       = first_colors,
                           tcl       = -0.4,
                           border    = NA,
                           space     = spaces_vec,
                           names.arg = rep.int("", nrow(colors_df)),
                           axes      = FALSE,
                           add       = TRUE,
                           plot      = create_plot
                           )[, 1]

  if (create_plot && (nrow(bars_mat) > 1)) {
    barplot(height    = bars_mat[2, ],
            offset    = bars_mat[1, ],
            col       = colors_df[["Dark"]],
            border    = NA,
            space     = spaces_vec,
            names.arg = rep.int("", nrow(colors_df)),
            axes      = FALSE,
            add       = TRUE
            )
  }

  results_list <- list(
    "bar_positions" = bar_positions,
    "y_limits"      = y_limits
  )
  return(invisible(results_list))
}





BarPlot_UniqueTwoGroups <- function(CRISPR_df,
                                    use_column,
                                    only_chosen  = FALSE,
                                    use_cutoff   = NULL,
                                    create_plot  = TRUE,
                                    show_title   = create_plot,
                                    use_y_limits = NULL
                                    ) {

  sum_up_SNPs <- grepl("^Expected_.+_SNP_AF_(sum|max)_", use_column)
  if (sum_up_SNPs) {
    use_column <- sub("^Expected_", "", use_column)
  }

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), use_column)

  ## Prepare the data
  if (!(sum_up_SNPs)) {
    categorical_list <- MakeCategoricalList(CRISPR_df, use_column, use_cutoff)
  }
  colors_df <- SubgroupColorsDf()
  colors_df <- colors_df[(colors_df[["Color_name"]] == "Grey") & !(colors_df[["Lightened"]]), ]

  if (only_chosen) {
    are_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list)
    CRISPR_df <- CRISPR_df[are_4sg, ]
    if (!(sum_up_SNPs)) {
      categorical_list[["logical_vec"]] <- categorical_list[["logical_vec"]][are_4sg]
    }
    num_subgroups <- 1L
    x_spaces <- NULL
  } else {
    num_subgroups <- 2L
    x_spaces <- 0.5
  }

  plot_df <- data.frame(
    "Chosen" = Are4sg(CRISPR_df, sublibraries_all_entrezs_list),
    stringsAsFactors = FALSE
  )
  if (sum_up_SNPs) {
    plot_df[[use_column]] <- CRISPR_df[[use_column]]
    counts_mat <- SumUpSNPs(plot_df, use_column, "Chosen")
    proportions_mat <- t(counts_mat[1, ] / counts_mat[2, ])
  } else {
    plot_df[["Data"]] <- categorical_list[["logical_vec"]]
    counts_mat <- as.matrix(table(factor(plot_df[["Data"]], levels = c(FALSE, TRUE)), plot_df[["Chosen"]]))
    proportions_mat <- prop.table(counts_mat, margin = 2)
  }

  ## Prepare the plot region
  if (only_chosen) {
    x_limits <- c(-0.34, 1.74)
  } else {
    x_limits <- c(0.1, 3.4)
  }
  old_mar <- par(mar = c(4, 4, 5, 0.2) + 0.1)

  if (sum_up_SNPs) {
    title_text <- aggregate_sum_column_labels[[use_column]]
    title_cex <- 0.9
  } else {
    title_text <- categorical_list[["title_text"]]
    title_cex <- categorical_list[["title_cex"]]
  }

  if (show_title) {
    title(title_text, cex.main = title_cex, line = 3.625)
  }
  bar_results <- PlotBars(proportions_mat, x_spaces, x_limits, colors_df,
                          create_plot = create_plot, use_y_limits = use_y_limits
                          )

  if (create_plot) {
    ## Draw the group labels
    plot_height <- par("usr")[[4]] - par("usr")[[3]]
    old_lheight <- par("lheight" = 1.1)
    if (only_chosen) {
      group_labels <- "All chosen guides"
    } else {
      group_labels <- c("Rest", "4sg")
    }
    text(x      = bar_results[["bar_positions"]],
         y      = par("usr")[[3]] - (plot_height * 0.055),
         labels = group_labels,
         adj    = c(0.5, 1),
         cex    = 0.9,
         col    =  "black",
         font   = 2,
         xpd    = NA
         )
    par(old_lheight)
    AnnotateBars(bar_results[["bar_positions"]], counts_mat, proportions_mat, text_cex = 0.5)

    ## Final steps
    par(old_mar)
  }

  results_list <- list("plot_df"    = plot_df,
                       "title_text" = title_text,
                       "title_cex"  = title_cex,
                       "y_limits"   = bar_results[["y_limits"]]
                       )
  return(invisible(results_list))
}



SumUpSNPs <- function(plot_df, SNP_column, group_column) {
  sums_vec <- tapply(plot_df[[SNP_column]], plot_df[[group_column]], sum, na.rm = TRUE) / 100
  num_per_group_vec <- tapply(plot_df[[SNP_column]],
                              plot_df[[group_column]],
                              function(x) sum(!(is.na(x)))
                              )
  results_mat <- rbind(sums_vec, num_per_group_vec)
  return(results_mat)
}



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

  if (show_sublibraries && show_rest_v_4sg) {
    stop("The 'show_sublibraries' and 'show_rest_v_4sg' arguments may not both be TRUE!")
  }
  sum_up_SNPs <- grepl("^Expected_.+_SNP_AF_(sum|max)_", use_column)
  if (sum_up_SNPs) {
    use_column <- sub("^Expected_", "", use_column)
  }
  if (!(filter_top4) && (use_column == "Are_overlapping")) {
    warning("It makes little sense to graph the fraction of overlapping guides across libraries unless the top 4 guides are selected!")
  }
  if (!(filter_top4) && (use_column == "Have_homologies")) {
    stop("The 'filter_top4' argument must be TRUE if the 'Have_homologies' parameter is selected!")
  }

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), use_column)

  if (filter_top4) {
    is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df)
    if (!(is_CRISPRko)) {
      are_main_transcript <- (GetMainTSS(CRISPR_df) == CRISPR_df[["AltTSS_ID"]]) %in% TRUE
      CRISPR_df <- CRISPR_df[are_main_transcript, ]
    }
  }

  plot_df <- ExpandedSubgroupsDf(CRISPR_df,
                                 collapse_GPP = collapse_GPP,
                                 show_rest_v_4sg = show_rest_v_4sg
                                 )

  if (use_column == "Have_homologies") {
    categorical_list <- list(
      "logical_vec" = NULL,
      "use_column"  = use_column,
      "title_text"  = "Top 4 guides share identical sub-sequences (8 base pairs or longer)",
      "title_cex"   = 0.9
    )
  } else if (!(sum_up_SNPs)) {
    categorical_list <- MakeCategoricalList(plot_df, use_column, use_cutoff)
  }

  if (!(sum_up_SNPs)) {
    new_column <- categorical_list[["use_column"]]
    plot_df[[new_column]] <- categorical_list[["logical_vec"]]
  } else {
    new_column <- use_column
  }

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

  if (sum_up_SNPs) {
    counts_mat <- SumUpSNPs(plot_df, use_column, "Groups_factor")
    bars_mat <- t(counts_mat[1, ] / counts_mat[2, ])
  } else if (new_column %in% categorical_4sg_columns) {
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

  if (!(sum_up_SNPs)) {
    counts_mat <- as.matrix(table(factor(plot_df[[new_column]], levels = c(FALSE, TRUE)), plot_df[["Groups_factor"]]))
    bars_mat <- prop.table(counts_mat, margin = 2)
  }

  spaces_vec <- ifelse(colors_df[["Are_new_color"]], 1.3, 0.3)
  x_limits <- BarPlotXlimits(spaces_vec, side_space = max(spaces_vec))
  old_mar <- par(mar = c(4, 4, 6, 3) + 0.1)

  bar_positions <- PlotBars(bars_mat, spaces_vec, x_limits, colors_df)[["bar_positions"]]

  if (sum_up_SNPs) {
    title_text <- aggregate_sum_column_labels[[use_column]]
    title_cex <- 0.9
  } else {
    title_text <- categorical_list[["title_text"]]
    title_cex <- categorical_list[["title_cex"]]
  }
  title(title_text, cex.main = title_cex, line = 3.625)

  SourcesSubgroupLabels(plot_df, bar_positions, colors_df,
                        show_sublibraries = show_sublibraries,
                        show_rest_v_4sg = show_rest_v_4sg
                        )

  AnnotateBars(bar_positions, counts_mat, bars_mat, text_cex = 0.4, smaller_plot_factor  = TRUE)

  colnames(counts_mat) <- colors_df[, "Group_label"]
  results_list <- list(
    "plot_df" = plot_df,
    "counts_mat" = counts_mat
  )

  ## Final steps
  par(old_mar)
  return(invisible(results_list))
}




BarPlot_UniqueLibraries <- function(CRISPR_df,
                                    use_column,
                                    use_cutoff   = NULL,
                                    show_title   = TRUE,
                                    only_chosen  = FALSE,
                                    create_plot  = TRUE,
                                    use_y_limits = NULL
                                    ) {


  sum_up_SNPs <- grepl("^Expected_.+_SNP_AF_(sum|max)_", use_column)
  if (sum_up_SNPs) {
    use_column <- sub("^Expected_", "", use_column)
  }

  CRISPR_df <- FixNumericColumns(FilterCRISPRDf(CRISPR_df), use_column)

  ## Prepare the data
  if (!(sum_up_SNPs)) {
    categorical_list <- MakeCategoricalList(CRISPR_df, use_column, use_cutoff)
    CRISPR_df[[categorical_list[["use_column"]]]] <- categorical_list[["logical_vec"]]
    use_column <- categorical_list[["use_column"]]
  }

  if (only_chosen) {
    are_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list)
    CRISPR_df <- CRISPR_df[are_4sg, ]
  }
  plot_df <- MutuallyExclusiveSubgroupsDf(CRISPR_df, use_column)
  colors_df <- SubgroupColorsDf()
  if (only_chosen) {
    colors_df <- colors_df[!(colors_df[["Lightened"]]), -2]
    group_column <- "Group"
  } else {
    group_column <- "Subgroup"
  }
  if (sum_up_SNPs) {
    counts_mat <- SumUpSNPs(plot_df, "Data", group_column)
    proportions_mat <- t(counts_mat[1, ] / counts_mat[2, ])
  } else {
    counts_mat <- as.matrix(table(factor(plot_df[["Data"]], levels = c(FALSE, TRUE)), plot_df[[group_column]]))
    proportions_mat <- prop.table(counts_mat, margin = 2)
  }
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

  bar_results <- PlotBars(proportions_mat, spaces_vec, x_limits, colors_df,
                          create_plot = create_plot, use_y_limits = use_y_limits
                          )
  bar_positions <- bar_results[["bar_positions"]]

  if (sum_up_SNPs) {
    title_text <- aggregate_sum_column_labels[[use_column]]
    title_cex <- 0.9
  } else {
    title_text <- categorical_list[["title_text"]]
    title_cex <- categorical_list[["title_cex"]]
  }
  if (show_title) {
    title(title_text, cex.main = title_cex, line = 3.625)
  }

  if (create_plot) {
    UniqueSequencesSubgroupLabels(plot_df, bar_positions, colors_df, only_chosen = only_chosen)
    AnnotateBars(bar_positions, counts_mat, proportions_mat, text_cex = 0.4)
  }

  ## Final steps
  par(old_mar)
  results_list <- list("plot_df"    = plot_df,
                       "title_text" = title_text,
                       "title_cex"  = title_cex,
                       "y_limits"   = bar_results[["y_limits"]]
                       )
  return(invisible(results_list))
}






UniqueSequencesBarPlots <- function(CRISPR_df) {

  is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    categorical_columns <- setdiff(categorical_columns, TSS_columns)
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
          pdf(file = file.path(output_plots_directory, "Bar charts",
                               paste0(folder_name, ".pdf")
                               ),
              width = use_width, height = use_height
              )
        }
        categorical_seq <- seq_along(use_columns)
        file_numbers <- FormatFixedWidthInteger(categorical_seq)
        for (i in categorical_seq) {
          use_column <- use_columns[[i]]
          if (grepl("^Expected_", use_column)) {
            assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
            assign("delete_use_column", use_column, envir = globalenv())
            assign("delete_only_chosen", only_chosen, envir = globalenv())
            TwoGroups_results <- BarPlot_UniqueTwoGroups(CRISPR_df, use_column, show_title = FALSE, only_chosen = only_chosen,
                                                         create_plot = FALSE
                                                         )
            Libraries_results <- BarPlot_UniqueLibraries(CRISPR_df, use_column, show_title = FALSE, only_chosen = only_chosen,
                                                         create_plot = FALSE
                                                         )
            y_upper_limit <- max(c(TwoGroups_results[["y_limits"]], Libraries_results[["y_limits"]]))
            combined_y_limits <- c(0, y_upper_limit)
          } else {
            combined_y_limits <- NULL
          }
          if (make_PNG) {
            folder_path <- file.path(output_plots_directory, "Bar charts", folder_name)
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
          BarPlot_UniqueTwoGroups(CRISPR_df, use_column, show_title = FALSE,
                                  only_chosen = only_chosen, use_y_limits = combined_y_limits
                                  )
          StackedBarPlot_results <- BarPlot_UniqueLibraries(CRISPR_df,
                                                            use_column,
                                                            show_title = FALSE,
                                                            only_chosen = only_chosen,
                                                            use_y_limits = combined_y_limits
                                                            )
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

  is_CRISPRko <- "Entrez_source_Brunello" %in% names(CRISPR_df)
  if (is_CRISPRko) {
    categorical_columns <- setdiff(categorical_columns, TSS_columns)
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
          pdf(file = file.path(output_plots_directory, "Bar charts", PDF_file_name),
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
            folder_path <- file.path(output_plots_directory, "Bar charts",
                                     paste0("Bar charts - ", file_name)
                                     )
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
          do.call(BarPlot_Sources, c(list(CRISPR_df   = CRISPR_df,
                                          use_column  = use_column
                                          ),
                                     args_list[[file_name]]
                                     )
                  )
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



PlotNumGenesInLibrary <- function() {
  use_width <- pdf_width * 0.85
  use_height <- pdf_height * 1.425
  for (make_PDF in c(FALSE, TRUE)) {
    if (make_PDF) {
      PDF_file_name <- "Library coverage.pdf"
      pdf(file = file.path(output_plots_directory, "Whole library", PDF_file_name),
          width = use_width, height = use_height
          )
    }
    NumGenesInLibrary()
    if (make_PDF) {
      dev.off()
    }
  }
}



NumGenesInLibrary <- function() {

  stopifnot("num_genes_in_library" %in% ls(envir = globalenv()))

  bars_mat <- t(num_genes_in_library[names(num_genes_in_library) != "GPP"])
  colnames(bars_mat)[[which(colnames(bars_mat) == "GPP_4_or_more")]] <- "GPP"

  colors_df <- SubgroupColorsDf()
  are_selected_colors <- !(colors_df[["Lightened"]]) &
                         (colors_df[["Color_name"]] %in% c("Blue", "Green", "Red", "Purple"))
  colors_df <- colors_df[are_selected_colors, ]

  spaces_vec <- rep(1.3, 4)

  old_mar <- par(mar = c(4, 4, 6, 3) + 0.1)

  bar_results <- PlotBars(bars_mat,
                          spaces_vec = spaces_vec,
                          x_limits = BarPlotXlimits(spaces_vec, side_space = 1.3),
                          colors_df = colors_df,
                          use_y_limits = c(0, 20000)
                          )
  title("Number of genes in the library", cex.main = 1, line = 3.625)

  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  text(x      = bar_results[["bar_positions"]],
       y      = par("usr")[[4]] + (plot_height * 0.06),
       labels = bars_mat[1, ],
       adj    = c(0.5, 1),
       cex    = 0.7,
       col    = "gray40",
       font   = 2,
       xpd    = NA
       )

  text(x      = bar_results[["bar_positions"]],
       y      = par("usr")[[3]] - (plot_height * 0.06),
       labels = AdjustTextVec(colnames(bars_mat)),
       adj    = c(0.5, 0.5),
       cex    = 0.8,
       col    = colors_df[["Dark"]],
       font   = 2,
       xpd    = NA
       )

  par(old_mar)

  return(invisible(bars_mat))
}








# Functions for generating histograms -------------------------------------

DrawHistogram <- function(overview_df, column_name) {
  are_4sg <- overview_df[["In_4sg_library"]] %in% "Yes"
  numeric_vec <- overview_df[[column_name]][are_4sg]
  hist_results <- hist(numeric_vec, breaks = 100, plot = FALSE)
  y_limits <- GetFlexibleAxisLimits(hist_results[["density"]], subtract_fraction = 0.01)
  hist(numeric_vec,
       las      = 1,
       mgp      = c(2.2, 0.45, 0),
       tcl      = -0.3,
       col      = brewer.pal(9, "Blues")[[8]],
       border   = NA,
       breaks   = 100,
       xlim     = c(-0.01, 1),
       ylim     = y_limits,
       xlab     = "",
       xaxs     = "i",
       yaxs     = "i",
       freq     = FALSE,
       main     = ""
       )
  box(bty = "l")
  title(paste0("Aggregate ", numeric_column_labels[[column_name]]),
        cex.main = 1, line = 2.2
        )
  return(invisible(NULL))
}





GetFlexibleAxisLimits <- function(numeric_vec, start_at_zero = TRUE, subtract_fraction = 0, add_fraction = 0) {
  use_range <- range(numeric_vec, na.rm = TRUE)
  if (start_at_zero && (use_range[[1]] >= 0)) {
    use_range[[1]] <- 0
  }
  axis_positions <- pretty(use_range)
  use_limits <- c(axis_positions[[1]], axis_positions[[length(axis_positions)]])
  use_span <- use_limits[[2]] - use_limits[[1]]
  use_limits[[1]] <- use_limits[[1]] - (subtract_fraction * use_span)
  use_limits[[2]] <- use_limits[[2]] + (add_fraction * use_span)
  return(use_limits)
}




DrawDeletionHistogram <- function(overview_df,
                                  column_name  = "Max_deletion_size",
                                  use_title    = "Size of the expected deletion (between the first and fourth cut sites)",
                                  fill_color   = brewer.pal(9, "Blues")[[8]],
                                  x_axis_label = "Number of base pairs (logarithmic scale)",
                                  y_axis_label = "Count",
                                  x_label_line = 2.2,
                                  y_label_line = 3,
                                  use_breaks   = 100,
                                  box_type     = "l",
                                  use_y_max    = NULL,
                                  box_color    = "black"
                                  ) {

  are_4sg <- overview_df[["In_4sg_library"]] %in% "Yes"
  numeric_vec <- overview_df[[column_name]][are_4sg]
  hist_results <- hist(log10(numeric_vec), breaks = use_breaks, plot = FALSE)
  x_limits <- GetFlexibleAxisLimits(hist_results[["breaks"]],
                                    start_at_zero = FALSE,
                                    add_fraction = 0.01
                                    )

  for_y_lim <- hist_results[["counts"]]
  if (!(is.null(use_y_max))) {
    for_y_lim <- c(use_y_max, for_y_lim)
  }
  y_limits <- GetFlexibleAxisLimits(for_y_lim, subtract_fraction = 0.02)

  hist(log10(numeric_vec),
       col      = fill_color,
       border   = NA,
       breaks   = use_breaks,
       xlim     = x_limits,
       ylim     = y_limits,
       xaxs     = "i",
       yaxs     = "i",
       freq     = TRUE,
       ann      = FALSE,
       axes     = FALSE,
       main     = "",
       mgp      = c(y_label_line, 1, 0),
       )

  axis(2,
       las = 1,
       mgp = c(3, 0.38, 0),
       tcl = -0.3,
       lwd = par("lwd")
       )

  axis_ticks <- axTicks(1)
  old_scipen <- options(scipen = -1)
  axis_labels <- as.character(10^axis_ticks)
  axis_labels <- sub("e+0", "0^", axis_labels, fixed = TRUE)
  options(old_scipen)

  axis(1,
       at     = axis_ticks,
       labels = parse(text = axis_labels),
       mgp    = c(3, 0.49, 0),
       tcl    = -0.3,
       lwd    = par("lwd")
       )

  old_scipen <- options(scipen = -1)
  as.character(10^axis_ticks)
  options(old_scipen)

  box(bty = box_type, col = box_color)

  title(use_title,
        cex.main = 1,
        line = 2.2
        )

  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  text(x      = par("usr")[[1]] - diff(grconvertX(c(0, y_label_line), from = "lines", to = "user")),
       y      = grconvertY(0.5, from = "npc", to = "user"),
       labels = VerticalAdjust(y_axis_label),
       srt    = 90,
       xpd    = NA,
       adj    = c(0.5, 0)
       )

  mtext(x_axis_label, line = x_label_line, side = 1, cex = par("cex"))
  return(invisible(hist_results))
}




GetLongestSubsequence <- function(CRISPR_df) {
  are_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list)
  if ("AltTSS_ID" %in% names(CRISPR_df)) {
    ID_column <- "AltTSS_ID"
  } else {
    ID_column <- "Combined_ID"
  }
  longest_subsequence_vec <- tapply(CRISPR_df[["sgRNA_sequence"]][are_4sg],
                                    CRISPR_df[[ID_column]][are_4sg],
                                    LongestSharedSubsequence
                                    )
  return(longest_subsequence_vec)
}




SharedSubsequencesBarplot <- function(CRISPR_df) {
  longest_subsequence_vec <- GetLongestSubsequence(CRISPR_df)
  longest_subsequence_table <- table(factor(longest_subsequence_vec, levels = 1:20))
  bar_positions <- barplot(longest_subsequence_table,
                           las       = 1,
                           mgp       = c(3, 0.45, 0),
                           tcl       = -0.3,
                           col       = brewer.pal(9, "Blues")[[8]],
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

  title("Length of the longest shared subsequence",
        cex.main = 1, line = 2.2
        )
  return(invisible(NULL))
}



Plot4sgData <- function(overview_df, CRISPR_df) {
  for (make_PDF in c(TRUE, FALSE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "Whole library",
                           "Histograms - 4sg combination.pdf"
                           ),
          width = pdf_width, height = pdf_height * 1.3
          )
    }
    par(oma = c(0, 1, 0.5, 1))
    if ("Max_deletion_size" %in% names(overview_df)) {
      DrawDeletionHistogram(overview_df)
    }
    SharedSubsequencesBarplot(CRISPR_df)
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
      pdf(file = file.path(output_plots_directory, "Whole library", "Venn diagrams.pdf"),
          width = pdf_width * 1.1, height = pdf_height * 1.05
          )
    }
    CRISPR_df <- FilterCRISPRDf(CRISPR_df)
    are_chosen <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list, show_messages = FALSE)
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




PlotVennDiagram <- function(sources_fac,
                            draw_plot       = FALSE,
                            quantities_cex  = 0.4,
                            quantities_font = 2,
                            use_padding     = 0.25,
                            use_colors      = c(brewer.pal(9, "Blues")[[3]], brewer.pal(9, "Greens")[[3]], brewer.pal(9, "Reds")[[3]])
                            ) {
  sources_table <- table(sources_fac)
  sources_names <- gsub(", ", "&", names(sources_table), fixed = TRUE)
  sources_table <- as.integer(sources_table)
  names(sources_table) <- sources_names
  eulerr_options("padding" = grid::unit(use_padding, "lines"))
  euler_plot <- plot(eulerr::euler(sources_table, shape = "circle"),
                     fills = use_colors,
                     edges = FALSE,
                     labels = list(cex = 1),
                     quantities = list(font = quantities_font, cex = quantities_cex)
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
                        embed_PNG                    = FALSE,
                        only_top4                    = FALSE
                        ) {

  CRISPR_df <- FilterCRISPRDf(CRISPR_df)

  if (only_top4) {
    CRISPR_df <- CRISPR_df[Are4sg(CRISPR_df, sublibraries_all_entrezs_list, show_messages = FALSE), ]
  }

  if (convert_CRISPOR_to_GuideScan) {
    CRISPR_df[["CRISPOR_CFD_specificity"]] <- ConvertCFDScores(CRISPR_df[["CRISPOR_CFD_specificity"]])
  }
  if (convert_GuideScan_to_CRISPOR) {
    CRISPR_df[["GuideScan_specificity"]] <- (100 / (100 + ((1 / CRISPR_df[["GuideScan_specificity"]]) - 1))) * 100
  }

  if (make_PNG || embed_PNG) {
    if (only_top4) {
      sub_folder <- "Scatter plots - B) only 4sg"
    } else {
      sub_folder <- "Scatter plots - A) all guides"
    }
  }

  if (make_PNG) {
    if (!("current_number" %in% ls(envir = globalenv()))) {
      current_number <- 1L
    }
    number_string <- formatC(current_number, width = 2, flag = "0")
    file_name <- paste0("Scatter plots - ",
                        number_string, ") ", gsub(":", " - ", show_title), ".png"
                        )
    file_path <- file.path(output_plots_directory, "Scatter plots",
                           sub_folder, file_name
                           )
    png(file = file_path,
        width = 5.75, height = 5.75, units = "in", res = 600
        )
    current_number <- current_number + 1L
    assign("current_number", current_number, envir = globalenv())
  }

  old_par <- par(mar = rep.int(5, 4))
  x_vec <- CRISPR_df[[x_column]]
  y_vec <- CRISPR_df[[y_column]]

  if (add_jitter) {
    set.seed(1)
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

  plot(1,
       las  = 1,
       mgp  = c(2.7, 0.5, 0),
       tcl  = -0.4,
       type = "n",
       xlab = numeric_column_labels[[x_column]],
       ylab = numeric_column_labels[[y_column]],
       xlim = x_axis_limits,
       ylim = y_axis_limits,
       xaxs = "i",
       yaxs = "i",
       bty  = "n"
       )

  alpha_hex <- substr(rgb(1, 1, 1, point_alpha), 8, 9)

  if (embed_PNG) {
    current_device <- dev.cur()
    temp_file_path <- file.path(output_plots_directory, "Scatter plots",
                                sub_folder, "temp.png"
                                )
    png(file = temp_file_path,
        width = 4.75, height = 4.75, units = "in", res = 600
        )
    par(mar = rep(0, 4))
    plot(x_vec,
         y_vec,
         col  = paste0(brewer.pal(9, "Blues")[[7]], alpha_hex),
         pch  = 16,
         cex  = point_cex,
         axes = FALSE,
         ann  = FALSE,
         xlim = x_axis_limits,
         ylim = y_axis_limits,
         xaxs = "i",
         yaxs = "i",
         bty  = "n"
         )

  }

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

  if (embed_PNG) {
    dev.off()
    raster_array <- readPNG(temp_file_path)
    file.remove(temp_file_path)
    dev.set(current_device)
    rasterImage(raster_array,
                xleft = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                )
  } else {
    points(x_vec,
           y_vec,
           col = paste0(brewer.pal(9, "Blues")[[7]], alpha_hex),
           pch = 16,
           cex = point_cex
           )
  }

  box()
  if (!(is.null(show_title))) {
    title(show_title, cex.main = par("cex") * 0.9)
  }
  par(old_par)

  if (make_PNG) {
    dev.off()
  }
  return(invisible(NULL))
}





MakeScatterPlots <- function(CRISPR_df, embed_raster_within_PDFs = TRUE) {

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
          pdf(file = file.path(output_plots_directory, "Scatter plots",
                               paste0("Scatter plots", append_to_filename, ".pdf")
                               ),
              width = plot_dimensions, height = plot_dimensions
              )
          embed_PNG <- embed_raster_within_PDFs
        } else {
          embed_PNG <- FALSE
        }
        ScatterPlot(CRISPR_df, "GuideScan_specificity", "GuideScan_efficiency",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "Efficacy vs. specificity (GuideScan)",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )
        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_3MM_specificity", "CRISPOR_Doench_efficacy",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title = "Efficacy vs. specificity (CRISPOR)",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_Doench_efficacy",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title = "Efficacy vs. specificity (CRISPOR, up to 4MM)",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_Doench_efficacy", "GuideScan_efficiency",
                      point_alpha = 0.2, point_cex = 0.4, add_jitter = TRUE,
                      show_title = "Original vs. updated Doench efficacy scores",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
        }

        ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                    show_title = "Specificity score vs. number of off-targets (GuideScan)",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )

        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_Num_2or3MM", "CRISPOR_3MM_specificity",
                      point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                      show_title = "Specificity score vs. number of off-targets (CRISPOR)",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
          if (FALSE) {
            ScatterPlot(CRISPR_df, "GuideScan_specificity", "CRISPOR_CFD_specificity",
                        point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                        make_PNG = make_PNG, only_top4 = only_top4,
                        embed_PNG = embed_PNG
                        )
            ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_CFD_specificity", # This is just for confirmation
                        convert_GuideScan_to_CRISPOR = TRUE, point_alpha = 0.2, point_cex = 0.4,
                        make_PNG = make_PNG, only_top4 = only_top4,
                        embed_PNG = embed_PNG
                        )
          }
        }

        ScatterPlot(CRISPR_df, "CRISPOR_Num_2or3MM", "GuideScan_Num_2or3MM",
                    point_alpha = 0.5, point_cex = 0.3, identical_axes = TRUE,
                    show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )
        ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "CRISPOR_Num_2or3MM",
                    point_alpha = 0.5, point_cex = 0.3,
                    identical_axes = TRUE, custom_axis_limits = c(0, 1000),
                    show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in)",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )
        ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "CRISPOR_Num_2or3MM",
                    point_alpha = 0.5, point_cex = 0.3,
                    identical_axes = TRUE, custom_axis_limits = c(0, 200),
                    show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in more)",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )

        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_Num_2MM", "GuideScan_Num_2MM",
                      point_alpha = 0.5, point_cex = 0.3, identical_axes = TRUE,
                      show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
          ScatterPlot(CRISPR_df, "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
                      point_alpha = 0.5, point_cex = 0.3,
                      identical_axes = TRUE,
                      custom_axis_limits = c(0, 300),
                      show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in)",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
          ScatterPlot(CRISPR_df, "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
                      point_alpha = 0.5, point_cex = 0.3,
                      identical_axes = TRUE, custom_axis_limits = c(0, 50),
                      show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in more)",
                      add_jitter = TRUE,
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_CFD_specificity", "GuideScan_specificity",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title  = "Specificity \u2013 GuideScan score vs. original CRISPOR CFD score",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
        }

        ScatterPlot(CRISPR_df, "CRISPOR_3MM_specificity", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4, identical_axes = TRUE,
                    show_title = "Specificity score \u2013 GuideScan vs. CRISPOR",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )

        ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "Specificity score \u2013 GuideScan vs. CRISPOR (4MM)",
                    make_PNG = make_PNG, only_top4 = only_top4,
                    embed_PNG = embed_PNG
                    )

        if (!(only_selected)) {
          ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_3MM_specificity",
                      point_alpha = 0.2, point_cex = 0.4,
                      show_title = "CRISPOR specificity scores: 3MM vs. 4MM",
                      make_PNG = make_PNG, only_top4 = only_top4,
                      embed_PNG = embed_PNG
                      )
        }
        if (make_PDF) {
          dev.off()
        }
      }
    }
  }
}





# Functions for producing hybrid bar plots/doughnut plots -----------------

doughnut <- function(x,
                     labels = names(x),
                     edges = 200,
                     outer.radius = 0.8,
                     inner.radius = 0.6,
                     clockwise = TRUE,
                     init.angle = if (clockwise) 90 else 0,
                     density = NULL,
                     angle = 45,
                     col = NULL,
                     border = FALSE,
                     lty = NULL,
                     main = NULL,
                     radius_factor = 1,
                     x_origin = 0,
                     y_origin = 0,
                     add = TRUE,
                     ...
                     ) {
  # Adapted from https://www.r-graph-gallery.com/130-ring-or-donut-chart.html
    if (!is.numeric(x) || any(is.na(x) | x < 0))
        stop("'x' values must be positive.")
    if (is.null(labels))
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    if (!(add)) {
      plot.new()
      pin <- par("pin")
      xlim <- ylim <- c(-1, 1)
      if (pin[1L] > pin[2L])
        xlim <- (pin[1L]/pin[2L]) * xlim
      else ylim <- (pin[2L]/pin[1L]) * ylim
      plot.window(xlim, ylim, "", asp = 1)
    }

    if (is.null(col))
        col <- if (is.null(density))
          palette()
        else par("fg")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    if (!(is.null(lty))) {
      lty <- rep(lty, length.out = nx)
    }
    angle <- rep(angle, length.out = nx)
    if (!(is.null(density))) {
      density <- rep(density, length.out = nx)
    }
    twopi <- if (clockwise)
        -2 * pi
    else 2 * pi
    t2xy <- function(t, radius) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p),
             y = radius * sin(t2p))
    }
    for (i in seq_len(nx)) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
                  outer.radius * radius_factor)
        polygon(c(P$x, 0) + x_origin, c(P$y, 0) + y_origin, density = density[i],
                angle = angle[i], border = border[i],
                col = col[i], lty = lty[i], xpd = NA
                )
        Pout <- t2xy(mean(x[i + 0:1]), outer.radius * radius_factor)
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            lines(c(1, 1.05) * Pout$x + x_origin, c(1, 1.05) * Pout$y + y_origin)
            text(1.1 * Pout$x + x_origin, 1.1 * Pout$y + y_origin, labels[i],
                 xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
                 ...
                 )
        }
        ## Add white disc
        Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                  inner.radius * radius_factor)
        polygon(Pin$x + x_origin, Pin$y + y_origin, density = density[i],
                angle = angle[i], border = border[i],
                col = "white", lty = lty[i], xpd = NA
                )
    }

    title(main = main, ...)
    invisible(NULL)
}





DonutBars <- function(use_factor         = NULL,
                      use_colors,
                      space              = 0.5,
                      use_title          = NULL,
                      title_font         = 2,
                      bar_labels         = NULL,
                      donut_radius       = 0.15,
                      donut_label        = "",
                      donut_text_size    = 0.5,
                      donut_x_mid        = 0.8,
                      donut_y_mid        = 0.2,
                      counts_vec         = NULL,
                      use_labels         = NULL,
                      show_percentages   = TRUE,
                      show_axis          = FALSE,
                      side_text_size     = 0.9,
                      bar_text_size      = 0.7,
                      bar_text_font      = 2,
                      text_dark_color    = NULL,
                      use_mai            = c(0.2, 1.9, 0.1, 0.3),
                      use_omi            = c(0, 0, 0.45, 0),
                      x_axis_label       = NULL,
                      leave_spaces_for   = NULL,
                      donut_inner_radius = 0.4,
                      use_line_height    = 1.15,
                      bar_label_line     = 0.6,
                      draw_box           = FALSE,
                      y_axis_label       = NULL,
                      percent_max        = NULL,
                      bottom_axis        = FALSE
                      ) {

  if (is.null(counts_vec)) {
    stopifnot(length(use_colors) == nlevels(use_factor))
    counts_vec <- tabulate(use_factor)
    if (is.null(use_labels)) {
      use_labels <- levels(use_factor)
    }
  } else {
    if (is.null(use_labels)) {
      use_labels <- names(counts_vec)
    }
  }

  if (show_percentages) {
    percentages_vec <- counts_vec / sum(counts_vec) * 100
    percentages_strings <- rep(NA, length(percentages_vec))
    percentages_strings[percentages_vec < 1] <- signif(percentages_vec[percentages_vec < 1], digits = 1)
    are_single_digit <- (percentages_vec < 10) & (percentages_vec >= 1)
    percentages_strings[are_single_digit] <- format(percentages_vec[are_single_digit], digits = 2)
    percentages_strings[percentages_vec > 10] <- round(percentages_vec[percentages_vec > 10])
    count_labels <- paste0(counts_vec, " (", percentages_strings, "%)")
  } else {
    count_labels <- counts_vec
  }

  num_bars <- length(counts_vec)

  if (!(is.null(leave_spaces_for))) {
    num_positions <- leave_spaces_for
  } else {
    num_positions <- num_bars
  }

  side_space <- space * (2 / 3)
  num_units <- num_positions + ((num_positions - 1) * space) + (2 * side_space)
  unit_width <- 1 / num_units

  bar_mid_positions <- cumsum(rep(unit_width, num_positions)) +
                       cumsum(rep(unit_width * space, num_positions)) -
                       (unit_width / 2) - ((side_space * unit_width) / 2)
  bar_mid_positions <- rev(bar_mid_positions)[seq_len(num_bars)]

  x_max <- max(counts_vec)
  if (show_axis) {
    total_num <- sum(counts_vec)
    if (is.null(percent_max)) {
      percent_max <- GetFlexibleAxisLimits(c(0, x_max / total_num * 100))[[2]]
    }
    x_max <- (percent_max / 100) * total_num
    pretty_pos <- pretty(c(0, percent_max))
  }
  bar_lengths <- counts_vec / x_max


  plot.new()
  old_par <- par(mai     = use_mai,
                 omi     = use_omi,
                 lheight = use_line_height
                 )
  plot.window(c(0, 1), c(0, 1), "", asp = 1)


  if (!(is.null(use_title))) {
    title(use_title, outer = TRUE, cex.main = par("cex") * 0.8, line = 0,
          font.main = title_font
          )
  }

  x_range <- par("usr")[[2]] - par("usr")[[1]]
  if (draw_box) {
    x_space_fraction <- 0.015
  } else {
    x_space_fraction <- 0
  }

  x_space <- x_space_fraction * x_range

  bar_y_mids <- grconvertY(bar_mid_positions, from = "npc", to = "user")
  y_half_bar <- (par("usr")[[4]] - par("usr")[[3]]) * unit_width * 0.5

  corr_bar_lengths <- ((1 - x_space_fraction) * bar_lengths)

  for (i in seq_len(num_bars)) {
    rect(xleft   = grconvertX(x_space, from = "npc", to = "user"),
         xright  = grconvertX(x_space + corr_bar_lengths[[i]], from = "npc", to = "user"),
         ybottom = bar_y_mids[[i]] - y_half_bar,
         ytop    = bar_y_mids[[i]] + y_half_bar,
         border  = use_colors[[i]],
         col     = use_colors[[i]],
         lwd     = 0.5,
         xpd     = NA
         )
  }

  are_too_small <- bar_lengths < (par("usr")[[2]] * 0.2)

  text(x      = grconvertX(x_space + corr_bar_lengths[are_too_small] + (x_range * 0.02), from = "npc", to = "user"),
       y      = bar_y_mids[are_too_small],
       xpd    = NA,
       labels = count_labels[are_too_small],
       adj    = c(0, 0.5),
       cex    = bar_text_size,
       font   = bar_text_font,
       col    = if (is.null(text_dark_color)) "gray52" else text_dark_color
       )

  text(x      = grconvertX(x_space + corr_bar_lengths[!(are_too_small)] / 2, from = "npc", to = "user"),
       y      = bar_y_mids[!(are_too_small)],
       xpd    = NA,
       labels = count_labels[!(are_too_small)],
       adj    = 0.5,
       cex    = bar_text_size,
       font   = bar_text_font,
       col    = ifelse((colMeans(col2rgb(use_colors[!(are_too_small)])) / 255) < 0.5,
                       "gray92",
                       if (is.null(text_dark_color)) "gray48" else text_dark_color
                       )
       )

  text(x      = par("usr")[[1]] - (diff(grconvertX(c(0, bar_label_line), from = "lines", to = "user"))),
       y      = bar_y_mids,
       xpd    = NA,
       labels = use_labels,
       adj    = c(1, 0.5),
       cex    = side_text_size
       )


  donut_x_mid <- grconvertX(donut_x_mid, from = "npc", to = "user")
  donut_y_mid <- grconvertY(donut_y_mid, from = "npc", to = "user")

  doughnut(rev(counts_vec), labels = NA, col = rev(use_colors), add = TRUE,
           x_origin = donut_x_mid, y_origin = donut_y_mid,
           radius_factor  = donut_radius, inner.radius = donut_inner_radius
           )

  text(x      = donut_x_mid,
       y      = donut_y_mid,
       labels = donut_label,
       font   = bar_text_font,
       cex    = donut_text_size
       )

  if (show_axis) {
    actual_pos <- seq(from = par("usr")[[1]] + x_space,
                      to   = par("usr")[[2]],
                      by   = (par("usr")[[2]] - (par("usr")[[1]] + x_space)) / (length(pretty_pos) - 1)
                      )

    axis(if (bottom_axis) 1 else 3,
         at     = actual_pos,
         labels = paste0(pretty_pos, "%"),
         mgp    = c(3, 0.38, 0),
         tcl    = -0.3,
         lwd    = par("lwd")
         )

    if (!(is.null(x_axis_label))) {
      mtext(x_axis_label, side = if (bottom_axis) 1 else 3,
            line = 1.7, cex = par("cex")
            )
    }
  }

  if (!(is.null(y_axis_label))) {
    text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 2.53), from = "lines", to = "user")),
         y      = grconvertY(0.5, from = "npc", to = "user"),
         labels = VerticalAdjust(y_axis_label),
         srt    = 90,
         xpd    = NA,
         adj    = c(0.5, 0)
         )
  }

  if (draw_box) {
    box_gray <- "gray60"
    box(col = box_gray)
  }

  par(old_par)
  return(invisible(NULL))
}




ReverseList <- function(my_list) {
  my_entries <- unlist(my_list)
  my_names <- rep(names(my_list), times = lengths(my_list, use.names = FALSE))
  names(my_names) <- my_entries
  return(my_names)
}


DeletionsDonutBar <- function(deletions_summary_df, ...) {
  CategoriesDonutBar(as.character(deletions_summary_df[, "Gene_targets_summary"]), ...)
}

SummaryDonutBar <- function(CRISPR_df, targets_df, ...) {
  categories_vec <- Summarize4sgTargets(CRISPR_df, targets_df)
  CategoriesDonutBar(categories_vec, ...)
}


TSSDonutBar <- function(CRISPR_df,
                        use_title       = "Number of TSSs targeted for each gene in the 4sg library",
                        donut_label     = "4sg",
                        donut_text_size = 0.8,
                        donut_radius    = 0.26,
                        donut_x_mid     = 0.72,
                        donut_y_mid     = 0.28,
                        use_colors      = carto_pal(6, "Emrld"),
                        ...
                        ) {

  num_TSS_labels <- c("Only 1 TSS targeted",
                    paste0(seq_len(20 - 1) + 1,
                           " TSSs targeted"
                           )
                    )

  are_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list)
  num_TSS_vec <- tapply(CRISPR_df[["AltTSS_ID"]][are_4sg],
                        CRISPR_df[["Combined_ID"]][are_4sg],
                        function(x) length(unique(x))
                        )
  num_TSS_fac <- factor(num_TSS_vec)
  levels(num_TSS_fac) <- num_TSS_labels[seq_len(nlevels(num_TSS_fac))]

  stopifnot(nlevels(num_TSS_fac) <= 6)

  DonutBars(num_TSS_fac,
            use_colors      = use_colors,
            donut_label     = donut_label,
            donut_text_size = donut_text_size,
            donut_radius    = donut_radius,
            donut_x_mid     = donut_x_mid,
            donut_y_mid     = donut_y_mid,
            use_title       = use_title,
            ...
            )
}



manuscript_map_list <- list(
  "Affect only the\nintended gene"         = "Only intended",
  "Also affect on-site\nunintended genes"  = c("Intended and unintended (in the same locus)",
                                             "Intended and unintended (in the same loci)"
                                             ),
  "Also affect off-site\nunintended genes" = "Intended and unintended (in other loci)",
  "Unknown (no\ngene location)"            = "No location data for the intended gene",
  "Do not affect the\nintended gene"       = c("Only unintended (in the same locus)",
                                               "Only unintended (in the same loci)",
                                               "Only unintended (multiple loci)",
                                               "No targets (single locus)",
                                               "No targets (multiple loci)",
                                               "No hits in the reference genome"
                                               )
)



CategoriesDonutBar <- function(character_vec,
                               donut_label     = "4sg",
                               donut_radius    = 0.26,
                               donut_x_mid     = 0.72,
                               donut_y_mid     = 0.28,
                               donut_text_size = 0.8,
                               use_map_list    = NULL,
                               ...
                               ) {

  categories_map_list <- list(
    "Affect only the\nintended gene"             = "Only intended",
    "Affect unintended genes\nin the same locus" = c("Intended and unintended (in the same locus)",
                                                     "Intended and unintended (in the same loci)"
                                                     ),
    "Affect unintended genes\nat other loci"     = "Intended and unintended (in other loci)",
    "Unknown (no location\ndata for this gene)"  = "No location data for the intended gene",
    "Do not seem to affect\nthe intended gene"   = c("Only unintended (in the same locus)",
                                                     "Only unintended (in the same loci)",
                                                     "Only unintended (multiple loci)",
                                                     "No targets (single locus)",
                                                     "No targets (multiple loci)",
                                                     "No hits in the reference genome"
                                                     )
  )
  if (!(is.null(use_map_list))) {
    categories_map_list <- use_map_list
  }

  new_categories_map <- ReverseList(categories_map_list)

  stopifnot(all(character_vec %in% names(new_categories_map)))
  category_fac <- factor(new_categories_map[character_vec],
                         levels = names(categories_map_list)
                         )
  category_colors <- carto_pal(12, "Safe")[c(4, 5, 2, 11, 10)]

  DonutBars(category_fac,
            category_colors,
            donut_label     = donut_label,
            donut_radius    = donut_radius,
            donut_x_mid     = donut_x_mid,
            donut_y_mid     = donut_y_mid,
            donut_text_size = donut_text_size,
            ...
            )
}




Summarize4sgTargets <- function(CRISPR_df, targets_df) {

  stopifnot(nrow(CRISPR_df) == nrow(targets_df))

  are_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list)

  if ("AltTSS_ID" %in% names(CRISPR_df)) {
    ID_column <- "AltTSS_ID"
  } else {
    ID_column <- "Entrez_ID"
  }

  summary_splits <- split(as.character(targets_df[["Gene_targets_summary"]][are_4sg]),
                          CRISPR_df[[ID_column]][are_4sg]
                          )

  only_unintended_categories <- c("Only unintended (in the same locus)",
                                  "Only unintended (in the same loci)",
                                  "Only unintended (multiple loci)"
                                  )
  intended_categories <- c("Only intended",
                           "Intended and unintended (in other loci)",
                           "Intended and unintended (in the same loci)",
                           "Intended and unintended (in the same locus)"
                           )
  are_divergent <- vapply(summary_splits,
                          function(x) {
                            (any(only_unintended_categories %in% x) &&
                               any(intended_categories %in% x)
                            )
                          },
                          logical(1)
                          )

  target_other_loci <- vapply(summary_splits,
                              function(x) {
                                any(c("Only unintended (multiple loci)", "Intended and unintended (in other loci)") %in% x)
                              },
                              logical(1)
                              )

  summary_splits[are_divergent & !(target_other_loci)] <- list("Intended and unintended (in the same loci)")
  summary_splits[are_divergent & target_other_loci] <- list("Intended and unintended (in other loci)")

  category_order <- c(
    "Intended and unintended (in other loci)",
    "Intended and unintended (in the same locus)",
    "Intended and unintended (in the same loci)",
    "Only unintended (in the same locus)",
    "Only unintended (in the same loci)",
    "Only unintended (multiple loci)",
    "Only intended",
    "No targets (multiple loci)",
    "No targets (single locus)",
    "No location data for the intended gene",
    "No hits in the reference genome"
  )
  stopifnot(identical(sort(category_order), sort(levels(targets_df[["Gene_targets_summary"]]))))

  summary_vec <- vapply(summary_splits,
                        function(x) {
                          x[[which.min(match(x, category_order))]]
                        },
                        ""
                        )
  return(summary_vec)
}




CompareTSSDonutBars <- function(CRISPR_df) {

  stopifnot("num_TSSs_hCRISPR" %in% ls(envir = globalenv()))

  TSSDonutBar(CRISPR_df)

  num_TSS_labels <- c("Only 1 TSS targeted",
                    paste0(seq_len(20 - 1) + 1,
                           " TSSs targeted"
                           )
                    )

  library_name <- names(dimnames(num_TSSs_hCRISPR))

  DonutBars(counts_vec      = num_TSSs_hCRISPR,
            use_labels      = num_TSS_labels[seq_along(num_TSSs_hCRISPR)],
            use_colors      = carto_pal(6, "Emrld")[seq_along(num_TSSs_hCRISPR)],
            donut_label     = library_name,
            donut_text_size = 0.5,
            donut_radius    = 0.26,
            donut_x_mid     = 0.72,
            donut_y_mid     = 0.28,
            use_title       = paste0("Number of TSSs targeted for each gene in the ",
                                     library_name, " library"
                                     )
            )
}



DrawAllDonutBars <- function(CRISPR_df = NULL) {

  for (use_PDF in c(FALSE, TRUE)) {

    if (use_PDF) {
      PDF_file_name <- "Doughnut charts.pdf"
      pdf(file = file.path(output_plots_directory, "Whole library", PDF_file_name),
          width = 6, height = 4.45
          )
    }

    if ("num_TSSs_hCRISPR" %in% ls(envir = globalenv())) {

      SummaryDonutBar(CRISPR_df,
                      TSS_targets_df,
                      use_title = "Genes (including non-coding RNAs) affected by each 4-guide combination"
                      )
      SummaryDonutBar(CRISPR_df,
                      TSS_protein_targets_df,
                      use_title = "Protein-coding genes affected by each 4-guide combination"
                      )
      SummaryDonutBar(CRISPR_df,
                      main_TSS_targets_df,
                      use_title = "Genes (main TSS only) affected by each 4-guide combination"
                      )
      SummaryDonutBar(CRISPR_df,
                      main_TSS_protein_targets_df,
                      use_title = "Protein-coding genes (main TSS only) affected by each 4-guide combination"
                      )

      CompareTSSDonutBars(CRISPR_df)

    } else {

      DeletionsDonutBar(deletions_CDS_df,
                        use_title = "Genes (including non-coding RNAs) affected by each 4-guide combination"
                        )
      DeletionsDonutBar(deletions_CDS_protein_df,
                        use_title = "Protein-coding genes affected by each 4-guide combination"
                        )

    }

    if (use_PDF) {
      dev.off()
    }
  }
}


manuscript_donut_args <- list(
  use_title          = NULL,
  show_axis          = TRUE,
  bar_text_size      = 1,
  side_text_size     = 1,
  show_percentages   = FALSE,
  bar_text_font      = 1,
  donut_radius       = 0.28,
  donut_text_size    = 0.9,
  donut_label        = "4sg",
  use_mai            = c(0.4, 1, 0.4, 0.15),
  use_omi            = rep(0, 4),
  donut_inner_radius = 0.35,
  use_line_height    = 0.95,
  donut_x_mid        = 0.80,
  donut_y_mid        = 0.26,
  bar_label_line     = 0.5,
  draw_box           = FALSE,
  bottom_axis        = FALSE
)








# Functions for generating multi-plot layouts -----------------------------

MakeEmptyPlot <- function() {
  plot(1, xlim = c(0, 1), ylim = c(0, 1), type = "n",
       xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE
       )
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







# Functions for formatting figures for the manuscript ---------------------


## Define plot dimensions for the manuscript
horizontal_width <- 2.4
horizontal_height <- 1.7

horizontal_mai <- c(0.43, 0.61, 0.16, 0.38)
vertical_mai   <- c(0.30, 0.47, 0.30, 0.32)

manuscript_wide_width <- 3.57
manuscript_wide_height <- 1.85
manuscript_wide_mai <- c(0.30, 1, 0.30, 0.32)

manuscript_cex <- 0.7
manuscript_lwd <- 0.8


## Define plot color schemes for the manuscript
manuscript_CRISPRo_colors <- brewer.pal(11, "RdBu")[c(7, 10)]
dark_blue <- Palify("#1A59B4", fraction_pale = 0.05)
manuscript_CRISPRo_colors[[2]] <- colorRampPalette(c(manuscript_CRISPRo_colors[[2]], dark_blue))(6)[[4]]

manuscript_CRISPRa_colors <- brewer.pal(11, "PiYG")[c(5, 2)]
# manuscript_CRISPRa_colors[[2]] <- colorRampPalette(c("#FF007F", manuscript_CRISPRa_colors[[2]]))(3)[[2]]
# manuscript_CRISPRa_colors[[2]] <- Palify(manuscript_CRISPRa_colors[[2]], fraction_pale = 0.15)
manuscript_CRISPRa_colors[[2]] <- "#E44499"




## Define dimensions and color scheme to PNGs

PNG_increment <- 0.05
PNG_multiplier <- 2
# PNG_height <- (horizontal_height + PNG_increment) * PNG_multiplier
PNG_height <- 3.75
PNG_width <- horizontal_width * PNG_multiplier
PNG_wide_width <- manuscript_wide_width * PNG_multiplier
PNG_wide_height <- manuscript_wide_height * PNG_multiplier
PNG_cex <- manuscript_cex * PNG_multiplier * 0.9
PNG_lwd <- manuscript_lwd * PNG_multiplier
PNG_vertical_mai <- vertical_mai
PNG_vertical_mai[[1]] <- PNG_vertical_mai[[1]] + PNG_increment
PNG_vertical_mai <- PNG_vertical_mai * PNG_multiplier
PNG_wide_mai <- manuscript_wide_mai * PNG_multiplier
PNG_CRISPRa_colors <- manuscript_CRISPRa_colors
PNG_CRISPRa_colors[[2]] <- "#D45498"






GetBarPlotStart <- function(num_bars, space, use_factor = 0.04) {
  space - (((space * (num_bars - 1)) + num_bars) * use_factor)
}
GetBarplotEnd <- function(num_bars, space, use_factor = 0.04) {
  space + (((space * (num_bars - 1)) + num_bars) * (1 + use_factor))
}



PlotBarplotMat <- function(barplot_mat,
                           colors_vec,
                           positions_vec = seq_len(ncol(barplot_mat)),
                           width = 2/3,
                           horizontal = FALSE
                           ) {

  num_categories <- nrow(barplot_mat)
  stopifnot(length(colors_vec) == num_categories)

  lower_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    if (x == 1) {
      rep(0, ncol(barplot_mat))
    } else {
      colSums(barplot_mat[seq_len(x - 1), , drop = FALSE])
    }
  })
  upper_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    colSums(barplot_mat[seq_len(x), , drop = FALSE])
  })

  final_width <- width * ((max(positions_vec) - min(positions_vec)) / (ncol(barplot_mat) - 1))

  for (i in seq_len(ncol(barplot_mat))) {
    for (j in seq_len(num_categories)) {
      if (lower_borders_vec_list[[j]][[i]] != upper_borders_vec_list[[j]][[i]]) { # Prevent zero-width rectangles from showing up as tiny lines on the resulting PDF
        if (horizontal) {
          rect(xleft   = lower_borders_vec_list[[j]][[i]],
               xright  = upper_borders_vec_list[[j]][[i]],
               ybottom = positions_vec[[i]] - (final_width / 2),
               ytop    = positions_vec[[i]] + (final_width / 2),
               col     = colors_vec[[j]],
               border  = NA,
               xpd     = NA
               )
        } else {
          rect(ybottom = lower_borders_vec_list[[j]][[i]],
               ytop    = upper_borders_vec_list[[j]][[i]],
               xleft   = positions_vec[[i]] - (final_width / 2),
               xright  = positions_vec[[i]] + (final_width / 2),
               col     = colors_vec[[j]],
               border  = NA,
               xpd     = NA
               )
        }
      }
    }
  }
  return(invisible(NULL))
}




ManuscriptAnnotate <- function(group_positions,
                               group_names,
                               modality_label,
                               axis_label,
                               horizontal         = TRUE,
                               use_colors         = "#000000",
                               colored_modality   = FALSE,
                               modality_on_top    = FALSE,
                               modality_on_side   = FALSE,
                               modality_on_bottom = FALSE
                               ) {

  if (horizontal) {
    mtext(group_names, line = 0.35, at = group_positions, side = 2, las = 2,
          cex = par("cex")
          )
  } else {
    text(x      = group_positions,
         y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.44), from = "lines", to = "user")),
         labels = group_names,
         adj    = c(0.5, 1),
         xpd    = NA
         )
  }

  if (colored_modality) {
    modality_color <- colorRampPalette(use_colors)(3)[[2]]
  } else {
    modality_color <- NULL
  }

  if (horizontal) {
    mtext(modality_label, line = 0.05, col = modality_color, cex = par("cex"))
  } else {
    if (modality_on_top) {
      text(x      = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.5),
           y      = par("usr")[[4]] + diff(grconvertY(c(0, 0.78), from = "lines", to = "user")),
           labels = modality_label,
           col    = modality_color,
           xpd    = NA
           )
    } else if (modality_on_side) {
      text(x      = par("usr")[[2]] + diff(grconvertX(c(0, 0.75), from = "lines", to = "user")), # Formerly 0.32 lines
           y      = grconvertY(0.5, from = "npc", to = "user"),
           labels = modality_label,
           srt    = 270,
           adj    = c(0.5, 0),
           col    = modality_color,
           xpd    = NA
           )
    } else if (modality_on_bottom) {
      text(x      = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.5),
           y      = par("usr")[[3]] - diff(grconvertY(c(0, 2), from = "lines", to = "user")),
           labels = modality_label,
           col    = modality_color,
           xpd    = NA
           )
    }
  }

  if (horizontal) {
    mtext(axis_label, line = 1.7, side = 1, cex = par("cex"))
  } else {
    text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 2.53), from = "lines", to = "user")),
         y      = grconvertY(0.5, from = "npc", to = "user"),
         labels = VerticalAdjust(axis_label),
         srt    = 90,
         xpd    = NA,
         adj    = c(0.5, 0)
         )
  }
}




ManuscriptGrid <- function(horizontal = TRUE) {
  axis_ticks <- axTicks(if (horizontal) 1 else 2)
  grid_positions <- seq(axis_ticks[[1]],
                        axis_ticks[[length(axis_ticks)]],
                        (axis_ticks[[length(axis_ticks)]] - axis_ticks[[1]]) / (length(axis_ticks) - 1) / 2
                        )
  grid_colors <- ifelse(grid_positions %in% axis_ticks, "gray85", "gray95")
  if (horizontal) {
    segments(y0   = par("usr")[[3]],
             y1   = par("usr")[[4]],
             x0   = grid_positions,
             col  = grid_colors,
             lend = "butt",
             xpd  = NA
             )
  } else {
    segments(x0   = par("usr")[[1]],
             x1   = par("usr")[[2]],
             y0   = grid_positions,
             col  = grid_colors,
             lend = "butt",
             xpd  = NA
             )
  }
  return(axis_ticks)
}




ManuscriptBars <- function(counts_mat,
                           axis_label           = "",
                           same_space_on_edge   = FALSE,
                           space_like_boxplot   = TRUE,
                           numeric_limits       = NULL,
                           lollipop             = FALSE,
                           horizontal           = FALSE,
                           use_space            = if (horizontal) 0.6 else 0.8,
                           abbreviate_thousands = TRUE,
                           modality_on_side     = FALSE,
                           modality_on_bottom   = FALSE,
                           abbreviate_libraries = TRUE,
                           use_cex              = manuscript_cex,
                           use_lwd              = manuscript_lwd,
                           use_mai              = NULL,
                           expected_SNP_percent = TRUE,
                           CRISPRa_colors       = manuscript_CRISPRa_colors,
                           CRISPRo_colors       = manuscript_CRISPRo_colors
                           ) {

  ## Prepare colors, percentages and labels
  group_names <- colnames(counts_mat)
  if ("Brunello" %in% colnames(counts_mat)) {
    modality_label <- "CRISPRo"
    use_colors <- manuscript_CRISPRo_colors
  } else if ("Calabrese" %in% colnames(counts_mat)) {
    modality_label <- "CRISPRa"
    use_colors <- CRISPRa_colors
    if (horizontal) {
      group_names[group_names == "hCRISPRa-v2"] <- "hCRISPRa\n-v2"
    } else {
      group_names[group_names == "hCRISPRa-v2"] <- "hCa-v2"
      if (abbreviate_libraries) {
        group_names[group_names == "Calabrese"] <- "Calab"
      }
    }
  } else {
    modality_label <- ""
    use_colors <- NULL
  }
  one_row <- nrow(counts_mat) == 1
  is_proportion <- !(one_row) && is.integer(counts_mat)
  if (is_proportion) {
    bar_colors <- use_colors
  } else {
    bar_colors <- use_colors[[2]]
  }

  if (one_row) {
    bars_mat <- counts_mat
  } else if (is_proportion) {
    bars_mat <- prop.table(counts_mat, margin = 2) * 100
  } else {
    if (expected_SNP_percent) {
      bars_mat <- t(counts_mat[1, ] / counts_mat[2, ]) * 100
    } else {
      bars_mat <- counts_mat[1, , drop = FALSE]
    }
  }

  ## Prepare the data axis
  if (!(is.null(numeric_limits))) {
    numeric_limits <- numeric_limits
    numeric_axis_pos <- pretty(numeric_limits)
  } else if (is_proportion) {
    numeric_axis_pos <- seq(0, 100, by = 20)
    numeric_limits <- c(0, 100)
  } else {
    use_numeric_limits <- c(0, max(bars_mat[1, ] * 1.00))
    numeric_axis_pos <- pretty(use_numeric_limits)
    numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  }
  if (horizontal) {
    numeric_axis_labels <- numeric_axis_pos
  } else {
    numeric_axis_labels <- format(numeric_axis_pos)
  }
  if (!(one_row) && !(!is_proportion && !expected_SNP_percent)) {
    numeric_axis_labels <- paste0(numeric_axis_labels, "%")
  }

  ## Prepare the groups axis
  num_groups <- ncol(counts_mat)
  if (space_like_boxplot) {
    ## Prepare the groups axis
    spaces_vec <- rep(1.15, num_groups - 1)
    group_positions <- cumsum(c(1, spaces_vec))
    if (horizontal) {
      group_positions <- rev(group_positions)
    }
    group_limits <- c((min(group_positions) - 0.5) - (num_groups * 0.04),
                      max(group_positions) + 0.5 + (num_groups * 0.04)
                      )

  } else {
    if (same_space_on_edge) {
      group_limits <- c(0, num_groups + (use_space * (num_groups + 1)))
    } else {
      use_factor <- 0.06
      group_limits <- c(GetBarPlotStart(num_groups, use_space, use_factor),
                        GetBarplotEnd(num_groups, use_space, use_factor)
                        )
    }
  }

  ## Draw the bar plot
  if (is.null(use_mai)) {
    if (horizontal) {
      use_mai <- horizontal_mai
    } else {
      use_mai <- vertical_mai
    }
  }
  old_par <- par(mai = use_mai, cex = use_cex, lwd = use_lwd)

  plot(1,
       xlim = if (horizontal) numeric_limits else group_limits,
       ylim = if (horizontal) group_limits else numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  axis_ticks <- ManuscriptGrid(horizontal)

  if (space_like_boxplot) {
    if (!(lollipop)) {
      PlotBarplotMat(bars_mat,
                     colors_vec = bar_colors,
                     positions_vec = group_positions,
                     width = 1 / (1 + use_space),
                     horizontal = horizontal
                     )
    }
  } else {
    group_positions <- barplot(bars_mat[, rev(seq_len(ncol(bars_mat)))],
                               horiz     = horizontal,
                               space     = use_space,
                               border    = NA,
                               col       = bar_colors,
                               add       = TRUE,
                               axes      = FALSE,
                               ann       = FALSE,
                               names.arg = rep("", ncol(bars_mat)),
                               plot      = !(lollipop)
                               )
    group_positions <- rev(group_positions)
  }

  if (lollipop) {
    segments(x0   = if (horizontal) par("usr")[[1]] else group_positions,
             x1   = if (horizontal) bars_mat[1, ] else group_positions,
             y0   = if (horizontal) group_positions else par("usr")[[3]],
             y1   = if (horizontal) group_positions else bars_mat[1, ],
             col  = colorRampPalette(use_colors)(10)[[4]],
             lwd  = 2 * par("lwd"),
             xpd  = NA,
             lend = "butt"
             )
    points(x   = if (horizontal) bars_mat[1, ] else group_positions,
           y   = if (horizontal) group_positions else bars_mat[1, ],
           col = use_colors[[2]],
           cex = 1.5,
           pch = 16,
           xpd = NA
           )
  }

  if (!(horizontal) && !(is_proportion) && abbreviate_thousands) {
    num_smaller_than_1000 <- sum(numeric_axis_pos[numeric_axis_pos != 0] < 1000)
    abbreviate_thousands <- num_smaller_than_1000 <= 1
    if (abbreviate_thousands) {
      numeric_axis_labels <- paste0(format(numeric_axis_pos / 1000), "k")
      abbreviated_thousands <- TRUE
    }
  }

  axis(if (horizontal) 1 else 2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, if (horizontal) 0.47 else 0.38, 0),
       gap.axis = 0,
       tcl      = if (horizontal) -0.35 else -0.3,
       las      = 1,
       lwd      = par("lwd")
       )

  ManuscriptAnnotate(group_positions,
                     group_names,
                     modality_label,
                     axis_label,
                     horizontal,
                     use_colors,
                     modality_on_side = modality_on_side,
                     modality_on_bottom = modality_on_bottom
                     )
  box(bty = "l")

  ## Annotate the bars with fractions (counts)
  if (one_row) {
    if (horizontal) {
      text(y      = group_positions,
           x      = par("usr")[[2]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.045),
           labels = bars_mat[1, ],
           adj    = c(0, 0.5),
           xpd    = NA
           )
    } else {
      text(x      = group_positions,
           y      = par("usr")[[4]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.1),
           labels = bars_mat[1, ],
           xpd    = NA
           )
    }

  } else {
    if (is_proportion) {
      dividends <- counts_mat[2, ]
      divisors <- colSums(counts_mat)
    } else {
      dividends <- RoundSmallPercentages(counts_mat[1, ])
      divisors <- counts_mat[2, ]
    }

    line_widths <- strwidth(divisors)
    if (horizontal) {
      x_position <- par("usr")[[2]] + ((par("usr")[[2]] - par("usr")[[1]]) * 0.15)
      y_fraction <- (par("usr")[[4]] - par("usr")[[3]]) * 0.055
      text(y      = group_positions - y_fraction,
           x      = x_position,
           labels = divisors,
           adj    = c(0.5, 0.5),
           xpd    = NA
           )
      segments(x0  = x_position - (line_widths / 2),
               x1  = x_position + (line_widths / 2),
               y0  = group_positions,
               xpd = NA
               )
      text(y      = group_positions + y_fraction,
           x      = x_position,
           labels = dividends,
           adj    = c(0.5, 0.5),
           xpd    = NA
           )
    } else {

      y_position <- par("usr")[[4]] + diff(grconvertY(c(0, 1.25), from = "lines", to = "user"))
      y_fraction <- diff(grconvertY(c(0, 0.456), from = "lines", to = "user"))
      text(x      = group_positions,
           y      = y_position - y_fraction,
           labels = divisors,
           adj    = c(0.5, 0.5),
           xpd    = NA
           )
      segments(x0  = group_positions - (line_widths / 2),
               x1  = group_positions + (line_widths / 2),
               y0  = y_position,
               xpd = NA
               )
      text(x      = group_positions,
           y      = y_position + y_fraction,
           labels = dividends,
           adj    = c(0.5, 0.5),
           xpd    = NA
           )
    }
  }

  par(old_par)
  return(invisible(NULL))
}




ManuscriptViolinBox <- function(plot_df,
                                axis_label           = "",
                                use_width            = horizontal_width,
                                use_height           = horizontal_height,
                                horizontal           = TRUE,
                                use_cex              = manuscript_cex,
                                use_lwd              = manuscript_lwd,
                                use_mai              = NULL,
                                abbreviate_libraries = TRUE,
                                modality_on_bottom   = FALSE,
                                CRISPRa_colors       = manuscript_CRISPRa_colors,
                                CRISPRo_colors       = manuscript_CRISPRo_colors
                                ) {

  ## Prepare colors and labels
  num_groups <- nlevels(plot_df[["Groups_factor"]])
  group_names <- levels(plot_df[["Group"]])
  if ("Brunello" %in% plot_df[["Group"]]) {
    modality_label <- "CRISPRo"
    use_colors <- CRISPRo_colors
  } else if ("Calabrese" %in% plot_df[["Group"]]) {
    modality_label <- "CRISPRa"
    use_colors <- CRISPRa_colors
    if (horizontal) {
      group_names[group_names == "hCRISPRa-v2"] <- "hCRISPRa\n-v2"
    } else {
      group_names[group_names == "hCRISPRa-v2"] <- "hCa-v2"
      if (abbreviate_libraries) {
        group_names[group_names == "Calabrese"] <- "Calab"
      }
    }
  } else {
    modality_label <- ""
    use_colors <- NULL
  }


  ## Prepare the data axis
  numeric_limits <- GetAxisLimits(plot_df[["Numeric_data"]],
                                  column_name = names(plot_df)[[7]],
                                  provide_other_limits = TRUE
                                  )
  if (names(plot_df)[[7]] %in% c("CRISPOR_Doench_efficacy", "GuideScan_efficiency")) {
    numeric_limits[[1]] <- 0
  }


  ## Prepare the groups axis
  spaces_vec <- rep(1.15, num_groups - 1)
  group_positions <- cumsum(c(1, spaces_vec))
  if (horizontal) {
    group_positions <- rev(group_positions)
  }
  group_limits <- c((min(group_positions) - 0.5) - (num_groups * 0.04),
                    max(group_positions) + 0.5 + (num_groups * 0.04)
                    )


  ## Prepare the data points
  set.seed(1)
  jittered_vec <- group_positions[as.integer(plot_df[["Groups_factor"]])] +
                  rnorm(n = nrow(plot_df), mean = 0, sd = 0.03)
  points_alpha <- 0.1
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)


  ## Prepare the raster graphics device
  if (is.null(use_mai)) {
    if (horizontal) {
      use_mai <- horizontal_mai
    } else {
      use_mai <- vertical_mai
    }
  }

  old_par <- par(mai = use_mai, cex = use_cex, lwd = use_lwd)
  main_folder_path <- file.path(output_plots_directory, "Manuscript")

  PDF_device <- dev.cur()
  temp_path <- file.path(main_folder_path, "temp.png")
  temp_width <- use_width - sum(use_mai[c(2, 4)])
  temp_height <- use_height - sum(use_mai[c(1, 3)])

  png(file   = temp_path,
      width  = temp_width,
      height = temp_height,
      units  = "in",
      res    = 900,
      bg     = "transparent"
      )


  ## Prepare the plot for the raster device
  old_par <- par(mai = rep(0, 4), cex = use_cex)
  plot(1,
       xlim = if (horizontal) numeric_limits else group_limits,
       ylim = if (horizontal) group_limits else numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )


  ## Draw the grid
  axis_ticks <- ManuscriptGrid(horizontal)

  ## Draw the violin plots
  vioplot(plot_df[["Numeric_data"]] ~ plot_df[["Groups_factor"]],
          add        = TRUE,
          at         = group_positions,
          pchMed     = NA,
          drawRect   = FALSE,
          col        = colorRampPalette(use_colors)(9)[[2]],
          border     = NA,
          wex        = 1.1,
          axes       = FALSE,
          horizontal = horizontal
          )


  ## Draw the jittered points
  points(x   = if (horizontal) plot_df[["Numeric_data"]] else jittered_vec,
         y   = if (horizontal) jittered_vec else plot_df[["Numeric_data"]],
         cex = 0.4,
         col = paste0(colorRampPalette(use_colors)(9)[[8]], alpha_hex),
         pch = 16
         )


  ## Capture and remove the temporary PNG file
  dev.off()
  raster_array <- readPNG(temp_path)
  file.remove(temp_path)
  dev.set(PDF_device)


  ## Prepare the final plot
  plot(1,
       xlim = if (horizontal) numeric_limits else group_limits,
       ylim = if (horizontal) group_limits else numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )

  rasterImage(raster_array,
              xleft = par("usr")[[1]], xright = par("usr")[[2]],
              ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
              )

  ## Draw the superimposed boxplots
  boxplot(plot_df[["Numeric_data"]] ~ plot_df[["Groups_factor"]],
          add        = TRUE,
          at         = group_positions,
          boxwex     = 0.35,
          outline    = FALSE,
          names      = rep.int("", length(group_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * 3,
          col        = Palify(use_colors[[1]]),
          border     = use_colors[[2]],
          axes       = FALSE,
          lwd        = 1,
          horizontal = horizontal
          )

  ## Draw the axes and main annotations

  axis(if (horizontal) 1 else 2,
       at       = axis_ticks,
       labels   = if (horizontal) axis_ticks else format(axis_ticks),
       mgp      = c(3, 0.38, 0),
       gap.axis = 0,
       tcl      = if (horizontal) -0.35 else -0.3,
       las      = 1,
       lwd      = par("lwd")
       )

  group_sizes <- tabulate(plot_df[["Groups_factor"]])
  one_group <- length(unique(group_sizes)) == 1

  ManuscriptAnnotate(group_positions,
                     group_names,
                     modality_label,
                     axis_label,
                     horizontal,
                     use_colors,
                     modality_on_top = !(modality_on_bottom) && one_group,
                     modality_on_side = !(modality_on_bottom) && !(one_group),
                     modality_on_bottom = modality_on_bottom
                     )

  text(x      = if (one_group) par("usr")[[2]] else group_positions,
       y      = par("usr")[[4]] + diff(grconvertY(c(0, 0.78), from = "lines", to = "user")),
       labels = if (one_group) paste0("(", group_sizes[[1]], ")") else group_sizes,
       adj    = if (one_group) c(1, 0.5) else 0.5,
       xpd    = NA
       )
  box()
  par(old_par)
}






manuscript_barplot_vars <- c(
  "Are_overlapping",
  "all22_SNP_AF_max_Kaviar",
  "Expected_all22_SNP_AF_max_Kaviar",
  "Affects_any_unintended_gene",
  "Affects_any_genes_at_other_loci",
  "Have_homologies"
)

manuscript_violin_vars <- c(
  "GuideScan_specificity",
  "CRISPOR_Doench_efficacy",
  "CRISPOR_3MM_specificity",
  "CRISPOR_4MM_specificity"
)




PrepareManuscriptPlots <- function(CRISPR_df) {

  mat_list_filtered <- sapply(
    manuscript_barplot_vars,
    function(x) {
      BarPlot_Sources(CRISPR_df,
                      x,
                      filter_top4           = TRUE,
                      show_sublibraries     = FALSE,
                      filter_complete_genes = TRUE
                      )[["counts_mat"]]
    }, simplify = FALSE
  )
  mat_list_unfiltered <- sapply(
    manuscript_barplot_vars,
    function(x) {
      BarPlot_Sources(CRISPR_df,
                      x,
                      filter_top4           = TRUE,
                      show_sublibraries     = FALSE,
                      filter_complete_genes = FALSE
                      )[["counts_mat"]]
    }, simplify = FALSE
  )

  df_list_filtered <- sapply(
    manuscript_violin_vars,
    function(x) {
      ViolinBox_Sources(CRISPR_df,
                        x,
                        filter_top4           = TRUE,
                        show_sublibraries     = FALSE,
                        filter_complete_genes = TRUE
                        )
    },
    simplify = FALSE
  )
  df_list_unfiltered <- sapply(
    manuscript_violin_vars,
    function(x) {
      ViolinBox_Sources(CRISPR_df,
                        x,
                        filter_top4           = TRUE,
                        show_sublibraries     = FALSE,
                        filter_complete_genes = FALSE
                        )
    },
    simplify = FALSE
  )

  num_genes_mat <- NumGenesInLibrary()
  mat_list_filtered   <- c(list("Num_genes" = num_genes_mat), mat_list_filtered)
  mat_list_unfiltered <- c(list("Num_genes" = num_genes_mat), mat_list_unfiltered)

  results_list <- list(
    "mat_list_filtered"   = mat_list_filtered,
    "mat_list_unfiltered" = mat_list_unfiltered,
    "df_list_filtered"    = df_list_filtered,
    "df_list_unfiltered"  = df_list_unfiltered
  )
  return(results_list)
}




DrawAllManuscriptPlots <- function(df_mat_list, make_PNGs = FALSE) {

  labels_list <- list(
    "Num_genes"                        = "Number of genes in library",
    "Are_overlapping"                  = expression("Spacing" < "50 bp"),
    "all22_SNP_AF_max_Kaviar"          = "Target polymorphism > 0.1%",
    "Expected_all22_SNP_AF_max_Kaviar" = "Expected to hit alternate allele",
    "Affects_any_unintended_gene"      = "Affect unintended gene",
    "Affects_any_genes_at_other_loci"  = "Affect off-site gene",
    "GuideScan_specificity"            = "GuideScan specificity score",
    "CRISPOR_Doench_efficacy"          = "Efficacy score (Rule Set 2)",
    "Have_homologies"                  = expression("Share subsequences" >= "8 bp"),
    "CRISPOR_3MM_specificity"          = "CRISPOR 3MM specificity",
    "CRISPOR_4MM_specificity"          = "CRISPOR 4MM specificity"
  )

  use_folder <- file.path(output_plots_directory, "Manuscript")

  if (make_PNGs) {
    use_folder <- file.path(use_folder, "PNGs")
    horizontal_height <- PNG_height
    horizontal_width <- PNG_width
    manuscript_wide_width <- PNG_wide_width
    manuscript_wide_height <- PNG_wide_height
    manuscript_cex <- PNG_cex
    manuscript_lwd <- PNG_lwd
    vertical_mai <- PNG_vertical_mai
    wide_mai <- PNG_wide_mai
    manuscript_CRISPRa_colors <- PNG_CRISPRa_colors
  }

  print(make_PNGs)
  print(pdf_width)
  original_PDF_height <- pdf_height

  for (var_name in names(labels_list)) {
    for (SNP_as_percent in c(TRUE, FALSE)) {
      if (SNP_as_percent && (var_name != "Expected_all22_SNP_AF_max_Kaviar")) {
        next
      }
      for (filter_genes in c(TRUE, FALSE)) {

        file_number <- formatC(match(var_name, names(labels_list)),
                               width = 2, flag = "0"
                               )
        file_name <- paste0(file_number, ") ", var_name)

        if (var_name == "Expected_all22_SNP_AF_max_Kaviar") {
          if (SNP_as_percent) {
            file_name <- paste0(file_name, " - percentage")
            use_numeric_limits <- c(0, 2)
          } else {
            file_name <- paste0(file_name, " - absolute number")
            if (filter_genes) {
              use_numeric_limits <- c(0, 1200)
            } else {
              next
            }
          }
        } else {
          use_numeric_limits <- NULL
        }

        make_wide <- (var_name == "Have_homologies") && filter_genes

        if (make_wide) {
          pdf_width <- manuscript_wide_width
          pdf_height <- manuscript_wide_height
        } else {
          pdf_width <- horizontal_width
          pdf_height <- horizontal_height
        }
        sub_folder <- paste0("Comparison - ", if (filter_genes) "filtered genes" else "all genes")
        if (make_PNGs) {
          png(file.path(use_folder,
                        paste0("Comparison - ", if (filter_genes) "filtered genes" else "all genes"),
                        paste0(file_name, ".png")
                        ),
              width  = pdf_width,
              height = pdf_height,
              units  = "in",
              res    = 900
              )
        } else {
          pdf(file.path(use_folder,
                        paste0("Comparison - ", if (filter_genes) "filtered genes" else "all genes"),
                        paste0(file_name, ".pdf")
                        ),
              width = pdf_width,
              height = pdf_height
              )
        }

        if (filter_genes) {
          mat_list <- df_mat_list[["mat_list_filtered"]]
          df_list <- df_mat_list[["df_list_filtered"]]
        } else {
          mat_list <- df_mat_list[["mat_list_unfiltered"]]
          df_list <- df_mat_list[["df_list_unfiltered"]]
        }
        if (make_wide) {
          use_mai <- manuscript_wide_mai
        } else if (make_PNGs) {
          use_mai <- vertical_mai
          use_mai[[1]] <- use_mai[[1]] + PNG_increment
        } else {
          use_mai <- NULL
        }

        if (var_name %in% names(mat_list)) {
          ManuscriptBars(mat_list[[var_name]],
                         axis_label           = labels_list[[var_name]],
                         numeric_limits       = use_numeric_limits,
                         lollipop             = var_name == "Num_genes",
                         horizontal           = FALSE,
                         modality_on_side     = make_wide && !(make_PNGs),
                         use_mai              = use_mai,
                         abbreviate_libraries = !(make_wide),
                         expected_SNP_percent = SNP_as_percent,
                         modality_on_bottom   = make_PNGs,
                         use_cex              = manuscript_cex,
                         use_lwd              = manuscript_lwd,
                         CRISPRa_colors       = manuscript_CRISPRa_colors,
                         CRISPRo_colors       = manuscript_CRISPRo_colors
                         )
        } else if (var_name %in% names(df_list)) {
          ManuscriptViolinBox(df_list[[var_name]],
                              axis_label           = labels_list[[var_name]],
                              horizontal           = FALSE,
                              use_mai              = use_mai,
                              use_width            = pdf_width,
                              use_height           = pdf_height,
                              abbreviate_libraries = !(make_wide),
                              modality_on_bottom   = make_PNGs,
                              use_cex              = manuscript_cex,
                              use_lwd              = manuscript_lwd,
                              CRISPRa_colors       = manuscript_CRISPRa_colors,
                              CRISPRo_colors       = manuscript_CRISPRo_colors
                              )
        }
        dev.off()
      }
    }
  }
}






# Functions for plotting the position relative to the TSS -----------------

AreGenesOutsideTSSRange <- function(input_df, x_range = c(-1000, 1000)) {
  stopifnot(identical(input_df[, "Combined_ID"], input_df[, "Entrez_ID"]))
  distance_vec <- input_df[["Distance_from_TSS"]]
  are_too_low <- distance_vec < x_range[[1]]
  are_too_high <- distance_vec > x_range[[2]]
  are_outside <- are_too_low | are_too_high
  outside_genes <- unique(input_df[["Entrez_ID"]][are_outside])
  are_outside_genes <- input_df[["Entrez_ID"]] %in% outside_genes
  return(are_outside_genes)
}



TSSHistogram <- function(distances_vec,
                         hist_color           = brewer.pal(9, "Blues")[[8]],
                         use_breaks           = 1000,
                         x_range              = c(-1000, 1000),
                         omit_outside_x_range = FALSE,
                         highlight_range      = c(-500, 500),
                         highlight_color      = brewer.pal(9, "Purples")[[2]],
                         modality_text        = "range for CRISPRoff",
                         make_plot            = TRUE,
                         use_y_max            = NULL,
                         use_title            = NULL
                         ) {

  are_NA <- is.na(distances_vec)
  if (any(are_NA)) {
    message(paste0("Out of ", length(distances_vec), " sgRNAs, ",
                   sum(are_NA), " had 'NA' values for the distance from the ",
                   "TSS and were excluded."
                   ))
    distances_vec <- distances_vec[!(are_NA)]
  }

  are_too_low <- distances_vec < x_range[[1]]
  are_too_high <- distances_vec > x_range[[2]]
  are_outside <- are_too_low | are_too_high
  if (any(are_outside)) {
    if (omit_outside_x_range) {
      message(paste0("Out of ", length(distances_vec), " sgRNAs, ",
                     sum(are_outside), " lay outside the range and were excluded."
                     ))
      distances_vec <- distances_vec[!(are_outside)]
    } else {
      distances_vec[are_too_low] <- x_range[[1]]
      distances_vec[are_too_high] <- x_range[[2]]
    }
  }
  message(paste0(length(distances_vec), " sgRNAs (corresponding to ",
                 length(distances_vec) / 4, " genes) are shown in the ",
                 "histogram.\n"
                 ))

  are_highlighted <- (distances_vec >= highlight_range[[1]]) &
                     (distances_vec <= highlight_range[[2]])
  fraction_highlighted <- sum(are_highlighted) / length(distances_vec)
  highlight_label <- as.expression(bquote(plain(.(as.character(round(fraction_highlighted * 100, digits = 1))) *
                                            "% of gRNAs lie within the " * .(modality_text)
                                          )))

  hist_results <- hist(distances_vec,
                       breaks = use_breaks,
                       plot   = FALSE
                       )

  if (make_plot) {

    if (is.null(use_y_max)) {
      y_span <- max(hist_results[["counts"]])
      use_y_max <- max(pretty(c(0, y_span)))
    }

    plot(1,
         xlim = range(x_range),
         ylim = c(-0.02 * use_y_max, use_y_max),
         xaxs = "i",
         yaxs = "i",
         type = "n",
         axes = FALSE,
         mgp  = c(2.5, 1, 0),
         ann  = FALSE
         )

    if (!(is.null(use_title))) {
      title(use_title, cex.main = 1, line = 2.5)
    }

    x_axis_pos <- seq(-1000, 1000, by = 200)
    x_axis_labels <- sapply(as.character(x_axis_pos), as.expression)
    if (!(omit_outside_x_range)) {
      if (x_axis_pos[[1]] == x_range[[1]]) {
        # x_axis_labels[[1]] <- bquote("" <= "\u2212" * .(as.character(abs(x_range[[1]]))))
        x_axis_labels[[1]] <- bquote("" <= .(as.character(x_range[[1]])))
      }
      if (x_axis_pos[[length(x_axis_pos)]] == x_range[[2]]) {
        x_axis_labels[[length(x_axis_labels)]] <- bquote("" >= .(as.character(x_range[[2]])))
      }
    }

    rect(xleft   = highlight_range[[1]],
         xright  = highlight_range[[2]],
         ybottom = par("usr")[[3]],
         ytop    = par("usr")[[4]],
         xpd     = NA,
         col     = highlight_color,
         border  = NA
         )
    text(x      = highlight_range[[1]] + ((highlight_range[[2]] - highlight_range[[1]]) / 2),
         y      = par("usr")[[4]] + diff(grconvertY(c(0, 0.8), from = "lines", to = "user")),
         labels = highlight_label,
         cex    = 0.8,
         xpd    = NA
         )

    segments(y0   = par("usr")[[3]],
             y1   = par("usr")[[4]],
             x0   = x_axis_pos,
             col  = "gray85",
             lend = "butt",
             xpd  = NA
             )

    between_distance <- (x_axis_pos[[2]] - x_axis_pos[[1]])
    between_pos <- seq(from = x_axis_pos[[1]] + (between_distance / 2),
                       to   = x_axis_pos[[length(x_axis_pos)]] - (between_distance / 2),
                       by   = between_distance
                       )
    segments(y0   = par("usr")[[3]],
             y1   = par("usr")[[4]],
             x0   = between_pos,
             col  = "gray85",
             lend = "butt",
             xpd  = NA
             )

    axis(1, at = x_axis_pos, labels = x_axis_labels, mgp = c(3, 0.7, 0), tcl = -0.45)
    axis(2, mgp = c(3, 0.6, 0), tcl = -0.45, las = 1)
    box(bty = "l")
    mtext("Position relative to the TSS", side = 1, line = 2.4)
    mtext("Count (gRNAs)", side = 2, line = 3)

    hist(distances_vec,
         breaks = use_breaks,
         col    = hist_color,
         border = NA,
         freq   = TRUE,
         add    = TRUE,
         axes   = FALSE,
         ylab   = "",
         xpd    = NA
         )

  }

  return(invisible(hist_results))
}




TSS_modalities_list <- list(
  "CRISPRoff" = list(
    "title" = "range for CRISPRoff",
    "range" = c(-500, 500),
    "color" = brewer.pal(9, "Purples")[[2]]
  ),
  "CRISPRoff (narrow)" = list(
    "title" = "narrow range for CRISPRoff",
    "range" = c(-200, 250),
    "color" = brewer.pal(9, "Purples")[[2]]
  ),
  "CRISPRa" = list(
    "title" = "range for CRISPRa",
    "range" = c(-400, 0),
    "color" = brewer.pal(9, "Greens")[[2]]
  ),
  "CRISPRa (narrow)" = list(
    "title" = "narrow range for CRISPRa",
    "range" = c(-150, -75),
    "color" = brewer.pal(9, "Greens")[[2]]
  ),
  "CRISPRi" = list(
    "title" = "range for CRISPRi",
    "range" = c(-50, 300),
    "color" = brewer.pal(9, "Blues")[[2]]
  ),
  "CRISPRi (narrow)" = list(
    "title" = "narrow range for CRISPRi",
    "range" = c(25, 75),
    "color" = brewer.pal(9, "Blues")[[2]]
  )
)



FilterDistanceByGroup <- function(input_df, use_library) {
  input_df[["Distance_from_TSS"]][input_df[["Group"]] %in% use_library]
}



TSSHistogramsForModality <- function(CRISPR_df,
                                     show_modality = "CRISPRoff",
                                     omit_outside_x_range = TRUE
                                     ) {

  stopifnot("TSS_modalities_list" %in% ls(envir = globalenv()))

  use_highlight_range <- TSS_modalities_list[[show_modality]][["range"]]
  title_postfix <- TSS_modalities_list[[show_modality]][["title"]]
  use_color <- TSS_modalities_list[[show_modality]][["color"]]

  are_4sg <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list, show_messages = FALSE)
  are_main_TSS <- (GetMainTSS(CRISPR_df) == CRISPR_df[["AltTSS_ID"]]) %in% TRUE
  sel_distances <- CRISPR_df[["Distance_from_TSS"]][are_4sg & are_main_TSS]
  TSSHistogram(sel_distances,
               use_breaks = 200,
               use_title = expression(phantom("gh") * bold("4sg library") * phantom("gh")),
               modality_text = title_postfix,
               highlight_range = use_highlight_range,
               highlight_color = use_color,
               omit_outside_x_range = omit_outside_x_range
               )

  TSS_dist_df <- ViolinBox_Sources(CRISPR_df,
                                   "Restricted_distance_from_TSS",
                                   show_sublibraries = FALSE,
                                   filter_top4 = TRUE,
                                   draw_plot = FALSE
                                   )

  shared_TSS_dist_df <- TSS_dist_df[!(AreGenesOutsideTSSRange(TSS_dist_df)), ]
  num_shared_genes <- length(unique(shared_TSS_dist_df[["Entrez_ID"]]))

  libraries_vec <- levels(shared_TSS_dist_df[["Group"]])
  library_titles <- lapply(libraries_vec, function(x) {
    as.expression(bquote(phantom("gh") * bold(.(x)) * " " * scriptstyle(" ") *
                        "(" * .(as.character(num_shared_genes)) *
                        " shared genes)" * phantom("gh")
    ))})

  for (i in 1:4) {
    TSSHistogram(distances_vec = FilterDistanceByGroup(shared_TSS_dist_df, libraries_vec[[i]]),
                 use_breaks = 200,
                 use_title = library_titles[[i]],
                 modality_text = title_postfix,
                 highlight_range = use_highlight_range,
                 highlight_color = use_color,
                 omit_outside_x_range = omit_outside_x_range
                 )
  }

  return(invisible(NULL))
}




CreateAllTSSHistograms <- function(CRISPR_df) {
  stopifnot(all(c("TSS_modalities_list", "output_plots_directory") %in% ls(envir = globalenv())))
  for (broad_range in c(TRUE, FALSE)) {
    for (omit_outside_range in c(TRUE, FALSE)) {
      are_narrow <- grepl("narrow", names(TSS_modalities_list), fixed = TRUE)
      if (broad_range) {
        use_modalities <- names(TSS_modalities_list)[!(are_narrow)]
        folder_prefix <- "Broad ranges"
      } else {
        use_modalities <- names(TSS_modalities_list)[are_narrow]
        folder_prefix <- "Narrow ranges"
      }
      if (omit_outside_range) {
        folder_postfix <- ""
      } else {
        folder_postfix <- " - outliers indicated"
      }
      folder_name <- paste0(folder_prefix, folder_postfix)
      for (use_modality in use_modalities) {
        file_name <- paste0("gRNA positions relative to the TSS - ",
                            sub(" (narrow)", "", use_modality, fixed = TRUE),
                            " - ",
                            if (broad_range) "broad range" else "narrow range"
                            )
        use_file_path <- file.path(output_plots_directory,
                                   "Positions relative to the TSS",
                                   folder_name,
                                   paste0(file_name, ".pdf")
                                   )
        pdf(file = use_file_path, width = 7, height = 5.5)
        par("oma" = c(0, 1, 0, 1))
        TSSHistogramsForModality(CRISPR_df,
                                 use_modality,
                                 omit_outside_x_range = omit_outside_range
                                 )
        dev.off()
      }
    }
  }
}






