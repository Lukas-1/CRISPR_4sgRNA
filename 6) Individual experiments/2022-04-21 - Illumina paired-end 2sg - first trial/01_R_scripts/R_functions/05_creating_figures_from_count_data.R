## 2022-05-23


# Import packages and source code -----------------------------------------

library("png")
library("pBrackets")




# Functions for data normalization ----------------------------------------

MedianRatioNormalization <- function(all_counts_mat, use_columns = NULL) {

  if (is.null(use_columns)) {
    use_counts_mat <- all_counts_mat
  } else {
    use_counts_mat <- all_counts_mat[, use_columns]
  }
  are_non_zero <- rowSums(use_counts_mat == 0) == 0
  non_zero_counts_mat <- use_counts_mat[are_non_zero, ]
  means_vec <- exp(rowMeans(log(non_zero_counts_mat)))

  ratios_mat <- apply(all_counts_mat[are_non_zero, ], 2, function(x) x / means_vec)
  ratios_vec <- apply(ratios_mat, 2, median)
  results_mat <- do.call(cbind,
                         lapply(seq_len(ncol(all_counts_mat)),
                                function(x) all_counts_mat[, x] / ratios_vec[[x]]
                                )
                         )
  colnames(results_mat) <- colnames(all_counts_mat)
  return(results_mat)
}



# Functions for accessing and processing counts data ----------------------

GetLog2FC <- function(input_df,
                      allow_switch         = TRUE,
                      allow_1MM            = TRUE,
                      choose_rep           = NULL,
                      baseline_indices     = 3:4,
                      intervention_indices = 5:6
                      ) {
  counts_mat <- GetCountsMat(input_df,
                             allow_switch = allow_switch,
                             allow_1MM = allow_1MM,
                             normalization_columns = c(baseline_indices, intervention_indices)
                             )
  if (is.null(choose_rep)) {
    baseline_counts <- rowMeans(counts_mat[, baseline_indices])
    intervention_counts <- rowMeans(counts_mat[, intervention_indices])
  } else {
    baseline_counts <- counts_mat[, baseline_indices[[choose_rep]]]
    intervention_counts <- counts_mat[, intervention_indices[[choose_rep]]]
  }
  log2fc_vec <- log2(intervention_counts / baseline_counts)
  results_df <- data.frame(
    "Log2FC"             = log2fc_vec,
    "Count_baseline"     = baseline_counts,
    "Count_intervention" = intervention_counts,
    stringsAsFactors = FALSE
  )
  return(results_df)
}



GetCountsMat <- function(use_counts_df,
                         allow_switch = TRUE,
                         allow_1MM = TRUE,
                         normalize = TRUE,
                         normalization_columns = NULL
                         ) {
  column_prefix <- paste0(if (allow_switch) "MaySwitch" else "NoSwitch", "_",
                          if (allow_1MM) "xMM" else "0MM", "_"
                          )
  column_names <- grep(paste0("^", column_prefix), names(use_counts_df), value = TRUE)
  counts_mat <- as.matrix(use_counts_df[, column_names])
  if (normalize) {
    counts_mat <- MedianRatioNormalization(counts_mat, use_columns = normalization_columns)
  }
  return(counts_mat)
}



GetAvailableGenes <- function(entrezs_vec, count_column = "NoSwitch_xMM", min_count = 200) {

  stopifnot(all(c("counts_df", "CRISPRoff_df") %in% ls(envir = globalenv())))
  stopifnot(is.numeric(entrezs_vec))
  stopifnot(!(any(duplicated(entrezs_vec))))
  stopifnot(identical(CRISPRoff_df[, "Entrez_ID"], counts_df[, "Entrez_ID"]))

  matches_vec <- match(entrezs_vec, CRISPRoff_df[, "Entrez_ID"])
  message(paste0("\nOf the ", length(entrezs_vec), " Entrez gene IDs, ",
                 sum(is.na(matches_vec)), " were not available in the library.\n"
                 ))
  entrezs_vec <- entrezs_vec[!(is.na(matches_vec))]
  matches_vec <- matches_vec[!(is.na(matches_vec))]

  have_multiple_plasmids <- CRISPRoff_df[, "Num_plasmids_for_Entrez"][matches_vec] >= 2
  message(paste0(sum(have_multiple_plasmids),
                 " genes were represented by multiple ",
                 "plasmids (targeting multiple TSSs).\nOnly the first plasmid ",
                 "was chosen for each gene.\n"
                 ))

  # sample_names <- c("T0_R1", "T0_R2")
  sample_names <- c("Tbefore_R1", "Tbefore_R2")
  column_names <- paste0(count_column, "_", sample_names)
  counts_mat <- as.matrix(counts_df[, column_names])
  total_counts_vec <- rowSums(counts_mat[matches_vec, ])

  if ((!(is.null(min_count))) && (min_count > 0)) {
    too_few_counts <- total_counts_vec < min_count
    message(paste0(sum(too_few_counts), " plasmids were represented by fewer ",
                   "than ", min_count, " reads across both replicates at the ",
                   "T0 timepoint,\nand were excluded from the analysis. ",
                   sum(!(too_few_counts)), " genes remained.\n"
                   ))
    entrezs_vec <- entrezs_vec[!(too_few_counts)]
  }
  return(entrezs_vec)
}



CountsMatForGenes <- function(entrezs_vec, ...) {
  stopifnot("counts_df" %in% ls(envir = globalenv()))
  matches_vec <- match(entrezs_vec, as.integer(counts_df[, "Entrez_ID"]))
  stopifnot(!(anyNA(matches_vec)))
  results_mat <- GetCountsMat(counts_df, ...)[matches_vec, ]
  return(results_mat)
}



Log2FCForGenes <- function(entrezs_vec, ...) {
  stopifnot("counts_df" %in% ls(envir = globalenv()))
  log2fc_df <- GetLog2FC(counts_df, ...)
  matches_vec <- match(entrezs_vec, as.integer(counts_df[, "Entrez_ID"]))
  stopifnot(!(anyNA(matches_vec)))
  results_df <- log2fc_df[matches_vec, ]
  row.names(results_df) <- NULL
  return(results_df)
}



AllSamplesLog2FC <- function(entrezs_vec, normalize_to_reps = 1:2, ...) {

  counts_mat <- CountsMatForGenes(entrezs_vec, ...)
  are_zero_mat <- counts_mat == 0
  are_all_zero <- rowSums(are_zero_mat) == 6
  are_both_zero <- rowSums(are_zero_mat[, normalize_to_reps, drop = FALSE]) == length(normalize_to_reps)
  are_additional <- are_both_zero & !(are_all_zero)
  message(paste0(sum(are_all_zero), " genes had zero reads in all samples ",
                 "and were excluded."
                 ))
  if (any(are_additional)) {
    message(paste0("An additional ", sum(are_additional), " genes had zero ",
                   "reads in the normalization replicates and were also excluded."
                   ))
  }

  counts_mat <- counts_mat[!(are_both_zero), ]
  baseline_means <- rowMeans(counts_mat[, normalize_to_reps, drop = FALSE])
  counts_mat <- do.call(rbind,
                        lapply(seq_len(nrow(counts_mat)),
                               function(x) counts_mat[x, ] / baseline_means[[x]]
                               ))
  results_mat <- log2(counts_mat + 0.001)
  return(results_mat)
}



# Functions for drawing violin + swarm plots for essential genes ----------

RepEssentialViolins <- function(baseline_indices      = 3:4,
                                intervention_indices  = 5:6,
                                use_title             = "CRISPRoff",
                                min_count_at_baseline = 0L,
                                allow_switch          = TRUE,
                                allow_1MM             = TRUE,
                                show_phenotype_score  = TRUE,
                                num_cell_divisions    = 10L,
                                ...
                                ) {

  stopifnot(all(c("essential_entrezs", "non_essential_entrezs") %in% ls(envir = globalenv())))
  use_args <- list(essential_genes       = essential_entrezs,
                   non_essential_genes   = non_essential_entrezs,
                   baseline_indices      = baseline_indices,
                   intervention_indices  = intervention_indices,
                   min_count_at_baseline = min_count_at_baseline,
                   allow_switch          = allow_switch,
                   allow_1MM             = allow_1MM
                   )
  R1_df <- do.call(GetEssentialROCDf, c(use_args, list(choose_rep = 1L)))
  R2_df <- do.call(GetEssentialROCDf, c(use_args, list(choose_rep = 2L)))
  rep_list <- c(split(R1_df[, "Mean_log2FC"], !(R1_df[, "Is_essential"])),
                split(R2_df[, "Mean_log2FC"], !(R2_df[, "Is_essential"]))
                )[c(1, 3, 2, 4)]

  old_mar <- par(mar = c(2.75, 4.3, 5.7, 2))

  if (show_phenotype_score) {
    rep_list <- lapply(rep_list, function(x) x / num_cell_divisions)
  }

  x_positions <- BeeViolinPlot(rep_list, point_cex = 0.3, use_spacing = 0.5, wex = 0.85,
                               violin_colors = rep(c(brewer.pal(9, "Purples")[[3]], "#c7e7c0"), each = 2),
                               point_colors  = rep(c("#7c7198", "#5b8669"), each = 2),
                               ...
                               )
  if (show_phenotype_score) {
    mtext(expression("Phenotype (" * gamma * ") = (log"[2] ~ "fold change) / 10"), side = 2, line = 2.75)
  } else {
    mtext(expression("Log"[2] ~ "fold change"), side = 2, line = 2.75)
  }

  segments(x0  = x_positions[c(1, 3)] - 0.25,
           x1  = x_positions[c(2, 4)] + 0.25,
           y0  = par("usr")[[4]] + diff(grconvertY(c(0, 0.45), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(c("Essential\ngenes", "Non-essential\ngenes"),
        at = c(mean(x_positions[1:2]), mean(x_positions[3:4])),
        line = 0.65, padj = 0, cex = par("cex")
        )

  mtext(text = rep(c("R1", "R2"), 2), at = x_positions, side = 1, line = 0.7, cex = par("cex"))

  title(as.expression(use_title), cex.main = par("cex"), line = 4.2)
  par(old_mar)
  return(invisible(NULL))
}



AllTimesLogFCViolins <- function(counts_mat,
                                 point_cex   = 0.4,
                                 use_spacing = 0.225,
                                 wex         = 0.85,
                                 upper_bound = 4,
                                 lower_bound = -8,
                                 y_limits    = NULL,
                                 use_title   = NULL,
                                 ...
                                 ) {

  stopifnot(ncol(counts_mat) == 5)

  if (is.null(y_limits) && (!(is.null(upper_bound))) && !(is.null(lower_bound))) {
    y_space <- (upper_bound - lower_bound) * 0.015
    y_limits <- c(lower_bound - y_space, upper_bound)
  }
  old_mar <- par(mar = c(4.1, 4.3, 4, 2.1))
  x_positions <- BeeViolinPlot(as.list(data.frame(counts_mat)),
                               point_cex     = point_cex,
                               use_spacing   = use_spacing,
                               groups_vec    = c(1, 2, 2, 3, 3),
                               wex           = wex,
                               lower_bound   = lower_bound,
                               upper_bound   = upper_bound,
                               y_limits      = y_limits,
                               ...
                               )

  mtext(text = c("R2", "R1", "R2", "R1", "R2"),
        at = x_positions, side = 1, line = 0.7, cex = par("cex")
        )

  segments(x0  = x_positions[c(2, 4)],
           x1  = x_positions[c(3, 5)],
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.85), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(c("Tbefore", "T0", "T12"), side = 1,
        at = c(1, mean(x_positions[2:3]), mean(x_positions[4:5])),
        line = 2.15, padj = 0, cex = par("cex")
        )

  mtext(expression("Log"[2] ~ "fold change (relative to R1 / Tbefore)"),
        side = 2, line = 2.5, cex = par("cex")
        )
  if (!(is.null(use_title))) {
    title(use_title, cex.main = par("cex"))
  }
  par(old_mar)
  return(invisible(NULL))
}



RawCountsViolins <- function(counts_mat,
                             point_cex       = 0.4,
                             use_spacing     = 0.225,
                             wex             = 0.85,
                             upper_bound     = 6000,
                             use_title       = NULL,
                             filter_all_zero = TRUE,
                             ...
                             ) {

  stopifnot(ncol(counts_mat) == 6)
  if (filter_all_zero) {
    are_all_zero <- rowSums(counts_mat == 0) == 6
    counts_mat <- counts_mat[!(are_all_zero), ]
    message(paste0(sum(are_all_zero), " genes had zero reads in all samples ",
                   "and were omitted."
                   )
            )
  }

  old_mar <- par(mar = c(4.1, 4.3, 4, 2.1))
  x_positions <- BeeViolinPlot(as.list(data.frame(counts_mat)),
                               point_cex     = point_cex,
                               use_spacing   = use_spacing,
                               groups_vec    = c(1, 1, 2, 2, 3, 3),
                               wex           = wex,
                               upper_bound   = upper_bound,
                               ...
                               )

  mtext(text = rep(c("R1", "R2"), times = 3),
        at = x_positions, side = 1, line = 0.7, cex = par("cex")
        )

  segments(x0  = x_positions[c(1, 3, 5)],
           x1  = x_positions[c(2, 4, 6)],
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.85), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(c("Tbefore", "T0", "T12"), side = 1,
        at = c(mean(x_positions[1:2]), mean(x_positions[3:4]), mean(x_positions[5:6])),
        line = 2.15, padj = 0, cex = par("cex")
        )

  mtext("Normalized count", side = 2, line = 2.75, cex = par("cex"))
  if (!(is.null(use_title))) {
    title(use_title, cex.main = par("cex"))
  }
  par(old_mar)
  return(invisible(NULL))
}



ViolinPlotEssentialDf <- function(essent_df,
                                  use_title            = "CRISPRoff",
                                  show_phenotype_score = TRUE,
                                  num_cell_divisions   = 10L,
                                  ...
                                  ) {
  old_mar <- par(mar = c(4, 4.3, 4, 2.1))
  show_list <- split(essent_df[, "Mean_log2FC"], !(essent_df[, "Is_essential"]))
  if (show_phenotype_score) {
    show_list <- lapply(show_list, function(x) x / num_cell_divisions)
  }
  BeeViolinPlot(show_list,
                violin_colors = c(brewer.pal(9, "Purples")[[3]], "#c7e7c0"),
                point_colors  = c("#645a7c", "#3b6d4c"),
                point_cex     = 0.35,
                use_spacing   = 0.6,
                wex           = 0.95,
                ...
                )
  if (show_phenotype_score) {
    mtext(expression("Phenotype (" * gamma * ") = (log"[2] ~ "fold change) / 10"),
          side = 2, line = 2.75
          )
  } else {
    mtext(expression("Log"[2] ~ "fold change"), side = 2, line = 2.75)
  }
  mtext(c("Essential\ngenes", "Non-essential\ngenes"),
        side = 1, line = 1.65, at = 1:2
        )
  title(use_title, cex.main = par("cex"))
  par(old_mar)
}



# Functions for plotting read counts for individual plasmids --------------

PlotCountsForPlasmid <- function(plasmid_ID, ...) {
  stopifnot("counts_df" %in% ls(envir = globalenv()))
  plasmid_index <- match(plasmid_ID, counts_df[, "Plasmid_ID"])
  if (is.na(plasmid_index)) {
    stop("A plasmid with the ID '", plasmid_ID, "' was not found!")
  }
  counts_mat <- GetCountsMat(counts_df, ...)
  group_names <- c("Pre-screen", "T0", "T12")
  groups_fac <- factor(rep(group_names, each = 2), levels = group_names)
  counts_vec <- counts_mat[plasmid_index, ]
  PlotCounts(counts_vec, groups_fac)
  title(plasmid_ID, font.main = 4)
}



PlotCounts <- function(numeric_vec,
                       groups_fac,
                       group_colors = c("Greys", "Blues", "Purples"),
                       group_labels = NULL
                       ) {

  group_means <- tapply(numeric_vec, groups_fac, mean)

  if (is.null(group_labels)) {
    group_labels <- levels(groups_fac)
  }

  ## Prepare the colors
  bar_colors   <- vapply(group_colors, function(x) brewer.pal(9, x)[[5]], "")
  point_colors <- vapply(group_colors, function(x) brewer.pal(9, x)[[9]], "")

  ## Set up the plot canvas
  SetUpBoxPlot(nlevels(groups_fac), range(numeric_vec, na.rm = TRUE))

  group_positions <- seq_len(nlevels(groups_fac))
  mtext(group_labels, at = group_positions, side = 1, line = 0.5, cex = par("cex"))

  ## Draw the bars
  bar_width <- 0.6
  rect(xleft   = group_positions - (bar_width / 2),
       xright  = group_positions + (bar_width / 2),
       ybottom = 0,
       ytop    = group_means,
       col     = bar_colors,
       border  = NA,
       xpd     = NA
       )

  ## Draw the points
  beeswarm_df <- beeswarm(numeric_vec ~ groups_fac,
                          at       = group_positions,
                          col      = point_colors,
                          priority = "density",
                          do.plot  = FALSE,
                          cex      = par("cex")
                          )
  points(beeswarm_df[, "x"],
         beeswarm_df[, "y"],
         pch = 16,
         col = beeswarm_df[, "col"],
         cex = par("cex") * 1.2,
         xpd = NA
         )

  return(invisible(NULL))
}




# Functions for plotting ROC curves (for gene essentiality) ---------------

GetEssentialROCDf <- function(essential_genes,
                              non_essential_genes,
                              min_count_at_baseline = 0L,
                              ...
                              ) {

  ROC_df <- data.frame(
    "Entrez_ID"    = c(essential_genes, non_essential_genes),
    "Is_essential" = c(rep(TRUE, length(essential_genes)),
                       rep(FALSE, length(non_essential_genes))
                       ),
    stringsAsFactors = FALSE
  )
  log2fc_df <- Log2FCForGenes(ROC_df[, "Entrez_ID"], ...)
  ROC_df[, "Mean_log2FC"] <- log2fc_df[, "Log2FC"]

  are_nan <- is.nan(ROC_df[, "Mean_log2FC"])
  message(paste0(sum(are_nan), " genes had a read count of zero in both ",
                     "conditions and had to be excluded."
                 ))
  are_included <- !(are_nan)
  if (min_count_at_baseline > 0) {
    are_too_sparse <- log2fc_df[, "Count_baseline"] < min_count_at_baseline
    num_additional <- sum(are_too_sparse & are_included)
    if (num_additional > 0) {
      message(paste0("An additional ", num_additional, " genes had an average ",
                     "normalized count below ", min_count_at_baseline, " at ",
                     "baseline and were also excluded."
                     ))
    }
    are_included[are_too_sparse] <- FALSE
  }
  message(paste0(sum(are_included), " genes remained."))

  ROC_df <- ROC_df[are_included, ]
  new_order <- order(ROC_df[, "Mean_log2FC"])
  ROC_df <- ROC_df[new_order, ]
  row.names(ROC_df) <- NULL

  ROC_df <- MakeROCDf(ROC_df, "Mean_log2FC")
  return(ROC_df)
}



PlotEssentialROCDf <- function(use_ROC_df, use_title = "Gene essentiality with CRISPRoff") {

  numbers_text <- paste0(sum(use_ROC_df[, "Is_essential"]), " essential genes and ",
                         sum(!(use_ROC_df[, "Is_essential"])), " non-essential genes"
                         )

  old_mar <- par(mar = c(3.7, 3.95, 4.1, 1.75))

  PlotROCDf(use_ROC_df, use_title = use_title, flip = TRUE)

  title(as.expression(use_title), cex.main = par("cex"), line = 2.5)
  mtext(numbers_text, line = 0.75, cex = par("cex") * 0.75, adj = 0.5)

  par(old_mar)
  return(invisible(NULL))
}




# Helper functions (required by ReplicateScatterPlot) ---------------------

ConcatenateExpressions <- function(expression_list, my_sep = "  \u2013  ") {
  literal_strings <- vapply(expression_list, StripExpression, "")
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}


VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(phantom("gh")), use_expression, expression(phantom("gh")))
  return(ConcatenateExpressions(my_list, my_sep = ""))
}


StripExpression <- function(my_expression) {
  if (is.character(my_expression)) {
    literal_string <- paste0("\"", capture.output(cat(my_expression)), "\"")
  } else {
    literal_string <- paste0(capture.output(my_expression), collapse = "")
    if (substr(literal_string, 1L, 11L) == "expression(") {
      literal_string <- substr(literal_string, 12L, nchar(literal_string) - 1L)
    }
  }
  return(literal_string)
}



DrawSideLegend <- function(labels_list,
                           use_colors,
                           border_colors        = NULL,
                           use_pch              = 16,
                           use_point_size       = 1.2,
                           lines_x_start        = 0.75,
                           y_mid                = 0.5,
                           small_gap_size       = 1.25,
                           large_gap_multiplier = 1.75,
                           point_x_start        = 0.15
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing legend
  small_gap <- diff(grconvertY(c(0, small_gap_size), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_multiplier

  if (all(lengths(labels_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(labels_list))
    are_first <- rep(TRUE, length(labels_list))
  } else {
    are_first <- unlist(lapply(labels_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }))
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")

  x_text  <- 1 + diff(grconvertX(c(0, lines_x_start), from = "lines", to = "npc"))
  x_point <- 1 + diff(grconvertX(c(0, lines_x_start + point_x_start), from = "lines", to = "npc"))

  ## Draw legend
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(labels_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )
  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  points(x   = rep(grconvertX(x = x_point, from = "npc", to = "user"), length(labels_list)),
         y   = tapply(y_pos, groups_vec, mean),
         cex = use_point_size,
         pch = use_pch,
         col = if (!(is.null(border_colors))) border_colors else use_colors,
         bg  = use_colors,
         xpd = NA
         )

  return(invisible(NULL))
}




# Functions for creating replicate scatter plots --------------------------

ReplicateScatterPlot <- function(input_df,
                                 show_phenotype_score = FALSE,
                                 lower_bound          = -8 * (if (show_phenotype_score) 0.1 else 1),
                                 upper_bound          = 8  * (if (show_phenotype_score) 0.1 else 1),
                                 axis_limits          = c(lower_bound, upper_bound),
                                 highlight_essential  = TRUE,
                                 use_blomen_hart      = TRUE,
                                 highlight_NT         = TRUE,
                                 use_title            = NULL,
                                 embed_PNG            = FALSE
                                 ) {

  required_objects <- c("essential_df", "essentials_2020Q2_df",
                        "non_essentials_2020Q2_df"
                        )
  stopifnot(all(required_objects %in% ls(envir = globalenv())))
  xy_mat <- as.matrix(input_df[, c("Rep1_data", "Rep2_data")])

  ## Remove NaN values (0 divided by 0)
  have_NaN <- rowSums(is.nan(xy_mat)) != 0
  xy_mat <- xy_mat[!(have_NaN), ]
  highlight_genes <- highlight_NT || highlight_essential
  if (highlight_genes) {
    entrezs_vec <- input_df[, "Entrez_ID"][!(have_NaN)]
    are_NT <- input_df[, "Is_NT"][!(have_NaN)]
  }
  message(paste0(sum(have_NaN), " plasmids had NaN Log2FC values (i.e. 0 divided",
                 " by 0) in one or both replicate experiments and were excluded."
                 ))

  ## Adjust axis limits
  xy_list <- lapply(1:2, function(x) {
    BringWithinLimits(xy_mat[, x], lower_bound = lower_bound, upper_bound = upper_bound)
  })
  xy_mat <- cbind(xy_list[[1]][["curtailed_vec"]], xy_list[[2]][["curtailed_vec"]])
  xy_range <- range(xy_mat[is.finite(xy_mat)])
  if (is.null(axis_limits)) {
    axis_limits <- xy_range
  }
  xy_space <- (axis_limits[[2]] - axis_limits[[1]]) * 0.025
  xy_lim <- c(axis_limits[[1]] - xy_space, axis_limits[[2]] + xy_space)

  ## Prepare axis tick labels
  tick_locations <- pretty(axis_limits, n = 6)
  tick_labels_list <- lapply(1:2, function(i) {
    CurtailedAxisLabels(tick_locations,
                        lower_bound          = lower_bound,
                        upper_bound          = upper_bound,
                        lower_bound_enforced = xy_list[[i]][["lower_bound_enforced"]],
                        upper_bound_enforced = xy_list[[i]][["lower_bound_enforced"]]
                        )
  })
  ## Prepare axis labels
  if (show_phenotype_score) {
    axis_labels_list <- list(
      expression("Phenotype score (" * gamma * ") \u2013 replicate 1"),
      expression("Phenotype score (" * gamma * ") \u2013 replicate 1")
    )
  } else {
    axis_labels_list <- list(
      expression("Log"[2] * " fold change \u2013 replicate 1"),
      expression("Log"[2] * " fold change \u2013 replicate 2")
    )
  }

  ## Prepare points that should be highlighted
  if (highlight_genes) {
    AddSum <- function(logical_vec) paste0("(", sum(logical_vec), ")")
    if (highlight_essential) {
      if (use_blomen_hart) {
        highlight_mat <- cbind(
          "is_essential"     = entrezs_vec %in% essentials_2020Q2_df[, "Entrez_ID"],
          "is_non_essential" = entrezs_vec %in% non_essentials_2020Q2_df[, "Entrez_ID"]
        )
        highlight_colors <- c("#6810c6", "#18a008")
        labels_list <- list(
          "essential" = c("Essential", "genes", AddSum(highlight_mat[, "is_essential"])),
          "non-essential" = c("Non-essential", "genes", AddSum(highlight_mat[, "is_non_essential"]))
        )
      } else {
        matches_vec <- match(entrezs_vec, essential_df[, "Entrez_ID"])
        are_essential <- essential_df[, "CRISPR_common"][matches_vec] %in% "Essential"
        highlight_mat <- cbind("is_essential" = are_essential)
        highlight_colors <- "#6810c6"
        labels_list <- list("essential" = c("Essential", "genes", AddSum(are_essential)))
      }
    } else {
      labels_list <- c()
      highlight_colors <- c()
    }
    if (highlight_NT) {
      if (highlight_essential) {
        highlight_mat <- cbind(highlight_mat, "is_NT" = are_NT)
      } else {
        highlight_mat <- cbind("is_NT" = are_NT)
      }
      highlight_colors <- c(highlight_colors, "#004ec2")
      labels_list <- c(labels_list, list("NT" = c("Non-targeting", "controls", AddSum(are_NT))))
    }
    stopifnot(!(any(rowSums(highlight_mat) > 1)))
    are_highlighted <- rowSums(highlight_mat) == 1
    set.seed(1)
    scrambled_indices <- sample(which(are_highlighted))
    are_highlighted_mat <- highlight_mat[scrambled_indices, , drop = FALSE]
    highlighted_xy_mat <- xy_mat[scrambled_indices, ]
    point_colors <- adjustcolor(highlight_colors, alpha.f = 0.7)
    colors_vec <- vapply(seq_len(nrow(are_highlighted_mat)),
                         function(x) point_colors[[which(are_highlighted_mat[x, ])]],
                         ""
                         )
  } else {
    are_highlighted <- rep(FALSE, nrow(xy_mat))
  }

  ## Prepare plot margins
  if (highlight_genes) {
    use_mar <- c(3.75, 3.75, 3.75, 7.5)
  } else {
    use_mar <- c(3.75, 3.75, 3.75, 2.1)
  }
  old_mar <- par(mar = use_mar)

  if (embed_PNG) {
    PDF_mar <- par("mar")
    PDF_device <- dev.cur()
    temp_path <- file.path(figures_dir, "temp.png")
    temp_width  <- par("pin")[[1]]
    temp_height <- par("pin")[[2]]
    current_par <- par(no.readonly = TRUE)
    png(filename = temp_path,
        width    = temp_width,
        height   = temp_height,
        units    = "in",
        res      = 900,
        bg       = "transparent"
        )
    par(lwd = current_par[["lwd"]])
    par(cex = current_par[["cex"]])
    par(mar = rep(0, 4))
  }

  ## Set up plot canvas
  plot(1, xlim = xy_lim, ylim = xy_lim, xaxs = "i", yaxs = "i", type = "n",
       axes = FALSE, ann = FALSE
       )
  abline(v = 0, h = 0, col = "gray85", lend = "butt")

  ## Draw points
  points(xy_mat[!(are_highlighted), ],
         pch = 16,
         cex = 0.3,
         col = adjustcolor("black", alpha.f = 0.35),
         xpd = NA
         )
  if (highlight_genes) {
    points(highlighted_xy_mat,
           pch = 16,
           cex = 0.3,
           col = colors_vec,
           xpd = NA
           )
  }

  if (embed_PNG) {
    dev.off()
    raster_array <- png::readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    par(PDF_mar)
    plot(1, xlim = xy_lim, ylim = xy_lim, xaxs = "i", yaxs = "i", type = "n",
         axes = FALSE, ann = FALSE
         )
    rasterImage(raster_array,
                xleft   = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop   = par("usr")[[4]]
                )
  }

  ## Annotate plot
  mtext(axis_labels_list[[1]], side = 1, line = 2)
  mtext(axis_labels_list[[1]], side = 2, line = 2.3)
  if (!(is.null(use_title))) {
    title(use_title, cex.main = par("cex"))
  }

  ## Draw axes
  tick_locations <- pretty(axis_limits, n = 6)
  for (i in 1:2) {
    axis(i,
         at       = tick_locations,
         labels   = tick_labels_list[[i]],
         mgp      = c(3, if (i == 1) 0.4 else 0.55, 0),
         tcl      = -0.375,
         las      = 1,
         lwd      = par("lwd"),
         cex.axis = par("cex") / 0.7
         )
  }
  box()
  if (highlight_genes) {
    DrawSideLegend(labels_list,
                   use_colors = adjustcolor(highlight_colors, alpha.f = 0.85),
                   use_point_size = 1, point_x_start = 0.2, lines_x_start = 0.6
                   )
  }
  par(old_mar)
  return(invisible(NULL))
}



Log2FCScatterPlot <- function(allow_switch         = FALSE,
                              allow_1MM            = TRUE,
                              baseline_indices     = 3:4,
                              intervention_indices = 5:6,
                              show_phenotype_score = FALSE,
                              num_cell_divisions   = 10L,
                              ...
                              ) {

  stopifnot(all(c("counts_df", "CRISPRoff_df") %in% ls(envir = globalenv())))

  use_args <- list(input_df             = counts_df,
                   baseline_indices     = baseline_indices,
                   intervention_indices = intervention_indices,
                   allow_switch         = allow_switch,
                   allow_1MM            = allow_1MM
                   )
  R1_df <- do.call(GetLog2FC, c(use_args, list(choose_rep = 1L)))
  R2_df <- do.call(GetLog2FC, c(use_args, list(choose_rep = 2L)))

  stopifnot(identical(counts_df[, "Plasmid_ID"], CRISPRoff_df[, "Plasmid_ID"]))

  replicates_df <- data.frame(
    counts_df[, c("Plasmid_ID", "Entrez_ID", "Gene_symbol")],
    "Is_NT" = CRISPRoff_df[, "gene"] == "negative_control",
    "Rep1_data" = if (show_phenotype_score) (R1_df[, "Log2FC"] / num_cell_divisions) else R1_df[, "Log2FC"],
    "Rep2_data" = if (show_phenotype_score) (R2_df[, "Log2FC"] / num_cell_divisions) else R2_df[, "Log2FC"],
    stringsAsFactors = FALSE
  )

  ReplicateScatterPlot(replicates_df, show_phenotype_score = show_phenotype_score, ...)
}




# Functions for creating bar charts for missing plasmids ------------------

FormatPercentages <- function(percentages_vec) {
  percentages_strings <- rep(NA, length(percentages_vec))
  percentages_strings[percentages_vec < 1] <- signif(percentages_vec[percentages_vec < 1], digits = 2)
  are_single_digit <- (percentages_vec < 10) & (percentages_vec >= 1)
  percentages_strings[are_single_digit] <- format(percentages_vec[are_single_digit], digits = 2)
  percentages_strings[percentages_vec > 10] <- round(percentages_vec[percentages_vec > 10], digits = 1)
  percentages_strings[percentages_vec > 99] <- round(percentages_vec[percentages_vec > 99], digits = 2)
  return(percentages_strings)
}


FourBars <- function(bar_values,
                     min_count    = 20L,
                     use_y_limits = NULL,
                     library_size = NULL,
                     title_text   = NULL
                     ) {

  stopifnot("columns_list" %in% ls(envir = globalenv()))
  stopifnot(length(bar_values) == 4)

  ## Prepare for drawing bar plots
  color_scheme <- c(brewer.pal(9, "Greys")[[5]], brewer.pal(9, "Blues")[[7]])
  bar_colors <- color_scheme[c(1, 2, 1, 2)]

  top_labels <- bar_values
  if (!(is.null(library_size))) {
    percentages_vec <- FormatPercentages((bar_values / library_size) * 100)
    top_labels <- paste0(top_labels, " (", percentages_vec, "%)")
  }
  reads_labels <- c(
    expression("" > "0 reads"),
    as.expression(bquote("" >= .(as.character(min_count)) * " reads"))
  )
  reads_labels <- sapply(bottom_labels, VerticalAdjust)
  x_positions <- RepositionByGroups(c(1, 1, 2, 2))

  ## Set up plot canvas
  old_mar <- par(mar = c(5.6, 4.1, 3.7, 2.1))
  SetUpBoxPlot(num_groups = 4L,
               data_range = range(bar_values, na.rm = TRUE),
               draw_box = FALSE,
               use_y_limits = use_y_limits
               )
  mtext(expression("Number of absent" * scriptstyle(" ") * "/" *
                   scriptscriptstyle(" ") * "scarce plasmids"
                   ),
        side = 2, line = 2.4
        )
  if (!(is.null(title_text))) {
    title(title_text, cex.main = 1)
  }

  ## Draw bars
  bar_width <- 0.6
  rect(xleft   = x_positions - (bar_width / 2),
       xright  = x_positions + (bar_width / 2),
       ybottom = 0,
       ytop    = bar_values,
       col     = bar_colors,
       border  = NA,
       xpd     = NA
       )
  box(bty = "l")

  ## Annotate x axis
  mtext("Both sgRNAs required", side = 1, line = 1.4)
  mtext(rep(c("No", "Yes"), 2), at = x_positions, side = 1, line = 0.15)
  line_y <- par("usr")[[3]] - diff(grconvertY(c(0, 3), from = "lines", to = "user"))
  for (i in 1:2) {
    brackets(x1 = x_positions[[c(1, 3)[[i]]]] - (bar_width / 2), y1 = line_y,
             x2 = x_positions[[c(2, 4)[[i]]]] + (bar_width / 2), y2 = line_y,
             h = -(diff(grconvertY(c(0, 0.5), from = "lines", to = "user"))),
             curvature = 0.3, xpd = NA
             )
  }
  mtext(reads_labels, at = c(mean(x_positions[1:2]), mean(x_positions[3:4])),
        side = 1, line = 3.8
        )

  ## Annotate bars
  text(x = x_positions,
       y = bar_values + diff(grconvertY(c(0, 0.3), from = "lines", to = "user")),
       labels = top_labels, cex = 0.4, col = "gray80", xpd = NA
       )

  par(old_mar)
  return(invisible(NULL))
}


