## 2022-05-23



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
  show_list <- split(essent_df[, "Mean_log2FC"], essent_df[, "Is_essential"])
  if (show_phenotype_score) {
    show_list <- lapply(show_list, function(x) x / num_cell_divisions)
  }
  BeeViolinPlot(show_list,
                violin_colors = c("#c7e7c0", brewer.pal(9, "Purples")[[3]]),
                point_colors  = c("#3b6d4c", "#645a7c"),
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
  SetUpEmptyPlot(nlevels(groups_fac), range(numeric_vec, na.rm = TRUE))

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
