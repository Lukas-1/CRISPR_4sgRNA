## 2022-05-23


# Import packages and source code -----------------------------------------

library("png")
library("pBrackets")
library("DescTools")
library("RColorBrewer")



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



GetAvailableGenes <- function(entrezs_vec,
                              count_column = "NoSwitch_xMM",
                              min_count = 200,
                              verbose = TRUE
                              ) {

  stopifnot(all("CRISPRoff_df" %in% ls(envir = globalenv())))
  stopifnot(is.numeric(entrezs_vec))
  stopifnot(!(any(duplicated(entrezs_vec))))

  matches_vec <- match(entrezs_vec, CRISPRoff_df[, "Entrez_ID"])
  if (verbose) {
    message(paste0("\nOf the ", length(entrezs_vec), " Entrez gene IDs, ",
                   sum(is.na(matches_vec)), " were not available in the library.\n"
                   ))
  }

  entrezs_vec <- entrezs_vec[!(is.na(matches_vec))]
  matches_vec <- matches_vec[!(is.na(matches_vec))]

  have_multiple_plasmids <- CRISPRoff_df[, "Num_plasmids_for_Entrez"][matches_vec] >= 2
  if (verbose) {
    message(paste0(sum(have_multiple_plasmids),
                   " genes were represented by multiple ",
                   "plasmids (targeting multiple TSSs).\nOnly the first plasmid ",
                   "was chosen for each gene.\n"
                   ))
  }

  if ((!(is.null(min_count))) && (min_count > 0)) {
    stopifnot(identical(CRISPRoff_df[, "Entrez_ID"], counts_df[, "Entrez_ID"]))
    sample_names <- c("Tbefore_R1", "Tbefore_R2")
    column_names <- paste0(count_column, "_", sample_names)
    counts_mat <- as.matrix(counts_df[, column_names])
    total_counts_vec <- rowSums(counts_mat[matches_vec, ])
    stopifnot("counts_mat" %in% ls(envir = globalenv()))
    too_few_counts <- total_counts_vec < min_count
    if (verbose) {
      message(paste0(sum(too_few_counts), " plasmids were represented by fewer ",
                     "than ", min_count, " reads across both replicates at the ",
                     "T0 timepoint,\nand were excluded from the analysis. ",
                     sum(!(too_few_counts)), " genes remained.\n"
                     ))
    }
    entrezs_vec <- entrezs_vec[!(too_few_counts)]
  }
  return(entrezs_vec)
}



GetEssentialGenes <- function(use_blomen_hart = TRUE) {
  required_objects <- c("essential_df", "essentials_2020Q2_df",
                        "non_essentials_2020Q2_df"
                        )
  stopifnot(all(required_objects %in% ls(envir = globalenv())))
  if (use_blomen_hart) {
    essential_entrezs <- GetAvailableGenes(
      essentials_2020Q2_df[, "Entrez_ID"], min_count = 0, verbose = FALSE
    )
    non_essential_entrezs <- GetAvailableGenes(
      non_essentials_2020Q2_df[, "Entrez_ID"], min_count = 0, verbose = FALSE
    )
  } else {
    essential_entrezs <- GetAvailableGenes(
      essential_df[, "Entrez_ID"][essential_df[, "Three_categories"] %in% "Essential"],
      min_count = 0, verbose = FALSE
    )
    non_essential_entrezs <- GetAvailableGenes(
      essential_df[, "Entrez_ID"][essential_df[, "Three_categories"] %in% "Non-essential"],
      min_count = 0, verbose = FALSE
    )
  }
  results_list <- list("essential_entrezs" = essential_entrezs,
                       "non_essential_entrezs" = non_essential_entrezs
                       )
  return(results_list)
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



StandardLog2FCDf <- function(input_counts_df) {
  stopifnot(identical(input_counts_df[, "Entrez_ID"], CRISPRoff_df[, "Entrez_ID"]))
  results_df <- GetLog2FC(input_counts_df,
                          baseline_indices = 1:2,
                          allow_switch = FALSE
                          )[, c("Count_baseline", "Count_intervention", "Log2FC")]
  names(results_df) <- paste0("Mean_",
                              substr(tolower(names(results_df)), 1, 1),
                              substr(names(results_df), 2, nchar(names(results_df)))
                              )
  results_df <- data.frame(
    input_counts_df[, c("Plasmid_ID", "Gene_symbol", "Entrez_ID")],
    "Is_NT" = CRISPRoff_df[, "gene"] == "negative_control",
    results_df,
    "Log2FC_rep1" = GetLog2FC(input_counts_df,
                              baseline_indices = 1:2,
                              allow_switch = FALSE,
                              choose_rep = 1
                              )[, "Log2FC"],
    "Log2FC_rep2" = GetLog2FC(input_counts_df,
                              baseline_indices = 1:2,
                              allow_switch = FALSE,
                              choose_rep = 2
                              )[, "Log2FC"]
  )
  return(results_df)
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
                                use_mar               = c(2.75, 4.3, 5.7, 2),
                                y_axis_label          = NULL,
                                y_label_line          = 2.75,
                                rep_label_line        = 0.7,
                                genes_label_line      = 0.65,
                                title_line            = 4.2,
                                point_cex             = 0.3,
                                write_rep             = FALSE,
                                wex                   = 0.85,
                                use_blomen_hart       = TRUE,
                                bracket_color         = "black",
                                ...
                                ) {

  essential_list <- GetEssentialGenes(use_blomen_hart)

  use_args <- list(essential_genes       = essential_list[["essential_entrezs"]],
                   non_essential_genes   = essential_list[["non_essential_entrezs"]],
                   baseline_indices      = baseline_indices,
                   intervention_indices  = intervention_indices,
                   min_count_at_baseline = min_count_at_baseline,
                   allow_switch          = allow_switch,
                   allow_1MM             = allow_1MM
                   )
  R1_df <- do.call(GetEssentialROCDf, c(use_args, list(choose_rep = 1L)))
  R2_df <- do.call(GetEssentialROCDf, c(use_args, list(choose_rep = 2L)))
  rep_list <- c(split(R1_df[, "Log2FC"], !(R1_df[, "Is_essential"])),
                split(R2_df[, "Log2FC"], !(R2_df[, "Is_essential"]))
                )[c(1, 3, 2, 4)]
  names(rep_list) <- c("Essential R1", "Essential R2",
                       "Non-essential R1", "Non-essential R2"
                       )

  old_mar <- par(mar = use_mar)

  if (show_phenotype_score) {
    rep_list <- lapply(rep_list, function(x) x / num_cell_divisions)
  }

  x_positions <- BeeViolinPlot(rep_list, point_cex = point_cex, use_spacing = 0.5, wex = wex,
                               violin_colors = rep(c(brewer.pal(9, "Purples")[[3]], "#c7e7c0"), each = 2),
                               point_colors  = rep(c("#7c7198", "#5b8669"), each = 2),
                               border_colors = rep(c("#d1cddb", "#bfd4c6"), each = 2),
                               use_swarm = use_blomen_hart,
                               cloud_alpha = 0.2, cloud_sd = 0.04,
                               ...
                               )
  if (is.null(y_axis_label)) {
    if (show_phenotype_score) {
      y_axis_label <- expression("Phenotype (" * gamma * ") = (log"[2] ~ "fold change) / 10")
    } else {
      y_axis_label <- expression("Log"[2] ~ "fold change")
    }
  }
  mtext(y_axis_label, side = 2, line = y_label_line, cex = par("cex"))

  segments(x0  = x_positions[c(1, 3)] - 0.25,
           x1  = x_positions[c(2, 4)] + 0.25,
           y0  = par("usr")[[4]] + diff(grconvertY(c(0, 0.45), from = "lines", to = "user")),
           col = bracket_color,
           xpd = NA
           )
  mtext(c("Essential\ngenes", "Non-essential\ngenes"),
        at = c(mean(x_positions[1:2]), mean(x_positions[3:4])),
        line = genes_label_line, padj = 0, cex = par("cex")
        )

  mtext(text = rep(paste0(if (write_rep) "Rep" else "R", 1:2), 2),
        at = x_positions, side = 1, line = rep_label_line, cex = par("cex")
        )

  title(as.expression(use_title), cex.main = 1, line = title_line)
  par(old_mar)
  return(invisible(rep_list))
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
  show_list <- split(essent_df[, "Log2FC"], !(essent_df[, "Is_essential"]))
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




# Functions for analyzing bidirectional promoters with violin plots -------

FilterRepeatedEntrezs <- function(bidirect_df, max_distance = 10000L) {
  bidirect_df <- bidirect_df[bidirect_df[, "Distance"] <= max_distance, ]
  all_entrezs <- c(bidirect_df[, "Entrez_ID_1"], bidirect_df[, "Entrez_ID_2"])
  num_occurrences <- table(all_entrezs)[as.character(unique(all_entrezs))]
  repeated_entrezs <- as.integer(names(num_occurrences)[num_occurrences >= 2])
  are_valid <- (!(bidirect_df[, "Entrez_ID_1"] %in% repeated_entrezs)) &
               (!(bidirect_df[, "Entrez_ID_2"] %in% repeated_entrezs))
  bidirect_df <- bidirect_df[are_valid, ]
  return(bidirect_df)
}


BidirectionalSets <- function(bidirect_df, max_distance = 10000L) {
  bidirect_df <- FilterRepeatedEntrezs(bidirect_df, max_distance)
  are_discordant <- bidirect_df[, "Combination"] %in% "Discordant"
  are_concordant <- bidirect_df[, "Combination"] %in% "Both non-essential"
  essential_entrezs <- c(
    bidirect_df[, "Entrez_ID_1"][are_discordant & (bidirect_df[, "Essentiality_1"] %in% "Essential")],
    bidirect_df[, "Entrez_ID_2"][are_discordant & (bidirect_df[, "Essentiality_2"] %in% "Essential")]
  )
  discordant_entrezs <- c(
    bidirect_df[, "Entrez_ID_1"][are_discordant & (bidirect_df[, "Essentiality_1"] %in% "Non-essential")],
    bidirect_df[, "Entrez_ID_2"][are_discordant & (bidirect_df[, "Essentiality_2"] %in% "Non-essential")]
  )
  concordant_entrezs <- bidirect_df[, "Entrez_ID_1"][are_concordant]
  results_list <- list(
    "essential next to non-essential"     = essential_entrezs,
    "non-essential next to essential"     = discordant_entrezs,
    "non-essential next to non-essential" = concordant_entrezs
  )
  return(results_list)
}



IntegrateData <- function(entrezs_vec,
                          bidirect_df,
                          logfc_df,
                          max_distance = 10000L,
                          only_complete_pairs = TRUE,
                          choose_rep = NULL
                          ) {
  bidirect_df <- FilterRepeatedEntrezs(bidirect_df, max_distance)
  matches_vec <- match(entrezs_vec, bidirect_df[, "Entrez_ID_1"])
  matches_vec <- ifelse(is.na(matches_vec),
                        match(entrezs_vec, bidirect_df[, "Entrez_ID_2"]),
                        matches_vec
                        )
  results_df <- data.frame(
    "Entrez_ID" = entrezs_vec,
    bidirect_df[matches_vec, ],
    stringsAsFactors = FALSE
  )
  matches_vec <- match(results_df[, "Entrez_ID"], logfc_df[, "Entrez_ID"])
  if (is.null(choose_rep)) {
    logfc_column <- "Mean_log2FC"
  } else {
    logfc_column <- paste0("Log2FC_rep", choose_rep)
  }
  logfc_vec <- logfc_df[, logfc_column][matches_vec]
  if (only_complete_pairs) {
    logfc_1_vec <- logfc_df[, logfc_column][match(results_df[, "Entrez_ID_1"], logfc_df[, "Entrez_ID"])]
    logfc_2_vec <- logfc_df[, logfc_column][match(results_df[, "Entrez_ID_2"], logfc_df[, "Entrez_ID"])]
    are_valid <- (!(is.na(logfc_1_vec))) & (!(is.na(logfc_2_vec)))
  } else {
    are_valid <- !(is.na(logfc_vec))
  }

  results_df <- data.frame(
    results_df[are_valid, ],
    "Log2FC" = logfc_vec[are_valid],
    row.names = NULL
  )
  new_order <- order(results_df[, "Entrez_ID_1"])
  results_df <- results_df[new_order, ]
  row.names(results_df) <- NULL
  return(results_df)
}



SplitByDistance <- function(pairs_df, distance_cutoff = 1000L) {
  split(pairs_df[, "Log2FC"], pairs_df[, "Distance"] > distance_cutoff)
}



PrettyScientific <- function(input_number, digits = 0, scientific = 4) {
  scientific_split <- strsplit(formatC(input_number, format = "e", digits = digits),
                               "e", fixed = TRUE
                               )[[1]]
  power_of_10 <- as.integer(scientific_split[[2]])
  if (abs(power_of_10) <= scientific) {
    result_text <- as.expression(formatC(input_number, digits = max(digits, 1) , format = "fg"))
  } else {
    result_text <- as.expression(bquote(.(scientific_split[[1]]) %*% 10^.(power_of_10)))
  }
  return(result_text)
}



IndicatePValue <- function(x_1,
                           x_2,
                           vec_1,
                           vec_2,
                           paired      = FALSE,
                           y_line_adj  = 0,
                           y_start_adj = 0,
                           p_value_cex = 0.3
                           ) {

  ## Compute statistics
  p_val <- t.test(vec_1, vec_2, paired = paired)[["p.value"]]
  p_value_label <- PrettyScientific(p_val)
  p_value_label <- ConcatenateExpressions(list(expression(italic("p")), p_value_label), my_sep = " = ")

  ## Determine x and y positions
  y_pos <- par("usr")[[4]] + diff(grconvertY(c(0, 0.5 + y_start_adj), from = "lines", to = "user"))
  y_pos_line <- y_pos + diff(grconvertY(c(0, y_line_adj), from = "lines", to = "user"))
  y_text_pos <- y_pos + diff(grconvertY(c(0, 0.5), from = "lines", to = "user"))
  y_bottom_line <-  y_pos_line - diff(grconvertY(c(0, 0.25), from = "lines", to = "user"))

  ## Draw a white background rectangle over the grid
  rect(xleft   = x_1 - diff(grconvertX(c(0, 1.5), from = "lines", to = "user")),
       xright  = x_2 + diff(grconvertX(c(0, 1.5), from = "lines", to = "user")),
       ybottom = y_bottom_line,
       ytop    = y_text_pos + diff(grconvertY(c(0, 0.5), from = "lines", to = "user")),
       col     = "white",
       border  = NA,
       xpd     = NA
       )

  ## Draw brackets to indicate between-group comparison
  segments(x0  = x_1,
           x1  = x_2,
           y0  = y_pos_line,
           col = "gray60",
           xpd = NA
           )
  segments(x0  = c(x_1, x_2),
           y0  = y_pos_line,
           y1  = y_bottom_line,
           col = "gray60",
           xpd = NA
           )

  ## Show p values
  text(x      = mean(c(x_1, x_2)),
       y      = y_text_pos,
       labels = VerticalAdjust(p_value_label),
       cex    = p_value_cex,
       xpd    = NA
       )
  return(invisible(NULL))
}


BidirectionalViolins <- function(bidirect_df,
                                 logfc_df,
                                 max_distance         = 10000L,
                                 distance_cutoff      = 1000L,
                                 show_phenotype_score = TRUE,
                                 num_cell_divisions   = 10L,
                                 y_limits             = if (show_phenotype_score) c(-0.55, 0.2) else c(-5.5, 0.2),
                                 lower_bound          = if (show_phenotype_score) -0.5 else -5,
                                 num_controls         = NULL,
                                 choose_rep           = NULL,
                                 compare_across       = TRUE,
                                 point_cex            = 0.5,
                                 annotation_cex       = 0.7,
                                 draw_groups_n        = TRUE,
                                 use_spacing          = 0.9,
                                 label_points         = FALSE,
                                 draw_group_labels    = TRUE,
                                 gap_ratio            = 1.35,
                                 y_start_adj          = 0,
                                 y_line_adj           = 0,
                                 p_value_cex          = 0.7,
                                 show_title           = NULL,
                                 ...
                                 ) {

  ## Assemble data
  entrezs_vec_list <- BidirectionalSets(bidirect_df, max_distance)
  pairs_df_list <- lapply(entrezs_vec_list,
                          function(x) IntegrateData(x,
                                                    bidirect_df,
                                                    logfc_df,
                                                    max_distance = max_distance,
                                                    choose_rep = choose_rep
                                                    )
                          )


  ## Split data into near and far subgroups
  set.seed(1)
  if (is.null(num_controls)) {
    random_indices <- sample(seq_len(nrow(pairs_df_list[[3]])), nrow(pairs_df_list[[1]]))
  } else {
    are_near <- pairs_df_list[[3]][, "Distance"] < distance_cutoff
    near_indices <- which(are_near)
    far_indices <- which(!(are_near))
    if (sum(!(are_near)) >= num_controls) {
      far_indices <- sample(far_indices, num_controls)
    }
    if (sum(are_near) >= num_controls) {
      near_indices <- sample(near_indices, num_controls)
    }
    random_indices <- c(near_indices, far_indices)
  }
  pairs_df_list[[3]] <- pairs_df_list[[3]][random_indices, ]
  pairs_df_list <- lapply(pairs_df_list, function(x) {
    x[order(x[, "Distance"] >= distance_cutoff), ]
  })
  numeric_list <- unlist(lapply(pairs_df_list, SplitByDistance), recursive = FALSE)
  if (show_phenotype_score) {
    numeric_list <- lapply(numeric_list, function(x) x / num_cell_divisions)
  }
  groups_vec <- rep(1:3, each = 2)

  if (label_points) {
    text_vec_list <- lapply(1:3, function(x) {
      vec_1 <- pairs_df_list[[x]][, "Gene_symbol_1"]
      vec_2 <- pairs_df_list[[x]][, "Gene_symbol_2"]
      are_essential <- pairs_df_list[[x]][, "Essentiality_1"] == "Essential"
      if (x == 1) {
        first_vec <- ifelse(are_essential, vec_1, vec_2)
        second_vec <- ifelse(are_essential, vec_2, vec_1)
      } else if (x == 2) {
        first_vec <- ifelse(are_essential, vec_2, vec_1)
        second_vec <- ifelse(are_essential, vec_1, vec_2)
      } else {
        first_vec <- vec_1
        second_vec <- vec_2
      }
      paste0(first_vec, "\n", second_vec)
    })
    text_vec <- unlist(text_vec_list)
  } else {
    text_vec <- NULL
  }

  ## Draw violin plots
  old_mar <- par("mar" = c(6, 4, 3, 1.5))
  x_positions <- BeeViolinPlot(numeric_list,
                               groups_vec,
                               y_limits       = y_limits,
                               lower_bound    = lower_bound,
                               gap_ratio      = gap_ratio,
                               point_cex      = point_cex,
                               use_spacing    = use_spacing,
                               violin_colors  = rep(c(brewer.pal(9, "Purples")[[3]], "#cbdde7", "#c7e7c0"), each = 2),
                               point_colors   = rep(c("#8f83af", "#7690ad", "#689c79"), each = 2),
                               line_colors    = rep(c("#7c7198", "#617b98", "#5b8669"), each = 2),
                               border_colors  = rep(c("#d1cddb", "#c5d3dd", "#bfd4c6"), each = 2),
                               draw_border    = TRUE,
                               draw_groups_n  = FALSE,
                               text_vec       = text_vec,
                               ...
                               )

  ## Draw y axis label
  if (is.null(choose_rep)) {
    ylab_prefix <- "Mean"
  } else {
    ylab_prefix <- paste0("Replicate ", choose_rep)
  }
  if (show_phenotype_score) {
    y_axis_label <- bquote(.(ylab_prefix) ~ "phenotype (" * gamma * ")")
  } else {
    y_axis_label <- bquote(.(ylab_prefix) ~ "log"[2] ~ "fold change")
  }
  mtext(y_axis_label, side = 2, line = 2.5, cex = par("cex"))

  ## Indicate the number of observations
  if (draw_groups_n) {
    mtext(lengths(numeric_list), at = x_positions, side = 1, line = -0.97,
          cex = par("cex") * 0.4, col = "gray60", padj = 0
          )
  }

  ## Draw x axis labels
  if (draw_group_labels) {
    label_positions <- tapply(x_positions, groups_vec, mean)
    more_space <- annotation_cex > 0.7
    distance_labels <- c(
      as.expression(bquote("" <= .(distance_cutoff / 1000) * scriptscriptstyle(" ") * "kb")),
      as.expression(bquote("1\u2013" * .(max_distance / 1000) * scriptscriptstyle(" ") * "kb"))
    )
    distance_labels <- sapply(distance_labels, VerticalAdjust)
    mtext(rep(distance_labels, 3),
          side = 1, at = x_positions, line = 0.75, cex = par("cex") * annotation_cex
          )
    segments(x0  = label_positions - 0.7,
             x1  = label_positions + 0.7,
             y0  = par("usr")[[3]] - diff(grconvertY(c(0, if (more_space) 1.65 else 1.3), from = "lines", to = "user")),
             col = "gray60",
             xpd = NA
             )
    old_lheight <- par("lheight" = 1.15)
    mtext(c("Essential genes\nwith a non-essential\ngene nearby",
            "Non-essential genes\nwith an essential\ngene nearby",
            "Non-essential genes\nwith a non-essential\ngene nearby"
            ),
          side = 1, at = label_positions, line = if (more_space) 1.45 else 0.9, padj = 1,
          cex = par("cex") * annotation_cex
          )
    par(old_lheight)
  }

  ## Draw brackets to indicate between-group comparisons
  numeric_list <- lapply(numeric_list, function(x) BringWithinLimits(x, lower_bound = lower_bound)[["curtailed_vec"]])
  x_space <- 0.05
  if (compare_across) {
    IndicatePValue(x_positions[[1]],
                   x_positions[[3]] - x_space,
                   numeric_list[[1]],
                   numeric_list[[3]],
                   paired = TRUE,
                   y_line_adj  = y_line_adj,
                   y_start_adj = y_start_adj,
                   p_value_cex = p_value_cex
                   )
    IndicatePValue(x_positions[[3]] + x_space,
                   x_positions[[5]],
                   numeric_list[[3]],
                   numeric_list[[5]],
                   y_line_adj  = y_line_adj,
                   y_start_adj = y_start_adj,
                   p_value_cex = p_value_cex
                   )
  } else {
    IndicatePValue(x_positions[[3]],
                   x_positions[[4]],
                   numeric_list[[3]],
                   numeric_list[[4]],
                   y_line_adj  = y_line_adj - 0.1,
                   y_start_adj = y_start_adj,
                   p_value_cex = p_value_cex
                   )
  }

  mtext(VerticalAdjust(show_title),
        at = grconvertX(0.015, from = "npc", to = "user"),
        adj = 0, line = -0.05, cex = par("cex")
        )
  par(old_mar)
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
  ROC_df[, "Log2FC"] <- log2fc_df[, "Log2FC"]

  are_nan <- is.nan(ROC_df[, "Log2FC"])
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
  new_order <- order(ROC_df[, "Log2FC"])
  ROC_df <- ROC_df[new_order, ]
  row.names(ROC_df) <- NULL

  ROC_df <- MakeROCDf(ROC_df, "Log2FC")
  return(ROC_df)
}



PlotEssentialROCDf <- function(use_ROC_df, use_title = "Gene essentiality with CRISPRoff") {

  numbers_text <- paste0(sum(use_ROC_df[, "Is_essential"]), " essential genes and ",
                         sum(!(use_ROC_df[, "Is_essential"])), " non-essential genes"
                         )

  old_mar <- par(mar = c(3.7, 3.95, 4.1, 1.75))

  PlotROCDf(use_ROC_df, flip = TRUE)

  title(as.expression(use_title), cex.main = par("cex"), line = 2.5)
  mtext(numbers_text, line = 0.75, cex = par("cex") * 0.75, adj = 0.5)

  par(old_mar)
  return(invisible(NULL))
}



MakeROCDfListList <- function(choose_rep = NULL, use_blomen_hart = TRUE) {
  essential_list <- GetEssentialGenes(use_blomen_hart)
  args_list <- list(essential_genes       = essential_list[["essential_entrezs"]],
                    non_essential_genes   = essential_list[["non_essential_entrezs"]],
                    min_count_at_baseline = 0L,
                    choose_rep            = choose_rep
                    )
  ROC_df_list_list <- lapply(c(FALSE, TRUE), function(allow_switch) {
    list(
      ROC_T0vT12 = do.call(GetEssentialROCDf, c(args_list, list(
        baseline_indices      = 3:4,
        intervention_indices  = 5:6,
        allow_switch          = allow_switch
      ))),
      ROC_BvT12  = do.call(GetEssentialROCDf, c(args_list, list(
        baseline_indices      = 1:2,
        intervention_indices  = 5:6,
        allow_switch          = allow_switch
      ))),
      ROC_BvT0 = do.call(GetEssentialROCDf, c(args_list, list(
        baseline_indices      = 1:2,
        intervention_indices  = 3:4,
        allow_switch          = allow_switch
      )))
    )
  })
  names(ROC_df_list_list) <- c("Template switch excluded", "Template switch allowed")
  return(ROC_df_list_list)
}



ROCInputDf <- function(input_df, essential_entrezs, non_essential_entrezs) {
  are_essential <- ifelse(input_df[, "Entrez_ID"] %in% essential_entrezs,
                          TRUE,
                          ifelse(input_df[, "Entrez_ID"] %in% non_essential_entrezs,
                                 FALSE,
                                 NA
                                 )
                          )
  input_df[, "Is_essential"] <- are_essential
  results_df <- input_df[!(is.na(are_essential)), ]
  new_order <- order(results_df[, "Is_essential"], decreasing = TRUE)
  results_df <- results_df[new_order, ]
  row.names(results_df) <- NULL
  return(results_df)
}



ROCDfForColumn <- function(input_df, use_column) {
  new_order <- order(input_df[, use_column])
  results_df <- input_df[new_order, ]
  results_df <- MakeROCDf(results_df, numeric_column = use_column)
  return(results_df)
}



# Helper functions (required by ReplicateScatterPlot) ---------------------

ConcatenateExpressions <- function(expression_list, my_sep = "  \u2013  ") {
  literal_strings <- vapply(expression_list, StripExpression, "")
  combined_string <- paste0(literal_strings, collapse = paste0(" * \"", my_sep, "\" * "))
  results_expression <- parse(text = combined_string)
  return(results_expression)
}


VerticalAdjust <- function(use_expression) {
  my_list <- list(expression(""^phantom("gh")), use_expression, expression(""[phantom("gh")]))
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
                           draw_lines           = FALSE,
                           use_lwd              = 1.25,
                           line_x_distance      = -0.2,
                           length_in_lines      = 0.5,
                           border_colors        = NULL,
                           use_pch              = 16,
                           use_point_size       = 1.2,
                           lines_x_start        = 0.75,
                           y_mid                = 0.5,
                           small_gap_size       = 1.25,
                           large_gap_multiplier = 1.75,
                           point_x_start        = 0.15,
                           title_vec            = NULL,
                           title_x_start        = NULL
                           ) {

  ## Perform checks
  stopifnot(identical(length(labels_list), length(use_colors)))

  ## Prepare for drawing legend
  small_gap <- diff(grconvertY(c(0, small_gap_size), from = "char", to = "npc"))
  medium_gap <- small_gap * 1.25
  large_gap <- small_gap * large_gap_multiplier

  if (is.null(title_vec)) {
    text_list <- labels_list
  } else {
    text_list <- c(list(title_vec), labels_list)
  }
  if (all(lengths(text_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(text_list))
    are_first <- rep(TRUE, length(text_list))
  } else {
    are_first <- unlist(lapply(text_list, function(x) {
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

  if (!(is.null(title_vec))) {
    if (is.null(title_x_start)) {
      title_x <- grconvertX(x = x_point, from = "npc", to = "user") - (strwidth(VerticalAdjust("")) / 2) - (strwidth("o") / 2)
    } else {
      x_title <- 1 + diff(grconvertX(c(0, title_x_start), from = "lines", to = "npc"))
      title_x <-  grconvertX(x = x_title, from = "npc", to = "user")
    }
    text(x      = title_x,
         y      = y_pos[seq_along(title_vec)],
         cex    = 1,
         labels = sapply(title_vec, VerticalAdjust),
         adj    = c(0, 0.5),
         xpd    = NA
         )
    y_pos <- y_pos[seq_along(unlist(labels_list)) + length(title_vec)]
  }

  ## Draw legend
  text(x      = grconvertX(x = x_text, from = "npc", to = "user"),
       y      = y_pos,
       labels = sapply(unlist(labels_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )
  groups_vec <- rep(seq_along(labels_list), lengths(labels_list))
  if (draw_lines) {
    x_user <- grconvertX(x = x_text, from = "npc", to = "user")
    line_x_start <- x_user + diff(grconvertX(c(0, line_x_distance), from = "lines", to = "user"))
    segments(x0  = line_x_start,
             x1  = line_x_start + diff(grconvertX(c(0, length_in_lines), from = "lines", to = "user")),
             y0  = tapply(y_pos, groups_vec, mean),
             col = use_colors,
             lwd = par("lwd") * use_lwd,
             xpd = NA
             )
  } else {
    points(x   = rep(grconvertX(x = x_point, from = "npc", to = "user"), length(labels_list)),
           y   = tapply(y_pos, groups_vec, mean),
           cex = use_point_size,
           pch = use_pch,
           col = if (!(is.null(border_colors))) border_colors else use_colors,
           bg  = use_colors,
           xpd = NA
           )
  }

  return(invisible(NULL))
}



CurtailedAxisLabels <- function(tick_positions, upper_bound, lower_bound,
                                lower_bound_enforced, upper_bound_enforced
                                ) {
  plain_axis_labels <- format(tick_positions, trim = TRUE)
  axis_labels <- lapply(plain_axis_labels, function(x) as.expression(bquote(""[.(x)])))
  axis_labels <- sapply(axis_labels, function(x) x)
  if (lower_bound_enforced && (abs(lower_bound - tick_positions[[1]]) < 1e-15)) {
    axis_labels[1] <- as.expression(bquote(""["" <= scriptscriptstyle(.(if (tick_positions[[1]] < 0) "" else " "))
                                              * .(plain_axis_labels[[1]])]
    ))
  }
  if (upper_bound_enforced) {
    num_ticks <- length(tick_positions)
    if (abs(upper_bound - tick_positions[[num_ticks]]) < 1e-15) {
      axis_labels[num_ticks] <- as.expression(bquote(""["" >= scriptscriptstyle(" ") *
                                                          .(plain_axis_labels[[num_ticks]])]
      ))
    }
  }
  return(axis_labels)
}



# Functions for creating replicate scatter plots --------------------------

Palify <- function(colors_vec, fraction_pale = 0.5) {
  adjustcolor(colors_vec,
              offset    = c(rep(fraction_pale, 3), 0),
              transform = diag(c(rep(1 - fraction_pale, 3), 1))
              )

}


ReplicateScatterPlot <- function(input_df,
                                 show_phenotype_score = FALSE,
                                 lower_bound          = -8 * (if (show_phenotype_score) 0.1 else 1),
                                 upper_bound          = 8  * (if (show_phenotype_score) 0.1 else 1),
                                 axis_limits          = c(lower_bound, upper_bound),
                                 highlight_essential  = TRUE,
                                 use_blomen_hart      = TRUE,
                                 highlight_NT         = TRUE,
                                 use_title            = NULL,
                                 title_font           = 2,
                                 embed_PNG            = FALSE,
                                 use_mar              = NULL,
                                 axis_labels_list     = NULL,
                                 essential_colors     = c("#6810c6", "#18a008"),
                                 point_cex            = 0.4,
                                 x_axis_label_line    = 2,
                                 y_axis_label_line    = 2.3,
                                 sparse_x_axis_labels = FALSE,
                                 use_tcl              = 0.375,
                                 x_axis_mgp           = 0.4,
                                 y_axis_mgp           = 0.55,
                                 legend_point_x_start = 0,
                                 legend_lines_x_start = 0.8,
                                 abbreviate_NT        = FALSE,
                                 show_y_labels        = TRUE
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
  if (is.null(axis_labels_list)) {
    if (show_phenotype_score) {
      axis_labels_list <- list(
        expression("Replicate 1 phenotype (" * gamma * ")"),
        expression("Replicate 2 phenotype (" * gamma * ")")
      )
    } else {
      axis_labels_list <- list(
        expression("Replicate 1 log"[2] * " fold change"),
        expression("Replicate 2 log"[2] * " fold change")
      )
    }
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
        highlight_colors <- essential_colors
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
      label_colors <- highlight_colors
    } else {
      labels_list <- list("all" = c("All", "genes", paste0("(", length(entrezs_vec), ")")))
      label_colors <- Palify("black", fraction_pale = 0.3)
      highlight_colors <- c()
    }
    if (highlight_NT) {
      if (highlight_essential) {
        highlight_mat <- cbind(highlight_mat, "is_NT" = are_NT)
      } else {
        highlight_mat <- cbind("is_NT" = are_NT)
      }
      NT_color <- "#297eff"
      highlight_colors <- c(highlight_colors, NT_color)
      label_colors <- c(label_colors, NT_color)
      labels_list <- c(labels_list,
                       list("NT" = c(if (abbreviate_NT) "NT" else "Non-targeting",
                                     "controls", AddSum(are_NT)
                                     )
                            )
                       )
    }
    stopifnot(!(any(rowSums(highlight_mat) > 1)))
    are_highlighted <- rowSums(highlight_mat) == 1
    set.seed(1)
    scrambled_indices <- sample(which(are_highlighted))
    are_highlighted_mat <- highlight_mat[scrambled_indices, , drop = FALSE]
    highlighted_xy_mat <- xy_mat[scrambled_indices, ]
    point_colors <- adjustcolor(highlight_colors, alpha.f = if (highlight_essential) 0.7 else 0.85)
    colors_vec <- vapply(seq_len(nrow(are_highlighted_mat)),
                         function(x) point_colors[[which(are_highlighted_mat[x, ])]],
                         ""
                         )
  } else {
    are_highlighted <- rep(FALSE, nrow(xy_mat))
  }

  ## Prepare plot margins
  if (is.null(use_mar)) {
    if (highlight_genes) {
      use_mar <- c(3.75, 3.75, 3.75, 7.5)
    } else {
      use_mar <- c(3.75, 3.75, 3.75, 2.1)
    }
  }
  old_mar <- par(mar = use_mar)

  if (embed_PNG) {
    current_device <- StartEmbedPNG(figures_dir)
  }

  ## Set up plot canvas
  plot(NA, xlim = xy_lim, ylim = xy_lim, xaxs = "i", yaxs = "i",
       axes = FALSE, ann = FALSE
       )
  abline(v = 0, h = 0,  col = "gray85", lend = "butt")
  abline(a = 0, b = 1,  col = "gray85", lend = "butt")
  abline(a = 0, b = -1, col = "gray85", lend = "butt")

  ## Draw points
  points(xy_mat[!(are_highlighted), ],
         pch = 16,
         cex = point_cex,
         col = adjustcolor("black", alpha.f = 0.35),
         xpd = NA
         )
  if (highlight_genes) {
    points(highlighted_xy_mat,
           pch = 16,
           cex = point_cex,
           col = colors_vec,
           xpd = NA
           )
  }

  if (embed_PNG) {
    StopEmbedPNG(current_device, figures_dir)
  }

  ## Annotate plot
  mtext(axis_labels_list[[1]], side = 1, line = x_axis_label_line, cex = par("cex"))
  if (show_y_labels) {
    mtext(axis_labels_list[[2]], side = 2, line = y_axis_label_line, cex = par("cex"))
  }
  if (!(is.null(use_title))) {
    title(use_title, cex.main = 1, font.main = title_font)
  }

  ## Draw axes
  for (i in 1:2) {
    tick_labels <- tick_labels_list[[i]]
    if ((i == 1) && sparse_x_axis_labels) {
      tick_labels <- ifelse(rep(c(TRUE, FALSE), length.out = length(tick_labels)),
                            tick_labels,
                            ""
                            )
    }
    axis(i,
         at       = tick_locations,
         labels   = if ((i == 2) && (!(show_y_labels))) NA else tick_labels,
         mgp      = c(3, if (i == 1) x_axis_mgp else y_axis_mgp, 0),
         tcl      = -(use_tcl),
         las      = 1,
         lwd      = par("lwd"),
         cex.axis = 1 / 0.7
         )
  }
  box()
  if (highlight_genes) {
    DrawSideLegend(labels_list,
                   use_colors = vapply(label_colors, Palify, fraction_pale = 0.15, ""),
                   use_point_size = 1, point_x_start = legend_point_x_start,
                   lines_x_start = legend_lines_x_start
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




# Functions for analyzing the proximity to essential genes ----------------

DistanceToEssential <- function(logfc_df,
                                use_TSS_df,
                                use_data_for_essential = FALSE,
                                essential_logfc_cutoff = -2
                                ) {

  ## Check assumptions
  stopifnot(!(anyNA(logfc_df[, "Gene_symbol"])))
  stopifnot(!(any(duplicated(logfc_df[, "Gene_symbol"]))))
  stopifnot(!(anyNA(use_TSS_df[, "Gene_symbol"])))
  stopifnot(!(any(duplicated(use_TSS_df[, "Gene_symbol"]))))
  stopifnot(!(anyNA(use_TSS_df[, "Chromosome"])))

  ## Prepare the TSS positions of essential genes
  if (use_data_for_essential) {
    matches_vec <- match(use_TSS_df[, "Entrez_ID"], logfc_df[, "Entrez_ID"])
    are_essential <- (logfc_df[, "Mean_log2FC"][matches_vec] <= essential_logfc_cutoff) %in% TRUE
  } else {
    are_essential <- use_TSS_df[, "Essentiality"] %in% "Essential"
  }
  positions_vec <- use_TSS_df[, "TSS"]
  names(positions_vec) <- use_TSS_df[, "Gene_symbol"]
  essential_positions_list <- split(positions_vec[are_essential],
                                    use_TSS_df[, "Chromosome"][are_essential]
                                    )

  ## Find the nearest essential genes
  df_list <- lapply(names(essential_positions_list), function(chromosome) {
    essential_positions_vec <- essential_positions_list[[chromosome]]
    are_this_chromosome <- use_TSS_df[, "Chromosome"] == chromosome
    nearest_list <- lapply(which(are_this_chromosome), function(x) {
      this_symbol <- use_TSS_df[, "Gene_symbol"][[x]]
      use_positions_vec <- essential_positions_vec[names(essential_positions_vec) != this_symbol]
      distances_vec <- use_TSS_df[, "TSS"][[x]] - use_positions_vec
      nearest_index <- which.min(abs(distances_vec))
      results_list <- list(
        "Gene_symbol"    = this_symbol,
        "Nearest_symbol" = names(distances_vec)[[nearest_index]],
        "Distance"       = distances_vec[[nearest_index]]
      )
      return(results_list)
    })
    results_df <- do.call(rbind.data.frame,
                          c(nearest_list, list(stringsAsFactors = FALSE, make.row.names = FALSE))
                          )
    return(results_df)
  })
  results_df <- do.call(rbind.data.frame,
                        c(df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE))
                        )

  ## Integrate with TSS data
  symbol_matches <- match(results_df[, "Gene_symbol"], use_TSS_df[, "Gene_symbol"])
  nearest_matches <- match(results_df[, "Nearest_symbol"], use_TSS_df[, "Gene_symbol"])
  results_df <- data.frame(
    use_TSS_df[symbol_matches, "Chromosome", drop = FALSE],
    results_df["Gene_symbol"],
    "Entrez_ID" = use_TSS_df[, "Entrez_ID"][symbol_matches],
    "DepMap_essentiality" = use_TSS_df[, "Essentiality"][symbol_matches],
    results_df["Nearest_symbol"],
    "Nearest_entrez" = use_TSS_df[, "Entrez_ID"][nearest_matches],
    results_df["Distance"],
    stringsAsFactors = FALSE
  )

  ## Integrate with log2 fold change data
  matches_vec <- match(logfc_df[, "Entrez_ID"], results_df[, "Entrez_ID"])
  are_numeric <- vapply(logfc_df, is.numeric, logical(1)) & (!(vapply(logfc_df, is.integer, logical(1))))
  results_df <- data.frame(
    logfc_df[, !(are_numeric)],
    logfc_df["Mean_log2FC"],
    "Nearest_mean_log2FC" = NA,
    results_df[matches_vec, !(names(results_df) %in% c("Entrez_ID", "Gene_symbol"))],
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = NULL
  )
  matches_vec <- match(results_df[, "Nearest_entrez"], results_df[, "Entrez_ID"])
  results_df[, "Nearest_mean_log2FC"] <- results_df[, "Mean_log2FC"][matches_vec]
  return(results_df)
}



RollingQuantiles <- function(x_vec,
                             y_vec,
                             use_window = 20L,
                             all_quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95)
                             ) {
  use_order <- order(x_vec)
  x_vec <- x_vec[use_order]
  y_vec <- y_vec[use_order]
  fill_vec <- rep(NA, use_window - 1L)
  padded_y_vec <- c(fill_vec, y_vec, fill_vec)
  start_vec <- seq(from = 1, to = length(padded_y_vec) - use_window + 1L)
  stop_vec <- seq(from = use_window, to = length(padded_y_vec))
  indices_list <- mapply(function(x, y) seq(x, y), start_vec, stop_vec, SIMPLIFY = FALSE)
  data_list <- lapply(indices_list, function(x) padded_y_vec[x])
  quantile_mat <- t(vapply(data_list, quantile, probs = all_quantiles, na.rm = TRUE, numeric(5)))
  use_indices <- seq(from = use_window / 2, to = (use_window / 2) + length(y_vec) - 1)
  quantile_mat <- quantile_mat[use_indices, ]
  results_mat <- cbind(
    "x" = x_vec,
    "y" = y_vec,
    quantile_mat
  )
  return(results_mat)
}



DistanceScatter <- function(distance_vec,
                            logfc_vec       = NULL,
                            lower_bound     = -0.6,
                            upper_bound     = 0.2,
                            num_ticks       = 5L,
                            x_lower_bound   = 10L,
                            y_axis_label    = NULL,
                            x_limit         = NULL,
                            y_limits        = c(lower_bound, upper_bound),
                            loess_span      = 0.5,
                            use_hue         = "#275dd3", #2777d3
                            point_cex       = 0.4,
                            point_color     = "black",
                            show_x_labels   = TRUE,
                            x_grid          = FALSE,
                            grid_lwd        = 1,
                            grid_color      = "gray85",
                            grid_highlight  = grid_color,
                            median_lwd      = 1.5,
                            median_color    = NULL,
                            embed_PNG       = FALSE,
                            png_res         = 900,
                            png_padding     = 0.1
                            ) {

  ## Prepare data
  if (is.null(logfc_vec) && ("data.frame" %in% class(distance_vec))) {
    logfc_vec <- distance_vec[, "Mean_log2FC"]
    distance_vec <- distance_vec[, "Distance"]
    if (is.null(y_axis_label)) {
      y_axis_label <- expression("Mean phenotype (" * gamma * ")")
    }
  } else {
    y_axis_label <- expression("Phenotype (" * gamma * ")")
  }
  x_vec <- abs(distance_vec)
  y_vec <- logfc_vec / 10

  are_NA <- is.na(x_vec) | (is.na(y_vec))
  x_vec <- x_vec[!(are_NA)]
  y_vec <- y_vec[!(are_NA)]
  use_order <- order(x_vec)
  x_vec <- x_vec[use_order]
  y_vec <- y_vec[use_order]

  ## Truncate very low distances (if there is only a single outlier)
  are_outliers <- x_vec < x_lower_bound
  truncated_x_vec <- x_vec
  if (length(unique(x_vec[are_outliers])) == 1) {
    curtailed_x <- TRUE
    truncated_x_vec[are_outliers] <- x_lower_bound
  } else {
    curtailed_x <- FALSE
  }
  x_vec <- log10(x_vec)
  truncated_x_vec <- log10(truncated_x_vec)

  ## Compute smoothened quantile estimates
  ## Inspired by: https://www.r-statistics.com/2010/04/quantile-loess-combining-a-moving-quantile-window-with-loess-r-function/
  finite_y_vec <- ifelse(y_vec == -Inf, -0.6, y_vec)
  quant_mat <- RollingQuantiles(x_vec, finite_y_vec, use_window = 20)
  loess_mat <- sapply(setdiff(colnames(quant_mat), c("x", "y")), function(x) {
    loess(quant_mat[, x] ~ quant_mat[, "x"], family = "symmetric", span = loess_span)[["fitted"]]
  })

  ## Ensure that the quantile estimates do not intersect
  loess_mat[, "25%"] <- ifelse(loess_mat[, "50%"] < loess_mat[, "25%"], loess_mat[, "50%"], loess_mat[, "25%"])
  loess_mat[, "5%"]  <- ifelse(loess_mat[, "25%"] < loess_mat[, "5%"],  loess_mat[, "25%"], loess_mat[, "5%"])
  loess_mat[, "75%"] <- ifelse(loess_mat[, "50%"] > loess_mat[, "75%"],  loess_mat[, "50%"], loess_mat[, "75%"])
  loess_mat[, "95%"] <- ifelse(loess_mat[, "75%"] > loess_mat[, "95%"],  loess_mat[, "75%"], loess_mat[, "95%"])

  ## Determine axis limits
  curtailed_y_list <- BringWithinLimits(y_vec, lower_bound = lower_bound, upper_bound = upper_bound)
  truncated_y_vec <- curtailed_y_list[["curtailed_vec"]]
  x_range <- range(truncated_x_vec)
  if (x_range[[1]] > 1) {
    x_range[[1]] <- 1
  }
  if ((!(is.null(x_limit))) && (x_range[[2]] < x_limit)) {
    x_range[[2]] <- x_limit
  }
  y_range <- range(truncated_y_vec)
  if (!(is.null(y_limits))) {
    if (y_range[[1]] > y_limits[[1]]) {
      y_range[[1]] <- y_limits[[1]]
    }
    if (y_range[[2]] < y_limits[[2]]) {
      y_range[[2]] <- y_limits[[2]]
    }
  }

  ## Set up the plot region and grid
  if (embed_PNG) {
    current_device <- StartEmbedPNG(figures_dir, png_res = png_res, add_padding = TRUE, padding_in_inches = png_padding)
  }
  MakeEmptyPlot(x_limits = x_range + (diff(x_range) * c(-0.015, if (is.null(x_limit)) 0.015 else 0)),
                y_limits = y_range + (diff(y_range) * c(if (x_grid) 0 else -0.025, 0))
                )
  x_ticks <- pretty(axTicks(1), n = 6)
  y_ticks <- pretty(axTicks(2), n = num_ticks)
  if (x_grid) {
    x_grid_bottom <- par("usr")[[3]]
    if (show_x_labels) {
      x_grid_bottom <- x_grid_bottom - diff(grconvertY(c(0, 0.275), from = "lines", to = "user"))
    }
    segments(x0  = x_ticks,
             y0  = x_grid_bottom,
             y1  = par("usr")[[4]],
             col = ifelse(x_ticks == log10(1000), grid_color, grid_color),
             lwd = grid_lwd * par("lwd"),
             xpd = NA
             )
  }
  segments(x0  = par("usr")[[1]],
           x1  = par("usr")[[2]],
           y0  = y_ticks,
           col = ifelse(y_ticks == 0, grid_highlight, grid_color),
           lwd = grid_lwd * par("lwd"),
           xpd = NA
           )

  ## Draw shading for the quantile estimates
  polygon(x      = c(quant_mat[, "x"], rev(quant_mat[, "x"])),
          y      = c(loess_mat[, "5%"], rev(loess_mat[, "95%"])),
          col    = adjustcolor(use_hue, alpha.f = 0.2),
          border = NA
          )
  polygon(x      = c(quant_mat[, "x"], rev(quant_mat[, "x"])),
          y      = c(loess_mat[, "25%"], rev(loess_mat[, "75%"])),
          col    =  adjustcolor(use_hue, alpha.f = 0.25),
          border = NA
          )

  ## Draw points and mean line
  points(x   = truncated_x_vec,
         y   = truncated_y_vec,
         pch = 16,
         cex = point_cex,
         col = adjustcolor(point_color, alpha.f = 0.4),
         xpd = NA
         )
  lines(x    = quant_mat[, "x"],
        y    = loess_mat[, "50%"],
        col  = if (is.null(median_color)) adjustcolor(use_hue, alpha.f = 0.85) else median_color,
        lwd  = median_lwd * par("lwd"),
        lend = "butt"
        )
  if (embed_PNG) {
    StopEmbedPNG(current_device, output_dir, add_padding = TRUE, padding_in_inches = png_padding)
  }

  ## Prepare axis labels
  old_scipen <- options(scipen = -1)
  x_axis_labels <- as.character(10^x_ticks)
  x_axis_labels <- sub("e+0", "0^", x_axis_labels, fixed = TRUE)
  if ((x_axis_labels[[1]] == "10") && (x_axis_labels[[length(x_axis_labels)]] == "10^7")) {
    x_axis_labels <- c("10 bp", "100 bp", "1 kb", "10 kb", "100 kb",
                       "1 Mb", "10 Mb"
                       )
    if (curtailed_x && (x_lower_bound == 10)) {
      x_axis_labels <- sapply(x_axis_labels, function(x) VerticalAdjust(x))
      x_axis_labels[1] <- VerticalAdjust(expression("" <= "10 bp"))
    }
  } else {
    x_axis_labels <- parse(text = x_axis_labels)
  }
  options(old_scipen)
  y_axis_labels <- CurtailedAxisLabels(
    y_ticks,
    lower_bound          = lower_bound,
    upper_bound          = upper_bound,
    lower_bound_enforced = curtailed_y_list[["lower_bound_enforced"]],
    upper_bound_enforced = curtailed_y_list[["upper_bound_enforced"]]
  )

  ## Draw axes
  if (show_x_labels) {
    axis(1,
         at     = x_ticks,
         labels = x_axis_labels,
         mgp    = c(3, if (is.expression(x_axis_labels)) 0.75 else 0.45, 0),
         tcl    = -0.35,
         lwd    = par("lwd"),
         tick   = !(x_grid)
         )
    mtext("Distance to nearest essential gene", side = 1, line = 1.45,
          cex = par("cex"), padj = 1
          )
  }
  axis(2,
       at       = y_ticks,
       labels   = y_axis_labels,
       mgp      = c(3, 0.5, 0),
       tcl      = -0.35,
       las      = 1,
       lwd      = par("lwd"),
       cex.axis = 1 / 0.7
       )
  mtext(y_axis_label, side = 2, line = 2.5, cex = par("cex"))
  if (!(x_grid)) {
    box(bty = "l")
  }
  return(invisible(NULL))
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
  reads_labels <- sapply(reads_labels, VerticalAdjust)
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



# Helper functions for plotting QC data -----------------------------------

PlotBarplotMat <- function(barplot_mat,
                           colors_vec,
                           positions_vec = seq_len(ncol(barplot_mat)),
                           bar_width = 2/3
                           ) {

  num_stacks <- nrow(barplot_mat)
  stopifnot(length(colors_vec) == num_stacks)

  barplot_mat <- barplot_mat[rev(seq_len(nrow(barplot_mat))), , drop = FALSE] # The first row should be on top

  lower_borders_vec_list <- lapply(seq_len(num_stacks), function(x) {
    if (x == 1) {
      rep(0, ncol(barplot_mat))
    } else {
      colSums(barplot_mat[seq_len(x - 1), , drop = FALSE])
    }
  })

  upper_borders_vec_list <- lapply(seq_len(num_stacks), function(x) {
    colSums(barplot_mat[seq_len(x), , drop = FALSE])
  })

  for (i in seq_len(ncol(barplot_mat))) {
    for (j in seq_len(num_stacks)) {
      rect(xleft   = positions_vec[[i]] - (bar_width / 2),
           xright  = positions_vec[[i]] + (bar_width / 2),
           ybottom = lower_borders_vec_list[[j]][[i]],
           ytop    = upper_borders_vec_list[[j]][[i]],
           col     = colors_vec[[j]],
           border  = NA,
           xpd     = NA
           )
    }
  }
  return(invisible(NULL))
}



DrawBottomLabels <- function(x_positions, groups_vec, are_included) {
  mtext(text = rep(paste0("R", 1:2), times = 3)[are_included],
        at = x_positions, side = 1, line = 0.3, cex = par("cex")
        )
  are_rep1 <- rep(c(TRUE, FALSE), times = 3)[are_included]
  segments(x0  = x_positions[are_rep1] - 0.25,
           x1  = x_positions[!(are_rep1)] + 0.25,
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.425), from = "lines", to = "user")),
           col = "gray50",
           xpd = NA
           )
  mtext(text = c("Baseline", "T0", "Endpoint")[are_included[are_rep1]],
        at   = tapply(x_positions, groups_vec[are_included], mean),
        side = 1,
        line = 1.7,
        cex  = par("cex")
        )
  return(invisible(NULL))
}



# Functions for displaying count-level QC data ----------------------------

CountBarPlot <- function(use_counts_mat,
                         include_timepoints = 1:3,
                         gini_index = FALSE,
                         lollipop   = FALSE,
                         bar_color  = brewer.pal(9, "Blues")[[8]],
                         stem_color = brewer.pal(9, "Blues")[[2]],
                         use_title  = NULL,
                         show_title = TRUE
                         ) {

  ## Prepare the timepoints included
  all_timepoints_vec <- rep(1:3, each = 2)
  are_included <- all_timepoints_vec %in% include_timepoints
  use_counts_mat <- use_counts_mat[, are_included]

  ## Determine bar positions
  bar_positions <- RepositionByGroups(all_timepoints_vec[are_included], gap_ratio = 1.5)
  num_bars <- length(bar_positions)
  bar_width <- 0.55
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - 0.5) - (num_bars * 0.04),
                     max(bar_positions) + 0.5  + (num_bars * 0.04)
                    )

  ## Prepare the data axis
  if (gini_index) {
    bars_vec <- vapply(seq_len(ncol(use_counts_mat)),
                       function(x) DescTools::Gini(use_counts_mat[, x]),
                       numeric(1)
                       )
    if (show_title && is.null(use_title)) {
      use_title <- "Plasmid count inequality"
    }
    if (show_title && isTRUE(grepl("count", use_title, ignore.case = TRUE))) {
      y_axis_label <- "Gini index"
    } else {
      y_axis_label <- "Gini index (count)"
    }
  } else {
    bars_vec <- colSums(use_counts_mat == 0)
    use_title <- "Absent plasmids"
    y_axis_label <- "Number of missed plasmids"
  }

  use_numeric_limits <- c(0, max(bars_vec) * 1.00)
  numeric_axis_pos <- pretty(use_numeric_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])

  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)

  if (lollipop) {
    grid_pos <- pretty(numeric_limits, n = 30)
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[2]],
             y0  = grid_pos,
             col = ifelse(((grid_pos * 100) %% 10) < (10^-12),
                          "gray88", "gray95"
                          ),
             xpd = NA
             )
    segments(x0  = bar_positions,
             y0  = 0,
             y1  = par("usr")[[4]],
             col = stem_color,
             xpd = NA,
             lwd = par("lwd") * 2
             )
    points(x   = bar_positions,
           y   = bars_vec,
           col = bar_color,
           cex = par("cex") * 1.5,
           pch = 16,
           xpd = NA
           )
  } else {
    PlotBarplotMat(t(bars_vec),
                   colors_vec    = bar_color,
                   positions_vec = bar_positions,
                   bar_width     = bar_width
                   )
  }

  ## Draw the y axis
  axis(2,
       las = 2,
       mgp = c(3, 0.55, 0),
       tcl = -0.375,
       lwd = par("lwd")
       )
  mtext(VerticalAdjust(y_axis_label),
        side = 2,
        line = 2.35,
        cex  = par("cex")
        )

  ## Draw the bar labels
  DrawBottomLabels(bar_positions, all_timepoints_vec, are_included)

  ## Final steps
  if (show_title) {
    title(use_title, font.main = 1, cex.main = 1)
  }
  box(bty = "l")

  names(bars_vec) <- colnames(use_counts_mat)[are_included]
  return(invisible(bars_vec))
}



CountBoxPlot <- function(use_counts_mat,
                         include_timepoints = 1:3,
                         embed_PNG = FALSE,
                         use_title = "Plasmid count"
                         ) {

  ## Prepare the timepoints included
  all_timepoints_vec <- rep(1:3, each = 2)
  are_included <- all_timepoints_vec %in% include_timepoints
  group_positions <- RepositionByGroups(all_timepoints_vec[are_included])
  numeric_mat <- log10(use_counts_mat[, are_included] + 1)

  BeeViolinPlot(as.list(data.frame(numeric_mat)),
                groups_vec    = all_timepoints_vec[are_included],
                use_swarm     = FALSE,
                cloud_alpha   = 0.08,
                cloud_sd      = 0.04,
                draw_border   = TRUE,
                violin_colors = brewer.pal(9, "Blues")[[3]],
                line_colors   = brewer.pal(9, "Blues")[[8]],
                border_colors = "#5e91c5",
                point_colors  = "#518dc2",
                # upper_bound   = 4,
                embed_PNG     = embed_PNG,
                draw_groups_n = FALSE
                )

  mtext(VerticalAdjust(expression("Log"[10] ~ "normalized count")),
        side = 2,
        line = 1.8,
        cex = par("cex")
        )
  DrawBottomLabels(group_positions, all_timepoints_vec, are_included)
  if (!(is.null(use_title))) {
    title(use_title, font.main = 1, cex.main = 1)
  }
  return(invisible(NULL))
}



GammaBoxPlot <- function(use_counts_df,
                         baseline_indices     = 1:2,
                         intervention_indices = 5:6,
                         both_timepoints      = TRUE,
                         num_cell_divisions   = 10L,
                         use_title            = "All genes",
                         cloud_alpha          = 0.15,
                         cloud_sd             = 0.025,
                         y_label_line         = 2.35,
                         show_y_axis          = TRUE,
                         use_swarm            = FALSE,
                         violin_colors        = brewer.pal(9, "Blues")[[3]],
                         line_colors          = brewer.pal(9, "Blues")[[8]],
                         border_colors        = brewer.pal(9, "Blues")[[8]],
                         point_colors         = "#518dc2",
                         ...
                         ) {

  if (both_timepoints) {
    logfc_1_vec <- GetLog2FC(use_counts_df,
                             choose_rep = 1,
                             baseline_indices = 1:2,
                             intervention_indices = 5:6,
                             allow_switch = FALSE
                             )[, "Log2FC"]
    logfc_2_vec <- GetLog2FC(use_counts_df,
                             choose_rep = 2,
                             baseline_indices = 1:2,
                             intervention_indices = 5:6,
                             allow_switch = FALSE
                             )[, "Log2FC"]
    logfc_3_vec <- GetLog2FC(use_counts_df,
                             choose_rep = 1,
                             baseline_indices = 3:4,
                             intervention_indices = 5:6,
                             allow_switch = FALSE
                             )[, "Log2FC"]
    logfc_4_vec <- GetLog2FC(use_counts_df,
                             choose_rep = 2,
                             baseline_indices = 3:4,
                             intervention_indices = 5:6,
                             allow_switch = FALSE
                             )[, "Log2FC"]
    logfc_list <- list(logfc_1_vec, logfc_2_vec, logfc_3_vec, logfc_4_vec)
    groups_vec <- c(1, 1, 2, 2)
    group_positions <- RepositionByGroups(groups_vec)
  } else {
    logfc_1_vec <- GetLog2FC(use_counts_df,
                             choose_rep = 1,
                             baseline_indices = baseline_indices,
                             intervention_indices = intervention_indices,
                             allow_switch = FALSE
                             )[, "Log2FC"]
    logfc_2_vec <- GetLog2FC(use_counts_df,
                             choose_rep = 2,
                             baseline_indices = baseline_indices,
                             intervention_indices = intervention_indices,
                             allow_switch = FALSE
                             )[, "Log2FC"]
    logfc_list <- list(logfc_1_vec, logfc_2_vec)
    groups_vec <- 1:2
  }

  logfc_list <- lapply(logfc_list, function(x) x[!(is.nan(x))])
  gamma_list <- lapply(logfc_list, function(x) x / num_cell_divisions)

  BeeViolinPlot(gamma_list,
                groups_vec,
                lower_bound   = -0.8,
                upper_bound   = 0.4,
                use_swarm     = use_swarm,
                cloud_alpha   = cloud_alpha,
                cloud_sd      = cloud_sd,
                draw_border   = TRUE,
                violin_colors = violin_colors,
                line_colors   = line_colors,
                border_colors = border_colors,
                point_colors  = point_colors,
                draw_groups_n = FALSE,
                show_y_axis   = show_y_axis,
                ...
                )

  if (show_y_axis) {
    mtext(VerticalAdjust(expression("Phenotype (" * gamma * ")")),
          side = 2,
          line = y_label_line,
          cex = par("cex")
          )
  }

  if (both_timepoints) {
    mtext(text = paste0("R", rep(1:2, 2)),
          at = group_positions, side = 1, line = 0.3, cex = par("cex")
          )
    segments(x0  = group_positions[c(1, 3)] - 0.25,
             x1  = group_positions[c(2, 4)] + 0.25,
             y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.45), from = "lines", to = "user")),
             col = "black",
             xpd = NA
             )
    mtext(text = c("Baseline vs. endpoint", "T0 vs. endpoint"),
          at   = c(mean(group_positions[1:2]), mean(group_positions[3:4])),
          side = 1,
          line = 1.7,
          cex  = par("cex")
          )
  } else {
    mtext(text = paste0("R", 1:2),
          at = 1:2, side = 1, line = 0.3, cex = par("cex")
          )
  }

  if (!(is.null(use_title))) {
    title(use_title, font.main = 1, cex.main = 1)
  }

  return(invisible(NULL))
}



# Functions for displaying read-level QC data -----------------------------

MappedReadsBarPlot <- function(num_reads_mat,
                               include_timepoints    = 1:3,
                               use_colors            = NULL,
                               use_title             = TRUE,
                               set_mar               = TRUE,
                               y_axis_mgp            = 0.5,
                               y_axis_tcl            = 0.375,
                               y_axis_label_line     = 2.35,
                               show_legend           = TRUE,
                               show_y_axis           = TRUE,
                               legend_point_x_start  = 0,
                               legend_lines_x_start  = 0.8,
                               y_upper_limit         = NULL,
                               bar_width             = 2/3,
                               gap_ratio             = 1.25,
                               side_gap              = 0.5,
                               unit_in_axis          = TRUE,
                               show_percentage       = FALSE,
                               large_gap_multiplier  = 1.5,
                               ...
                               ) {

  if (isTRUE(use_title)) {
    if (show_percentage) {
       use_title <- "Mapped reads"
    } else {
      use_title <- "All sequencing reads"
    }
  }

  if (is.null(use_colors)) {
    use_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:8])(nrow(num_reads_mat))
    if (nrow(num_reads_mat) == 4) {
      use_colors <- use_colors[c(1, 3, 2, 4)]
    }
  }

  ## Prepare the timepoints included
  all_timepoints_vec <- rep(1:3, each = 2)
  are_included <- all_timepoints_vec %in% include_timepoints
  numeric_mat <- num_reads_mat[, are_included]

  ## Determine bar positions
  bar_positions <- RepositionByGroups(all_timepoints_vec[are_included], gap_ratio = gap_ratio)
  num_bars <- length(bar_positions)
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - side_gap) - (num_bars * 0.04),
                     max(bar_positions) + side_gap  + (num_bars * 0.04)
                    )

  ## Prepare the data axis
  use_numeric_limits <- c(0, max(colSums(numeric_mat)) * 1.00)
  numeric_axis_pos <- pretty(use_numeric_limits)
  if (is.null(y_upper_limit)) {
    y_upper_limit <- numeric_axis_pos[[length(numeric_axis_pos)]]
  }
  numeric_limits <- c(numeric_axis_pos[[1]], y_upper_limit)

  ## Draw the barplot
  if (set_mar) {
    old_mar <- par(mar = c(4, 3.75, 3.75, 6.5))
  }
  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)
  if (show_percentage) {
    grid_pos <- pretty(use_numeric_limits, n = 30)
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[2]],
             y0  = grid_pos,
             col = ifelse(rep_len(c(TRUE, FALSE), length.out = length(grid_pos)),
                          "gray88", "gray95"
                          ),
             xpd = NA
             )
  }
  PlotBarplotMat(numeric_mat,
                 colors_vec    = use_colors,
                 positions_vec = bar_positions,
                 bar_width     = bar_width
                 )

  ## Draw the y axis
  if (show_y_axis) {
    tick_pos <- axTicks(2)
    if (show_percentage) {
      tick_labels <- paste0(tick_pos * 100, if (unit_in_axis) "%" else "")
    } else {
      tick_labels <- paste0(tick_pos / 10^6, if (unit_in_axis) "M" else "")
    }
    axis(2,
         at     = tick_pos,
         labels = tick_labels,
         las    = 2,
         mgp    = c(3, y_axis_mgp, 0),
         tcl    = -(y_axis_tcl),
         lwd    = par("lwd")
         )
    if (show_percentage) {
      y_axis_label <- "Percentage of reads"
    } else {
      if (unit_in_axis) {
        y_axis_label <- "Number of reads"
      } else {
        y_axis_label <- "Number of reads (millions)"
      }
    }
    mtext(VerticalAdjust(y_axis_label),
          side = 2,
          line = y_axis_label_line,
          cex  = par("cex")
          )
  }

  ## Draw the bar labels
  DrawBottomLabels(bar_positions, all_timepoints_vec, are_included)

  ## Draw the legend
  if (show_legend) {
    if (show_percentage) {
      if (length(use_colors) == 2) {
        title_vec <- c("Template", "switch")
        labels_list <- list("Yes", "No")
      } else {
        title_vec <- c("One", "mismatch", "tolerance")
        labels_list <- list(c("Both reads"),
                            c("Read 1"),
                            c("Read 2"),
                            c("None")
                            )
      }
    } else {
      labels_list <- list(c("Both reads", "unmapped"),
                          c("Read 1", "unmapped"),
                          c("Read 2", "unmapped"),
                          c("Both reads", "mapped")
                          )
      title_vec <- NULL
    }
    DrawSideLegend(labels_list,
                   use_colors     = rev(use_colors),
                   border_colors  = "gray50",
                   use_pch        = 22,
                   point_x_start  = legend_point_x_start,
                   lines_x_start  = legend_lines_x_start,
                   title_vec      = title_vec,
                   large_gap_multiplier = large_gap_multiplier,
                   ...
                   )
  }

  ## Final steps
  if (!(is.null(use_title))) {
    title(use_title, font.main = 1, cex.main = 1)
  }
  if (show_y_axis) {
    box(bty = "l")
  } else {
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[2]],
             y0  = par("usr")[[3]],
             xpd = NA
             )
  }
  if (set_mar) {
    par(mar = old_mar)
  }
  return(invisible(NULL))
}



PercentageBarPlot <- function(percentages_vec,
                              include_timepoints = 1:3,
                              dark_color         = brewer.pal(9, "Blues")[[8]],
                              light_color        = brewer.pal(9, "Blues")[[3]],
                              use_title          = "Template switch"
                              ) {

  ## Prepare the timepoints included
  all_timepoints_vec <- rep(1:3, each = 2)
  are_included <- all_timepoints_vec %in% include_timepoints
  percentages_vec <- percentages_vec[are_included]

  ## Determine bar positions
  bar_positions <- RepositionByGroups(all_timepoints_vec[are_included])
  num_bars <- length(bar_positions)
  bar_width <- 2/3
  final_width <- bar_width * ((max(bar_positions) - min(bar_positions)) / (num_bars - 1))
  group_limits <- c((min(bar_positions) - 0.5) - (num_bars * 0.04),
                     max(bar_positions) + 0.5  + (num_bars * 0.04)
                    )

  ## Prepare the data axis
  if (isTRUE(grepl("template switch", use_title, ignore.case = TRUE))) {
    y_axis_label <- "% mapped reads"
  } else {
    y_axis_label <- "Template switch"
  }
  y_ticks_pos <- seq(0, 1, by = 0.2)

  plot(NA,
       xlim = group_limits,
       ylim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       ann  = FALSE
       )

  PlotBarplotMat(rbind(1 - percentages_vec, percentages_vec),
                 colors_vec    = c(dark_color, light_color),
                 positions_vec = bar_positions,
                 bar_width     = bar_width
                 )

  ## Draw the y axis
  axis(2,
       at     = y_ticks_pos,
       labels = paste0(y_ticks_pos * 100, "%"),
       las    = 2,
       mgp    = c(3, 0.5, 0),
       tcl    = -0.375,
       lwd    = par("lwd")
       )
  mtext(VerticalAdjust(y_axis_label),
        side = 2,
        line = 2.5,
        cex  = par("cex")
        )

  ## Draw the bar labels
  DrawBottomLabels(bar_positions, all_timepoints_vec, are_included)

  ## Final steps
  if (!(is.null(use_title))) {
    title(use_title, font.main = 1, cex.main = 1)
  }
  box(bty = "l")

  return(invisible(NULL))
}



MakeEmptyPlot <- function(x_limits = c(0, 1), y_limits = c(0, 1)) {
  plot(NA,
       xlim = x_limits,
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       ann  = FALSE
       )
}



GetLibraryDensities <- function(sg_mat) {
  sg_mat <- toupper(sg_mat)
  char_vec_1 <- strsplit(sg_mat[, 1], "")
  char_vec_2 <- strsplit(sg_mat[, 2], "")
  GC_vec_1 <- vapply(char_vec_1, function(x) sum(x %in% c("C", "G")), integer(1))
  GC_vec_2 <- vapply(char_vec_2, function(x) sum(x %in% c("C", "G")), integer(1))
  density_1 <- density(GC_vec_1 / 20, adj = 3, from = 0, to = 1)
  density_2 <- density(GC_vec_2 / 20, adj = 3, from = 0, to = 1)
  density_list <- list(density_1, density_2)
  return(density_list)
}



TwoDensities <- function(show_GC               = TRUE,
                         include_timepoints    = 1:3,
                         semitransparent_lines = TRUE,
                         show_title            = TRUE,
                         use_title             = NULL,
                         include_zero          = TRUE,
                         label_read_on_y_axis  = TRUE,
                         show_y_axis_label     = TRUE,
                         show_legend           = TRUE,
                         omit_zero_label       = FALSE,
                         x_axis_mgp            = 0.375,
                         x_axis_label_line     = 2.1,
                         y_axis_label_line     = 0.5,
                         legend_x_lines        = 1,
                         legend_y_lines        = 1,
                         embed_PNG             = FALSE,
                         grid_lwd              = 1,
                         title_y_pos           = 0.5,
                         darker_box            = FALSE,
                         broad_margins         = FALSE
                         ) {

  if (show_GC) {
    density_object <- "GC_content_densities"
  } else {
    density_object <- "sequence_qual_densities"
    stopifnot("library_densities" %in% ls(envir = globalenv()))
  }
  stopifnot(density_object %in% ls(envir = globalenv()))
  density_list <- get(density_object)[c("Read 1", "Read 2")]

  ## Prepare the timepoints included
  timepoint_labels <- c("Baseline", "T0", "Endpoint")[include_timepoints]
  all_timepoints_vec <- rep(1:3, each = 2)
  are_included <- all_timepoints_vec %in% include_timepoints
  original_colors <- c(
    brewer.pal(9, "Blues")[[8]],
    brewer.pal(9, "Purples")[[8]],
    brewer.pal(9, "Reds")[[8]]
  )[include_timepoints]
  if (semitransparent_lines) {
    timepoint_colors <- adjustcolor(original_colors, alpha.f = 0.6)
  } else {
    timepoint_colors <- original_colors
  }
  colors_vec <- rep(timepoint_colors, each = 2)
  density_list <- lapply(density_list, function(x) x[are_included])

  ## Prepare the x axis
  x_range <- range(unlist(lapply(density_list, function(x) sapply(x, "[[", "x"))))
  if (include_zero) {
    x_range[[1]] <- min(0, x_range[[1]])
  }

  x_ticks <- pretty(x_range)
  x_limits <- range(x_ticks)
  x_grid <- pretty(x_limits, n = 10)
  x_grid <- x_grid[x_grid != 0]
  if (show_GC) {
    x_axis_label <- "GC content"
    x_tick_labels <- paste0(x_ticks * 100, "%")
    if (omit_zero_label) {
      x_tick_labels[x_tick_labels == "0%"] <- NA
    }
    if (is.null(use_title)) {
      use_title <- "Density (all reads)"
    }
  } else {
    x_axis_label <- "Mean sequence quality"
    x_tick_labels <- x_ticks
    if (is.null(use_title)) {
      use_title <- "Density (all reads)"
    }
  }

  ## Prepare the y axis
  y_max <- max(unlist(lapply(density_list, function(x) sapply(x, "[[", "y"))))
  y_limits <- c(0, y_max * 1.04)

  ## Additional preparations
  if (show_GC) {
    grid_color <- "gray86"
  } else {
    grid_color <- "gray80"
  }
  layout_mat <- rbind(
    c(3, 1, 4),
    c(3, 6, 4),
    c(3, 5, 4),
    c(3, 7, 4),
    c(3, 2, 4)
  )

  ## Set up the multi-plot layout
  original_cex <- par("cex")
  if (broad_margins) {
    use_widths <- c(2.25, 7, 0.75)
    use_heights <- c(1.2, 3.27, 0.21, 3.27, 2)
  } else {
    use_widths <- c(1, 8, 1)
    use_heights <- c(1.2, 3.4, 0.25, 3.4, 1.75)
  }
  layout(layout_mat,
         widths  = use_widths,
         heights = use_heights
         )
  old_mar <- par(mar = rep(0, 4), cex = original_cex * par("cex") * (1 / 0.66))
  MakeEmptyPlot()
  if (show_title) {
    text(x = 0.5, y = title_y_pos, labels = use_title)
  }
  for (i in 1:4) {
    MakeEmptyPlot()
  }

  ## Draw the two density plots
  for (i in 1:2) {

    MakeEmptyPlot(x_limits, y_limits)

    if (embed_PNG) {
      current_device <- StartEmbedPNG(figures_dir, use_cairo = TRUE, add_padding = TRUE)
      MakeEmptyPlot(x_limits, y_limits)
    }

    if (!(show_GC)) {
      ## Draw the background
      rect(xleft   = 28,
           xright  = 40,
           ybottom = par("usr")[[3]],
           ytop    = par("usr")[[4]],
           col     = adjustcolor("#C6E8BF", alpha.f = 0.4),
           border  = NA
           )
      rect(xleft   = 20,
           xright  = 28,
           ybottom = par("usr")[[3]],
           ytop    = par("usr")[[4]],
           col     = adjustcolor("#FFEDA0", alpha.f = 0.4),
           border  = NA
           )
      rect(xleft   = 0,
           xright  = 20,
           ybottom = par("usr")[[3]],
           ytop    = par("usr")[[4]],
           col     = adjustcolor("#FCC9B3", alpha.f = 0.4),
           border  = NA
           )
    }

    ## Draw the grid
    if (show_legend) {
      if (label_read_on_y_axis) {
        are_behind_read <- rep(FALSE, length(x_grid))
      } else {
        are_behind_read <- seq_along(x_grid) %in% 1:2
      }
      are_behind_legend <- rep(FALSE, length(x_grid))
      if (i == 1) {
        ## Calculate how much white space to leave for the legend
        legend_x_start <- par("usr")[[1]] + diff(grconvertX(c(0, legend_x_lines), from = "lines", to = "user"))
        legend_text_x <- legend_x_start + diff(grconvertX(c(0, 0.9), from = "lines", to = "user"))
        legend_y_start <- par("usr")[[4]] - diff(grconvertY(c(0, legend_y_lines), from = "lines", to = "user"))
        timepoints_seq <- seq_along(timepoint_colors) - 1L
        legend_y_vec <- legend_y_start - diff(grconvertY(c(0, 1.2), from = "lines", to = "user")) * timepoints_seq

        legend_x_end <- legend_text_x + max(strwidth(timepoint_labels)) + strwidth("OO")
        are_behind_legend <- x_grid < legend_x_end
        legend_y_end <- min(legend_y_vec) - (strheight(timepoint_labels[[length(timepoint_labels)]])) * 1.2
        if ((length(include_timepoints) == 2) && (grconvertY(legend_y_end, from = "user", to = "npc") > 0.5)) {
          grid_top <- 0.5
        } else {
          grid_top <- 1/3
        }
        segments(x0   = x_grid[are_behind_legend & !(are_behind_read)],
                 y0   = par("usr")[[3]],
                 y1   = grconvertY(grid_top, from = "npc", to = "user"),
                 col  = grid_color,
                 lend = "butt",
                 lwd  = par("lwd") * grid_lwd
                 )
      } else {
        are_behind_legend <- rep(FALSE, length(x_grid))
        if (!(label_read_on_y_axis)) {
          segments(x0   = x_grid[are_behind_read],
                   y0   = grconvertY(1/3, from = "npc", to = "user"),
                   y1   = par("usr")[[4]],
                   col  = grid_color,
                   lend = "butt",
                   lwd  = par("lwd") * grid_lwd
                   )
        }
      }
    } else {
      are_behind_legend <- rep(FALSE, length(x_grid))
      are_behind_read <- rep(FALSE, length(x_grid))
    }
    abline(v = x_grid[!(are_behind_legend | are_behind_read)],
           col = grid_color, lwd = par("lwd") * grid_lwd
           )

    box(col = if (darker_box) "gray50" else grid_color,
        lwd = par("lwd") * grid_lwd,
        xpd = NA
        )

    ## Draw the density lines
    if (show_GC) {
      ref_x_vec <- library_densities[[i]][["x"]]
      polygon(x      = c(ref_x_vec[[1]], ref_x_vec, ref_x_vec[[length(ref_x_vec)]]),
              y      = c(0, library_densities[[i]][["y"]], 0),
              col    = adjustcolor("gray80", alpha.f = 0.25),
              border = NA,
              xpd    = NA
              )
      lines(library_densities[[i]], col = "gray80", lwd = par("lwd"))
    }

    for (j in seq_along(density_list[[i]])) {
      lines(density_list[[i]][[j]],
            col  = colors_vec[[j]],
            lwd  = par("lwd") * 2,
            xpd  = NA
            )
    }

    if (embed_PNG) {
      StopEmbedPNG(current_device, figures_dir, add_padding = TRUE,
                   make_empty_plot = FALSE
                   )
    }

    ## Draw the y axis
    if (show_y_axis_label) {
      mtext(if (label_read_on_y_axis) paste0("Read ", i) else "Density",
            side = 2,
            line = y_axis_label_line,
            cex  = par("cex")
            )
    }

    if (show_legend && (i == 1)) {
      ## Draw the legend
      segments(x0  = legend_x_start,
               x1  = legend_x_start + diff(grconvertX(c(0, 0.45), from = "lines", to = "user")),
               y0  = legend_y_vec,
               col = if (embed_PNG) original_colors else timepoint_colors,
               lwd = par("lwd") * 2,
               xpd = NA
               )
      text(x      = legend_text_x,
           y      = legend_y_vec,
           labels = timepoint_labels,
           adj    = c(0, 0.5),
           xpd    = NA
           )
    } else if (i == 2) {
      ## Draw the x axis
      axis(1,
           at       = x_ticks,
           labels   = x_tick_labels,
           mgp      = c(3, x_axis_mgp, 0.3),
           tcl      = -0.375,
           gap.axis = 0.5,
           lwd      = par("lwd")
           )
      mtext(x_axis_label,
            side = 1,
            line = x_axis_label_line,
            cex  = par("cex")
            )
    }
    if (show_legend && !(label_read_on_y_axis)) {
      mtext(paste0("Read ", i), side = 1, line = -1.5,
            at = grconvertX(0.15, from = "npc", to = "user"), cex = par("cex")
            )
    }
  }

  ## Final steps
  par(old_mar)
  layout(1)
  return(invisible(NULL))
}



PerBaseQuality <- function(qual_mat,
                           include_timepoints    = 1:3,
                           semitransparent_lines = TRUE,
                           use_title             = "Per-base quality (all reads)",
                           show_title            = TRUE,
                           show_y_axis           = TRUE,
                           show_legend           = TRUE,
                           omit_zero_label       = FALSE,
                           x_axis_mgp            = 0.375,
                           y_axis_mgp            = 0.55,
                           x_axis_label_line     = 1.9,
                           y_axis_label_line     = 2.2,
                           title_y_pos           = 0.5,
                           separate_x_labels     = FALSE,
                           legend_x_lines        = 2.5,
                           legend_y_lines        = 3,
                           embed_PNG             = FALSE,
                           small_middle_gap      = FALSE,
                           broad_margins         = FALSE,
                           x_axis_tcl            = 0.375
                           ) {

  ## Prepare the timepoints included
  timepoint_labels <- c("Baseline", "T0", "Endpoint")[include_timepoints]
  all_timepoints_vec <- rep(1:3, each = 2)
  are_included <- all_timepoints_vec %in% include_timepoints
  original_colors <- c(
    brewer.pal(9, "Blues")[[8]],
    brewer.pal(9, "Purples")[[8]],
    brewer.pal(9, "Reds")[[8]]
  )[include_timepoints]
  colors_vec <- rep(original_colors, each = 2)
  if (semitransparent_lines) {
    colors_vec <- adjustcolor(colors_vec, alpha.f = 0.6)
    timepoint_colors <- adjustcolor(original_colors, alpha.f = 0.9) # For the legend
  } else {
    timepoint_colors <- original_colors
  }

  ## Set up the multi-plot layout
  layout_mat <- rbind(
    c(3, 1, 1, 1, 4),
    c(3, 6, 5, 7, 4),
    c(3, 2, 2, 2, 4)
  )
  original_cex <- par("cex")
  if (broad_margins) {
    if (small_middle_gap) {
      use_widths <- c(2.25, 3.425, 0.15, 3.425, 0.75)
    } else {
      use_widths <- c(2.25, 3.25, 0.5, 3.25, 0.75)
    }
    use_heights <- c(1.2, 6.8, 2)
  } else {
    if (small_middle_gap) {
      use_widths <- c(1.5, 3.925, 0.15, 3.925, 0.5)
    } else {
      use_widths <- c(1.5, 3.725, 0.55, 3.725, 0.5)
    }
    use_heights <- c(1.2, 7.05, 1.75)
  }
  layout(layout_mat,
         widths  = use_widths,
         heights = c(1.2, 7.05, 1.75)
         )
  old_mar <- par(mar = rep(0, 4), cex = original_cex * par("cex") * (1 / 0.66))
  MakeEmptyPlot()
  if (show_title) {
    text(x = 0.5, y = title_y_pos, labels = use_title)
  }
  for (i in 1:4) {
    MakeEmptyPlot()
  }

  ## Draw the two line graphs
  for (i in 1:2) {
    MakeEmptyPlot(x_limits = c(0, 21), y_limits = c(0, 40))

    ## Determine the sample indices
    are_this_read <- grepl(paste0("_read", i), colnames(qual_mat), fixed = TRUE)
    use_indices <- which(are_this_read)[are_included]

    if (embed_PNG) {
      current_device <- StartEmbedPNG(figures_dir, use_cairo = TRUE)
      MakeEmptyPlot(x_limits = c(0, 21), y_limits = c(0, 40))
    }

    ## Draw the background
    rect(xleft   = par("usr")[[1]],
         xright  = if (i == 1) 20 else par("usr")[[2]],
         ybottom = 28,
         ytop    = 40,
         col     = Palify("#C6E8BF", fraction_pale = 0.6),
         border  = NA
         )
    rect(xleft   = par("usr")[[1]],
         xright  = if (i == 1) 20 else par("usr")[[2]],
         ybottom = 20,
         ytop    = 28,
         col     = Palify("#FFEDA0", fraction_pale = 0.6),
         border  = NA
         )
    rect(xleft   = par("usr")[[1]],
         xright  = if (i == 1) 20 else par("usr")[[2]],
         ybottom = 0,
         ytop    = 20,
         col     = Palify("#FCC9B3", fraction_pale = 0.6),
         border  = NA
         )

    ## Draw the lines
    for (j in seq_along(use_indices)) {
      lines(x    = seq_len(nrow(qual_mat)),
            y    = qual_mat[, use_indices[[j]]],
            col  = colors_vec[[j]],
            lwd  = par("lwd") * 2,
            lend = "butt",
            xpd  = NA
            )
    }

    if (embed_PNG) {
      StopEmbedPNG(current_device, figures_dir, make_empty_plot = FALSE)
    }

    ## Draw the x axis
    x_ticks <- axTicks(1)
    x_tick_labels <- x_ticks
    if (omit_zero_label && (i == 2)) {
      x_tick_labels[x_tick_labels == 0] <- NA
    }
    axis(1,
         at       = x_ticks,
         labels   = x_tick_labels,
         mgp      = c(3, x_axis_mgp, 0),
         tcl      = -(x_axis_tcl),
         gap.axis = 0.5,
         lwd      = par("lwd")
         )
    mtext(if (separate_x_labels) "Base" else paste0("Read ", i, " (base)"),
          side = 1,
          line = x_axis_label_line,
          cex  = par("cex")
          )
    if (separate_x_labels) {
      text(x      = grconvertX(0.5, from = "npc", to = "user"),
           y      = grconvertY(0.6, from = "npc", to = "user"),
           labels = paste0("Read ", i),
           xpd     = NA
           )
    }

    if (i == 1) {
      ## Draw the y axis
      if (show_y_axis) {
        axis(2,
             las = 2,
             mgp = c(3, y_axis_mgp, 0),
             tcl = -0.375,
             lwd = par("lwd")
             )
        mtext(VerticalAdjust("Mean quality"),
              side = 2,
              line = y_axis_label_line,
              cex  = par("cex")
              )
      } else {
        segments(x0 = par("usr")[[1]], y0 = par("usr")[[3]],
                 y1 = par("usr")[[4]], xpd = NA
                 )
      }
      segments(x0 = par("usr")[[1]], x1 = 20,
               y0 = par("usr")[[4]], xpd = NA
               )
    } else if (i == 2) {
      box(bty = "]")
      if (show_legend) {
        ## Draw the legend
        x_start <- par("usr")[[2]] -
                   max(strwidth(timepoint_labels)) -
                   diff(grconvertX(c(0, legend_x_lines), from = "lines", to = "user"))
        y_start <- par("usr")[[3]] + diff(grconvertY(c(0, legend_y_lines), from = "lines", to = "user"))
        timepoints_seq <- seq_along(timepoint_colors) - 1L
        y_vec <- rev(y_start + diff(grconvertY(c(0, 1.2), from = "lines", to = "user")) * timepoints_seq)
        segments(x0  = x_start,
                 x1  = x_start + diff(grconvertX(c(0, 0.45), from = "lines", to = "user")),
                 y0  = y_vec,
                 col = if (embed_PNG) original_colors else timepoint_colors,
                 lwd = par("lwd") * 2,
                 xpd = NA
                 )
        text(x      = x_start + diff(grconvertX(c(0, 0.9), from = "lines", to = "user")),
             y      = y_vec,
             labels = timepoint_labels,
             adj    = c(0, 0.5),
             xpd    = NA
             )
      }
    }
  }

  ## Final steps
  par(old_mar)
  layout(1)
  return(invisible(NULL))
}



# Functions for assessing the separation between groups -------------------

FilterNonFinite <- function(vec1, vec2, lower_bound = NULL) {
  are_finite_vec1 <- is.finite(vec1)
  are_finite_vec2 <- is.finite(vec2)
  if (!(is.null(lower_bound))) {
    are_neginf_vec1 <- (!(are_finite_vec1)) & (vec1 < 0)
    are_neginf_vec2 <- (!(are_finite_vec2)) & (vec2 < 0)
    vec1[are_neginf_vec1] <- lower_bound
    vec2[are_neginf_vec2] <- lower_bound
    if ((any(are_neginf_vec1) || any(are_neginf_vec2))) {
      message(paste0(sum(are_neginf_vec1) + sum(are_neginf_vec2),
                     " -Inf values were replaced with the lower limit (i.e. ",
                     lower_bound, ")."
                     ))
    }
  }
  are_finite_vec1 <- is.finite(vec1)
  are_finite_vec2 <- is.finite(vec2)
  if (!(all(are_finite_vec1) && all(are_finite_vec2))) {
    message(paste0(sum(!(are_finite_vec1)) + sum(!(are_finite_vec2)),
                   " non-finite values were excluded."
                   ))
  }
  vec1 <- vec1[are_finite_vec1]
  vec2 <- vec2[are_finite_vec2]
  results_list <- list(vec1, vec2)
  return(results_list)
}


ZPrimeFactor <- function(neg_vec, pos_vec, lower_bound = NULL) {
  filtered_list <- FilterNonFinite(neg_vec, pos_vec, lower_bound = lower_bound)
  neg_vec <- filtered_list[[1]]
  pos_vec <- filtered_list[[2]]
  z_prime <- (1 - (3 * (sd(pos_vec) + sd(neg_vec))) / abs(mean(pos_vec) - mean(neg_vec)))
  return(z_prime)
}


RobustZPrimeFactor <- function(neg_vec, pos_vec, lower_bound = NULL) {
  filtered_list <- FilterNonFinite(neg_vec, pos_vec, lower_bound = lower_bound)
  neg_vec <- filtered_list[[1]]
  pos_vec <- filtered_list[[2]]
  z_prime <- (1 - (3 * (mad(pos_vec) + mad(neg_vec))) / abs(median(pos_vec) - median(neg_vec)))
  return(z_prime)
}


SSMDControls <- function(neg_vec, pos_vec, lower_bound = NULL) {
  filtered_list <- FilterNonFinite(neg_vec, pos_vec, lower_bound = lower_bound)
  neg_vec <- filtered_list[[1]]
  pos_vec <- filtered_list[[2]]
  SSMD <- (mean(pos_vec) - mean(neg_vec)) / (sqrt(var(pos_vec) + var(neg_vec)))
  return(SSMD)
}


RobustSSMDControls <- function(neg_vec, pos_vec, lower_bound = NULL) {
  filtered_list <- FilterNonFinite(neg_vec, pos_vec, lower_bound = lower_bound)
  neg_vec <- filtered_list[[1]]
  pos_vec <- filtered_list[[2]]
  SSMD <- (median(pos_vec) - median(neg_vec)) / (1.4826 * sqrt((mad(pos_vec))^2 + (mad(neg_vec))^2))
  return(SSMD)
}


SeparationMetrics <- function(input_list) {
  all_methods <- c("ZPrimeFactor", "RobustZPrimeFactor",
                   "SSMDControls", "RobustSSMDControls"
                   )
  vec_list_all <- lapply(1:2, function(rep_index) vapply(all_methods, function(x) {
    get(x)(neg_vec = input_list[[paste0("Non-essential R", rep_index)]],
           pos_vec = input_list[[paste0("Essential R", rep_index)]],
           lower_bound = -0.6
           )
  }, numeric(1)))
  vec_list_finite <- lapply(1:2, function(rep_index) vapply(all_methods, function(x) {
    get(x)(neg_vec = input_list[[paste0("Non-essential R", rep_index)]],
           pos_vec = input_list[[paste0("Essential R", rep_index)]],
           lower_bound = NULL
           )
  }, numeric(1)))
  results_mat <- rbind(do.call(cbind, vec_list_all), do.call(cbind, vec_list_finite))
  dimnames(results_mat) <- list(c("Z' factor", "Robust z'", "SSMD", "Robust SSMD",
                                  "Z' finite", "Robust z' finite",
                                  "SSMD finite", "Robust SSMD finite"
                                  ),
                                c("R1", "R2")
                                )
  return(results_mat)
}




