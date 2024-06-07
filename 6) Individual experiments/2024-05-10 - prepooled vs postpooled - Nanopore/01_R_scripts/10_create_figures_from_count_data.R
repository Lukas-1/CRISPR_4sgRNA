### 2023-11-28


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

project_dir     <- file.path(experiments_directory, "2023-10-05 - prepooled vs postpooled - Nanopore")
rdata_dir       <- file.path(project_dir, "03_R_objects")
figures_dir     <- file.path(project_dir, "05_output", "Figures")
PDFs_dir        <- file.path(figures_dir, "PDFs")
manuscript_dir  <- file.path(PDFs_dir, "Manuscript-style")
first_rdata_dir <- file.path(first_illumina_trial_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "04_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "07_produce_read_counts_and_metrics.RData"))
# load(file.path(rdata_dir, "04_assign_sgRNAs_to_plasmids_num_mapped_reads.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))

load(file.path(CRISPR_root_directory, "3) RData files", "1) General",
               "24) Enumerate pairs of genes at bidirectional promoters.RData"
               ))



# Prepare R objects -------------------------------------------------------

CRISPRoff_df <- sg_sequences_df

CRISPRoff_df[, "Num_plasmids_for_Entrez"] <- table(CRISPRoff_df[, "Entrez_ID"])[CRISPRoff_df[, "Entrez_ID"]]

are_obsolete <- CRISPRoff_df[, "Is_obsolete"] %in% "Yes"
new_order <- order(
  as.integer(CRISPRoff_df[, "Entrez_ID"]),
  are_obsolete,
  match(CRISPRoff_df[, "Is_main_TSS"], c("Main", "Single", "Other"))
)
CRISPRoff_df <- CRISPRoff_df[new_order, ]
row.names(CRISPRoff_df) <- NULL

CRISPRoff_df[, "gene"] <- ifelse(startsWith(CRISPRoff_df[, "Gene_symbol"], "Control_"),
                                 "negative_control", ""
                                 )

counts_df <- counts_df[new_order, ]
row.names(counts_df) <- NULL

original_counts_df <- counts_df


use_count_columns <- grep("^Count_all4_", names(counts_df), value = TRUE)
other_columns <- grep("^Count_", names(counts_df), value = TRUE, invert = TRUE)
renamed_count_columns <- sub("^Count_all4_p", "NoSwitch_xMM_P", use_count_columns)
renamed_count_columns <- sub("_rep", "_R", renamed_count_columns, fixed = TRUE)
counts_df <- counts_df[, c(other_columns, use_count_columns)]
names(counts_df) <- c(other_columns, renamed_count_columns)




# Examine only the genes used in Nunez et al. -----------------------------
## (intersection of Blomen and Hart)

are_missing <- rowSums(counts_df[, c("NoSwitch_xMM_Prepool_T0_R1", "NoSwitch_xMM_Prepool_T0_R2")]) == 0
are_present <- rowSums(counts_df[, c("NoSwitch_xMM_Prepool_T0_R1", "NoSwitch_xMM_Prepool_T0_R2")]) != 0
missing_entrezs <- setdiff(counts_df[, "Entrez_ID"][are_missing], NA)
present_entrezs <- setdiff(counts_df[, "Entrez_ID"][are_present], NA)

common_entrezs <- intersect(missing_entrezs, present_entrezs)
missing_entrezs <- setdiff(missing_entrezs, common_entrezs)
present_entrezs <- setdiff(present_entrezs, common_entrezs)



# Plot the mean essentiality scores for present vs. absent sgRNAs ---------

missing_df <- data.frame(
  "Entrez_ID"  = as.integer(c(missing_entrezs, present_entrezs)),
  "Is_missing" = c(rep(TRUE,  length(missing_entrezs)),
                   rep(FALSE, length(present_entrezs))
                   )
)

matches_vec <- match(missing_df[, "Entrez_ID"], essential_df[, "Entrez_ID"])
missing_df[, "CRISPR_mean_probability"] <- essential_df[, "CRISPR_mean_probability"][matches_vec]

probs_list <- split(missing_df[, "CRISPR_mean_probability"],
                    missing_df[, "Is_missing"]
                    )
probs_list <- lapply(probs_list, function(x) x[!(is.na(x))])
set.seed(1)
probs_list[[1]] <- probs_list[[1]][seq_along(probs_list[[2]])]

for (create_PDF in c(FALSE)) {

  if (create_PDF) {
    pdf(file.path(PDFs_dir, "Essentiality of missing genes.pdf"),
        width = 4.5, height = 4.5
        )
  }

  old_mar <- par(mar = c(4, 4.3, 4, 2.1))
  BeeViolinPlot(probs_list, use_swarm = FALSE, point_cex = 0.35, wex = 0.9,
                y_limits = c(0, 1), sina_plot = TRUE, draw_violins = TRUE
                )
  mtext("Mean probability of essentiality", side = 2, line = 2.6)
  mtext(c("Plasmids\npresent", "Plasmids\nabsent"),
        side = 1, line = 1.65, at = 1:2
        )
  title("Essentiality of missing genes", cex.main = par("cex"))
  par(old_mar)

  if (create_PDF) {
    dev.off()
  }
}



# Draw ROC curves for the identification of essential genes ---------------

MakeROCDfList <- function(choose_rep = NULL, use_blomen_hart = TRUE) {
  essential_list <- GetEssentialGenes(use_blomen_hart)
  args_list <- list(essential_genes       = essential_list[["essential_entrezs"]],
                    non_essential_genes   = essential_list[["non_essential_entrezs"]],
                    min_count_at_baseline = 0L,
                    choose_rep            = choose_rep
                    )
  ROC_df_list <- list(
      ROC_prepool = do.call(GetEssentialROCDf, c(args_list, list(
        baseline_indices      = 1:2,
        intervention_indices  = 3:4,
        allow_switch          = FALSE
      ))),
      ROC_postpool = do.call(GetEssentialROCDf, c(args_list, list(
        baseline_indices      = 5:6,
        intervention_indices  = 7:8,
        allow_switch          = FALSE
      )))
  )
  return(ROC_df_list)
}


prepooled_title <- expression("CRISPRoff with T.gonfio:" ~ bold("pre-pooled"))
postpooled_title <- expression("CRISPRoff with T.gonfio:" ~ bold("post-pooled"))

both_reps_ROC_df_list <- MakeROCDfList()
rep1_ROC_df_list <- MakeROCDfList(choose_rep = 1)
rep2_ROC_df_list <- MakeROCDfList(choose_rep = 2)

both_reps_more_genes_ROC_df_list <- MakeROCDfList(use_blomen_hart = FALSE)
rep1_ROC_more_genes_df_list <- MakeROCDfList(choose_rep = 1, use_blomen_hart = FALSE)
rep2_ROC_more_genes_df_list <- MakeROCDfList(choose_rep = 2, use_blomen_hart = FALSE)


old_oma <- par(oma = rep(0.3, 4))
for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(PDFs_dir, "ROC curves - both replicates.pdf"),
        width = 3 + 1.14, height = 3 + 1.56
        )
  }
  AUC_prepooled_both_reps <- PlotEssentialROCDf(both_reps_ROC_df_list[["ROC_prepool"]], use_title = prepooled_title)
  AUC_postpooled_both_reps <- PlotEssentialROCDf(both_reps_ROC_df_list[["ROC_postpool"]],  use_title = postpooled_title)
  PlotEssentialROCDf(both_reps_more_genes_ROC_df_list[["ROC_prepool"]], use_title = prepooled_title)
  PlotEssentialROCDf(both_reps_more_genes_ROC_df_list[["ROC_postpool"]],  use_title = postpooled_title)
  if (create_PDF) {
    dev.off()
    pdf(file.path(PDFs_dir, "ROC curves - replicate 1.pdf"),
        width = 3 + 1.14, height = 3 + 1.56
    )
  }
  AUC_prepooled_rep1 <- PlotEssentialROCDf(rep1_ROC_df_list[["ROC_prepool"]],   use_title = expression("Replicate 1 only:" ~ bold("pre-pooled")))
  AUC_postpooled_rep1 <- PlotEssentialROCDf(rep1_ROC_df_list[["ROC_postpool"]],  use_title = expression("Replicate 1 only:" ~ bold("post-pooled")))
  PlotEssentialROCDf(rep1_ROC_more_genes_df_list[["ROC_prepool"]],   use_title = expression("Replicate 1 only:" ~ bold("pre-pooled")))
  PlotEssentialROCDf(rep1_ROC_more_genes_df_list[["ROC_postpool"]],  use_title = expression("Replicate 1 only:" ~ bold("post-pooled")))
  if (create_PDF) {
    dev.off()
    pdf(file.path(PDFs_dir, "ROC curves - replicate 2.pdf"),
        width = 3 + 1.14, height = 3 + 1.56
    )
  }
  AUC_prepooled_rep2 <- PlotEssentialROCDf(rep2_ROC_df_list[["ROC_prepool"]],   use_title = expression("Replicate 2 only:" ~ bold("pre-pooled")))
  AUC_postpooled_rep2 <- PlotEssentialROCDf(rep2_ROC_df_list[["ROC_postpool"]],  use_title = expression("Replicate 2 only:" ~ bold("post-pooled")))
  PlotEssentialROCDf(rep2_ROC_more_genes_df_list[["ROC_prepool"]],   use_title = expression("Replicate 2 only:" ~ bold("pre-pooled")))
  PlotEssentialROCDf(rep2_ROC_more_genes_df_list[["ROC_postpool"]],  use_title = expression("Replicate 2 only:" ~ bold("post-pooled")))
  if (create_PDF) {
    dev.off()
  }
}
par(old_oma)


use_df <- both_reps_ROC_df_list[["ROC_prepool"]]
pdf(file.path(manuscript_dir, "ROC curve - pre-pooled.pdf"),
    width = 2, height = 2
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
PlotROCDf(use_df, flip = TRUE, xlab_line = 1.6, ylab_line = 2.1, ROC_lwd = 1.5)
title("Pre-pooled T.gonfio library", cex.main = 1, font.main = 1, line = 0.7)
par(old_par)
dev.off()

use_df <- both_reps_ROC_df_list[["ROC_postpool"]]
pdf(file.path(manuscript_dir, "ROC curve - post-pooled.pdf"),
    width = 2, height = 2
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
PlotROCDf(use_df, flip = TRUE, xlab_line = 1.6, ylab_line = 2.1, ROC_lwd = 1.5)
title("Post-pooled T.gonfio library", cex.main = 1, font.main = 1, line = 0.7)
par(old_par)
dev.off()



# Draw violin plots showing the log2FC ------------------------------------

y_span <- 0.7 * 0.02
custom_y_limits <- c(-0.6 - y_span, 0.2)
custom_args <- list(lower_bound = -0.6, upper_bound = 0.2, y_limits = custom_y_limits)

for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(PDFs_dir, "Violin plots - mean of replicates.pdf"),
        width = 4.5, height = 4.5
        )
  }
  do.call(ViolinPlotEssentialDf, c(list(both_reps_ROC_df_list[["ROC_prepool"]], use_title = prepooled_title), custom_args))
  do.call(ViolinPlotEssentialDf, c(list(both_reps_ROC_df_list[["ROC_postpool"]],  use_title = postpooled_title), custom_args))
  if (create_PDF) {
    dev.off()
  }
}



# Draw violin plots showing the log2FC for both replicates ----------------

RepEssentialViolins(1:2, 3:4, use_title = prepooled_title, lower_bound = -0.6, upper_bound = 0.2, min_count_at_baseline = 0, y_limits = c(-0.61, 0.2),
                    allow_switch = FALSE
                    )
RepEssentialViolins(5:6, 7:8, use_title = postpooled_title, lower_bound = -0.6, upper_bound = 0.2, min_count_at_baseline = 0, y_limits = c(-0.61, 0.2),
                    allow_switch = FALSE
                    )


custom_y_limits <- c(-0.6 - (0.86 * 0.02), 0.2)
for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(PDFs_dir, "Violin plots - R1 and R2.pdf"), width = 5, height = 4.5)
  }
  for (use_blomen_hart in c(TRUE, FALSE)) {
    RepEssentialViolins(1:2, 3:4, use_title = prepooled_title,
                        lower_bound = -0.6, upper_bound = 0.2, y_limits = custom_y_limits,
                        allow_switch = FALSE, use_blomen_hart = use_blomen_hart
                        )
    RepEssentialViolins(5:6, 7:8, use_title = postpooled_title,
                        lower_bound = -0.6, upper_bound = 0.2, y_limits = custom_y_limits,
                        allow_switch = FALSE, use_blomen_hart = use_blomen_hart
                        )
  }
  if (create_PDF) {
    dev.off()
  }
}



reps_list_list <- lapply(c(TRUE, FALSE), function(prepooled) {

  if (prepooled) {
    file_name_suffix <- "prepool"
    use_title <- "Prepool T.gonfio"
  } else {
    file_name_suffix <- "postpool"
    use_title <- "Postpool T.gonfio"
  }
  file_name <- paste0("Violin plots - ", file_name_suffix, ".pdf")

  pdf(file.path(manuscript_dir, file_name),
      width = 2, height = 2
      )
  old_par <- par(cex = 0.6, lwd = 0.8, lheight = 0.9)

  reps_list <- RepEssentialViolins(
    if (prepooled) 1:2 else 5:6,
    if (prepooled) 3:4 else 7:8,
    use_title        = as.expression(bquote(bold(.(use_title)))),
    lower_bound      = -0.6,
    upper_bound      = 0.2,
    y_limits         = custom_y_limits,
    show_truncation  = FALSE,
    allow_switch     = FALSE,
    use_mar          = c(3, 4, 4, 1),
    y_axis_label     = expression("Phenotype (" * gamma * ")"),
    y_label_line     = 2.1,
    rep_label_line   = 0.4,
    genes_label_line = 0.6,
    draw_groups_n    = FALSE,
    point_cex        = 0.175,
    title_line       = 3.3,
    draw_border      = TRUE,
    wex              = 0.88,
    essential_labels = gsub("-", "\uad", c("essential\ngenes", "non-essential\ngenes"), fixed = TRUE)
  )
  par(old_par)
  dev.off()

  return(reps_list)
})

separation_prepooled_mat <- SeparationMetrics(reps_list_list[[1]])
separation_postpooled_mat <- SeparationMetrics(reps_list_list[[2]])



# Draw violin plots showing normalized count data -------------------------

RawCountsViolins <- function(counts_mat,
                             point_cex       = 0.4,
                             use_spacing     = 0.225,
                             wex             = 0.85,
                             upper_bound     = 6000,
                             use_title       = NULL,
                             filter_all_zero = TRUE,
                             ...
                             ) {

  stopifnot(ncol(counts_mat) == 8)

  if (filter_all_zero) {
    are_all_zero <- rowSums(counts_mat == 0) == ncol(counts_mat)
    counts_mat <- counts_mat[!(are_all_zero), ]
    message(paste0(sum(are_all_zero), " genes had zero reads in all samples ",
                   "and were omitted."
                   )
            )
  }

  counts_mat <- counts_mat[, c(1:2, 5:6, 3:4, 7:8)]


  old_mar <- par(mar = c(5.7, 4.3, 4, 2.1))
  x_positions <- BeeViolinPlot(as.list(data.frame(counts_mat)),
                               point_cex     = point_cex,
                               use_spacing   = use_spacing,
                               groups_vec    = c(1, 1, 2, 2, 3, 3, 4, 4),
                               wex           = wex,
                               upper_bound   = upper_bound,
                               ...
                               )

  mtext(text = rep(c("R1", "R2"), times = 4),
        at = x_positions, side = 1, line = 0.7, cex = par("cex")
        )

  segments(x0  = x_positions[c(1, 3, 5, 7)],
           x1  = x_positions[c(2, 4, 6, 8)],
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.85), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(rep(c("pre-pool", "post-pool"), times = 2), side = 1,
        at = c(mean(x_positions[1:2]), mean(x_positions[3:4]), mean(x_positions[5:6]), mean(x_positions[7:8])),
        line = 2.1, padj = 0, cex = par("cex")
        )


  segments(x0  = x_positions[c(1, 5)],
           x1  = x_positions[c(4, 8)],
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 3.4), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(c("T0", "T12"), side = 1,
        at = c(mean(x_positions[1:4]), mean(x_positions[5:8])),
        line = 3.8, padj = 0, cex = par("cex")
        )

  mtext("Normalized count", side = 2, line = 2.75, cex = par("cex"))
  if (!(is.null(use_title))) {
    title(use_title, cex.main = par("cex"))
  }
  par(old_mar)
  return(invisible(NULL))
}



essential_entrezs     <- GetAvailableGenes(essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)
non_essential_entrezs <- GetAvailableGenes(non_essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)

essential_mat     <- CountsMatForGenes(essential_entrezs, allow_switch = FALSE)
non_essential_mat <- CountsMatForGenes(non_essential_entrezs, allow_switch = FALSE)


for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(PDFs_dir, "Violin plots - normalized counts - all timepoints.pdf"),
        width = 6.2, height = 4.5
        )
  }
  for (draw_points in c(TRUE)) {
    RawCountsViolins(essential_mat,
                     violin_colors = brewer.pal(9, "Purples")[[3]],
                     point_colors  = "#7c7198",
                     use_title     = "Essential genes",
                     gap_ratio     = 1.2,
                     )
    RawCountsViolins(non_essential_mat,
                     violin_colors = "#c7e7c0",
                     point_colors  = "#5b8669",
                     use_title     = "Non-essential genes",
                     gap_ratio     = 1.2
    )
  }
  if (create_PDF) {
    dev.off()
  }
}



# Create scatter plots ----------------------------------------------------

Log2FCScatterPlot(baseline_indices = 1:2,
                  intervention_indices = 3:4,
                  allow_switch = FALSE,
                  highlight_NT = TRUE, highlight_essential = FALSE,
                  show_phenotype_score = TRUE
                  )

for (create_PDF in c(FALSE, TRUE)) {
  for (highlight_option in c("none", "essential", "NT")) {
    if (highlight_option == "none") {
      highlight_NT <- FALSE
      highlight_essential <- FALSE
      file_postfix <- "plain"
      use_width <- 4.37
    } else if (highlight_option == "essential") {
      highlight_NT <- FALSE
      highlight_essential <- TRUE
      use_width <- 5.45
      file_postfix <- "essential genes highlighted"
    } else if (highlight_option == "NT") {
      highlight_NT <- TRUE
      highlight_essential <- FALSE
      use_width <- 5.45
      file_postfix <- "NT controls highlighted"
    }

    if (create_PDF) {
      pdf(file.path(PDFs_dir, paste0("Scatter plots - ", file_postfix, ".pdf")),
          width = use_width, height = 4.7
          )
    }
    args_list <- list(allow_switch         = FALSE,
                      highlight_NT         = highlight_NT,
                      highlight_essential  = highlight_essential,
                      embed_PNG            = create_PDF,
                      show_phenotype_score = TRUE
                      )
    do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 1:2, intervention_indices = 3:4, use_title = prepooled_title)))
    do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 5:6, intervention_indices = 7:8, use_title = postpooled_title)))
    if (create_PDF) {
      dev.off()
    }
  }
}



# Draw histograms ---------------------------------------------------------

columns_vec <- c(
  paste0("NoSwitch_xMM_Prepool_T0_R", 1:2),
  paste0("NoSwitch_xMM_Prepool_T12_R", 1:2),
  paste0("NoSwitch_xMM_Postpool_T0_R", 1:2),
  paste0("NoSwitch_xMM_Postpool_T12_R", 1:2)
)
columns_list <- list(
  "All pre-pooled samples"     = columns_vec[1:4],
  "All post-pooled samples"    = columns_vec[5:8],
  "Pre-pooled \u2013 T0"       = columns_vec[1:2],
  "Post-pooled \u2013 T0"      = columns_vec[5:6],
  "Pre-pooled \u2013 T12"      = columns_vec[3:4],
  "Post-pooled \u2013 T12"     = columns_vec[7:8],

  "Pre-pooled \u2013 T0 \u2013 replicate 1"   = columns_vec[[1]],
  "Pre-pooled \u2013 T0 \u2013 replicate 2"   = columns_vec[[2]],
  "Post-pooled \u2013 T0 \u2013 replicate 1"  = columns_vec[[5]],
  "Post-pooled \u2013 T0 \u2013 replicate 2"  = columns_vec[[6]],
  "Pre-pooled \u2013 T12 \u2013 replicate 1"  = columns_vec[[3]],
  "Pre-pooled \u2013 T12 \u2013 replicate 2"  = columns_vec[[4]],
  "Post-pooled \u2013 T12 \u2013 replicate 1" = columns_vec[[7]],
  "Post-pooled \u2013 T12 \u2013 replicate 2" = columns_vec[[8]]
)

DrawHistogram(rowMeans(counts_df[, columns_vec]),
              truncation_limit = 15000,
              num_breaks = 150,
              x_axis_label = "Mean read count",
              y_axis_label = "Number of plasmids"
              )

for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(PDFs_dir, paste0("Histograms - read counts.pdf")),
        width = 5, height = 4
        )
  }
  old_mar <- par(mar = c(4.25, 4.1, 3.25, 2.1))
  for (use_title in names(columns_list)) {
    use_columns <- columns_list[[use_title]]
    DrawHistogram(rowMeans(as.matrix(counts_df[, use_columns, drop = FALSE])),
                  truncation_limit   = 140L,
                  use_exact_breaks   = TRUE,
                  num_breaks         = 71,
                  breaks_range       = c(0, 140),
                  title_text         = use_title,
                  x_axis_label       = if (length(use_columns) == 1) "Read count" else "Mean read count",
                  y_axis_label       = "Number of plasmids",
                  y_axis_limits      = c(-50, 3500),
                  x_axis_upper_limit = 140,
                  x_axis_space       = 1 / 100
                  )
  }
  par(old_mar)
  if (create_PDF) {
    dev.off()
  }
}




# Try stuff ---------------------------------------------------------------

BeeViolinPlot(lapply(counts_df[, columns_vec], function(x) log10(x + 1)), brewer_pals = c(rep("Blues", 4), rep("Purples", 4)),
              sina_plot = TRUE
              )



# Display individual genes ------------------------------------------------

PlotCountsForPlasmid <- function(gene_symbol, ...) {
  stopifnot("counts_df" %in% ls(envir = globalenv()))
  plasmid_index <- match(gene_symbol, counts_df[, "Gene_symbol"])
  if (is.na(plasmid_index)) {
    stop("A plasmid with the ID '", gene_symbol, "' was not found!")
  }
  counts_mat <- GetCountsMat(counts_df, ...)
  group_names <- c("Pre-pooled \u2013 T0", "Pre-pooled \u2013 T12",
                   "Post-pooled \u2013 T0", "Post-pooled \u2013 T12"
                   )
  groups_fac <- factor(rep(group_names, each = 2), levels = group_names)
  counts_vec <- counts_mat[plasmid_index, ]
  PlotCounts(counts_vec, groups_fac, group_colors = c("Greys", "Blues", "Greys", "Purples"))
  title(gene_symbol, font.main = 4)
}

PlotCountsForPlasmid("HNRNPK", allow_switch = FALSE)




# Create QC plots ---------------------------------------------------------


## Prepare count data
are_obsolete <- CRISPRoff_df[, "Is_obsolete"] %in% "Yes"
prepooled_counts_mat <- GetCountsMat(counts_df[!(are_obsolete), ],
                                     allow_switch = FALSE,
                                     allow_1MM    = TRUE,
                                     normalization_columns = 1:4,
                                     normalize = TRUE
                                     )[, 1:4]
prepooled_raw_counts_mat <- GetCountsMat(counts_df[!(are_obsolete), ],
                                         allow_switch = FALSE,
                                         allow_1MM    = TRUE,
                                         normalization_columns = 1:4,
                                         normalize = FALSE
                                         )[, 1:4]

postpooled_counts_mat <- GetCountsMat(counts_df[!(are_obsolete), ],
                                      allow_switch = FALSE,
                                      allow_1MM    = TRUE,
                                      normalization_columns = 5:8,
                                      normalize = TRUE
                                      )[, 5:8]
postpooled_raw_counts_mat <- GetCountsMat(counts_df[!(are_obsolete), ],
                                          allow_switch = FALSE,
                                          allow_1MM    = TRUE,
                                          normalization_columns = 5:8,
                                          normalize = FALSE
                                          )[, 5:8]



for (prepooled in c(TRUE, FALSE)) {
  file_name <- "QC plots - "
  if (prepooled) {
    counts_mat <- prepooled_counts_mat
    raw_counts_mat <- prepooled_raw_counts_mat
    file_name <- paste0(file_name, "pre-pooled")
    four_indices <- c(NA, NA, 1:4)
  } else {
    counts_mat <- postpooled_counts_mat
    raw_counts_mat <- postpooled_raw_counts_mat
    file_name <- paste0(file_name, "post-pooled")
    four_indices <- c(NA, NA, 5:8)
  }

  for (make_PDF in c(FALSE, TRUE)) {

    if (make_PDF) {
      pdf(file.path(PDFs_dir, paste0(file_name, ".pdf")),
          width = use_width, height = 4.7
          )
    }

    ## Display count-level data
    RawCountsHistogram(raw_counts_mat,
                       y_axis_upper_limit = 4000,
                       fixed_y_upper_limit = TRUE,
                       title_text = "Plasmid count at baseline (both replicates)",
                       x_axis_upper_limit = 4.3
                       )
    RawCountsHistogram(raw_counts_mat,
                       y_axis_upper_limit = 2000,
                       fixed_y_upper_limit = TRUE,
                       show_replicates = TRUE,
                       semitransparent_lines = make_PDF,
                       x_axis_upper_limit = 4.3
                       )
    CountBoxPlot(counts_mat[, c(NA, NA, 1:4)], include_timepoints = 2:3, embed_PNG = TRUE,
                 upper_bound = 4.3, y_limits = c(-0.01, 4.3)
                 )
    gini_indices <- CountBarPlot(raw_counts_mat[, c(NA, NA, 1:4)], gini_index = TRUE, include_timepoints = 2:3,
                                 lollipop = TRUE, y_limit = 0.5
                                 )
    num_missing <- CountBarPlot(raw_counts_mat[, c(NA, NA, 1:4)], include_timepoints = 2:3,
                                y_limit = 1000
                                )
    par(old_par)

    GammaBoxPlot(counts_df,
                 embed_PNG            = TRUE,
                 both_timepoints      = FALSE,
                 baseline_indices     = if (prepooled) 1:2 else 5:6,
                 intervention_indices = if (prepooled) 3:4 else 7:8
                 )

    Log2FCScatterPlot(baseline_indices     = if (prepooled) 1:2 else 5:6,
                      intervention_indices = if (prepooled) 3:4 else 7:8,
                      allow_switch         = FALSE,
                      highlight_NT         = TRUE,
                      highlight_essential  = FALSE,
                      show_phenotype_score = TRUE,
                      use_title            = "Replicate scatter plot",
                      title_font           = 1,
                      use_mar              = c(4, 4, 3.5, 7.5),
                      embed_PNG            = TRUE
                      )

    if (make_PDF) {
      dev.off()
    }
  }
  if (prepooled) {
    gini_indices_prepooled <- gini_indices
    num_missing_prepooled <- num_missing
  } else {
    gini_indices_postpooled <- gini_indices
    num_missing_postpooled <- num_missing
  }
}



# Compile log2FC data -----------------------------------------------------

sg_CRISPRoff_df <- CRISPRoff_df
logfc_prepooled_df  <- StandardLog2FCDf(counts_df[, !(grepl("Postpool_", names(counts_df), fixed = TRUE))],
                                       intervention_indices = 3:4
                                       )
logfc_postpooled_df <- StandardLog2FCDf(counts_df[, !(grepl("Prepool_", names(counts_df), fixed = TRUE))],
                                       intervention_indices = 3:4
                                       )
ROC_prepooled_df <- both_reps_ROC_df_list[["ROC_prepool"]]
ROC_postpooled_df <- both_reps_ROC_df_list[["ROC_postpool"]]

AUC_prepooled_vec <- c(
  "both_reps" = AUC_prepooled_both_reps,
  "rep1"      = AUC_prepooled_rep1,
  "rep2"      = AUC_prepooled_rep2
)
AUC_postpooled_vec <- c(
  "both_reps" = AUC_postpooled_both_reps,
  "rep1"      = AUC_postpooled_rep1,
  "rep2"      = AUC_postpooled_rep2
)



# Examine bidirectional promoters -----------------------------------------

BidirectionalViolins(bidirectional_df, logfc_prepooled_df, max_distance = 10000)
BidirectionalViolins(bidirectional_df, logfc_prepooled_df, max_distance = 15000)
BidirectionalViolins(bidirectional_df, logfc_prepooled_df, max_distance = 20000,
                     num_controls = 30L
                     )

BidirectionalViolins(bidirectional_df, logfc_postpooled_df, max_distance = 10000)
BidirectionalViolins(bidirectional_df, logfc_postpooled_df, max_distance = 15000)
BidirectionalViolins(bidirectional_df, logfc_postpooled_df, max_distance = 20000,
                     num_controls = 30L
                     )



# Export data -------------------------------------------------------------

stopifnot(identical(counts_df[, "Plasmid_ID"], logfc_prepooled_df[, "Plasmid_ID"]))
stopifnot(identical(counts_df[, "Plasmid_ID"], logfc_postpooled_df[, "Plasmid_ID"]))

prepooled_count_columns <- c(
  "NoSwitch_xMM_Prepool_T0_R1", "NoSwitch_xMM_Prepool_T0_R2",
  "NoSwitch_xMM_Prepool_T12_R2", "NoSwitch_xMM_Prepool_T12_R2"
)
prepooled_counts_df <- counts_df[, prepooled_count_columns]
export_new_column_names <- c(
  "Count_baseline_rep1", "Count_baseline_rep2",
  "Count_endpoint_rep1", "Count_endpoint_rep2"
)
names(prepooled_counts_df) <- export_new_column_names


postpooled_count_columns <- sub("_Prepool_", "_Postpool_", prepooled_count_columns, fixed = TRUE)
postpooled_counts_df <- counts_df[, postpooled_count_columns]
names(postpooled_counts_df) <- export_new_column_names


tidy_tgonfio_df <- PrepareTgonfioForExport(CRISPRoff_df)
tables_dir <- file.path(project_dir, "05_output", "Tables")

ExportResultsDf(tidy_tgonfio_df, logfc_prepooled_df, prepooled_counts_df,
                file_path = file.path(tables_dir, "nanopore_prepooled_counts.csv")
                )
ExportResultsDf(tidy_tgonfio_df, logfc_postpooled_df, postpooled_counts_df,
                file_path = file.path(tables_dir, "nanopore_postpooled_counts.csv")
                )



# Save data ---------------------------------------------------------------

save(list = c("logfc_prepooled_df", "logfc_postpooled_df",
              "gini_indices_prepooled", "gini_indices_postpooled",
              "num_missing_prepooled", "num_missing_postpooled",
              "separation_prepooled_mat", "separation_postpooled_mat",
              "ROC_prepooled_df", "ROC_postpooled_df",
              "AUC_prepooled_vec", "AUC_postpooled_vec"
              ),
     file = file.path(rdata_dir, "09_create_figures_from_count_data.RData")
     )

