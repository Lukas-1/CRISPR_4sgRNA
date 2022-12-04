### 2022-06-23


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_nanopore_dir       <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

project_dir     <- file.path(experiments_directory, "2022-06-21 - Illumina paired-end 2sg - correct reference")
rdata_dir       <- file.path(project_dir, "03_R_objects")
figures_dir     <- file.path(project_dir, "04_output_data", "Figures")
PDFs_dir        <- file.path(figures_dir, "PDFs")
with_switch_dir <- file.path(PDFs_dir, "Including template switch")
no_switch_dir   <- file.path(PDFs_dir, "Excluding template switch")
first_rdata_dir <- file.path(first_illumina_trial_dir, "03_R_objects")
manuscript_dir  <- file.path(PDFs_dir, "Manuscript")
thesis_dir      <- file.path(PDFs_dir, "Thesis")



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))
load(file.path(rdata_dir, "03_disambiguate_dJR072_CRISPRoff_library.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__num_mapped_reads.RData"))
load(file.path(rdata_dir, "11_compute_read_quality_metrics.RData"))

load(file.path(CRISPR_root_directory, "3) RData files", "1) General",
               "24) Enumerate pairs of genes at bidirectional promoters.RData"
               ))



# Prepare R objects -------------------------------------------------------

sg_mat <- as.matrix(CRISPRoff_df[, c("protospacer_A", "protospacer_B")])
library_densities <- GetLibraryDensities(sg_mat)

are_missing <- rowSums(counts_df[, c("NoSwitch_0MM_Tbefore_R1", "NoSwitch_0MM_Tbefore_R2")]) == 0
are_present <- rowSums(counts_df[, c("NoSwitch_0MM_Tbefore_R1", "NoSwitch_0MM_Tbefore_R1")]) != 0
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

for (create_PDF in c(FALSE, TRUE)) {

  if (create_PDF) {
    pdf(file.path(PDFs_dir, "Essentiality of missing genes.pdf"),
        width = 4.5, height = 4.5
        )
  }

  old_mar <- par(mar = c(4, 4.3, 4, 2.1))
  BeeViolinPlot(probs_list, use_swarm = FALSE, point_cex = 0.35, wex = 0.9,
                y_limits = c(0, 1)
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

T0vT12_title <- expression(bold("CRISPRoff:" ~ bolditalic("T0")      ~ "vs." ~ bolditalic("T12")))
BvT12_title  <- expression(bold("CRISPRoff:" ~ bolditalic("Tbefore") ~ "vs." ~ bolditalic("T12")))
BvT0_title   <- expression(bold("CRISPRoff:" ~ bolditalic("Tbefore") ~ "vs." ~ bolditalic("T0")))

both_reps_ROC_df_list_list <- MakeROCDfListList()
rep1_ROC_df_list_list <- MakeROCDfListList(choose_rep = 1)
rep2_ROC_df_list_list <- MakeROCDfListList(choose_rep = 2)

both_reps_more_genes_ROC_df_list_list <- MakeROCDfListList(use_blomen_hart = FALSE)
rep1_ROC_more_genes_df_list_list <- MakeROCDfListList(choose_rep = 1, use_blomen_hart = FALSE)
rep2_ROC_more_genes_df_list_list <- MakeROCDfListList(choose_rep = 2, use_blomen_hart = FALSE)

for (allow_switch in c(FALSE, TRUE)) {

  if (allow_switch) {
    use_dir <- file.path(PDFs_dir, "Including template switch")
  } else {
    use_dir <- file.path(PDFs_dir, "Excluding template switch")
  }

  old_oma <- par(oma = rep(0.3, 4))
  for (create_PDF in c(FALSE, TRUE)) {
    if (create_PDF) {
      pdf(file.path(use_dir, "ROC curves - both replicates.pdf"),
          width = 3 + 1.14, height = 3 + 1.56
          )
    }
    PlotEssentialROCDf(both_reps_ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = T0vT12_title)
    PlotEssentialROCDf(both_reps_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = BvT12_title)
    PlotEssentialROCDf(both_reps_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = BvT0_title)
    PlotEssentialROCDf(both_reps_more_genes_ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = T0vT12_title)
    PlotEssentialROCDf(both_reps_more_genes_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = BvT12_title)
    PlotEssentialROCDf(both_reps_more_genes_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = BvT0_title)
    if (create_PDF) {
      dev.off()
      pdf(file.path(use_dir, "ROC curves - replicate 1.pdf"),
          width = 3 + 1.14, height = 3 + 1.56
          )
    }
    PlotEssentialROCDf(rep1_ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = expression(bold("T0 vs. T12 (replicate 1)")))
    PlotEssentialROCDf(rep1_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = expression(bold("Tbefore vs. T12 (replicate 1)")))
    PlotEssentialROCDf(rep1_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = expression(bold("Tbefore vs. T0 (replicate 1)")))
    PlotEssentialROCDf(rep1_ROC_more_genes_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = expression(bold("T0 vs. T12 (replicate 1)")))
    PlotEssentialROCDf(rep1_ROC_more_genes_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = expression(bold("Tbefore vs. T12 (replicate 1)")))
    PlotEssentialROCDf(rep1_ROC_more_genes_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = expression(bold("Tbefore vs. T0 (replicate 1)")))
    if (create_PDF) {
      dev.off()
      pdf(file.path(use_dir, "ROC curves - replicate 2.pdf"),
          width = 3 + 1.14, height = 3 + 1.56
          )
    }
    PlotEssentialROCDf(rep2_ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = expression(bold("T0 vs. T12 (replicate 2)")))
    PlotEssentialROCDf(rep2_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = expression(bold("Tbefore vs. T12 (replicate 2)")))
    PlotEssentialROCDf(rep2_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = expression(bold("Tbefore vs. T0 (replicate 2)")))
    PlotEssentialROCDf(rep2_ROC_more_genes_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = expression(bold("T0 vs. T12 (replicate 2)")))
    PlotEssentialROCDf(rep2_ROC_more_genes_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = expression(bold("Tbefore vs. T12 (replicate 2)")))
    PlotEssentialROCDf(rep2_ROC_more_genes_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = expression(bold("Tbefore vs. T0 (replicate 2)")))
    if (create_PDF) {
      dev.off()
    }
  }
  par(old_oma)
}


use_df <- both_reps_ROC_df_list_list[[1]][["ROC_BvT12"]]
pdf(file.path(manuscript_dir, "2sg ROC curve BvT12.pdf"),
    width = 2, height = 2
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
PlotROCDf(use_df, flip = TRUE, xlab_line = 1.6, ylab_line = 2.1, ROC_lwd = 1.5)
title("CRISPRoff library", cex.main = 1, font.main = 1, line = 0.7)
par(old_par)
dev.off()



# Draw violin plots showing the log2FC ------------------------------------

y_span <- 0.7 * 0.02
custom_y_limits <- c(-0.5 - y_span, 0.2)
custom_args <- list(lower_bound = -0.5, upper_bound = 0.2, y_limits = custom_y_limits)

for (allow_switch in c(FALSE, TRUE)) {

  if (allow_switch) {
    use_dir <- file.path(PDFs_dir, "Including template switch")
  } else {
    use_dir <- file.path(PDFs_dir, "Excluding template switch")
  }

  for (create_PDF in c(FALSE, TRUE)) {
    if (create_PDF) {
      pdf(file.path(use_dir, "Violin plots - mean of replicates.pdf"),
          width = 4.5, height = 4.5
          )
    }
    for (draw_points in c(TRUE)) {
      do.call(ViolinPlotEssentialDf, c(list(both_reps_ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12"]], use_title = T0vT12_title, draw_points = draw_points), custom_args))
      do.call(ViolinPlotEssentialDf, c(list(both_reps_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12"]],  use_title = BvT12_title,  draw_points = draw_points), custom_args))
      do.call(ViolinPlotEssentialDf, c(list(both_reps_ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0"]],   use_title = BvT0_title,   draw_points = draw_points), custom_args))
    }
    if (create_PDF) {
      dev.off()
    }
  }
}



# Draw violin plots showing the log2FC for both replicates ----------------

RepEssentialViolins(3:4, 5:6, use_title = T0vT12_title, lower_bound = -0.6, upper_bound = 0.4, min_count_at_baseline = 0)
RepEssentialViolins(1:2, 5:6, use_title = BvT12_title,  lower_bound = -0.6, upper_bound = 0.4, min_count_at_baseline = 0)
RepEssentialViolins(1:2, 3:4, use_title = BvT0_title,   lower_bound = -0.6, upper_bound = 0.4, min_count_at_baseline = 0)


RepEssentialViolins(1:2, 5:6, use_title = BvT12_title,  lower_bound = -0.6, upper_bound = 0.4, min_count_at_baseline = 0,
                    use_blomen_hart = FALSE
                    )


custom_y_limits <- c(-0.6 - (0.86 * 0.02), 0.25)
for (allow_switch in c(FALSE, TRUE)) {
  for (create_PDF in c(FALSE, TRUE)) {
    if (allow_switch) {
      use_dir <- file.path(PDFs_dir, "Including template switch")
    } else {
      use_dir <- file.path(PDFs_dir, "Excluding template switch")
    }
    if (create_PDF) {
      pdf(file.path(use_dir, "Violin plots - R1 and R2.pdf"), width = 5, height = 4.5)
    }
    for (draw_points in c(TRUE)) {
      for (use_blomen_hart in c(TRUE, FALSE)) {
        RepEssentialViolins(3:4, 5:6, use_title = T0vT12_title,draw_points = draw_points,
                            lower_bound = -0.6, upper_bound = 0.25, y_limits = custom_y_limits,
                            allow_switch = allow_switch, use_blomen_hart = use_blomen_hart
                            )
        RepEssentialViolins(1:2, 5:6, use_title = BvT12_title, draw_points = draw_points,
                            lower_bound = -0.6, upper_bound = 0.25, y_limits = custom_y_limits,
                            allow_switch = allow_switch, use_blomen_hart = use_blomen_hart
                            )
        RepEssentialViolins(1:2, 3:4, use_title = BvT0_title, draw_points = draw_points,
                            lower_bound = -0.6, upper_bound = 0.25, y_limits = custom_y_limits,
                            allow_switch = allow_switch, use_blomen_hart = use_blomen_hart
                            )
      }
    }
    if (create_PDF) {
      dev.off()
    }
  }
}



pdf(file.path(manuscript_dir, "2sg violin plots BvT12.pdf"),
    width = 2, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.8, lheight = 0.9)

reps_list <- RepEssentialViolins(
  1:2, 5:6,
  use_title        = expression(bold("CRISPRoff library")),
  lower_bound      = -0.6,
  upper_bound      = 0.25,
  y_limits         = custom_y_limits,
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
  wex              = 0.88
)
par(old_par)
dev.off()

separation_mat <- SeparationMetrics(reps_list)





# Draw violin plots showing all samples -----------------------------------

essential_entrezs     <- GetAvailableGenes(essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)
non_essential_entrezs <- GetAvailableGenes(non_essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)

# rdata_dir_4sg <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing", "03_R_objects")
# load(file.path(rdata_dir_4sg, "10_identify_genes_present_in_both_libraries.RData"))
# essential_entrezs     <- intersect_essential_entrezs
# non_essential_entrezs <- intersect_non_essential_entrezs

for (allow_switch in c(FALSE, TRUE)) {

  essential_allsamples_mat     <- AllSamplesLog2FC(essential_entrezs, normalize_to_reps = 1, allow_switch = allow_switch)
  non_essential_allsamples_mat <- AllSamplesLog2FC(non_essential_entrezs, normalize_to_reps = 1, allow_switch = allow_switch)

  for (create_PDF in c(FALSE, TRUE)) {
    if (allow_switch) {
      use_dir <- file.path(PDFs_dir, "Including template switch")
    } else {
      use_dir <- file.path(PDFs_dir, "Excluding template switch")
    }
    if (create_PDF) {
      pdf(file.path(use_dir, "Violin plots - log2FC - all timepoints.pdf"),
          width = 5.25, height = 4.5
          )
    }
    for (draw_points in c(TRUE)) {
      AllTimesLogFCViolins(essential_allsamples_mat[, -1],
                           violin_colors = brewer.pal(9, "Purples")[[3]],
                           point_colors  = "#7c7198",
                           use_title     = "Essential genes",
                           draw_points   = draw_points,
                           gap_ratio     = 1.2
                           )
      AllTimesLogFCViolins(non_essential_allsamples_mat[, -1],
                           violin_colors = "#c7e7c0",
                           point_colors  = "#5b8669",
                           use_title     = "Non-essential genes",
                           draw_points   = draw_points,
                           gap_ratio     = 1.2
                           )
    }
    if (create_PDF) {
      dev.off()
    }
  }
}



# Draw violin plots showing normalized count data -------------------------

for (allow_switch in c(FALSE, TRUE)) {

  essential_mat     <- CountsMatForGenes(essential_entrezs, allow_switch = allow_switch)
  non_essential_mat <- CountsMatForGenes(non_essential_entrezs, allow_switch = allow_switch)

  for (create_PDF in c(FALSE, TRUE)) {
    if (allow_switch) {
      use_dir <- file.path(PDFs_dir, "Including template switch")
    } else {
      use_dir <- file.path(PDFs_dir, "Excluding template switch")
    }
    if (create_PDF) {
      pdf(file.path(use_dir, "Violin plots - normalized counts - all timepoints.pdf"),
          width = 6.2, height = 4.5
          )
    }
    for (draw_points in c(TRUE)) {
      RawCountsViolins(essential_mat,
                       violin_colors = brewer.pal(9, "Purples")[[3]],
                       point_colors  = "#7c7198",
                       use_title     = "Essential genes",
                       draw_points   = draw_points,
                       gap_ratio     = 1.2,
                       )
      RawCountsViolins(non_essential_mat,
                       violin_colors = "#c7e7c0",
                       point_colors  = "#5b8669",
                       use_title     = "Non-essential genes",
                       draw_points   = draw_points,
                       gap_ratio     = 1.2
                       )
    }
    if (create_PDF) {
      dev.off()
    }
  }
}



# Create scatter plots ----------------------------------------------------

Log2FCScatterPlot(baseline_indices = 1:2, intervention_indices = 5:6,
                  allow_switch = FALSE,
                  highlight_NT = TRUE, highlight_essential = FALSE,
                  show_phenotype_score = TRUE
                  )

for (allow_switch in c(FALSE, TRUE)) {
  if (allow_switch) {
    use_dir <- file.path(PDFs_dir, "Including template switch")
  } else {
    use_dir <- file.path(PDFs_dir, "Excluding template switch")
  }
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
        pdf(file.path(use_dir, paste0("Scatter plots - ", file_postfix, ".pdf")),
            width = use_width, height = 4.7
            )
      }
      args_list <- list(allow_switch        = allow_switch,
                        highlight_NT        = highlight_NT,
                        highlight_essential = highlight_essential,
                        embed_PNG           = create_PDF
                        )
      do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 3:4, intervention_indices = 5:6, use_title = T0vT12_title, show_phenotype_score = TRUE)))
      do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 1:2, intervention_indices = 5:6, use_title = BvT12_title, show_phenotype_score = TRUE)))
      do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 1:2, intervention_indices = 3:4, use_title = BvT0_title, show_phenotype_score = TRUE)))
      if (create_PDF) {
        dev.off()
      }
    }
  }
}



# Draw histograms ---------------------------------------------------------

columns_vec <- c(
  paste0("NoSwitch_xMM_Tbefore_R", 1:2),
  paste0("NoSwitch_xMM_T0_R", 1:2),
  paste0("NoSwitch_xMM_T12_R", 1:2)
)
columns_list <- list(
  "All samples and timepoints" = columns_vec,
  "Tbefore (both replicates)"  = columns_vec[1:2],
  "T0 (both replicates)"       = columns_vec[3:4],
  "T12 (both replicates)"      = columns_vec[5:6],
  "Tbefore replicate 1"        = columns_vec[[1]],
  "Tbefore replicate 2"        = columns_vec[[2]],
  "T0 replicate 1"             = columns_vec[[3]],
  "T0 replicate 2"             = columns_vec[[4]],
  "T12 replicate 1"            = columns_vec[[5]],
  "T12 replicate 2"            = columns_vec[[6]]
)

DrawHistogram(counts_df[, "Sum_MaySwitch_xMM"] / 6,
              truncation_limit = 5000,
              num_breaks = 150,
              x_axis_label = "Mean read count",
              y_axis_label = "Number of plasmids"
              )

for (allow_switch in c(FALSE, TRUE)) {
  if (allow_switch) {
    use_dir <- file.path(PDFs_dir, "Including template switch")
  } else {
    use_dir <- file.path(PDFs_dir, "Excluding template switch")
  }
  for (create_PDF in c(FALSE, TRUE)) {
    if (create_PDF) {
      pdf(file.path(use_dir, paste0("Histograms - read counts.pdf")),
          width = 5, height = 4
          )
    }
    old_mar <- par(mar = c(4.25, 4.1, 3.25, 2.1))
    for (use_title in names(columns_list)) {
      use_columns <- columns_list[[use_title]]
      if (allow_switch) {
        use_columns <- sub("^NoSwitch_", "MaySwitch_", use_columns)
      }
      DrawHistogram(rowMeans(as.matrix(counts_df[, use_columns, drop = FALSE])),
                    truncation_limit = 5000,
                    num_breaks       = 150L,
                    title_text       = use_title,
                    x_axis_label     = if (length(use_columns) == 1) "Read count" else "Mean read count",
                    y_axis_label     = "Number of plasmids",
                    y_axis_limits    = c(-20, 1600),
                    x_axis_space     = 1 / 100
                    )
    }
    par(old_mar)
    if (create_PDF) {
      dev.off()
    }
  }
}



# Create bar plots --------------------------------------------------------

use_min_count <- 20L
four_metrics_list <- lapply(columns_list, function(noswitch_columns) {
  mayswitch_columns <- sub("^NoSwitch_", "MaySwitch_", noswitch_columns)
  c("mayswitch_zero_reads" = sum(rowMeans(as.matrix(counts_df[, mayswitch_columns])) == 0),
    "noswitch_zero_reads"  = sum(rowMeans(as.matrix(counts_df[, noswitch_columns])) == 0),
    "mayswitch_below_min"  = sum(rowMeans(as.matrix(counts_df[, mayswitch_columns])) < use_min_count),
    "noswitch_below_min"   = sum(rowMeans(as.matrix(counts_df[, noswitch_columns])) < use_min_count)
    )
})

FourBars(four_metrics_list[[1]], library_size = nrow(counts_df), title_text = "All samples")

for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(PDFs_dir, paste0("Bar charts - missing plasmids.pdf")),
        width = 4, height = 4.5
        )
  }
  for (use_title in names(columns_list)) {
    FourBars(four_metrics_list[[use_title]],
             use_y_limits = c(0, 300),
             title_text = sub(" (both replicates)", " (mean of both replicates)", use_title, fixed = TRUE),
             library_size = nrow(counts_df)
             )
  }
  par(old_mar)
  if (create_PDF) {
    dev.off()
  }
}



# Display individual genes ------------------------------------------------

PlotCountsForPlasmid("HNRNPK")



# Create QC plots ---------------------------------------------------------

## Prepare read-level statistics
num_reads_detailed_mat <- rbind(
  "Unmapped_both" = samples_df[, "Num_reads"] - num_mapped_df[, "Num_either_read_mapped"],
  "Unmapped_sg1"  = num_mapped_df[, "Num_unmapped_read1_only"],
  "Unmapped_sg2"  = num_mapped_df[, "Num_unmapped_read2_only"],
  "Mapped"        = num_mapped_df[, "Num_both_reads_mapped"]
)
percent_switch_vec <- num_mapped_df[, "Num_template_switch"] / num_mapped_df[, "Num_both_reads_mapped"]
percent_1MM_mat <- rbind(
  "Tolerates1MM_both" = num_mapped_df[, "Num_1MM_both_reads"],
  "Tolerates1MM_sg1"  = num_mapped_df[, "Num_1MM_read1_only"],
  "Tolerates1MM_sg2"  = num_mapped_df[, "Num_1MM_read2_only"],
  "Without_mismatch"  = num_mapped_df[, "Num_both_reads_mapped"] -
                        (num_mapped_df[, "Num_1MM_both_reads"] +
                         num_mapped_df[, "Num_1MM_read1_only"] +
                         num_mapped_df[, "Num_1MM_read2_only"]
                         )
)
percent_1MM_mat <- prop.table(percent_1MM_mat, margin = 2)



## Prepare count data
counts_mat <- GetCountsMat(counts_df,
                           allow_switch = FALSE,
                           allow_1MM    = TRUE,
                           normalization_columns = c(1:2, 5:6),
                           normalize = TRUE
                           )
raw_counts_mat <- GetCountsMat(counts_df,
                               allow_switch = FALSE,
                               allow_1MM    = TRUE,
                               normalization_columns = c(1:2, 5:6),
                               normalize = FALSE
                               )


for (all_timepoints in c(FALSE, TRUE)) {

  file_name <- "QC plots - "
  if (all_timepoints) {
    use_timepoints <- 1:3
    file_name <- paste0(file_name, "all timepoints")
  } else {
    use_timepoints <- c(1, 3)
    file_name <- paste0(file_name, "without T0")
  }

  for (make_PDF in c(FALSE, TRUE)) {

    if ((!(make_PDF) && (!(all_timepoints)))) {
      next
    }

    if (make_PDF) {
      pdf(file.path(no_switch_dir, paste0(file_name, ".pdf")),
          width = use_width, height = 4.7
          )
    }

    ## Display read-level data
    TwoDensities(show_GC = TRUE, semitransparent_lines = make_PDF, include_timepoints = use_timepoints)
    PerBaseQuality(base_qual_mat, semitransparent_lines = make_PDF, include_timepoints = use_timepoints)
    TwoDensities(show_GC = FALSE, semitransparent_lines = make_PDF, include_timepoints = use_timepoints)
    MappedReadsBarPlot(num_reads_detailed_mat, include_timepoints = use_timepoints)
    MappedReadsBarPlot(percent_1MM_mat, include_timepoints = use_timepoints, show_percentage = TRUE)
    old_mar <- par(mar = c(4, 4, 3.75, 2.1))
    PercentageBarPlot(percent_switch_vec, include_timepoints = use_timepoints)

    ## Display count-level data
    RawCountsHistogram(raw_counts_mat,
                       y_axis_upper_limit = 3000,
                       fixed_y_upper_limit = TRUE
                       )
    RawCountsHistogram(raw_counts_mat,
                       y_axis_upper_limit = 1600,
                       fixed_y_upper_limit = TRUE,
                       show_replicates = TRUE,
                       semitransparent_lines = make_PDF
                       )
    CountBoxPlot(counts_mat, include_timepoints = use_timepoints, embed_PNG = TRUE)
    gini_indices <- CountBarPlot(raw_counts_mat, gini_index = TRUE, include_timepoints = use_timepoints)
    num_missing <- CountBarPlot(raw_counts_mat, include_timepoints = use_timepoints)
    par(old_par)

    GammaBoxPlot(counts_df, embed_PNG = TRUE)

    Log2FCScatterPlot(baseline_indices     = 1:2,
                      intervention_indices = 5:6,
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
}



# Export QC plots for the manuscript --------------------------------------

pdf(file.path(manuscript_dir, "2sg - Fig. S7A - count histograms.pdf"),
    width = 2.1, height = 1.75
    )
ManuscriptRawCountsHistogram(raw_counts_mat, "CRISPRoff library")
dev.off()


use_mai <- c(0.7, 0.8, 0.4, 1.4)
use_cex <- 0.6
base_height <- 1.15
pdf(file.path(manuscript_dir, "2sg - Fig. S7C - scatter plot.pdf"),
    width  = base_height + (sum(use_mai[c(2, 4)] * use_cex)),
    height = base_height + (sum(use_mai[c(1, 3)]) * use_cex)
    )
old_par <- par(cex = 0.6, lwd = 0.8)
Log2FCScatterPlot(baseline_indices     = 1:2,
                  intervention_indices = 5:6,
                  allow_switch         = FALSE,
                  highlight_NT         = TRUE,
                  highlight_essential  = FALSE,
                  show_phenotype_score = TRUE,
                  use_title            = "CRISPRoff library",
                  title_font           = 1,
                  use_mar              = use_mai * 5,
                  embed_PNG            = TRUE,
                  x_axis_label_line    = 1.8,
                  y_axis_label_line    = 2.5,
                  sparse_x_axis_labels = TRUE,
                  use_tcl              = 0.3,
                  x_axis_mgp           = 0.35,
                  y_axis_mgp           = 0.5,
                  point_cex            = 0.45,
                  legend_lines_x_start = 0.65,
                  legend_point_x_start = 0.05,
                  abbreviate_NT        = TRUE
                  )
dev.off()


pdf(file.path(manuscript_dir, "2sg - Fig. S7D - phenotype violin plots.pdf"),
    width = 1.4, height = 1.75
    )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
GammaBoxPlot(counts_df, embed_PNG = TRUE, both_timepoints = FALSE,
             use_title = "CRISPRoff library",
             cloud_alpha = 0.2, cloud_sd = 0.015, point_cex = 0.2,
             use_lwd = 0.6,
             zero_lty = "solid", zero_lwd = 0.6, zero_color = "gray70",
             y_label_line = 2.1,
             png_res = 1200, wex = 0.85, side_gap = 0.525
             )
dev.off()




# Export read-level QC plots for the thesis -------------------------------

devEMF::emf(file.path(thesis_dir, "1A - GC content.emf"),
            width  = 2.15,
            height = 2.27, emfPlus = FALSE
            )
par("cex" = 0.7, "lwd" = 0.8)
TwoDensities(show_GC = TRUE, include_timepoints = c(1, 3),
             semitransparent_lines = TRUE, use_title = "CRISPRoff library",
             y_axis_label_line = 0.4, legend_x_lines = 0.7, legend_y_lines = 0.8,
             label_read_on_y_axis = FALSE, embed_PNG = TRUE, grid_lwd = 0.75,
             show_y_axis_label = TRUE, broad_margins = TRUE,
             title_y_pos = 0.4, x_axis_label_line = 1.8, x_axis_mgp = 0.39
             )
dev.off()



devEMF::emf(file.path(thesis_dir, "1B - Per-base quality.emf"),
            width  = 2.15,
            height = 2.27, emfPlus = FALSE
            )
par("cex" = 0.7, "lwd" = 0.8)
PerBaseQuality(base_qual_mat, include_timepoints = c(1, 3),
               semitransparent_lines = FALSE, use_title = "CRISPRoff library",
               y_axis_label_line = 1.75, x_axis_label_line = 1.5,
               broad_margins = TRUE, x_axis_tcl = 0.35,
               title_y_pos = 0.4, separate_x_labels = TRUE, x_axis_mgp = 0.39,
               embed_PNG = TRUE, small_middle_gap = TRUE, omit_zero_label = TRUE,
               show_legend = FALSE, y_axis_mgp = 0.525
               )
dev.off()



devEMF::emf(file.path(thesis_dir, "1C - Mean sequence quality.emf"),
            width  = 2.15,
            height = 2.27, emfPlus = FALSE
            )
par("cex" = 0.7, "lwd" = 0.8)
TwoDensities(show_GC = FALSE, include_timepoints = c(1, 3),
             semitransparent_lines = TRUE, use_title = "CRISPRoff library",
             y_axis_label_line = 0.4, legend_x_lines = 0.7, legend_y_lines = 0.8,
             label_read_on_y_axis = FALSE, embed_PNG = TRUE, grid_lwd = 0.75,
             darker_box = TRUE, broad_margins = TRUE,
             show_y_axis_label = TRUE,
             title_y_pos = 0.4, x_axis_label_line = 1.8, x_axis_mgp = 0.39
             )
dev.off()



devEMF::emf(file.path(thesis_dir, "1D - Percentage of mapped reads.emf"),
            width  = 2.4,
            height = 2.27, emfPlus = FALSE
            )
par("cex" = 0.7, "lwd" = 0.8, mai = c(0.412, 0.48375, 0.2724, 0.68))
MappedReadsBarPlot(num_reads_detailed_mat, include_timepoints = c(1, 3),
                   set_mar = FALSE, use_title = NA,
                   y_axis_mgp = 0.525, y_axis_label_line = 1.75,
                   y_axis_tcl = 0.375, show_legend = FALSE,
                   y_upper_limit = 80 * 10^6, bar_width = 0.6, gap_ratio = 1.4,
                   side_gap = 0.6, unit_in_axis = FALSE
                   )
text(x      = grconvertX(0.5, from = "npc", to = "user"),
     y      = par("usr")[[4]] + (diff(grconvertY(c(0, par("mai")[[3]]), from = "inches", to = "user")) * 0.4),
     labels = "CRISPRoff library",
     xpd    = NA
     )
dev.off()



devEMF::emf(file.path(thesis_dir, "1E - Percentage of 1MM reads.emf"),
            width  = 2.4,
            height = 2.27, emfPlus = FALSE
            )
par("cex" = 0.7, "lwd" = 0.8, mai = c(0.412, 0.48375, 0.2724, 0.68))
MappedReadsBarPlot(percent_1MM_mat, include_timepoints = c(1, 3),
                   set_mar = FALSE, use_title = NA,
                   y_axis_mgp = 0.525, y_axis_label_line = 1.75,
                   y_axis_tcl = 0.375, show_legend = FALSE,
                   y_upper_limit = 1, bar_width = 0.6, gap_ratio = 1.4,
                   side_gap = 0.6,
                   show_percentage = TRUE, unit_in_axis = FALSE
                   )
text(x      = grconvertX(0.5, from = "npc", to = "user"),
     y      = par("usr")[[4]] + (diff(grconvertY(c(0, par("mai")[[3]]), from = "inches", to = "user")) * 0.4),
     labels = "CRISPRoff library",
     xpd    = NA
     )
dev.off()



devEMF::emf(file.path(thesis_dir, "1F - Template switch.emf"),
            width  = 2.4,
            height = 2.27, emfPlus = FALSE
            )
par("cex" = 0.7, "lwd" = 0.8, mai = c(0.412, 0.48375, 0.2724, 0.68))
MappedReadsBarPlot(rbind(percent_switch_vec, 1 - percent_switch_vec),
                   include_timepoints = c(1, 3),
                   set_mar = FALSE, use_title = NA,
                   y_axis_mgp = 0.525, y_axis_label_line = 1.75,
                   y_axis_tcl = 0.375, show_legend = FALSE,
                   y_upper_limit = 1, bar_width = 0.6, gap_ratio = 1.4,
                   side_gap = 0.6,
                   show_percentage = TRUE, unit_in_axis = FALSE,
                   use_colors = brewer.pal(9, "Blues")[c(3, 6)]
                   )
text(x      = grconvertX(0.5, from = "npc", to = "user"),
     y      = par("usr")[[4]] + (diff(grconvertY(c(0, par("mai")[[3]]), from = "inches", to = "user")) * 0.4),
     labels = "CRISPRoff library",
     xpd    = NA
     )
dev.off()




# Export count-level QC plots for the thesis ------------------------------

devEMF::emf(file.path(thesis_dir, "2A - Count histograms.emf"),
            width = 2.47, height = 1.75, emfPlus = FALSE
            )
ManuscriptRawCountsHistogram(raw_counts_mat, "CRISPRoff library",
                             use_mai = c(0.6, 0.8, 0.4, 1.4) * 0.6,
                             semitransparent_lines = FALSE,
                             x_axis_upper_limit = 4.2
                             )
dev.off()


use_mai <- c(0.7, 0.8, 0.4, 1.4)
use_cex <- 0.6
base_height <- 1.15
devEMF::emf(file.path(thesis_dir, "2D - Scatter plot.emf"),
            width  = base_height + (sum(use_mai[c(2, 4)] * use_cex)),
            height = base_height + (sum(use_mai[c(1, 3)]) * use_cex),
            emfPlus = FALSE
            )
old_par <- par(cex = 0.6, lwd = 0.8)
Log2FCScatterPlot(baseline_indices     = 1:2,
                  intervention_indices = 5:6,
                  allow_switch         = FALSE,
                  highlight_NT         = TRUE,
                  highlight_essential  = FALSE,
                  show_phenotype_score = TRUE,
                  use_title            = "CRISPRoff library",
                  title_font           = 1,
                  use_mar              = use_mai * 5,
                  embed_PNG            = TRUE,
                  x_axis_label_line    = 1.8,
                  y_axis_label_line    = 2.5,
                  sparse_x_axis_labels = TRUE,
                  use_tcl              = 0.3,
                  x_axis_mgp           = 0.35,
                  y_axis_mgp           = 0.5,
                  point_cex            = 0.45,
                  legend_lines_x_start = 0.65,
                  legend_point_x_start = 0.05,
                  abbreviate_NT        = TRUE
                  )
dev.off()


devEMF::emf(file.path(thesis_dir, "2E - Phenotype violin plots.emf"),
            width = 1.4, height = 1.75, emfPlus = FALSE
            )
old_par <- par(mar = c(3, 4, 2, 1), cex = 0.6, lwd = 0.8)
GammaBoxPlot(counts_df, embed_PNG = TRUE, both_timepoints = FALSE,
             use_title = "CRISPRoff library",
             cloud_alpha = 0.2, cloud_sd = 0.015, point_cex = 0.2,
             use_lwd = 0.6,
             use_swarm = "sina", sina_wex_factor = 0.93,
             violin_colors = brewer.pal(9, "Blues")[[7]], line_colors = brewer.pal(9, "Blues")[[4]],
             zero_lty = "solid", zero_lwd = 0.6, zero_color = "gray70",
             y_label_line = 2.1,
             png_res = 1200, wex = 0.85, side_gap = 0.525, right_gap = 0.45,
             )
dev.off()




# Export main figures for the thesis --------------------------------------

devEMF::emf(file.path(thesis_dir, "3A) Violin plot - ii) CRISPRoff.emf"),
            width = 2, height = 2, emfPlus = FALSE
            )
old_par <- par(cex = 0.6, lwd = 0.8, lheight = 0.9)

reps_list <- RepEssentialViolins(
  1:2, 5:6,
  use_title        = expression(bold("CRISPRoff library")),
  lower_bound      = -0.6,
  upper_bound      = 0.25,
  y_limits         = custom_y_limits,
  allow_switch     = FALSE,
  use_mar          = c(3, 4, 4, 1),
  y_axis_label     = NA,
  y_label_line     = 2.1,
  rep_label_line   = 0.2,
  genes_label_line = 0.6,
  draw_groups_n    = FALSE,
  point_cex        = 0.175,
  title_line       = 3.3,
  draw_border      = TRUE,
  wex              = 0.88,
  quantiles_lty    = c("dashed", "longdash", "dashed"), # compatibility with emf device
  right_gap        = 0.4,
  bracket_color    = "gray50",
  draw_grid        = TRUE,
  indicate_zero    = FALSE,
  show_x_axis      = FALSE,
  show_y_axis      = FALSE
)
par(old_par)
dev.off()




# Compile log2FC data -----------------------------------------------------

logfc_CRISPRoff_df <- StandardLog2FCDf(counts_df)
logfc_CRISPRoff_df <- data.frame("sgID" = CRISPRoff_df[, "sgID_AB"],
                                 logfc_CRISPRoff_df, stringsAsFactors = FALSE
                                 )
logfc_CRISPRoff_df[, "Entrez_ID"] <- as.integer(logfc_CRISPRoff_df[, "Entrez_ID"])



# Examine bidirectional promoters -----------------------------------------

BidirectionalViolins(bidirectional_df, logfc_CRISPRoff_df, max_distance = 10000)
BidirectionalViolins(bidirectional_df, logfc_CRISPRoff_df, max_distance = 15000)
BidirectionalViolins(bidirectional_df, logfc_CRISPRoff_df, max_distance = 20000,
                     num_controls = 30L
                     )



# Save data ---------------------------------------------------------------

num_missing_CRISPRoff <- num_missing
gini_indices_CRISPRoff <- gini_indices
separation_CRISPRoff_mat <- separation_mat
ROC_CRISPRoff_df <- both_reps_ROC_df_list_list[[1]][["ROC_BvT12"]]

save(list = c("logfc_CRISPRoff_df", "gini_indices_CRISPRoff",
              "separation_CRISPRoff_mat", "num_missing_CRISPRoff",
              "ROC_CRISPRoff_df"
              ),
     file = file.path(rdata_dir, "09_create_figures_from_count_data.RData")
     )



