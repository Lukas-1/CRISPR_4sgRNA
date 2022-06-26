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



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_disambiguate_dJR072_CRISPRoff_library.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))



## Examine only the genes used in Nunez et al. (intersection of Blomen and Hart):

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

essential_entrezs     <- GetAvailableGenes(essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)
non_essential_entrezs <- GetAvailableGenes(non_essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)

T0vT12_title <- expression(bold("CRISPRoff:" ~ bolditalic("T0")      ~ "vs." ~ bolditalic("T12")))
BvT12_title  <- expression(bold("CRISPRoff:" ~ bolditalic("Tbefore") ~ "vs." ~ bolditalic("T12")))
BvT0_title   <- expression(bold("CRISPRoff:" ~ bolditalic("Tbefore") ~ "vs." ~ bolditalic("T0")))


ROC_df_list_list <- lapply(c(FALSE, TRUE), function(x) {
  list(
    ROC_T0vT12_df = GetEssentialROCDf(essential_entrezs,
                                      non_essential_entrezs,
                                      baseline_indices      = 3:4,
                                      intervention_indices  = 5:6,
                                      min_count_at_baseline = 0L,
                                      allow_switch          = x
                                      ),
    ROC_BvT12_df  = GetEssentialROCDf(essential_entrezs,
                                      non_essential_entrezs,
                                      baseline_indices      = 1:2,
                                      intervention_indices  = 5:6,
                                      min_count_at_baseline = 0L,
                                      allow_switch          = x
                                      ),
    ROC_BvT0_df   = GetEssentialROCDf(essential_entrezs,
                                      non_essential_entrezs,
                                      baseline_indices      = 1:2,
                                      intervention_indices  = 3:4,
                                      min_count_at_baseline = 0L,
                                      allow_switch          = x
                                      )
  )
})


for (allow_switch in c(FALSE, TRUE)) {

  if (allow_switch) {
    use_dir <- file.path(PDFs_dir, "Including template switch")
  } else {
    use_dir <- file.path(PDFs_dir, "Excluding template switch")
  }

  old_oma <- par(oma = rep(0.3, 4))
  for (create_PDF in c(FALSE, TRUE)) {
    if (create_PDF) {
      pdf(file.path(use_dir, "ROC curves.pdf"),
          width = 3 + 1.14, height = 3 + 1.56
          )
    }
    PlotEssentialROCDf(ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12_df"]], use_title = T0vT12_title)
    PlotEssentialROCDf(ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12_df"]],  use_title = BvT12_title)
    PlotEssentialROCDf(ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0_df"]],   use_title = BvT0_title)
    if (create_PDF) {
      dev.off()
    }
  }
  par(old_oma)
}



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
      do.call(ViolinPlotEssentialDf, c(list(ROC_df_list_list[[allow_switch + 1]][["ROC_T0vT12_df"]], use_title = T0vT12_title, draw_points = draw_points), custom_args))
      do.call(ViolinPlotEssentialDf, c(list(ROC_df_list_list[[allow_switch + 1]][["ROC_BvT12_df"]],  use_title = BvT12_title,  draw_points = draw_points), custom_args))
      do.call(ViolinPlotEssentialDf, c(list(ROC_df_list_list[[allow_switch + 1]][["ROC_BvT0_df"]],   use_title = BvT0_title,   draw_points = draw_points), custom_args))
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
      RepEssentialViolins(3:4, 5:6, use_title = T0vT12_title,draw_points = draw_points,
                          lower_bound = -0.6, upper_bound = 0.25, y_limits = custom_y_limits,
                          allow_switch = allow_switch
                          )
      RepEssentialViolins(1:2, 5:6, use_title = BvT12_title, draw_points = draw_points,
                          lower_bound = -0.6, upper_bound = 0.25, y_limits = custom_y_limits,
                          allow_switch = allow_switch
                          )
      RepEssentialViolins(1:2, 3:4, use_title = BvT0_title, draw_points = draw_points,
                          lower_bound = -0.6, upper_bound = 0.25, y_limits = custom_y_limits,
                          allow_switch = allow_switch
                          )
    }
    if (create_PDF) {
      dev.off()
    }
  }
}



# Draw violin plots showing all samples -----------------------------------

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
      do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 3:4, intervention_indices = 5:6, use_title = T0vT12_title)))
      do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 1:2, intervention_indices = 5:6, use_title = BvT12_title)))
      do.call(Log2FCScatterPlot, c(args_list, list(baseline_indices = 1:2, intervention_indices = 3:4, use_title = BvT0_title)))
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
      DrawHistogram(numeric_vec      = rowMeans(as.matrix(counts_df[, use_columns, drop = FALSE])),
                    truncation_limit = 5000L,
                    num_breaks       = 150L,
                    title_text       = use_title,
                    x_axis_label     = if (length(use_columns) == 1) "Read count" else "Mean read count",
                    y_axis_label     = "Number of plasmids",
                    y_axis_limits    = c(-20, 1600)
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



