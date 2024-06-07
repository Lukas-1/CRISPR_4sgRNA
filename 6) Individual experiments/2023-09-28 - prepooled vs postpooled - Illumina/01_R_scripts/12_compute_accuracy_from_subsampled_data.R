### 2024-05-05


# Load packages and source code -------------------------------------------

library("scales") # for CombinedTemplateSwitch

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")
project_dir              <- file.path(experiments_directory, "2023-09-28 - prepooled vs postpooled - Illumina")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))
source(file.path(project_dir, "01_R_scripts", "R_functions", "02_annotating_plots.R"))



# Define paths ------------------------------------------------------------

first_rdata_dir <- file.path(first_illumina_trial_dir, "03_R_objects")
rdata_dir       <- file.path(project_dir, "03_R_objects")
figures_dir     <- file.path(project_dir, "04_output", "Figures")
PDFs_dir        <- file.path(figures_dir, "Subsampling")


# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "11_produce_subsampled_read_counts.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))



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


samples_in_order <- c(
  "Prepool_T0_R1", "Prepool_T0_R2", "Prepool_T12_R1", "Prepool_T12_R2",
  "Postpool_T0_R1", "Postpool_T0_R2", "Postpool_T12_R1", "Postpool_T12_R2"
)

subsampled_counts_mat_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    y <- y[new_order, ]
    colnames(y) <- sub("-", "", colnames(y), fixed = TRUE)
    colnames(y) <- sub("_S[0-8]$", "", colnames(y))
    y <- y[, samples_in_order]
    return(y)
  })
})



# Compute ROC curves ------------------------------------------------------

essential_list <- GetEssentialGenes(use_blomen_hart = TRUE)

RocDfForMat <- function(input_mat, library_df, use_essential_list, postpool = FALSE, choose_rep = NULL) {

  if (postpool) {
    use_baseline_indices <- 5:6
    use_intervention_indices <- 7:8
  } else {
    use_baseline_indices <- 1:2
    use_intervention_indices <- 3:4
  }

  colnames(input_mat) <- paste0("NoSwitch_xMM_", colnames(input_mat))
  make_counts_df <- data.frame(
    CRISPRoff_df[, c("Plasmid_ID", "Gene_symbol", "Entrez_ID")],
    input_mat,
    stringsAsFactors = FALSE
  )
  assign("counts_df", make_counts_df, envir = globalenv())

  results_df <- GetEssentialROCDf(
    essential_genes       = use_essential_list[["essential_entrezs"]],
    non_essential_genes   = essential_list[["non_essential_entrezs"]],
    min_count_at_baseline = 0L,
    baseline_indices      = use_baseline_indices,
    intervention_indices  = use_intervention_indices,
    choose_rep            = choose_rep,
    allow_switch          = FALSE
  )

  rm(counts_df, envir = globalenv())
  return(results_df)
}


prepool_bothreps_AUC_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    ROC_df <- RocDfForMat(y,
                          CRISPRoff_df,
                          essential_list,
                          postpool = FALSE,
                          choose_rep = NULL
                          )
    PlotEssentialROCDf(ROC_df)
  })
})

prepool_rep1_AUC_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    ROC_df <- RocDfForMat(y,
                          CRISPRoff_df,
                          essential_list,
                          postpool = FALSE,
                          choose_rep = 1
                          )
    PlotEssentialROCDf(ROC_df)
  })
})

prepool_rep2_AUC_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    ROC_df <- RocDfForMat(y,
                          CRISPRoff_df,
                          essential_list,
                          postpool = FALSE,
                          choose_rep = 2
                          )
    PlotEssentialROCDf(ROC_df)
  })
})


postpool_bothreps_AUC_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    ROC_df <- RocDfForMat(y,
                          CRISPRoff_df,
                          essential_list,
                          postpool = TRUE,
                          choose_rep = NULL
                          )
    PlotEssentialROCDf(ROC_df)
  })
})

postpool_rep1_AUC_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    ROC_df <- RocDfForMat(y,
                          CRISPRoff_df,
                          essential_list,
                          postpool = TRUE,
                          choose_rep = 1
                          )
    PlotEssentialROCDf(ROC_df)
  })
})

postpool_rep2_AUC_list <- lapply(subsampled_counts_mat_list, function(x) {
  lapply(x, function(y) {
    ROC_df <- RocDfForMat(y,
                          CRISPRoff_df,
                          essential_list,
                          postpool = TRUE,
                          choose_rep = 2
                          )
    PlotEssentialROCDf(ROC_df)
  })
})



# Plot violins ------------------------------------------------------------

RepEssentialViolinsForMat <- function(input_mat, postpool = FALSE, use_title = "") {

  if (postpool) {
    use_baseline_indices <- 5:6
    use_intervention_indices <- 7:8
  } else {
    use_baseline_indices <- 1:2
    use_intervention_indices <- 3:4
  }

  colnames(input_mat) <- paste0("NoSwitch_xMM_", colnames(input_mat))
  make_counts_df <- data.frame(
    CRISPRoff_df[, c("Plasmid_ID", "Gene_symbol", "Entrez_ID")],
    input_mat,
    stringsAsFactors = FALSE
  )
  assign("counts_df", make_counts_df, envir = globalenv())

  custom_y_limits <- c(-0.6 - (0.86 * 0.02), 0.2)

  reps_list <- RepEssentialViolins(
    baseline_indices = use_baseline_indices,
    intervention_indices = use_intervention_indices,
    use_title        = use_title,
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

  rm(counts_df, envir = globalenv())
  return(invisible(reps_list))
}



pdf(file.path(PDFs_dir, "Violin plots - subsampling - prepooled.pdf"),
    width = 2, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.8, lheight = 0.9)
prepool_separation_mat_list <- lapply(names(subsampled_counts_mat_list), function(x) {
  lapply(names(subsampled_counts_mat_list[[x]]), function(y) {
    use_title <- paste0("Prepool \u2013 ", x, " \u2013 ", y)
    reps_list <- RepEssentialViolinsForMat(subsampled_counts_mat_list[[x]][[y]], postpool = FALSE, use_title = use_title)
    SeparationMetrics(reps_list)
  })
})
dev.off()



pdf(file.path(PDFs_dir, "Violin plots - subsampling - postpooled.pdf"),
    width = 2, height = 2
    )
old_par <- par(cex = 0.6, lwd = 0.8, lheight = 0.9)
postpool_separation_mat_list <- lapply(names(subsampled_counts_mat_list), function(x) {
  lapply(names(subsampled_counts_mat_list[[x]]), function(y) {
    use_title <- paste0("Postpool \u2013 ", x, " \u2013 ", y)
    reps_list <- RepEssentialViolinsForMat(subsampled_counts_mat_list[[x]][[y]], postpool = TRUE, use_title = use_title)
    SeparationMetrics(reps_list)
  })
})
dev.off()



prepool_SSMD_rep1_list <- lapply(prepool_separation_mat_list, function(x) {
  lapply(x, function(y) y["Robust SSMD finite", "R1"])
})
prepool_SSMD_rep2_list <- lapply(prepool_separation_mat_list, function(x) {
  lapply(x, function(y) y["Robust SSMD finite", "R2"])
})
postpool_SSMD_rep1_list <- lapply(postpool_separation_mat_list, function(x) {
  lapply(x, function(y) y["Robust SSMD finite", "R1"])
})
postpool_SSMD_rep2_list <- lapply(postpool_separation_mat_list, function(x) {
  lapply(x, function(y) y["Robust SSMD finite", "R2"])
})



# Compile results ---------------------------------------------------------

subsampling_df <- data.frame(
  "Fraction_sampled"   = rep(as.numeric(sub("% sampled", "", names(subsampled_counts_mat_list))),
                             lengths(subsampled_counts_mat_list)
                             ),
  "Resampling_rep"     = unlist(lapply(subsampled_counts_mat_list, seq_along)),
  "Prepool_AUC_both"   = unlist(prepool_bothreps_AUC_list),
  "Prepool_AUC_rep1"   = unlist(prepool_rep1_AUC_list),
  "Prepool_AUC_rep2"   = unlist(prepool_rep1_AUC_list),
  "Postpool_AUC_both"  = unlist(postpool_bothreps_AUC_list),
  "Postpool_AUC_rep1"  = unlist(postpool_rep1_AUC_list),
  "Postpool_AUC_rep2"  = unlist(postpool_rep1_AUC_list),
  "Prepool_SSMD_rep1"  = unlist(prepool_SSMD_rep1_list),
  "Prepool_SSMD_rep2"  = unlist(prepool_SSMD_rep2_list),
  "Postpool_SSMD_rep1" = unlist(postpool_SSMD_rep1_list),
  "Postpool_SSMD_rep2" = unlist(postpool_SSMD_rep2_list),
  row.names            = NULL
)



# Save data ---------------------------------------------------------------

save(subsampling_df,
     file = file.path(rdata_dir, "12_compute_accuracy_from_subsampled_data.RData")
     )

