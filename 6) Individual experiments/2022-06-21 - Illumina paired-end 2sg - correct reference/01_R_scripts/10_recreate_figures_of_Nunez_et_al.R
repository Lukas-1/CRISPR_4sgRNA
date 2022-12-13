## 2022-05-23



# Load packages and source code -------------------------------------------

library("readxl")

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "01_violin_swarm_plots.R"))
source(file.path(R_functions_dir, "02_ROC_curves.R"))
source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R")) # For scatter plots



# Define paths ------------------------------------------------------------

project_dir       <- file.path(experiments_directory, "2022-06-21 - Illumina paired-end 2sg - correct reference")
rdata_dir         <- file.path(project_dir, "03_R_objects")
figures_dir       <- file.path(project_dir, "04_output_data", "Figures")
PDFs_dir          <- file.path(figures_dir, "PDFs")
recreate_figs_dir <- file.path(PDFs_dir, "Re-create figures of Nunez et al")
library_path      <- file.path(first_illumina_trial_dir, "02_input_data", "2021 - Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing - Table S3.xlsx")
first_rdata_dir   <- file.path(first_illumina_trial_dir, "03_R_objects")
thesis_dir        <- file.path(figures_dir, "PDFs", "Thesis")


# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))



# Read in data ------------------------------------------------------------

Nunez_df <- data.frame(readxl::read_excel(library_path, skip = 3, sheet = 2),
                       check.names = FALSE, stringsAsFactors = FALSE
                       )



# Define functions --------------------------------------------------------

RocOffvMutant <- function(CRISPRoff_ROC_df,
                          mutant_ROC_df,
                          title_text = expression(bold("Original data of Nu\u00f1ez" ~ bolditalic("et al."))),
                          num_digits_AUC = 2
                          ) {

  line_colors <- c("#6daca0", "#c9944a")

  old_mar <- par(mar = c(3.7, 3.95, 3.6, 1.75))

  CRISPRoff_AUC <- PlotROCDf(CRISPRoff_ROC_df, line_color = "#6daca0",
                             show_AUC = FALSE, flip = TRUE
                             )
  mutant_AUC <- PlotROCDf(mutant_ROC_df, line_color = "#c9944a", add = TRUE,
                          show_AUC = FALSE, flip = TRUE
                          )

  legend_vec <- c(
    as.expression(bquote("CRISPRoff" ~ scriptstyle("(AUC" ~
                         .(format(round(CRISPRoff_AUC, digits = num_digits_AUC), nsmall = num_digits_AUC)) * ")"))),
    as.expression(bquote("CRISPRoff mutant" ~ scriptstyle("(AUC" ~
                         .(format(round(mutant_AUC, digits = num_digits_AUC), nsmall = num_digits_AUC)) * ")")))
  )

  y_start <- diff(grconvertY(c(0, 0.7), from = "lines", to = "user"))
  y_pos_vec <- y_start + c(diff(grconvertY(c(0, 1), from = "lines", to = "user")), 0)

  text(x      = grconvertX(0.3, from = "npc", to = "user"),
       y      = y_pos_vec,
       labels = legend_vec,
       adj    = c(0, 0.5),
       xpd    = NA
       )

  segments(x0   = grconvertX(0.245, from = "npc", to = "user"),
           x1   = grconvertX(0.285, from = "npc", to = "user"),
           y0   = y_pos_vec,
           col  = line_colors,
           lend = "butt",
           lwd  = 2
           )

  title(title_text, cex.main = par("cex"), line = 1.7)

  par(old_mar)
  return(invisible(NULL))
}



# Tidy data ---------------------------------------------------------------

names(Nunez_df) <- gsub("[- ]", "_", names(Nunez_df))



# Explore data ------------------------------------------------------------

range(Nunez_df[, "CRISPRoff_Rep1"]) * 10
range(Nunez_df[, "CRISPRoff_Rep2"]) * 10
range(Nunez_df[, "CRISPRoff_mut_Rep1"]) * 10
range(Nunez_df[, "CRISPRoff_mut_Rep2"]) * 10

table(Nunez_df[, "CpG"])
table(Nunez_df[, "hit_and_run_hit"])
table(Nunez_df[, "mutant_hit"])



# Add Entrez gene IDs -----------------------------------------------------

matches_vec <- match(Nunez_df[, "sgID"], CRISPRoff_df[, "sgID_AB"])
stopifnot(!(anyNA(matches_vec)))
for (column_name in c("Entrez_ID", "Gene_symbol", "Is_preferred_plasmid")) {
  Nunez_df[, column_name] <- CRISPRoff_df[, column_name][matches_vec]
}
new_order <- order(
  match(Nunez_df[, "Entrez_ID"], Nunez_df[, "Entrez_ID"]),
  !(Nunez_df[, "Is_preferred_plasmid"])
)
Nunez_df <- Nunez_df[new_order, ]
row.names(Nunez_df) <- NULL



# Check for the availability of essential and non-essential genes ---------

essential_entrezs     <- intersect(essentials_2020Q2_df[, "Entrez_ID"], Nunez_df[, "Entrez_ID"])
non_essential_entrezs <- intersect(non_essentials_2020Q2_df[, "Entrez_ID"], Nunez_df[, "Entrez_ID"])



# Prepare for re-creating the ROC curves of Nunez et al. ------------------

ROC_columns <- c(
  "Entrez_ID", "Gene_symbol",
  "CRISPRoff_average", "mutant_average",
  "CRISPRoff_Rep1", "CRISPRoff_Rep2", "CRISPRoff_mut_Rep1", "CRISPRoff_mut_Rep2"
)

ROC_input_df <- Nunez_df[, ROC_columns]
ROC_input_df <- ROCInputDf(ROC_input_df, essential_entrezs, non_essential_entrezs)



# Draw ROC curves ---------------------------------------------------------

CRISPRoff_both_reps_ROC_df <- ROCDfForColumn(ROC_input_df, "CRISPRoff_average")
mutant_both_reps_ROC_df <- ROCDfForColumn(ROC_input_df, "mutant_average")

CRISPRoff_rep1_ROC_df <- ROCDfForColumn(ROC_input_df, "CRISPRoff_Rep1")
mutant_rep1_ROC_df <- ROCDfForColumn(ROC_input_df, "CRISPRoff_mut_Rep1")

CRISPRoff_rep2_ROC_df <- ROCDfForColumn(ROC_input_df, "CRISPRoff_Rep2")
mutant_rep2_ROC_df <- ROCDfForColumn(ROC_input_df, "CRISPRoff_mut_Rep2")


for (create_PDF in c(FALSE, TRUE)) {

  if (create_PDF) {
    pdf(file.path(recreate_figs_dir, "Fig. 4F - ROC curves.pdf"),
        width = 3 + 1.14, height = 3 + 1.46
        )
  }

  RocOffvMutant(CRISPRoff_both_reps_ROC_df, mutant_both_reps_ROC_df)
  RocOffvMutant(CRISPRoff_rep1_ROC_df, mutant_rep1_ROC_df,
                title_text = expression(bold("Nu\u00f1ez" ~ bolditalic("et al.") ~ bold("\u2013 replicate 1"))),
                num_digits_AUC = 4
                )
  RocOffvMutant(CRISPRoff_rep2_ROC_df, mutant_rep2_ROC_df,
                title_text = expression(bold("Nu\u00f1ez" ~ bolditalic("et al.") ~ bold("\u2013 replicate 2"))),
                num_digits_AUC = 4
                )

  if (create_PDF) {
    dev.off()
  }
}





# Create violin plots -----------------------------------------------------

are_for_violins <- !(duplicated(ROC_input_df[, "Entrez_ID"]))
violin_columns <- c("Is_essential", "CRISPRoff_Rep1", "CRISPRoff_Rep2",
                    "CRISPRoff_mut_Rep1", "CRISPRoff_mut_Rep2"
                    )
violins_df <- ROC_input_df[are_for_violins, violin_columns]
row.names(violins_df) <- NULL

rep_list <- c(split(violins_df[, "CRISPRoff_Rep1"], !(violins_df[, "Is_essential"])),
              split(violins_df[, "CRISPRoff_Rep2"], !(violins_df[, "Is_essential"])),
              split(violins_df[, "CRISPRoff_mut_Rep1"], !(violins_df[, "Is_essential"])),
              split(violins_df[, "CRISPRoff_mut_Rep2"], !(violins_df[, "Is_essential"]))
              )[c(1, 3, 2, 4, 5, 7, 6, 8)]


for (create_PDF in c(FALSE, TRUE)) {

  if (create_PDF) {
    pdf(file.path(recreate_figs_dir, "Fig. 4E - violin plots.pdf"), width = 7, height = 4)
  }

  old_mar <- par(mar = c(2.75, 4.3, 5.7, 2))

  x_positions <- BeeViolinPlot(rep_list,
                               violin_colors = rep(c("#d6e6e3", "#edddc5"), each = 4),
                               point_colors  = rep(c("#4e7e76", "#ac7c35"), each = 4), #
                               point_cex     = 0.25,
                               use_spacing   = 0.45,
                               adjust        = 1,
                               lower_bound   = -0.6,
                               upper_bound   = 0.1,
                               y_limits      = c(-0.608, 0.1),
                               groups_vec    = c(rep(1, 4), rep(2, 4)),
                               gap_ratio     = 1.5
                               )

  mtext(expression("Phenotype (" * gamma * ")"), side = 2, line = 2.75)

  segments(x0  = x_positions[c(1, 3, 5, 7)] - 0.25,
           x1  = x_positions[c(2, 4, 6, 8)] + 0.25,
           y0  = par("usr")[[4]] + diff(grconvertY(c(0, 0.5), from = "lines", to = "user")),
           col = "black",
           xpd = NA
           )
  mtext(rep(c("Essential\ngenes", "Non-essential\ngenes"), 2),
        at = tapply(x_positions, rep(1:4, each = 2), mean),
        line = 0.65, padj = 0, cex = par("cex")
        )
  mtext(text = rep(c("R1", "R2"), 4), at = x_positions, side = 1, line = 0.7, cex = par("cex"))

  mtext(c("CRISPRoff", "CRISPRoff mutant"),
        at = tapply(x_positions, rep(1:2, each = 4), mean),
        line = 3.2, padj = 0, cex = par("cex"), font = 2
        )

  gap_pos <- mean(x_positions[4:5])
  gap_width <- (x_positions[[5]] - x_positions[[4]])
  rect(xleft   = gap_pos - 0.1,
       xright  = gap_pos + 0.1,
       ybottom = par("usr")[[3]] - GetHalfLineWidth() * 3,
       ytop    = par("usr")[[3]] + GetHalfLineWidth() * 3,
       col     = "white",
       border  = NA,
       xpd     = NA
       )

  if (create_PDF) {
    dev.off()
  }
}



# Create violin plots for CRISPRoff only ----------------------------------

devEMF::emf(file.path(thesis_dir, "3A) Violin plot - i) re-analysis.emf"),
            width = 2, height = 2, emfPlus = FALSE, coordDPI = 3000
            )
old_par <- par(cex = 0.6, lwd = 0.7, lheight = 0.9, mar = c(3, 4, 4, 1))
custom_y_limits <- c(-0.6 - (0.86 * 0.02), 0.25)
x_positions <- BeeViolinPlot(
  rep_list[1:4],
  point_cex        = 0.175,
  use_spacing      = 0.5,
  wex              = 0.88,
  violin_colors    = rep(c(brewer.pal(9, "Purples")[[3]], "#c7e7c0"), each = 2),
  point_colors     = rep(c("#7c7198", "#5b8669"), each = 2),
  border_colors    = rep(c("#d1cddb", "#bfd4c6"), each = 2),
  cloud_alpha      = 0.2,
  cloud_sd         = 0.04,
  lower_bound      = -0.6,
  upper_bound      = 0.25,
  y_limits         = custom_y_limits,
  draw_groups_n    = FALSE,
  draw_border      = TRUE,
  quantiles_lty    = c("dotted", "dashed", "dotted"),
  right_gap        = 0.4,
  show_x_axis      = FALSE,
  indicate_zero    = FALSE,
  draw_grid        = TRUE,
  grid_lwd         = 0.8
)

mtext(expression("Phenotype (" * gamma * ")"), side = 2, line = 2.1, cex = par("cex"))
segments(x0  = x_positions[c(1, 3)] - 0.25,
         x1  = x_positions[c(2, 4)] + 0.25,
         y0  = par("usr")[[4]] + diff(grconvertY(c(0, 0.45), from = "lines", to = "user")),
         col = "gray50",
         xpd = NA
         )
mtext(c("Essential\ngenes", "Non-essential\ngenes"),
      at = c(mean(x_positions[1:2]), mean(x_positions[3:4])),
      line = 0.6, padj = 0, cex = par("cex")
      )
mtext(text = rep(paste0("R", 1:2), 2),
      at = x_positions, side = 1, line = 0.2, cex = par("cex")
      )
title(expression(bold("Data reanalysis (Nu\u00F1ez et al.)")), cex.main = 1, line = 3.3)
dev.off()




# Create scatter plots ----------------------------------------------------

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
      pdf(file.path(recreate_figs_dir, paste0("Scatter plots - ", file_postfix, ".pdf")),
          width = use_width, height = 4.7
          )
    }
    scatter_df <- data.frame(
      Nunez_df[, c("Entrez_ID", "Gene_symbol")],
      "Is_NT" = Nunez_df[, "gene"] == "non-targeting",
      "Rep1_data" = Nunez_df[, "CRISPRoff_Rep1"],
      "Rep2_data" = Nunez_df[, "CRISPRoff_Rep2"],
      stringsAsFactors = FALSE
    )
    for (show_mutant in c(FALSE, TRUE)) {
      if (show_mutant) {
        scatter_df[, "Rep1_data"] <- Nunez_df[, "CRISPRoff_mut_Rep1"]
        scatter_df[, "Rep2_data"] <- Nunez_df[, "CRISPRoff_mut_Rep2"]
      }
      ReplicateScatterPlot(scatter_df,
                           show_phenotype_score = TRUE,
                           use_title            = if (show_mutant) "CRISPRoff mutant" else "CRISPRoff",
                           highlight_NT         = highlight_NT,
                           highlight_essential  = highlight_essential,
                           embed_PNG            = create_PDF
                           )
    }
    rm(scatter_df)
    if (create_PDF) {
      dev.off()
    }
  }
}



# Save data ---------------------------------------------------------------

separation_original_mat <- SeparationMetrics(
  setNames(rep_list[1:4],
           paste0(rep(c("Essential R", "Non-essential R"), each = 2), 1:2)
           )
)[1:4, ]

logfc_original_df <- data.frame(
  Nunez_df[, c("sgID", "Gene_symbol")],
  "Entrez_ID" = as.integer(Nunez_df[, "Entrez_ID"]),
  "Mean_log2FC" = rowMeans(Nunez_df[, c("CRISPRoff_Rep1", "CRISPRoff_Rep2")]) * 10,
  "Log2FC_rep1" = Nunez_df[, "CRISPRoff_Rep1"] * 10,
  "Log2FC_rep2" = Nunez_df[, "CRISPRoff_Rep2"] * 10,
  "Is_NT" = Nunez_df[, "gene"] == "non-targeting",
  stringsAsFactors = FALSE
)

ROC_original_df <- CRISPRoff_both_reps_ROC_df

save(list = c("logfc_original_df", "separation_original_mat", "ROC_original_df"),
     file = file.path(rdata_dir, "10_recreate_figures_of_Nunez_et_al.RData")
     )


