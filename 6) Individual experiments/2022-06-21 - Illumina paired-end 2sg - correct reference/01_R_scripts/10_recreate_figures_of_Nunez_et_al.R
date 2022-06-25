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



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))



# Read in data ------------------------------------------------------------

Nunez_df <- data.frame(readxl::read_excel(library_path, skip = 3, sheet = 2),
                       check.names = FALSE, stringsAsFactors = FALSE
                       )


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
Nunez_df[, "Entrez_ID"] <- CRISPRoff_df[, "Entrez_ID"][matches_vec]
Nunez_df[, "Gene_symbol"] <- CRISPRoff_df[, "Gene_symbol"][matches_vec]




# Check for the availability of essential and non-essential genes ---------

essential_entrezs     <- intersect(essentials_2020Q2_df[, "Entrez_ID"], Nunez_df[, "Entrez_ID"])
non_essential_entrezs <- intersect(non_essentials_2020Q2_df[, "Entrez_ID"], Nunez_df[, "Entrez_ID"])



# Prepare for re-creating the ROC curves of Nunez et al. ------------------

are_essential <- ifelse(Nunez_df[, "Entrez_ID"] %in% essential_entrezs,
                        TRUE,
                        ifelse(Nunez_df[, "Entrez_ID"] %in% non_essential_entrezs,
                               FALSE,
                               NA
                               )
                        )
Nunez_df[, "Is_essential"] <- are_essential

ROC_columns <- c(
  "Entrez_ID", "Gene_symbol", "Is_essential", "CRISPRoff_average", "mutant_average"
)

ROC_df <- Nunez_df[!(is.na(are_essential)), ROC_columns]
new_order <- order(ROC_df[, "Is_essential"], decreasing = TRUE)
ROC_df <- ROC_df[new_order, ]
row.names(ROC_df) <- NULL




# Draw ROC curves ---------------------------------------------------------

ROCDfForColumn <- function(input_df, use_column) {
  new_order <- order(input_df[, use_column])
  results_df <- input_df[new_order, ]
  results_df <- MakeROCDf(results_df, numeric_column = use_column)
  return(results_df)
}

CRISPRoff_ROC_df <- ROCDfForColumn(ROC_df, "CRISPRoff_average")
mutant_ROC_df <- ROCDfForColumn(ROC_df, "mutant_average")

line_colors <- c("#6daca0", "#c9944a")

for (create_PDF in c(FALSE, TRUE)) {

  if (create_PDF) {
    pdf(file.path(recreate_figs_dir, "Fig. 4F - ROC curves.pdf"),
        width = 3 + 1.14, height = 3 + 1.46
        )
  }

  old_mar <- par(mar = c(3.7, 3.95, 3.6, 1.75))

  CRISPRoff_AUC <- PlotROCDf(CRISPRoff_ROC_df, line_color = "#6daca0",
                             show_AUC = FALSE, flip = TRUE
                             )
  mutant_AUC <- PlotROCDf(mutant_ROC_df, line_color = "#c9944a", add = TRUE,
                          show_AUC = FALSE, flip = TRUE
                          )

  legend_vec <- c(
    as.expression(bquote("CRISPRoff" ~ scriptstyle("(AUC" ~ .(round(CRISPRoff_AUC, digits = 2)) * ")"))),
    as.expression(bquote("CRISPRoff mutant" ~ scriptstyle("(AUC" ~ .(round(mutant_AUC, digits = 2)) * ")")))
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

  title(expression(bold("Original data of Nu\u00f1ez" ~ bolditalic("et al."))),
        cex.main = par("cex"), line = 1.7
        )

  par(old_mar)

  if (create_PDF) {
    dev.off()
  }
}





# Create violin plots -----------------------------------------------------

are_for_violins <- (!(duplicated(Nunez_df[, "Entrez_ID"]))) &
                   (!(is.na(Nunez_df[, "Is_essential"])))
violin_columns <- c("Is_essential", "CRISPRoff_Rep1", "CRISPRoff_Rep2",
                    "CRISPRoff_mut_Rep1", "CRISPRoff_mut_Rep2"
                    )
violins_df <- Nunez_df[are_for_violins, violin_columns]
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



