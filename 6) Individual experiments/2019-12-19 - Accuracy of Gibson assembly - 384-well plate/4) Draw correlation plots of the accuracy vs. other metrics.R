### 19th December 2019 ###



# Import packages and source code -----------------------------------------

library("beeswarm")
library("RColorBrewer")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")
file_output_directory            <- file.path(file_directory, "4) Output")




# Load data ---------------------------------------------------------------

load(file.path(intermediate_R_objects_directory, "3) Import data from external tools.RData"))





# Functions for drawing plots ---------------------------------------------

y_labels <- c(
  "Fraction_correct_3sg" = "% correct (average of sg2, sg3, sg4)",
  "Fraction_correct_4sg" = "% correct (average of 4 guide sequences)"
)


MakeCorrTitle <- function(numeric_vec_1, numeric_vec_2) {
  corr_results <- cor.test(numeric_vec_1, numeric_vec_2)
  title_expression <- as.expression(bquote(
    bold("Pearson's " * bolditalic("r") * " = " *
           .(as.character(signif(corr_results[["estimate"]], digits = 2)))) *
      " ("  * italic("p") * " = " *
      .(as.character(signif(corr_results[["p.value"]], digits = 1))) * ")"
  ))
  text(labels = title_expression, x = par("usr")[[1]] + (par("usr")[[2]] - par("usr")[[1]]) / 2, y = 111, xpd = NA)
  return(invisible(corr_results))
}



PlotBySharedSubsequence <- function(use_assembly_df, only_sg2to4 = FALSE, exclude_outliers = TRUE) {

  if (only_sg2to4) {
    accuracy_column <- "Fraction_correct_3sg"
  } else {
    accuracy_column <- "Fraction_correct_4sg"
  }

  if (exclude_outliers) {
    use_assembly_df <- use_assembly_df[!(use_assembly_df[, "Are_accuracy_outliers"]), ]
  }

  light_color <- brewer.pal(9, "Blues")[[2]]
  dark_color <- brewer.pal(9, "Blues")[[7]]

  display_smallest <- 3
  display_biggest <- 19

  plot(1,
       xlim = c(display_smallest, display_biggest),
       ylim = c(0, 100),
       xaxs = "i",
       yaxs = "i",
       ylab = y_labels[[accuracy_column]],
       xlab = "Length of shared subsequence",
       type = "n",
       axes = FALSE,
       mgp  = c(2.7, 1, 0)
       )

  abline(h = seq(10, 100, by = 10), col = "gray90")#, lwd = 0.75)


  corr_results <- MakeCorrTitle(use_assembly_df[, accuracy_column], use_assembly_df[, "Longest_subsequence"])


  axis(2, las = 1, tcl = -0.4, mgp = c(3, 0.7, 0))
  box(bty = "l")

  accuracy_vec <- use_assembly_df[, accuracy_column] * 100
  longest_fac <- factor(use_assembly_df[, "Longest_subsequence"], levels = 1:20)

  boxplot(accuracy_vec ~ longest_fac,
          ylim      = c(0.5, 1),
          boxwex    = 0.7,
          outline   = FALSE,
          names     = rep("", 20),
          whisklty  = "blank",
          staplewex = 0,
          axes      = FALSE,
          whisklwd  = 0,
          staplelty = 0,
          col       = light_color,
          # lwd       = 0.75,
          add       = TRUE
          )

  beeswarm(accuracy_vec ~ longest_fac,
           pch     = 16,
           cex     = 0.35,
           spacing = 0.6,
           col     = dark_color,
           xpd     = NA,
           add     = TRUE
           )

  present_lengths <- seq(from = min(use_assembly_df[, "Longest_subsequence"]),
                         to = max(use_assembly_df[, "Longest_subsequence"])
                         )

  text(x      = present_lengths,
       y      = -6,
       labels = present_lengths,
       font   = 1,
       cex    = 0.8,
       xpd    = NA
       )

  return(invisible(corr_results))
}




x_labels = list(
  "Only_fwd_Tm_Celcius"           = "Melting temperature [?C] \u2013 only forward strands",
  "Fwd_and_rev_Tm_Celcius"        = "Melting temperature [?C] \u2013 both strands",
  "Only_fwd_enthalpy"             = "Enthalpy [kcal/mol] \u2013 only forward strands",
  "Fwd_and_rev_enthalpy"          = "Enthalpy [kcal/mol] \u2013 both strands",
  "Only_fwd_entropy"              = "Entropy [cal/mol/K] \u2013 only forward strands",
  "Fwd_and_rev_entropy"           = "Entropy [cal/mol/K] \u2013 both strands",
  "Only_fwd_Gibbs_free_energy"    = "Gibbs free energy [kcal/mol] \u2013 only forward strands",
  "Fwd_and_rev_Gibbs_free_energy" = "Gibbs free energy [kcal/mol] \u2013 both strands",

  "MinMFE_4sg"                    = "ViennaRNA minimum free energy (MFE) \u2013 minimum of sg1-4",
  "MinMFE_3sg"                    = "ViennaRNA minimum free energy (MFE) \u2013 minimum of sg2, sg3, s4"
)




AccuracyScatterPlot <- function(use_assembly_df, sequence_column, accuracy_column, exclude_accuracy_outliers = TRUE) {

  old_par <- par(mar = rep(6, 4))

  if (exclude_accuracy_outliers) {
    use_assembly_df <- use_assembly_df[!(use_assembly_df[, "Are_accuracy_outliers"]), ]
  }

  longest_subsequence_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(20)

  plot(use_assembly_df[, sequence_column],
       use_assembly_df[, accuracy_column] * 100,
       col  = paste0(longest_subsequence_colors[use_assembly_df[, "Longest_subsequence"]], "99"), #60% alpha
       ylim = c(0, 100),
       yaxs = "i",
       las  = 1,
       mgp  = c(2.7, 0.5, 0),
       tcl  = -0.4,
       lwd  = 1.25,
       xlab = x_labels[[sequence_column]],
       ylab = y_labels[[accuracy_column]]
       )

  longest_subseq_seq <- seq(min(use_assembly_df[, "Longest_subsequence"]), max(use_assembly_df[, "Longest_subsequence"]))
  legend_space <- paste0(rep(" ", 11), collapse = "")

  legend("right",
         title     = paste0(legend_space, "Longest\n", legend_space, "subsequence"),
         inset     = -0.31,
         legend    = as.character(longest_subseq_seq),
         col       = longest_subsequence_colors[longest_subseq_seq],
         pch       = 1,
         pt.lwd    = 2,
         xpd       = NA,
         cex       = 0.8,
         pt.cex    = 1,
         bty       = "n",
         adj       = c(0, 0),
         title.adj = c(0, 0),
         y.intersp = 1.2
         )

  corr_results <- MakeCorrTitle(use_assembly_df[, accuracy_column], use_assembly_df[, sequence_column])

  par(old_par)

  return(invisible(corr_results))
}





# Make plots of accuracy vs. the longest subsequence ----------------------

options(scipen = 6)

PlotBySharedSubsequence(assembly_df)
PlotBySharedSubsequence(assembly_df, only_sg2to4 = TRUE)

pdf(file = file.path(file_output_directory, "Accuracy vs. length of shared subsequences.pdf"),
    width = 8, height = 5.5
    )
PlotBySharedSubsequence(assembly_df)
PlotBySharedSubsequence(assembly_df, only_sg2to4 = TRUE)
dev.off()





# Plot data from the DINAMelt tool ----------------------------------------

plot(only_fwd_melting_df[, "enthalpy"], only_fwd_melting_df[, "entropy"])
plot(only_fwd_melting_df[, "enthalpy"], only_fwd_melting_df[, "Tm_Celcius"])

cor.test(assembly_df[, "Only_fwd_enthalpy"],          assembly_df[, "Fraction_correct_4sg"])
cor.test(assembly_df[, "Only_fwd_entropy"],           assembly_df[, "Fraction_correct_4sg"])
cor.test(assembly_df[, "Only_fwd_Gibbs_free_energy"], assembly_df[, "Fraction_correct_4sg"])
cor.test(assembly_df[, "Only_fwd_Tm_Celcius"],        assembly_df[, "Fraction_correct_4sg"])

cor.test(assembly_df[, "Only_fwd_enthalpy"],          assembly_df[, "Fraction_correct_3sg"])
cor.test(assembly_df[, "Only_fwd_entropy"],           assembly_df[, "Fraction_correct_3sg"])
cor.test(assembly_df[, "Only_fwd_Gibbs_free_energy"], assembly_df[, "Fraction_correct_3sg"])
cor.test(assembly_df[, "Only_fwd_Tm_Celcius"],        assembly_df[, "Fraction_correct_3sg"])

cor.test(assembly_df[, "Fwd_and_rev_enthalpy"],          assembly_df[, "Fraction_correct_4sg"])
cor.test(assembly_df[, "Fwd_and_rev_entropy"],           assembly_df[, "Fraction_correct_4sg"])
cor.test(assembly_df[, "Fwd_and_rev_Gibbs_free_energy"], assembly_df[, "Fraction_correct_4sg"])
cor.test(assembly_df[, "Fwd_and_rev_Tm_Celcius"],        assembly_df[, "Fraction_correct_4sg"])

cor.test(assembly_df[, "Fwd_and_rev_enthalpy"],          assembly_df[, "Fraction_correct_3sg"])
cor.test(assembly_df[, "Fwd_and_rev_entropy"],           assembly_df[, "Fraction_correct_3sg"])
cor.test(assembly_df[, "Fwd_and_rev_Gibbs_free_energy"], assembly_df[, "Fraction_correct_3sg"])
cor.test(assembly_df[, "Fwd_and_rev_Tm_Celcius"],        assembly_df[, "Fraction_correct_3sg"])


AccuracyScatterPlot(assembly_df, "Only_fwd_Tm_Celcius",           "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_Gibbs_free_energy",    "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_enthalpy",             "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_entropy",              "Fraction_correct_3sg")

AccuracyScatterPlot(assembly_df, "Only_fwd_Tm_Celcius",           "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_Gibbs_free_energy",    "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_enthalpy",             "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_entropy",              "Fraction_correct_4sg")

AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Tm_Celcius",        "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Gibbs_free_energy", "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_enthalpy",          "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_entropy",           "Fraction_correct_3sg")

AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Tm_Celcius",        "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Gibbs_free_energy", "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_enthalpy",          "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_entropy",           "Fraction_correct_4sg")




plot_dimensions <- 6

pdf(file = file.path(file_output_directory, "Accuracy vs. melting temperature.pdf"),
    width = plot_dimensions, height = plot_dimensions
    )

AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Tm_Celcius",        "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_Tm_Celcius",           "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Gibbs_free_energy", "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_Gibbs_free_energy",    "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_enthalpy",          "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_enthalpy",             "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_entropy",           "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_entropy",              "Fraction_correct_3sg")

AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Tm_Celcius",        "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_Tm_Celcius",           "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_Gibbs_free_energy", "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_Gibbs_free_energy",    "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_enthalpy",          "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_enthalpy",             "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Fwd_and_rev_entropy",           "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "Only_fwd_entropy",              "Fraction_correct_4sg")

dev.off()





# Plot MFE scores from ViennaRNA ------------------------------------------

AccuracyScatterPlot(assembly_df, "MinMFE_4sg", "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "MinMFE_3sg", "Fraction_correct_4sg")

AccuracyScatterPlot(assembly_df, "MinMFE_4sg", "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "MinMFE_3sg", "Fraction_correct_3sg")






# Plot the data for individual guides -------------------------------------

for (make_PDF in c(FALSE, TRUE)) {

  if (make_PDF) {
    plot_dimensions <- 6
    pdf(file = file.path(file_output_directory, "Accuracy vs. minimum free energy.pdf"),
        width = plot_dimensions, height = plot_dimensions
        )
  }
  old_par <- par(mar = rep(6, 4))

  plot(individual_guides_df[, "MFE"], individual_guides_df[, "Fraction_correct"] * 100,
       las  = 1,
       mgp  = c(2.7, 1, 0),
       tcl  = -0.4,
       ylab = "% correct (invididual guides)",
       xlab = "ViennaRNA minimum free energy (MFE)",
       type = "n"
       )

  abline(v = -6, h = 80, lty = "dotted", col = "gray80")

  points(individual_guides_df[, "MFE"], individual_guides_df[, "Fraction_correct"] * 100,
         lwd  = 1.5,
         col  = paste0(brewer.pal(9, "Blues")[[8]], "99"),
         )

  box()

  MakeCorrTitle(individual_guides_df[, "MFE"], individual_guides_df[, "Fraction_correct"])

  par(old_par)

  if (make_PDF) {
    dev.off()
  }

}



