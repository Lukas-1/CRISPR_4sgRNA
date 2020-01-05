### 19th December 2019 ###





# Import packages and source code -----------------------------------------

library("Biostrings")

library("readxl")
library("beeswarm")
library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
file_input_directory             <- file.path(file_directory, "1) Input")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")
exchange_with_tools_directory    <- file.path(file_directory, "3) Exchange with other tools")
file_output_directory            <- file.path(file_directory, "4) Output")





# Read in data ------------------------------------------------------------

accuracies_df <- read.table(file.path(file_input_directory, "identical.perc.txt"),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                            check.names = FALSE
                            )

guides_df <- data.frame(read_excel(file.path(file_input_directory, "1-384 oligo sequence.xlsx"),
                                            sheet = "1-384 primers", n_max = 384, col_names = FALSE
                                            ),
                                 stringsAsFactors = FALSE, check.names = FALSE
                                 )





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





# Functions for using DINAMelt --------------------------------------------

FwdAndRevCombinations <- function(four_sg_vec) {
  stopifnot(length(four_sg_vec) == 4)
  rev_sequences <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(four_sg_vec)))
  mat_list <- lapply(1:4, function(x) {
    second_sequence_vec <- c(four_sg_vec, rev_sequences[(1:4) != x])
    sub_mat <- cbind(rep(four_sg_vec[[x]], length(second_sequence_vec)), second_sequence_vec)
    return(sub_mat)
  })
  results_mat <- do.call(rbind, mat_list)
  colnames(results_mat) <- paste0("Sequence_", 1:2)
  rownames(results_mat) <- NULL
  return(results_mat)
}


OnlyFwdCombinations <- function(four_sg_vec) {
  stopifnot(length(four_sg_vec) == 4)
  results_mat <- as.matrix(expand.grid(four_sg_vec, four_sg_vec, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))[, 2:1]
  colnames(results_mat) <- paste0("Sequence_", 1:2)
  return(results_mat)
}


MatListToDf <- function(mat_list) {
  df_list <- lapply(1:384, function(x) data.frame("Gene_symbol" = assembly_df[x, "Gene_symbol"],
                                                  "Well_number" = x,
                                                  mat_list[[x]],
                                                  stringsAsFactors = FALSE
                                                  )
                    )
  results_df <- do.call(rbind.data.frame, c(df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}


ExportToTwoStateMelting <- function(sequence_vec, file_prefix) {

  num_sequences_per_file <- 1000L
  num_sequences <- length(sequence_vec)
  num_files <- ceiling(num_sequences / num_sequences_per_file)

  file_sequence <- rep(seq_len(num_files), each = num_sequences_per_file)
  file_sequence <- file_sequence[seq_len(num_sequences)]

  for (i in seq_len(num_files)) {
    sequence_sub_vec <- sequence_vec[file_sequence == i]
    sequence_sub_vec <- paste0(sequence_sub_vec, c(rep(";", times = length(sequence_sub_vec) - 1), ""))
    write.table(sequence_sub_vec,
                file = file.path(exchange_with_tools_directory, "1) Input for DINAMelt", paste0(file_prefix, " - file ", i, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )

  }
}







# Process guides_df -------------------------------------------------------

are_empty_columns <- vapply(seq_len(ncol(guides_df)), function(x) all(is.na(guides_df[, x])), logical(1))
guides_df <- guides_df[, !(are_empty_columns)]

symbols_mat <- do.call(cbind, unique(lapply(seq(1, ncol(guides_df), by = 2),
                                            function(x) sapply(strsplit(guides_df[, x], "[_ ]"), "[[", 1)
                                            )
                                     )
                       )

are_identical <- apply(symbols_mat, 1, function(x) length(unique(x)) == 1)
symbols_mat[!(are_identical), ] # MIR20A is very likely an error introduced by Excel




# Compile data ------------------------------------------------------------

sequences_mat <- as.matrix(guides_df[, seq(2, ncol(guides_df), by = 2)])
colnames(sequences_mat) <- paste0("sg_", 1:4)

accuracy_colnames <- paste0("sg", 1:4, "_cr", 1:4, "_identical")

accuracies_1sg_mat <- do.call(cbind, lapply(accuracy_colnames, function(x) accuracies_df[, x] / accuracies_df[, "tot"]))
colnames(accuracies_1sg_mat) <- paste0("Fraction_correct_sg", 1:4)

all_4_sums_vec <- rowSums(as.matrix(accuracies_df[, accuracy_colnames]))

assembly_df <- data.frame(
  "Gene_symbol"           = symbols_mat[, 1],
  "Num_reads"             = accuracies_df[, "tot"],
  "Fraction_correct_4sg"  = all_4_sums_vec / (accuracies_df[, "tot"] * 4),
  "Fraction_correct_3sg"  = rowSums(as.matrix(accuracies_df[, paste0("sg", 1:3, "_cr", 1:3, "_identical")])) /
                            (accuracies_df[, "tot"] * 3),
  accuracies_1sg_mat,
  "Longest_subsequence"   = NA,
  "Are_accuracy_outliers" = all_4_sums_vec == 0,
  sequences_mat,
  stringsAsFactors        = FALSE,
  row.names               = NULL
)

guides_list <- lapply(1:384, function(x) unname(sequences_mat[x, ]))

assembly_df[, "Longest_subsequence"] <- vapply(guides_list, LongestSharedSubsequence, integer(1))





# Re-order assembly_df ----------------------------------------------------

assembly_df <- assembly_df[order(assembly_df[, "Fraction_correct_4sg"], decreasing = FALSE), ]

rownames(assembly_df) <- NULL






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






# Generate input for the DINAMelt two-state melting tool ------------------
# http://unafold.rna.albany.edu/?q=DINAMelt/Two-state-melting

only_fwd_df <- MatListToDf(lapply(1:384, function(x) OnlyFwdCombinations(sequences_mat[x, ])))
fwd_and_rev_df <- MatListToDf(lapply(1:384, function(x) FwdAndRevCombinations(sequences_mat[x, ])))


ExportToTwoStateMelting(only_fwd_df[, "Sequence_1"], "1a) Only forward sequences")
ExportToTwoStateMelting(only_fwd_df[, "Sequence_2"], "1b) Only forward sequences")
ExportToTwoStateMelting(fwd_and_rev_df[, "Sequence_1"], "2a) Forward and reverse sequences")
ExportToTwoStateMelting(fwd_and_rev_df[, "Sequence_2"], "2b) Forward and reverse sequences")






# Process output from the DINAMelt two-state melting tool -----------------

ReadTwoStateMeltingOutput <- function(file_prefix) {

  all_files <- list.files(file.path(exchange_with_tools_directory, "2) Output from DINAMelt"))
  selected_files <- grep(file_prefix, all_files, fixed = TRUE, value = TRUE)

  file_numbers <- as.integer(sub("file ", "", sapply(strsplit(selected_files, " - ", fixed = TRUE), "[[", 2), fixed = TRUE))
  selected_files <- selected_files[order(file_numbers)]

  output_mat_list <- lapply(selected_files, function(x) {
    output_df <- read.table(file.path(exchange_with_tools_directory, "2) Output from DINAMelt", x),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                            )[, 2:5]

    for (i in 1:4) {
      output_df[, i] <- as.numeric(gsub("[^-+.0-9]", "", output_df[, i]))
    }
    output_mat <- as.matrix(output_df)
    return(output_mat)
  })
  results_mat <- do.call(rbind, output_mat_list)
  colnames(results_mat) <- c("Gibbs_free_energy", "enthalpy", "entropy", "Tm_Celcius")
  return(results_mat)
}


only_fwd_melting_mat <- ReadTwoStateMeltingOutput("1) Only forward sequences")
only_fwd_melting_df <- data.frame(only_fwd_df,
                                  only_fwd_melting_mat,
                                  stringsAsFactors = FALSE,
                                  row.names = NULL
                                  )


fwd_and_rev_melting_mat <- ReadTwoStateMeltingOutput("2) Forward and reverse sequences")
fwd_and_rev_melting_df <- data.frame(fwd_and_rev_df,
                                     fwd_and_rev_melting_mat,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL
                                     )



for (column_name in colnames(only_fwd_melting_mat)[1:3]) {
  assembly_df[, paste0("Only_fwd_", column_name)] <- tapply(only_fwd_melting_df[, column_name],
                                                            only_fwd_melting_df[, "Well_number"],
                                                            min
                                                            )
}
assembly_df[, "Only_fwd_Tm_Celcius"] <- tapply(only_fwd_melting_df[, "Tm_Celcius"],
                                               only_fwd_melting_df[, "Well_number"],
                                               max
                                               )


for (column_name in colnames(fwd_and_rev_melting_mat)[1:3]) {
  assembly_df[,  paste0("Fwd_and_rev_", column_name)] <- tapply(fwd_and_rev_melting_df[, column_name],
                                                                fwd_and_rev_melting_df[, "Well_number"],
                                                                min
                                                                )
}
assembly_df[, "Fwd_and_rev_Tm_Celcius"] <- tapply(fwd_and_rev_melting_df[, "Tm_Celcius"],
                                                  fwd_and_rev_melting_df[, "Well_number"],
                                                  max
                                                  )







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






# Export sequences for viennaRNA self-annealing prediction ----------------

single_sequences <- unlist(lapply(paste0("sg_", 1:4), function(x) assembly_df[, x]))

write.table(single_sequences,
            file = file.path(exchange_with_tools_directory, "Single sequences for self-annealing prediction.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )






# Process minimum free energy (MFE) scores from ViennaRNA -----------------

load(file.path(exchange_with_tools_directory, "ViennaRNA_minimum_free_energies.RData"))

ViennaRNA_mat <- matrix(minimum_free_energies_vec,
                        nrow = nrow(assembly_df),
                        dimnames = list(NULL, paste0("MFE_sg", 1:4))
                        )

ViennaRNA_mat <- cbind(ViennaRNA_mat,
                       "MinMFE_4sg" = apply(ViennaRNA_mat, 1, min),
                       "MinMFE_3sg" = apply(ViennaRNA_mat[, 2:4], 1, min)
                       )

assembly_df <- data.frame(assembly_df, ViennaRNA_mat, stringsAsFactors = FALSE, row.names = NULL)





# Plot MFE scores from ViennaRNA ------------------------------------------

AccuracyScatterPlot(assembly_df, "MinMFE_4sg", "Fraction_correct_4sg")
AccuracyScatterPlot(assembly_df, "MinMFE_3sg", "Fraction_correct_4sg")

AccuracyScatterPlot(assembly_df, "MinMFE_4sg", "Fraction_correct_3sg")
AccuracyScatterPlot(assembly_df, "MinMFE_3sg", "Fraction_correct_3sg")





# Build a data frame for individual guides --------------------------------

individual_guides_df <- data.frame(
  "Gene_symbol"      = rep(assembly_df[, "Gene_symbol"], 4),
  "Sequence"         = as.vector(sequences_mat),
  "Num_reads"        = rep(assembly_df[, "Num_reads"], 4),
  "Fraction_correct" = as.vector(accuracies_1sg_mat),
  "MFE"              = minimum_free_energies_vec,
  stringsAsFactors   = FALSE,
  row.names          = NULL
)




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




# Save data ---------------------------------------------------------------

save(list = c("assembly_df", "individual_guides_df"),
     file = file.path(intermediate_R_objects_directory,
                      "Correlate the accuracy of the Gibson assembly with the length of shared subsequences.RData"
                      )
     )






