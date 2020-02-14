### 4th February 2020 ###


# Import packages and source code -----------------------------------------

library("vioplot")
library("eulerr")
library("gridExtra")
library("ggplot2")
library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R"))





# Define maps -------------------------------------------------------------

core_numeric_column_labels <- c(
  "all22_SNP_AF_max_Kaviar" = "Most frequent polymorphism within the sgRNA/PAM region (Kaviar database)",
  "CRISPOR_Doench_efficacy" = "Cutting efficiency score (Doench et al.) from CRISPOR",
  "GuideScan_efficiency"    = "Cutting efficiency score (Doench et al.) from GuideScan",
  "GuideScan_specificity"   = "GuideScan specificity score",
  "CRISPOR_3MM_specificity" = "CRISPOR specificity score (up to 3 mismatches)",
  "CRISPOR_4MM_specificity" = "CRISPOR specificity score (up to 4 mismatches)",
  "CRISPOR_CFD_specificity" = "CRISPOR CFD specificity score",
  "CRISPOR_MIT_specificity" = "CRISPOR MIT specificity score"
)

numeric_column_labels <- c(
  core_numeric_column_labels,
  "GuideScan_Num_2MM"    = "GuideScan \u2013 number of 2MM sites",
  "GuideScan_Num_3MM"    = "GuideScan \u2013 number of 3MM sites",
  "GuideScan_Num_2or3MM" = "GuideScan \u2013 number of 2MM or 3MM sites",
  "CRISPOR_Num_2MM"      = "CRISPOR \u2013 number of 2MM sites",
  "CRISPOR_Num_3MM"      = "CRISPOR \u2013 number of 3MM sites",
  "CRISPOR_Num_2or3MM"   = "CRISPOR \u2013 number of 2MM or 3MM sites",
  "CRISPOR_Num_4MM"      = "CRISPOR \u2013 number of 4MM sites"
)




specificity_score_cutoffs <- c(
  "GuideScan_specificity"   = 0.2,
  "CRISPOR_3MM_specificity" = 0.2,
  "CRISPOR_4MM_specificity" = 0.04,
  "CRISPOR_CFD_specificity" = 50,
  "CRISPOR_MIT_specificity" = 50
)

categorical_columns <- c(
  "Do_not_meet_criteria",
  "Not_mapped",
  preferred_AF_max_column,
  "Poly_T",
  "Non_canonical_PAM",
  "CRISPOR_Graf_status",
  names(core_numeric_column_labels)
)

libraries_order <- c(
  "Calabrese",
  "hCRISPRa-v2",

  "Brunello",
  "TKOv3",
  "GPP"
)





# Settings for generating multi-plot layouts ------------------------------

use_layout_mat <- rbind(c(3, 3), c(1, 2), c(4, 4))
layout_widths <- c(2, 6)
layout_heights <- c(1, 19, 0.3)
pdf_width <- 8
pdf_height <- 4







# General helper functions ------------------------------------------------

FilterCRISPRDf <- function(CRISPR_df) {
  CRISPR_df[(CRISPR_df[["Is_control"]] == "No") & (CRISPR_df[["Source"]] != "Curated"), ]
}

palify_cache_101 <- list()
Palify <- function(myhex, fraction_pale = 0.5) {
  if (myhex %in% names(palify_cache_101)) {
    color_vec <- palify_cache_101[[myhex]]
  } else {
    color_vec <- colorRampPalette(c(myhex, "#FFFFFF"))(101)
    palify_cache_101[[myhex]] <- color_vec
    assign("palify_cache_101", palify_cache_101, envir = globalenv())
  }
  color_vec[[round(fraction_pale * 100) + 1]]
}






# Helper functions for plotting -------------------------------------------

SubgroupColorsDf <- function() {
  dark_colors <- vapply(c("Blues", "Greens", "Reds", "Purples", "Greys"),
                        function(x) colorRampPalette(brewer.pal(9, x))(18)[[17]],
                        ""
                        )
  names(dark_colors) <- substr(names(dark_colors), 1, nchar(names(dark_colors)) - 1)
  dark_colors <- c(
    dark_colors[1:3],
    "BlueGreen" = "#006969",
    dark_colors["Purple"],
    "Brown" = "#581200",
    dark_colors["Grey"]
  )
  colors_df <- data.frame(
    "Color_name" = names(dark_colors),
    "Lightened"  = FALSE,
    "Pale"       = vapply(dark_colors, Palify, fraction_pale = 0.9, ""),
    "Medium"     = vapply(dark_colors, Palify, fraction_pale = 0.775, ""),
    "Dark"       = dark_colors,
    stringsAsFactors = FALSE
  )
  colors_df <- colors_df[rep(seq_len(nrow(colors_df)), each = 2), ]
  rownames(colors_df) <- NULL
  colors_df[["Lightened"]] <- rep(c(TRUE, FALSE), length.out = nrow(colors_df))
  for (column_name in c("Pale", "Medium", "Dark")) {
    colors_df[[column_name]][colors_df[["Lightened"]]] <- vapply(colors_df[[column_name]][colors_df[["Lightened"]]], Palify, fraction_pale = 0.5, "")
  }
  return(colors_df)
}



ReformatSourceToFactor <- function(source_vec) {

  source_vec <- sub("Curated, ", "", source_vec, fixed = TRUE)

  unique_sources <- unique(source_vec)
  reordered_unique_sources <- vapply(strsplit(unique_sources, ", ", fixed = TRUE),
                                     function(x) paste0(x[order(x == "GPP")], collapse = ", "),
                                     ""
                                     )
  for (i in which(unique_sources != reordered_unique_sources)) {
    source_vec[source_vec == unique_sources[[i]]] <- reordered_unique_sources[[i]]
  }

  levels_present <- libraries_order[libraries_order %in% source_vec]
  stopifnot(length(levels_present) == 3)
  levels_expanded <- c(apply(combn(levels_present, 2), 2, paste0, collapse = ", "),
                       paste0(levels_present, collapse = ", ")
                       )

  source_fac <- factor(source_vec, levels = c(levels_present, levels_expanded))
  return(source_fac)
}




MakeSubgroupsPlotDf <- function(CRISPR_df, use_column) {

  ## Rename and re-order the groups
  groups_fac <- ReformatSourceToFactor(CRISPR_df[["Source"]])


  ## Split into subgroups
  chosen_fac <- factor(ifelse(CRISPR_df[["Rank"]] %in% 1:4, "4sg", "Rest"), levels = c("Rest", "4sg"))
  subgroups_fac <- interaction(groups_fac, chosen_fac, sep = "_", lex.order = TRUE)

  ## Build a data frame for plotting
  plot_df <- data.frame(
    "Data"     = CRISPR_df[[use_column]],
    "Group"    = groups_fac,
    "Subgroup" = subgroups_fac,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[order(plot_df[["Subgroup"]]), ]
  rownames(plot_df) <- NULL
  return(plot_df)
}



DrawSubgroupLabels <- function(plot_df, colors_df, x_positions, extra_space = FALSE, only_chosen = FALSE) {
  if (extra_space) {
    first_line_factor  <- 0.065
    second_line_factor <- 0.13
  } else {
    first_line_factor  <- 0.035
    second_line_factor <- 0.1
  }
  if (only_chosen) {
    second_line_factor <- second_line_factor <- 0.04
  }
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  num_subgroups <- nrow(colors_df)

  ## Draw the small subgroup labels (Rest/4sg)
  if (!(only_chosen)) {
    text(x      = x_positions,
         y      = par("usr")[[3]] - (plot_height * first_line_factor),
         labels = c("Rest", "4sg"),
         adj    = c(0.5, 1),
         cex    = 0.75,
         col    = rep(c(Palify("#000000", fraction_pale = 0.5), "#000000"),
                      length.out = num_subgroups
                      ),
         font   = 2,
         xpd    = NA
         )
  }

  ## Draw the group labels
  group_labels <- sub(", ", ",\n", levels(plot_df[["Group"]]))
  group_labels[[length(group_labels)]] <- "All 3\nsources"
  old_lheight <- par("lheight" = 1.15)
  colors_vec <- colors_df[["Dark"]]
  if (!(only_chosen)) {
    colors_vec <- colors_vec[!(colors_df[["Lightened"]])]
    x_positions <- tapply(x_positions, rep(seq_len(num_subgroups / 2), each = 2), mean)
  }
  text(x      = x_positions,
       y      = par("usr")[[3]] - (plot_height * second_line_factor),
       labels = group_labels,
       adj    = c(0.5, 1),
       cex    = 0.8,
       col    = colors_vec,
       font   = 2,
       xpd    = NA
       )
  par(old_lheight)
  return(invisible(NULL))
}






# Functions for plotting numerical data -----------------------------------


FixSNPColumn <- function(CRISPR_df, y_column) {
  if (grepl("SNP", y_column, fixed = TRUE)) {
    CRISPR_df[[y_column]] <- ifelse(is.na(CRISPR_df[[y_column]]) & !(is.na(CRISPR_df[["Start"]])),
                                    0,
                                    CRISPR_df[[y_column]]
                                    ) * 100
  }
  return(CRISPR_df)
}



GetAxisLimits <- function(numeric_vec, provide_other_limits = FALSE, check_limits = FALSE, extend_range_fraction = 0.02) {
  max_value <- max(numeric_vec, na.rm = TRUE)
  upper_limit <- 10^ceiling(log10(max_value))
  if (upper_limit == 1) {
    limits <- c(-0.02, 1.02)
  } else if (upper_limit == 100) {
    limits <- c(-2, 102)
    if (min(numeric_vec, na.rm = TRUE) > 2) {
      limits[[1]] <- 0
    }
    if (max_value < 98) {
      limits[[2]] <- 100
    }
  } else {
    if (check_limits) {
      stop("The numerical values to plot on the y axis did not fall into the expected ranges!")
    } else {
      if (provide_other_limits) {
        limits <- range(numeric_vec, na.rm = TRUE)
        span <- limits[[2]] - limits[[1]]
        limits[[1]] <- limits[[1]] - (span * 0.02)
        limits[[2]] <- limits[[2]] + (span * 0.02)
      } else {
        return(NULL)
      }
    }
  }
  return(limits)
}


DrawViolinGridAndAxes <- function(y_limits, y_column, show_title = TRUE) {
  y_limits <- par("usr")[c(3, 4)]
  abline(h = seq(0 + (y_limits[[2]] / 10), y_limits[[2]], by = y_limits[[2]] / 5),
         col = "gray95", lwd = 0.5
         )
  even_seq <- seq(0, y_limits[[2]], by = y_limits[[2]] / 5)
  if (y_column %in% c("CRISPOR_3MM_specificity", "GuideScan_specificity")) {
    even_seq <- even_seq[-2]
    abline(h = 0.2, col = "gray50", lwd = 0.5)
  }
  abline(h = even_seq, col = "gray90", lwd = 0.5)
  axis_ticks <- axTicks(2)
  tick_labels <- format(axis_ticks)
  if (grepl("SNP", y_column, fixed = TRUE)) {
    tick_labels <- paste0(axis_ticks, "%")
  }
  axis(2, labels = tick_labels, at = axis_ticks,
       las = 1, mgp = c(3, 0.45, 0), tcl = -0.3,
       cex.axis = 0.8, lwd = 0.75
       )
  title_text <- numeric_column_labels[[y_column]]
  if (show_title) {
    title(title_text, cex.main = 1)
  }
  return(invisible(title_text))
}



TwoViolinBox <- function(CRISPR_df, y_column, show_title = TRUE) {

  CRISPR_df <- FixSNPColumn(FilterCRISPRDf(CRISPR_df), y_column)

  ## Prepare the group colors
  colors_df <- SubgroupColorsDf()
  colors_df <- colors_df[colors_df[["Color_name"]] == "Grey", ]
  points_alpha <- 0.1
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)


  ## Build a data frame for plotting
  plot_df <- data.frame(
    "Numeric" = CRISPR_df[[y_column]],
    "Chosen"  = CRISPR_df[["Rank"]] %in% 1:4,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[order(plot_df[["Chosen"]]), ]
  rownames(plot_df) <- NULL
  plot_df[["Point_color"]] <- rep(colors_df[["Dark"]], as.integer(table(plot_df[["Chosen"]])))


  ## Prepare the x positions
  x_positions <- 1:2
  plot_df[["Jittered_x"]] <- x_positions[as.integer(plot_df[["Chosen"]]) + 1] +
                             rnorm(n = nrow(plot_df), mean = 0, sd = 0.03)


  ## Prepare the plot region
  old_mar <- par(mar = c(6, 4, 3, 0.2) + 0.1)
  y_limits <- GetAxisLimits(plot_df[["Numeric"]], provide_other_limits = TRUE)
  num_groups <- 2
  plot(1,
       xlim = c(0.5 - (num_groups * 0.04), num_groups + 0.5 + (num_groups * 0.04)),
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  title_text <- DrawViolinGridAndAxes(y_limits, y_column, show_title = show_title)


  ## Draw the violin plots
  vioplot(plot_df[["Numeric"]] ~ plot_df[["Chosen"]],
          add      = TRUE,
          at       = x_positions,
          pchMed   = NA,
          drawRect = FALSE,
          col      = colors_df[["Medium"]],
          border   = NA,
          axes     = FALSE
          )

  ## Draw the jittered points
  points(x   = plot_df[["Jittered_x"]],
         y   = plot_df[["Numeric"]],
         cex = 0.5,
         col = paste0(plot_df[["Point_color"]], alpha_hex),
         pch = 16
         )

  ## Draw the superimposed boxplots
  boxplot(plot_df[["Numeric"]] ~ plot_df[["Chosen"]],
          add       = TRUE,
          at        = x_positions,
          cex       = 0.2,
          boxwex    = 0.325,
          outline   = FALSE,
          names     = rep("", num_groups),
          whisklty  = "blank",
          staplewex = 0,
          whisklwd  = 0,
          staplelty = 0,
          col       = colors_df[["Pale"]],
          border    = colors_df[["Dark"]],
          axes      = FALSE,
          lwd       = 1
          )


  ## Draw the subgroup numbers
  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  assign("delete_plot_df", plot_df, envir = globalenv())

  text(x      = x_positions,
       y      = par("usr")[[3]] - (plot_height * 0.0275),
       labels = table(plot_df[["Chosen"]][!(is.na(plot_df[["Numeric"]]))]),
       adj    = c(0.5, 1),
       cex    = 0.7,
       col    = "gray40",
       xpd    = NA
       )


  ## Draw the group labels
  old_lheight <- par("lheight" = 1.15)
  text(x      = x_positions,
       y      = par("usr")[[3]] - (plot_height * 0.085),
       labels = c("Rest", "4sg"),
       adj    = c(0.5, 1),
       cex    = 0.9,
       col    =  "black",
       font   = 2,
       xpd    = NA
       )
  par(old_lheight)


  ## Final steps
  box(xpd = NA, lwd = 0.75)
  par(old_mar)
  results_list <- list("plot_df" = plot_df, "title_text" = title_text)
  return(invisible(results_list))
}





ViolinBox <- function(CRISPR_df, y_column, show_title = TRUE) {

  CRISPR_df <- FixSNPColumn(FilterCRISPRDf(CRISPR_df), y_column)

  plot_df <- MakeSubgroupsPlotDf(CRISPR_df, y_column)
  num_subgroups <- nlevels(plot_df[["Subgroup"]])

  ## Prepare the group colors
  colors_df <- SubgroupColorsDf()
  points_alpha <- 0.1
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)

  plot_df[["Point_color"]] <- rep(colors_df[["Dark"]], as.integer(table(plot_df[["Subgroup"]])))


  ## Prepare the x positions
  x_positions <- seq_len(num_subgroups) + (rep_len(c(1, -1), length.out = num_subgroups) * 0.18)
  plot_df[["Jittered_x"]] <- x_positions[as.integer(plot_df[["Subgroup"]])] +
                             rnorm(n = nrow(plot_df), mean = 0, sd = 0.03)


  ## Prepare the plot region
  old_mar <- par(mar = c(6, 4, 3, 3) + 0.1)

  y_limits <- GetAxisLimits(plot_df[["Data"]], provide_other_limits = TRUE)
  plot(1,
       xlim = c(0.85 - (num_subgroups * 0.04), num_subgroups + 0.15 + (num_subgroups * 0.04)),
       ylim = y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  title_text <- DrawViolinGridAndAxes(y_limits, y_column, show_title = show_title)

  ## Draw the violin plots
  vioplot(plot_df[["Data"]] ~ plot_df[["Subgroup"]],
          add      = TRUE,
          at       = x_positions,
          pchMed   = NA,
          drawRect = FALSE,
          col      = colors_df[["Medium"]],
          border   = NA,
          wex      = 0.75,
          axes     = FALSE
          )

  ## Draw the jittered points
  points(x   = plot_df[["Jittered_x"]],
         y   = plot_df[["Data"]],
         cex = 0.5,
         col = paste0(plot_df[["Point_color"]], alpha_hex),
         pch = 16
         )

  ## Draw the superimposed boxplots
  boxplot(plot_df[["Data"]] ~ plot_df[["Subgroup"]],
          add       = TRUE,
          at        = x_positions,
          cex       = 0.2,
          boxwex    = 0.275,
          outline   = FALSE,
          names     = rep("", num_subgroups),
          whisklty  = "blank",
          staplewex = 0,
          whisklwd  = 0,
          staplelty = 0,
          col       = colors_df[["Pale"]],
          border    = colors_df[["Dark"]],
          axes      = FALSE,
          lwd       = 0.75
          )


  ## Draw the subgroup numbers
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  text(x      = x_positions,
       y      = par("usr")[[3]] - (plot_height * 0.0125),
       labels = table(plot_df[["Subgroup"]][!(is.na(plot_df[["Data"]]))]),
       adj    = c(0.5, 1),
       cex    = 0.5,
       col    = "gray40",
       # font   = 2,
       xpd    = NA
       )
  DrawSubgroupLabels(plot_df, colors_df, x_positions, extra_space = TRUE)

  ## Final steps
  box(xpd = NA, lwd = 0.75)
  par(old_mar)
  results_list <- list("plot_df" = plot_df, "title_text" = title_text)
  return(invisible(results_list))
}




PlotNumericalData <- function(CRISPR_df) {
  for (make_PDF in c(FALSE, TRUE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "4sg library overview - boxplots.pdf"),
          width = pdf_width, height = pdf_height
          )
    }
    for (y_column in names(numeric_column_labels)) {
      layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
      TwoViolinBox(CRISPR_df, y_column, show_title = FALSE)
      ViolinBox_results <- ViolinBox(CRISPR_df, y_column, show_title = FALSE)
      OuterTitleForLayout(ViolinBox_results[["title_text"]], extra_space = TRUE)
    }
    if (make_PDF) {
      dev.off()
    }
  }
  numeric_seq <- seq_along(numeric_column_labels)
  file_numbers <- FormatFixedWidthInteger(numeric_seq)
  for (i in seq_along(numeric_column_labels)) {
    file_name <- paste0(file_numbers[[i]], ") ", names(numeric_column_labels)[[i]])
    png(file = file.path(output_plots_directory, paste0("4sg library overview - boxplots - ", file_name, ".png")),
        width = pdf_width, height = pdf_height, units = "in", res = 300
        )
    layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
    TwoViolinBox(CRISPR_df, names(numeric_column_labels)[[i]], show_title = FALSE)
    ViolinBox_results <- ViolinBox(CRISPR_df, names(numeric_column_labels)[[i]], show_title = FALSE)
    OuterTitleForLayout(ViolinBox_results[["title_text"]])
    dev.off()
  }
}







# Functions for plotting categorical data ---------------------------------

MakeCategoricalList <- function(CRISPR_df, use_column, use_cutoff) {
  assign("delete_use_column", use_column, envir = globalenv())
  assign("delete_CRISPR_df", CRISPR_df, envir = globalenv())
  assign("delete_use_cutoff", use_cutoff, envir = globalenv())
  title_cex <- 1
  if (use_column == "Are_overlapping") {
    logical_vec <- CRISPR_df[["Num_overlaps"]] > 0
    use_title <- "Overlap with other guides (fewer than 50 base pairs apart)"
  } else if (use_column == "Not_mapped") {
    logical_vec <- is.na(CRISPR_df[["Start"]])
    use_title <- "Cannot be mapped to a specific location in the genome"
  } else if (use_column == "Do_not_meet_criteria") {
    logical_vec <- !(MeetCriteria(CRISPR_df))
    use_title <- "Do not meet all criteria"
  } else if (use_column == "CRISPOR_Graf_status") {
    logical_vec <- CRISPR_df[["CRISPOR_Graf_status"]] != "GrafOK"
    use_title <- "Violate criteria of Graf et al."
    use_column <- "Graf_not_OK"
  } else if (use_column == "Poly_T") {
    logical_vec <- grepl("TTTT", CRISPR_df[["sgRNA_sequence"]], ignore.case = TRUE)
    use_title <- "Contain TTTT sequence"
  } else if (use_column == "Non_canonical_PAM") {
    logical_vec <- substr(CRISPR_df[["PAM"]], 2, 3) != "GG"
    use_title <- "Non-canonical PAM (i.e. not NGG)"
  } else if (use_column %in% grep("_SNP_AF_(max|sum)_", colnames(CRISPR_df), value = TRUE)) {
    if (is.null(use_cutoff)) {
      use_cutoff <- SNP_frequency_cutoff
    }
    overlap_vec <- CRISPR_df[[use_column]] > use_cutoff
    logical_vec <- ifelse(is.na(CRISPR_df[["Start"]]),
                          NA,
                          ifelse(is.na(overlap_vec), FALSE, overlap_vec)
                          )
    use_title <- paste0("Overlap with a genetic polymorphism (frequency >",
                        paste0(use_cutoff * 100, "%)")
                        )
    use_column <- "Overlap_with_SNP"
  } else if (use_column %in% names(specificity_score_cutoffs)) {
    if (is.null(use_cutoff)) {
      use_cutoff <- specificity_score_cutoffs[[use_column]]
    }
    logical_vec <- CRISPR_df[[use_column]] < use_cutoff
    use_title <- paste0("Unspecific (<", use_cutoff, ") according to the ", numeric_column_labels[[use_column]])
    title_cex <- 0.9
    use_column <- "Are_unspecific"
  } else if (use_column %in% c("GuideScan_efficiency", "CRISPOR_Doench_efficacy")) {
    if (is.null(use_cutoff)) {
      use_cutoff <- 50
    }
    logical_vec <- CRISPR_df[[use_column]] < use_cutoff
    if (use_column == "GuideScan_efficiency") {
      via_string  <- "GuideScan"
    } else if (use_column == "CRISPOR_Doench_efficacy") {
      via_string <- "CRISPOR"
    }
    use_title <- paste0("Suboptimal cutting efficiency (<", use_cutoff,
                        ") according to the Rule Set 2 score (via ", via_string, ")"
                        )
    title_cex <- 0.9
    use_column <- "Is_inefficient"
  } else {
    stop("Unknown column selected!")
  }
  results_list <- list(
    "logical_vec" = logical_vec,
    "use_column"  = use_column,
    "title_text"  = use_title,
    "title_cex"   = title_cex
  )
  return(results_list)
}




DrawBarplotGridAndAxis <- function() {
  segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]], y0 = seq(0, 1, by = 0.2),   col = "gray95", lwd = 0.5)
  segments(x0 = par("usr")[[1]], x1 = par("usr")[[2]], y0 = seq(0.1, 1, by = 0.2), col = "gray95", lwd = 0.5)
  axis(2, labels = paste0(seq(0, 100, by = 20), "%"), at = seq(0, 1, by = 0.2),
       mgp = c(3, 0.45, 0), las = 1, cex.axis = 0.8, lwd = 0.75, tcl = -0.3
       )
  return(invisible(NULL))
}



AnnotateBars <- function(bar_positions_vec, counts_mat, proportions_mat, text_cex) {

  plot_height <- par("usr")[[4]] - par("usr")[[3]]

  ## Annotate the bars with percentages
  proportions_vec <- ifelse(proportions_mat[2, ] < 0.01,
                            signif(proportions_mat[2, ] * 100, digits = 1),
                            signif(proportions_mat[2, ] * 100, digits = 2)
                            )
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * 0.03),
       labels = paste0(proportions_vec, "%"),
       adj    = c(0.5, 1),
       cex    = text_cex * 1.3,
       col    = "black",
       font   = 2,
       xpd    = NA
       )

  ## Annotate the bars with fractions (counts)
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * 0.08),
       labels = colSums(counts_mat),
       adj    = c(0.5, 0.5),
       cex    = text_cex,
       col    = "black",
       xpd    = NA
       )
  line_widths <- strwidth(colSums(counts_mat), cex = text_cex)
  segments(x0  = bar_positions_vec - (line_widths / 2),
           x1  = bar_positions_vec + (line_widths / 2),
           y0  = par("usr")[[4]] + (plot_height * 0.0925),
           col = "black",
           xpd = NA,
           lwd = 0.3
           )
  text(x      = bar_positions_vec,
       y      = par("usr")[[4]] + (plot_height * 0.105),
       labels = counts_mat[2, ],
       adj    = c(0.5, 0.5),
       cex    = text_cex,
       col    = "black",
       xpd    = NA
       )
  return(invisible(NULL))
}





TwoStackedBarPlot <- function(CRISPR_df, use_column, use_cutoff = NULL, use_dummy_plot = FALSE, show_title = TRUE, only_chosen = FALSE) {

  CRISPR_df <- FilterCRISPRDf(CRISPR_df)

  ## Prepare the data
  categorical_list <- MakeCategoricalList(CRISPR_df, use_column, use_cutoff)
  colors_df <- SubgroupColorsDf()
  colors_df <- colors_df[(colors_df[["Color_name"]] == "Grey") & !(colors_df[["Lightened"]]), ]

  if (only_chosen) {
    are_top_4 <- CRISPR_df[["Rank"]] %in% 1:4
    CRISPR_df <- CRISPR_df[are_top_4, ]
    categorical_list[["logical_vec"]] <- categorical_list[["logical_vec"]][are_top_4]
    num_subgroups <- 1L
    use_space <- NULL
  } else {
    num_subgroups <- 2L
    use_space <- 0.5
  }

  plot_df <- data.frame(
    "Data"   = categorical_list[["logical_vec"]],
    "Chosen" = CRISPR_df[["Rank"]] %in% 1:4,
    stringsAsFactors = FALSE
  )
  counts_mat <- as.matrix(table(factor(plot_df[["Data"]], levels = c(FALSE, TRUE)), plot_df[["Chosen"]]))
  proportions_mat <- prop.table(counts_mat, margin = 2)


  ## Prepare the plot region
  if (use_dummy_plot) {
    barplot(height = proportions_mat[1, ], space = use_space) # dummy plot, just for par("usr")
    x_coordinates <- par("usr")[c(1, 2)]
    assign("delete_x_coordinates_1", x_coordinates, envir = globalenv())
    x_coordinates <- c(x_coordinates[[1]] - 0.5, x_coordinates[[2]] + 0.5)
  } else {
    if (only_chosen) {
      x_coordinates <- c(-0.34, 1.74)
    } else {
      x_coordinates <- c(0.1, 3.4)
    }
  }
  old_mar <- par(mar = c(4, 4, 5, 0.2) + 0.1)
  plot(1,
       xlim = x_coordinates,
       ylim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       ann  = FALSE,
       type = "n"
       )
  if (show_title) {
    title(categorical_list[["title_text"]], cex.main = categorical_list[["title_cex"]], line = 2.8)
  }

  DrawBarplotGridAndAxis()

  ## Draw the bar plots
  bar_positions <- barplot(height    = proportions_mat[1, ],
                           col       = colors_df[["Pale"]],
                           ylim      = c(0, 1),
                           las       = 1,
                           tcl       = -0.4,
                           border    = NA,
                           space     = use_space,
                           names.arg = rep("", num_subgroups),
                           axes      = FALSE,
                           add       = TRUE,
                           )

  barplot(height    = proportions_mat[2, ],
          offset    = proportions_mat[1, ],
          col       = colors_df[["Dark"]],
          border    = NA,
          space     = use_space,
          names.arg = rep("", num_subgroups),
          axes      = FALSE,
          add       = TRUE
          )

  ## Draw the group labels
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  old_lheight <- par("lheight" = 1.15)
  if (only_chosen) {
    group_labels <- "All chosen guides"
  } else {
    group_labels <- c("Rest", "4sg")
  }
  text(x      = bar_positions[, 1],
       y      = par("usr")[[3]] - (plot_height * 0.055),
       labels = group_labels,
       adj    = c(0.5, 1),
       cex    = 0.9,
       col    =  "black",
       font   = 2,
       xpd    = NA
       )
  par(old_lheight)
  AnnotateBars(bar_positions[, 1], counts_mat, proportions_mat, text_cex = 0.5)

  ## Final steps
  par(old_mar)
  results_list <- list("plot_df" = plot_df, categorical_list[c("title_text", "title_cex")])
  return(invisible(results_list))
}




StackedBarPlot <- function(CRISPR_df, use_column, use_cutoff = NULL, use_dummy_plot = FALSE, show_title = TRUE, only_chosen = FALSE) {

  CRISPR_df <- FilterCRISPRDf(CRISPR_df)

  ## Prepare the data
  categorical_list <- MakeCategoricalList(CRISPR_df, use_column, use_cutoff)
  CRISPR_df[[categorical_list[["use_column"]]]] <- categorical_list[["logical_vec"]]
  if (only_chosen) {
    CRISPR_df <- CRISPR_df[CRISPR_df[["Rank"]] %in% 1:4, ]
  }
  plot_df <- MakeSubgroupsPlotDf(CRISPR_df, categorical_list[["use_column"]])
  colors_df <- SubgroupColorsDf()
  if (only_chosen) {
    colors_df <- colors_df[!(colors_df[["Lightened"]]), -2]
    group_column <- "Group"
  } else {
    group_column <- "Subgroup"
  }
  counts_mat <- as.matrix(table(factor(plot_df[["Data"]], levels = c(FALSE, TRUE)), plot_df[[group_column]]))
  proportions_mat <- prop.table(counts_mat, margin = 2)
  num_subgroups <- nlevels(plot_df[[group_column]])
  if (only_chosen) {
    spaces_vec <- 0.8
  } else {
    spaces_vec <- rep(0.8, num_subgroups) + (rep_len(c(1, -1), length.out = num_subgroups) * 0.5)
  }
  ## Prepare the plot region
  if (use_dummy_plot) {
    barplot(height = proportions_mat[1, ], space = spaces_vec) # dummy plot, just for par("usr")
    x_coordinates <- par("usr")[c(1, 2)]
    assign("delete_x_coordinates_3", x_coordinates, envir = globalenv())
  } else if (only_chosen && (num_subgroups == 7)) {
    x_coordinates <- c(0.328, 13.072)
  } else if (num_subgroups == 14) {
    x_coordinates <- c(0.344, 26.156)
  } else {
     stop("The barplot x axis limits were hard-coded, and are only applicable for 7 or 14 bars!")
  }
  if (only_chosen) {
    x_coordinates <- c(x_coordinates[[1]] - 0.3, x_coordinates[[2]] + 0.3)
  }
  old_mar <- par(mar = c(4, 4, 5, 3) + 0.1)
  plot(1,
       xlim = x_coordinates,
       ylim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       axes = FALSE,
       ann  = FALSE,
       type = "n"
       )
  if (show_title) {
    title(categorical_list[["title_text"]], cex.main = categorical_list[["title_cex"]], line = 2.8)
  }

  DrawBarplotGridAndAxis()

  ## Draw the bar plots
  bar_positions <- barplot(height    = proportions_mat[1, ],
                           col       = colors_df[["Pale"]],
                           ylim      = c(0, 1),
                           las       = 1,
                           tcl       = -0.4,
                           border    = NA,
                           space     = spaces_vec,
                           names.arg = rep("", num_subgroups),
                           axes      = FALSE,
                           add       = TRUE
                           )

  barplot(height    = proportions_mat[2, ],
          offset    = proportions_mat[1, ],
          col       = colors_df[["Dark"]],
          border    = NA,
          space     = spaces_vec,
          names.arg = rep("", num_subgroups),
          axes      = FALSE,
          add       = TRUE
          )

  DrawSubgroupLabels(plot_df, colors_df, bar_positions[, 1], only_chosen = only_chosen)
  AnnotateBars(bar_positions[, 1], counts_mat, proportions_mat, text_cex = 0.4)

  ## Final steps
  par(old_mar)
  results_list <- c(list("plot_df" = plot_df), categorical_list[c("title_text", "title_cex")])
  return(invisible(results_list))
}



PlotCategoricalData <- function(CRISPR_df) {
  file_name <- "4sg library overview - barplots"
  for (only_chosen in c(TRUE, FALSE)) {
    if (only_chosen) {
      use_file_name <- paste0(file_name, " - only 4sg")
      use_columns <- c("Are_overlapping", categorical_columns)
    } else {
      use_file_name <- file_name
      use_columns <- categorical_columns
    }
    for (make_PDF in c(FALSE, TRUE)) {
      if (make_PDF) {
        pdf(file = file.path(output_plots_directory, paste0(use_file_name, ".pdf")),
            width = pdf_width * 0.85, height = pdf_height
            )
      }
      for (y_column in use_columns) {
        layout(use_layout_mat, widths = layout_widths, heights = layout_heights)
        TwoStackedBarPlot(CRISPR_df, y_column, show_title = FALSE, only_chosen = only_chosen)
        StackedBarPlot_results <- StackedBarPlot(CRISPR_df, y_column, show_title = FALSE, only_chosen = only_chosen)
        OuterTitleForLayout(StackedBarPlot_results[["title_text"]], extra_space = FALSE)
      }
      if (make_PDF) {
        dev.off()
      }
    }
  }
}









# Functions for generating histograms -------------------------------------

DrawHistogram <- function(overview_df, column_name) {
  hist_results <- hist(overview_df[[column_name]], breaks = 100, plot = FALSE)
  y_max <- max(hist_results[["density"]])
  hist(overview_df[[column_name]],
       las      = 1,
       mgp      = c(2.2, 0.45, 0),
       tcl      = -0.3,
       col      = brewer.pal(9, "Blues")[[8]],
       border   = NA,
       breaks   = 100,
       xlim     = c(-0.02, 1),
       ylim     = c(0 - (y_max * 0.02), y_max + (y_max * 0.02)),
       main     = paste0("Aggregate ", numeric_column_labels[[column_name]]),
       cex.main = 1,
       xlab     = "",
       xaxs     = "i",
       yaxs     = "i",
       freq     = FALSE
       )
  box(bty = "l")
  return(invisible(NULL))
}


SharedSubsequencesBarplot <- function(overview_df) {
  longest_subsequence_table <- table(factor(sgRNAs_overview_df[["Longest_subsequence"]], levels = 1:20))
  bar_positions <- barplot(longest_subsequence_table,
                           las       = 1,
                           mgp       = c(3, 0.45, 0),
                           tcl       = -0.3,
                           col       = brewer.pal(9, "Blues")[[8]],
                           main      = "Length of the longest shared subsequence",
                           cex.main  = 1,
                           names.arg = "",
                           border    = NA,
                           space     = 0.4,
                           ylim      = c(0, 10000),
                           xlim      = c(0, 28.4),
                           xaxs      = "i",
                           ylab      = "Count"
                           )
  box(xpd = NA, bty = "l")
  plot_height <- par("usr")[[4]] - par("usr")[[3]]
  text(x      = bar_positions[, 1],
       y      = par("usr")[[3]] - (plot_height * 0.1),
       labels = as.character(1:20),
       cex    = 0.8,
       font   = 2,
       xpd    = NA
       )
  text(x      = bar_positions[, 1],
       y      = par("usr")[[3]] - (plot_height * 0.03),
       labels = longest_subsequence_table,
       cex    = 0.5,
       font   = 2,
       xpd    = NA,
       col    = "gray50"
       )
  text(x      = par("usr")[[1]] + ((par("usr")[[2]] - par("usr")[[1]]) / 2),
       y      = par("usr")[[3]] - (plot_height * 0.2),
       labels = "Number of base pairs",
       xpd    = NA
       )
  return(invisible(NULL))
}




Plot4sgData <- function(overview_df) {
  for (make_PDF in c(TRUE, FALSE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "4sg library overview - histograms - 4sg combination.pdf"),
          width = pdf_width, height = pdf_height * 1.3
          )
    }
    par("oma" = c(0, 1, 0, 1))
    SharedSubsequencesBarplot(overview_df)
    for (column_name in c("GuideScan_specificity", "CRISPOR_3MM_specificity", "CRISPOR_4MM_specificity")) {
      DrawHistogram(overview_df, column_name)
    }
    if (make_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}






# Functions for generating Venn diagrams ----------------------------------

PlotVennDiagrams <- function(CRISPR_df) {
  for (make_PDF in c(FALSE, TRUE)) {
    if (make_PDF) {
      pdf(file = file.path(output_plots_directory, "4sg library overview - Venn diagrams.pdf"),
          width = pdf_width * 1.1, height = pdf_height * 1.05
          )
    }
    CRISPR_df <- FilterCRISPRDf(CRISPR_df)
    are_chosen <- CRISPR_df[["Rank"]] %in% 1:4
    all_sources_fac <- ReformatSourceToFactor(CRISPR_df[["Source"]])

    euler_not_chosen <- PlotVennDiagram(all_sources_fac[!(are_chosen)])
    euler_4sg <- PlotVennDiagram(all_sources_fac[are_chosen])
    empty_plot <- ggplot() + theme_void()
    title_rest <- ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = 'bold("Rest")', parse = TRUE, size = 5)
    title_4sg  <- ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = 'bold("4sg")',  parse = TRUE, size = 5)
    grid.arrange(empty_plot, empty_plot, empty_plot, empty_plot, empty_plot,
                 empty_plot, title_rest, empty_plot, title_4sg, empty_plot,
                 empty_plot, euler_not_chosen, empty_plot, euler_4sg, empty_plot,
                 nrow = 3, widths = c(0.1, 1, 0.1, 1, 0.1), heights = c(0.075, 0.075, 1)
                 )
    if (make_PDF) {
      dev.off()
    }
  }
  return(invisible(NULL))
}


PlotVennDiagram <- function(sources_fac, draw_plot = FALSE) {
  sources_table <- table(sources_fac)
  sources_names <- gsub(", ", "&", names(sources_table), fixed = TRUE)
  sources_table <- as.integer(sources_table)
  names(sources_table) <- sources_names
  eulerr_options("padding" = grid::unit(0.25, "lines"))
  euler_plot <- plot(eulerr::euler(sources_table, shape = "circle"),
                     fill = c(brewer.pal(9, "Blues")[[3]], brewer.pal(9, "Greens")[[3]], brewer.pal(9, "Reds")[[3]]),
                     edges = FALSE,
                     labels = list(cex = 1),
                     quantities = list(font = 2, cex = 0.4)
                     )
  if (draw_plot) {
    print(euler_plot)
  }
  return(invisible(euler_plot))
}





# Functions for producing scatter plots -----------------------------------

ConvertCFDScores <- function(numeric_vec) {
  1 / (1 + (10000 / numeric_vec) - 100)
}


GaussianJitter <- function(numeric_vec) {
  numeric_vec +
  (rnorm(n = length(numeric_vec), mean = 0, sd = 0.05) *
  ((max(numeric_vec, na.rm = TRUE) - min(numeric_vec, na.rm = TRUE)) / 50))
}


ScatterPlot <- function(CRISPR_df,
                        x_column,
                        y_column,
                        identical_axes               = FALSE,
                        mark_diagonal                = identical_axes,
                        convert_CRISPOR_to_GuideScan = FALSE,
                        convert_GuideScan_to_CRISPOR = FALSE,
                        point_cex                    = 0.5,
                        point_alpha                  = 0.6,
                        custom_axis_limits           = NULL,
                        show_title                   = NULL,
                        add_jitter                   = FALSE,
                        make_PNG                     = FALSE
                        ) {

  if (convert_CRISPOR_to_GuideScan) {
    CRISPR_df[["CRISPOR_CFD_specificity"]] <- ConvertCFDScores(CRISPR_df[["CRISPOR_CFD_specificity"]])
  }
  if (convert_GuideScan_to_CRISPOR) {
    CRISPR_df[["GuideScan_specificity"]] <- (100 / (100 + ((1 / CRISPR_df[["GuideScan_specificity"]]) - 1))) * 100
  }

  if (make_PNG) {
    if (!("current_number" %in% ls(envir = globalenv()))) {
      current_number <- 1L
    }
    number_string <- formatC(current_number, width = 2, flag = "0")
    file_name <- paste0("4sg library overview - scatter plots - ",
                        number_string, ") ", gsub(":", " - ", show_title), ".png"
                        )
    png(file = file.path(output_plots_directory, file_name),
        width = 5.75, height = 5.75, units = "in", res = 300
        )
    current_number <- current_number + 1L
    assign("current_number", current_number, envir = globalenv())
  }

  old_par <- par(mar = rep.int(5, 4))
  x_vec <- CRISPR_df[[x_column]]
  y_vec <- CRISPR_df[[y_column]]

  if (add_jitter) {
    x_vec <- GaussianJitter(x_vec)
    y_vec <- GaussianJitter(y_vec)
  }

  if (identical_axes) {
    if (is.null(custom_axis_limits)) {
      x_axis_limits <- GetAxisLimits(x_vec, provide_other_limits = FALSE)
      stopifnot(identical(x_axis_limits, GetAxisLimits(y_vec, provide_other_limits = FALSE)))
      if (is.null(x_axis_limits)) {
        are_not_NA <- (!(is.na(x_vec))) & (!(is.na(y_vec)))
        x_axis_limits <- GetAxisLimits(c(x_vec[are_not_NA], y_vec[are_not_NA]), provide_other_limits = TRUE)
      }
    } else {
      x_axis_limits <- custom_axis_limits
    }
    y_axis_limits <- x_axis_limits
  } else {
    x_axis_limits <- GetAxisLimits(x_vec, provide_other_limits = TRUE)
    y_axis_limits <- GetAxisLimits(y_vec, provide_other_limits = TRUE)
  }

  plot(x_vec,
       y_vec,
       las  = 1,
       mgp  = c(2.7, 0.5, 0),
       tcl  = -0.4,
       type = "n",
       xlab = numeric_column_labels[[x_column]],
       ylab = numeric_column_labels[[y_column]],
       xlim = x_axis_limits,
       ylim = y_axis_limits,
       xaxs = "i",
       yaxs = "i"
       )

  assign("delete_mark_diagonal", mark_diagonal, envir = globalenv())
  if (mark_diagonal) {
    abline(a = 0, b = 1, col = "gray88", lwd = 0.5)
  }

  line_color <- "gray88"
  line_type <- "solid"
  line_width <- 0.75
  if (x_column %in% c("GuideScan_specificity", "CRISPOR_3MM_specificity")) {
    abline(v = 0.2, col = line_color, lty = line_type, lwd = line_width)
  }
  if (x_column %in% c("CRISPOR_CFD_specificity")) {
    abline(v = 80, col = line_color, lty = line_type, lwd = line_width)
  }
  if (y_column %in% c("GuideScan_specificity", "CRISPOR_3MM_specificity")) {
    abline(h = 0.2, col = line_color, lty = line_type, lwd = line_width)
  }
  if (y_column %in% c("CRISPOR_CFD_specificity")) {
    abline(h = 80, col = line_color, lty = line_type, lwd = line_width)
  }
  box()

  alpha_hex <- substr(rgb(1, 1, 1, point_alpha), 8, 9)
  points(x_vec,
         y_vec,
         col = paste0(brewer.pal(9, "Blues")[[7]], alpha_hex),
         pch = 16,
         cex = point_cex
         )
  box()
  if (!(is.null(show_title))) {
    title(show_title, cex.main = par("cex") * 0.9)
  }
  assign("delete_x_vec", x_vec, envir = globalenv())
  assign("delete_y_vec", y_vec, envir = globalenv())
  par(old_par)

  if (make_PNG) {
    dev.off()
  }

  return(invisible(NULL))
}





MakeScatterPlots <- function(CRISPR_df) {

  for (i in 1:3) {

    if (i == 1) {
      make_PDF <- FALSE
      make_PNG <- FALSE
    } else if (i == 2) {
      make_PDF <- FALSE
      make_PNG <- TRUE
    } else {
      make_PDF <- TRUE
      make_PNG <- FALSE
    }
    for (only_selected in c(TRUE, FALSE)) {

      if ((only_selected) && !(make_PDF)) {
        next
      }

      if (make_PDF) {
        if (only_selected) {
          append_to_filename <- " - selected"
        } else {
          append_to_filename <- ""
        }
        plot_dimensions <- 5.75
        pdf(file = file.path(output_plots_directory, paste0("4sg library overview - scatterplots", append_to_filename, ".pdf")),
            width = plot_dimensions, height = plot_dimensions
            )
      }

      ScatterPlot(CRISPR_df, "GuideScan_specificity", "GuideScan_efficiency",
                  point_alpha = 0.2, point_cex = 0.4,
                  show_title = "Efficacy vs. specificity (GuideScan)",
                  make_PNG = make_PNG
                  )
      if (!(only_selected)) {
        ScatterPlot(CRISPR_df, "CRISPOR_3MM_specificity", "CRISPOR_Doench_efficacy",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "Efficacy vs. specificity (CRISPOR)",
                    make_PNG = make_PNG
                    )
        ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_Doench_efficacy",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "Efficacy vs. specificity (CRISPOR, up to 4MM)",
                    make_PNG = make_PNG
                    )
        ScatterPlot(CRISPR_df, "CRISPOR_Doench_efficacy", "GuideScan_efficiency",
                    point_alpha = 0.2, point_cex = 0.4, add_jitter = TRUE,
                    show_title = "Original vs. updated Doench efficacy scores",
                    make_PNG = make_PNG
                    )
      }

      ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "GuideScan_specificity",
                  point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                  show_title = "Specificity score vs. number of off-targets (GuideScan)",
                  make_PNG = make_PNG
                  )

      if (!(only_selected)) {
        ScatterPlot(CRISPR_df, "CRISPOR_Num_2or3MM", "CRISPOR_3MM_specificity",
                    point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                    show_title = "Specificity score vs. number of off-targets (CRISPOR)",
                    make_PNG = make_PNG
                    )
        if (FALSE) {
          ScatterPlot(CRISPR_df, "GuideScan_specificity", "CRISPOR_CFD_specificity",
                      point_alpha = 0.2, point_cex = 0.4, convert_CRISPOR_to_GuideScan = TRUE,
                      make_PNG = make_PNG
                      )
          ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_CFD_specificity", # This is just for confirmation
                      convert_GuideScan_to_CRISPOR = TRUE, point_alpha = 0.2, point_cex = 0.4,
                      make_PNG = make_PNG
                      )
        }
      }

      ScatterPlot(CRISPR_df, "CRISPOR_Num_2or3MM", "GuideScan_Num_2or3MM",
                  point_alpha = 0.5, point_cex = 0.2, identical_axes = TRUE,
                  show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR",
                  make_PNG = make_PNG
                  )
      ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "CRISPOR_Num_2or3MM",
                  point_alpha = 0.5, point_cex = 0.2,
                  identical_axes = TRUE, custom_axis_limits = c(0, 1000),
                  show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in)",
                  make_PNG = make_PNG
                  )
      ScatterPlot(CRISPR_df, "GuideScan_Num_2or3MM", "CRISPOR_Num_2or3MM",
                  point_alpha = 0.5, point_cex = 0.2,
                  identical_axes = TRUE, custom_axis_limits = c(0, 200),
                  show_title = "Off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in more)",
                  make_PNG = make_PNG
                  )

      if (!(only_selected)) {
        ScatterPlot(CRISPR_df, "CRISPOR_Num_2MM", "GuideScan_Num_2MM",
                    point_alpha = 0.5, point_cex = 0.2, identical_axes = TRUE,
                    show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR",
                    make_PNG = make_PNG
                    )
        ScatterPlot(CRISPR_df, "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
                    point_alpha = 0.5, point_cex = 0.2,
                    identical_axes = TRUE,
                    custom_axis_limits = c(0, 300),
                    show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in)",
                    make_PNG = make_PNG
                    )
        ScatterPlot(CRISPR_df, "GuideScan_Num_2MM", "CRISPOR_Num_2MM",
                    point_alpha = 0.5, point_cex = 0.2,
                    identical_axes = TRUE, custom_axis_limits = c(0, 50),
                    show_title = "2MM off-target sites \u2013 GuideScan vs. CRISPOR (zoomed in more)",
                    add_jitter = TRUE,
                    make_PNG = make_PNG
                    )
        ScatterPlot(CRISPR_df, "CRISPOR_CFD_specificity", "GuideScan_specificity",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title  = "Specificity \u2013 GuideScan score vs. original CRISPOR CFD score",
                    make_PNG = make_PNG
                    )
      }

      ScatterPlot(CRISPR_df, "CRISPOR_3MM_specificity", "GuideScan_specificity",
                  point_alpha = 0.2, point_cex = 0.4, identical_axes = TRUE,
                  show_title = "Specificity score \u2013 GuideScan vs. CRISPOR",
                  make_PNG = make_PNG
                  )

      ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "GuideScan_specificity",
                  point_alpha = 0.2, point_cex = 0.4,
                  show_title = "Specificity score \u2013 GuideScan vs. CRISPOR (4MM)",
                  make_PNG = make_PNG
                  )

      if (!(only_selected)) {
        ScatterPlot(CRISPR_df, "CRISPOR_4MM_specificity", "CRISPOR_3MM_specificity",
                    point_alpha = 0.2, point_cex = 0.4,
                    show_title = "CRISPOR specificity scores: 3MM vs. 4MM",
                    make_PNG = make_PNG
                    )
      }
      if (make_PDF) {
        dev.off()
      }
    }
  }
}









# Functions for generating multi-plot layouts -----------------------------

MakeEmptyPlot <- function() {
  plot(1, xlim = c(0, 1), ylim = c(0, 1), type = "n", ann = FALSE, axes = FALSE)
}


OuterTitleForLayout <- function(title_text, title_cex = 1, extra_space = FALSE) {
  old_mar <- par("mar" = rep(0, 4))
  if (extra_space) {
    title_y <- -0.4
  } else {
    title_y <- -0.7
  }
  MakeEmptyPlot()
  text(x      = 0.5,
       y      = title_y,
       labels = title_text,
       cex    = title_cex,
       font   = 2,
       xpd    = NA
       )
  MakeEmptyPlot()
  par(old_mar)
}






