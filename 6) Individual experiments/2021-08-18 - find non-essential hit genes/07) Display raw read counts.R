### 2nd December 2021 ###



# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("beeswarm")
library("edgeR")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2021-08-18 - find non-essential hit genes")
file_input_directory  <- file.path(file_directory, "1) Input")
screen_data_directory <- file.path(file_input_directory, "Pooled screen", "RawReadsCounts")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")
plots_directory       <- file.path(file_output_directory, "Plots", "8) Raw counts")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "05) Rank hit genes.RData"))
load(file.path(R_objects_directory, "06) Show differential enrichment or depletion.RData"))




# Define functions --------------------------------------------------------

MakeEmptyPlot <- function() {
  plot(1, type = "n", axes = FALSE, ann = FALSE,
       xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i"
       )
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



Darken <- function(color, factor = 1.4) {
  # from https://gist.github.com/Jfortin1/72ef064469d1703c6b30
    col <- col2rgb(color)
    col <- col / factor
    col <- rgb(t(col), maxColorValue = 255 )
    col
}



DataForGene <- function(gene_symbol, counts_df = norm_by_NT_df) {
  stopifnot("raw_df" %in% ls(envir = globalenv()))
  are_this_gene <- counts_df[, "Gene_symbol"] == gene_symbol
  if (!(any(are_this_gene))) {
    stop(paste0("The gene symbol '", gene_symbol, "' was not found!"))
  }
  data_columns <- grep("File", names(counts_df), fixed = TRUE, value = TRUE)
  data_mat <- as.matrix(counts_df[are_this_gene, data_columns])
  row.names(data_mat) <- NULL
  colnames(data_mat) <- vapply(strsplit(data_columns, "_", fixed = TRUE),
                               function(x) paste0(x[[2]], "_", x[[3]]),
                               ""
                               )
  return(data_mat)
}



CountsForSg <- function(gene_symbol, sg_number, plot_inhibitor = TRUE) {

  ## Prepare the data
  numeric_vec <- DataForGene(gene_symbol)[sg_number, ]


  ## Define the color scheme
  baseline_color <- brewer.pal(9, "Greys")[[6]]
  DMSO_color     <- brewer.pal(9, "Blues")[[5]]
  YM_color       <- brewer.pal(9, "Reds")[[7]]
  NT_color       <- brewer.pal(9, "Greens")[[6]]
  PIKFYVE_color  <- brewer.pal(9, "Purples")[[8]]


  ## Prepare the groupings
  if (plot_inhibitor) {
    numeric_vec <- numeric_vec[1:8]
    groups_fac <- factor(rep(1:4, 2))
    timepoints <- c(1, 3, 5, 5)
    group_colors <- c(baseline_color, DMSO_color, DMSO_color, YM_color)
    group_labels <- c("B1", "D3", "D5", "YM5")
    have_replicates <- TRUE
  } else {
    numeric_vec <- numeric_vec[c(9, 10, 12, 11, 14, 13)]
    groups_fac <- factor(1:6)
    timepoints <- c(4, 8, 10, 10, 14, 14)
    group_colors <- c(baseline_color, baseline_color,
                      NT_color, PIKFYVE_color,
                      NT_color, PIKFYVE_color
                      )
    group_labels <- c("B4", "B8", "C10", "P10", "C14", "P14")
    have_replicates <- FALSE
  }

  groups_list <- split(levels(groups_fac), timepoints)
  group_means <- tapply(numeric_vec, groups_fac, mean)


  ## Prepare the final colors
  if (have_replicates) {
    use_fraction_pale <- 0.4
  } else {
    use_fraction_pale <- 0.4
  }
  bar_colors <- vapply(group_colors, Palify, "", fraction_pale = use_fraction_pale)


  ## Determine group positions
  small_gap <- 0.75
  medium_gap <- 1
  large_gap <- 1
  if (all(lengths(groups_list) == 1)) {
    gaps_vec <- rep(medium_gap, length(groups_list))
    are_first <- rep(TRUE, length(groups_list))
  } else {
    are_first <- unlist(lapply(groups_list, function(x) {
      c(TRUE, rep(FALSE, length(x) - 1))
    }), use.names = FALSE)
    gaps_vec <- ifelse(are_first, large_gap, small_gap)
  }
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_position <- 1
  group_positions <- start_position + cumsum(gaps_vec)

  width <- 0.6 #2 / 3
  final_width <- width * ((max(group_positions) - min(group_positions)) / (nlevels(groups_fac) - 1))
  side_gap <- final_width
  group_limits <- c(group_positions[[1]] - side_gap,
                    group_positions[[length(group_positions)]] + side_gap
                    )


  ## Prepare the data axis
  use_numeric_limits <- c(0, max(numeric_vec) * 1.00)
  numeric_axis_pos <- pretty(use_numeric_limits)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  numeric_axis_labels <- format(numeric_axis_pos)


  ## Draw the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )
  mtext(group_labels, at = group_positions, side = 1, line = 0.5, cex = par("cex"))


  ## Draw the bars
  rect(xleft   = group_positions - (final_width / 2),
       xright  = group_positions + (final_width / 2),
       ybottom = 0,
       ytop    = group_means,
       col     = bar_colors,
       border  = NA,
       xpd     = NA
       )


  ## Draw the points
  if (have_replicates) {
    beeswarm_df <- beeswarm(numeric_vec ~ groups_fac,
                            at       = group_positions,
                            col      = vapply(group_colors, Darken, "", factor = 1.2),
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
  }

  axis(2,
       at       = numeric_axis_pos,
       labels   = numeric_axis_labels,
       mgp      = c(3, 0.54, 0),
       gap.axis = 0,
       tcl      = -0.4,
       las      = 1,
       lwd      = par("lwd"),
       cex      = par("cex"),
       xpd      = NA
       )

  box(bty = "l", lwd = par("lwd"))

  return(invisible(NULL))
}



VerticalCountsForGene <- function(gene_symbol) {

  ## Prepare the plot layout
  layout_mat <- matrix(nrow = 9, ncol = 5)

  layout_mat[3, ] <- 5L
  layout_mat[5, ] <- 6L
  layout_mat[7, ] <- 7L

  layout_mat[, 1] <- 3L
  layout_mat[, 5] <- 4L
  layout_mat[1, ] <- 1L
  layout_mat[9, ] <- 2L

  for (row_number in seq_len(nrow(layout_mat))) {
    for (column_number in seq_len(ncol(layout_mat))) {
      if (is.na(layout_mat[row_number, column_number])) {
        layout_mat[row_number, column_number] <- max(layout_mat, na.rm = TRUE) + 1L
      }
    }
  }
  vertical_space <- 0.5
  layout(layout_mat,
         widths = c(0.4, 1, 0.35, 1, 0.15),
         heights = c(0.7, 1, vertical_space, 1, vertical_space, 1, vertical_space, 1, 0.3),
         )


  ## Activate the plot layout
  par(mar = rep(0, 4))
  title_text <- as.expression(bquote(bolditalic(.(gene_symbol)) * " \u2013 sgRNA counts"))
  MakeEmptyPlot()
  text(x = 0.5, y = 0.6, labels = title_text, cex = par("cex") * 2)

  for (i in 1:6) {
    MakeEmptyPlot()
  }

  ## Plot the data
  for (i in 1:4) {
    CountsForSg(gene_symbol, i)
    mtext(paste0("sgRNA ", i), side = 2, line = 3, cex = par("cex"))
    MakeEmptyPlot()
    CountsForSg(gene_symbol, i, plot_inhibitor = FALSE)
  }
  return(invisible(NULL))
}




HorizontalCountsForGene <- function(gene_symbol) {

  ## Prepare the plot layout
  layout_mat <- matrix(nrow = 5, ncol = 9)

  layout_mat[3, ] <- 5L
  layout_mat[, 1] <- 3L
  layout_mat[, 9] <- 4L
  layout_mat[1, ] <- 1L
  layout_mat[5, ] <- 2L

  for (row_number in seq_len(nrow(layout_mat))) {
    for (column_number in seq_len(ncol(layout_mat))) {
      if (is.na(layout_mat[row_number, column_number])) {
        layout_mat[row_number, column_number] <- max(layout_mat, na.rm = TRUE) + 1L
      }
    }
  }

  horizontal_space <- 0.4
  layout(layout_mat,
         widths  = c(0.4, 1, horizontal_space, 1, horizontal_space, 1, horizontal_space, 1, 0.15),
         heights = c(0.8, 1, 0.5, 1, 0.3),
         )


  ## Activate the plot layout
  par(mar = rep(0, 4))

  title_text <- as.expression(bquote(bolditalic(.(gene_symbol)) * " \u2013 sgRNA counts"))
  MakeEmptyPlot()
  text(x = 0.5, y = 0.65, labels = title_text, cex = par("cex") * 2)

  for (i in 1:4) {
    MakeEmptyPlot()
  }

  ## Plot the data
  for (show_inhibitor in c(TRUE, FALSE)) {
    for (i in 1:4) {
      CountsForSg(gene_symbol, i, plot_inhibitor = show_inhibitor)
      if (show_inhibitor) {
        mtext(paste0("sgRNA ", i), side = 3, line = 0.8, cex = par("cex"))
      }
      if (i == 1) {
        mtext("Read count", side = 2, line = 3, cex = par("cex"))

      }
      if (i != 4) {
        MakeEmptyPlot()
      }
    }
  }

  return(invisible(NULL))
}





ExportPNGsForGene <- function(gene_symbol,
                              gene_number = NULL,
                              sub_folder = "c) Additional genes"
                              ) {

  file_name <- gene_symbol
  if (!(is.null(gene_number))) {
    file_name <- paste0(gene_number, ") ", file_name)
  }

  png(file = file.path(plots_directory,
                       "Horizontal",
                       sub_folder,
                       paste0(file_name, ".png")
                       ),
      width = 9.5, height = 4.2,
      units = "in", res = 600
      )
  HorizontalCountsForGene(gene_symbol)
  dev.off()

  png(file = file.path(plots_directory,
                       "Vertical",
                       sub_folder,
                       paste0(file_name, ".png")
                       ),
      width = 4.7, height = 7,
      units = "in", res = 600
      )
  VerticalCountsForGene(gene_symbol)
  dev.off()

  return(invisible(NULL))
}





CountsHistogram <- function(use_counts_df,
                            file_name,
                            top_title = "",
                            only_NT = TRUE,
                            set_axis_limits = TRUE,
                            truncate_outliers = TRUE
                            ) {


  ## Prepare the data
  sample_columns <- grep("File", names(use_counts_df), value = TRUE, fixed = TRUE)
  use_counts_mat <- as.matrix(use_counts_df[, sample_columns])
  sample_names <- vapply(strsplit(sample_columns, "_", fixed = TRUE),
                         function(x) paste0(x[-1], collapse = "_"),
                         ""
                         )
  x_axis_limits <- NULL
  y_axis_limits <- NULL
  if (only_NT) {
    are_NT <- grepl("non-targeting", use_counts_df[, "Gene_symbol"], fixed = TRUE)
    use_counts_mat <- use_counts_mat[are_NT, ]
    rownames(use_counts_mat) <- NULL
    if (set_axis_limits) {
      x_axis_limits <- c(0, 4050)
      y_axis_limits <- c(-2, 130)
    }
  } else {
    if (set_axis_limits) {
      x_axis_limits <- c(0, 4050)
      y_axis_limits <- c(-100, 10000)
    }
  }

  if (truncate_outliers) {
    use_counts_mat[use_counts_mat > 4000] <- 4000
  }

  ## Set up the graphics device

  file_path <- file.path(file_output_directory,
                         "Plots", "7) Quality control",
                         paste0(file_name, ".png")
                         )
  png(file_path, height = 13, width = 6, res = 900, units = "in")


  ## Set up the plot layout
  layout_mat <- cbind(1:8, 9:16)
  layout_mat <- rbind(c(1, 1), layout_mat + 1L)
  layout(layout_mat,
         heights = c(0.5, rep(1, 8))
         )


  ## Display an overall title at the top
  par(mar = c(0, 0, 0, 0))
  MakeEmptyPlot()
  text(x      = 0.5,
       y      = 0.5,
       labels = top_title,
       font   = 2,
       cex    = par("cex") * 2,
       xpd    = NA
       )


  par(mar = c(3.3, 4, 2, 2.5))

  num_samples <- ncol(use_counts_mat)

  for (i in seq_len(num_samples)) {

    trimmed_mean <- mean(use_counts_mat[, i], trim = 0.05)
    legend_text <- paste0("95% trimmed mean\n= ",
                          round(trimmed_mean, digits = 1)
                          )

    if (is.null(x_axis_limits)) {
      use_breaks <- 200
    } else {
      use_breaks <- seq(from = x_axis_limits[[1]], to = x_axis_limits[[2]], by = 50)
    }

    hist_args <- list(
      x = use_counts_mat[, i],
      freq     = TRUE,
      xlab     = "",
      ylab     = "",
      breaks   = use_breaks,
      xlim     = x_axis_limits,
      ylim     = y_axis_limits,
      xaxs     = "i",
      yaxs     = "i",
      main     = "",
      mgp      = c(2.7, 1, 0),
      col      = brewer.pal(9, "Greys")[[7]],
      border   = NA,
      axes     = FALSE
    )
    if (is.null(x_axis_limits)) {
      hist_args[["xlim"]] <- NULL
    }
    hist_results <- do.call(hist, hist_args)

    axis(2,
         las = 1,
         mgp = c(3, 0.55, 0),
         tcl = -0.4,
         lwd = par("lwd"),
         cex = par("cex")
         )

    axis(1,
         las = 1,
         mgp = c(3, 0.55, 0),
         tcl = -0.4,
         lwd = par("lwd"),
         cex = par("cex"),
         gap.axis = 0.5
         )

    abline(v   = trimmed_mean,
           col = brewer.pal(9, "Oranges")[[5]],
           lwd = 1.5
           )

    text(x      = grconvertX(0.25, from = "npc", to = "user"),
         y      = grconvertY(1.2, from = "npc", to = "user"),
         adj    = c(0.5, 1),
         labels = sample_names[[i]],
         xpd    = NA,
         cex    = par("cex") * 1.6
         )
    text(x      = grconvertX(0.5, from = "npc", to = "user"),
         y      = grconvertY(0.8, from = "npc", to = "user"),
         adj    = c(0, 1),
         labels = legend_text,
         xpd    = NA,
         cex    = par("cex") * 1.2
         )
    box(bty = "l", lwd = par("lwd"))

    if (i %in% c(8, num_samples)) {
      mtext("Read count", side = 1, line = 2, cex = par("cex"), xpd = NA)
    }

  }

  for (i in seq_len(16 - num_samples)) {
    MakeEmptyPlot()
  }

  dev.off()

  return(invisible(NULL))
}





# Read in data ------------------------------------------------------------

count_files <- list.files(screen_data_directory)
df_list <- lapply(count_files, function(x) read.delim(file.path(screen_data_directory, x), stringsAsFactors = FALSE))




# Compile raw data --------------------------------------------------------

stopifnot(length(unique(lapply(df_list, function(x) x[, 1:5]))) == 1)

counts_mat <- do.call(cbind, lapply(df_list, function(x) x[, 6]))
file_names <- sub("-result.txt", "", count_files, fixed = TRUE)
file_names <- sub("o25448_1_", "File", file_names, fixed = TRUE)
file_names <- gsub("-", "_", file_names, fixed = TRUE)
colnames(counts_mat) <- file_names

raw_df <- data.frame(df_list[[1]][, 1:5], counts_mat, stringsAsFactors = FALSE)





# Make the gene names consistent ------------------------------------------

matches_vec <- match(raw_df[, "Sequence"], D5_v_YM5_df[, "Gene_Name"])
stopifnot(!(anyNA(matches_vec)))

raw_df[, "Gene_symbol"] <- D5_v_YM5_df[matches_vec, "Gene_symbol"]




# Compare read counts across samples --------------------------------------

are_NT <- grepl("non-targeting", raw_df[, "Gene_symbol"], fixed = TRUE)

cbind("NT"    = apply(counts_mat, 2, function(x) sum(x[are_NT])),
      "Total" = colSums(counts_mat)
      )

dgList <- edgeR::calcNormFactors(counts_mat, method = "TMM")
dgList
cpm(dgList)




# Normalize raw data using non-targeting control sgRNAs -------------------

norm_by_NT_df <- raw_df[, c("Sequence", "Gene_symbol")]
for (column_name in grep("File", names(raw_df), value = TRUE, fixed = TRUE)) {
  trimmed_mean_NT <- mean(raw_df[are_NT, column_name], trim = 0.05)
  norm_by_NT_df[[column_name]] <- raw_df[, column_name] / trimmed_mean_NT * 1000
}




# Plot counts data for individual genes -----------------------------------

HorizontalCountsForGene("Ostm1")


for (i in seq_along(top_10_HS)) {
  ExportPNGsForGene(gene_symbol = top_10_HS[[i]], gene_number = i,
                    sub_folder = "a) Top 10 - overall (by hit strength)"
                    )
}

for (i in seq_along(top_10_NE)) {
  ExportPNGsForGene(gene_symbol = top_10_NE[[i]], gene_number = i,
                    sub_folder = "b) Top 10 - filtered by essentiality"
                    )
}


for (additional_gene in c("Prnp", "Atp6v1b2", paste0("Atp1a", 1:2))) {
  ExportPNGsForGene(additional_gene)
}




# Compile normalized data -------------------------------------------------

stopifnot(length(unique(list(B1_v_D5_df[, "Gene_Name"], B1_v_YM5_df[, "Gene_Name"], D5_v_YM5_df[, "Gene_Name"]))) == 1)
stopifnot(length(unique(list(B1_v_D5_df[, "Gene_symbol"], B1_v_YM5_df[, "Gene_symbol"], D5_v_YM5_df[, "Gene_symbol"]))) == 1)


sample_names <- c("01-YM2-B1", "05-YM3-B1", "03-YM2-BD5", "07-YM3-BD5")
column_names <- paste0("o25448_1_", sample_names, " [normalized count]")
B1_v_D5_counts_mat <- B1_v_D5_df[, column_names]
name_splits <- strsplit(sample_names, "-", fixed = TRUE)
new_sample_names <- vapply(name_splits, function(x) paste0(x[-1], collapse = "_"), "")
new_sample_names <- paste0("File", sapply(name_splits, "[[", 1), "_",
                           "B1-v-D5_", new_sample_names
                           )
colnames(B1_v_D5_counts_mat) <- new_sample_names



sample_names <- c("01-YM2-B1", "05-YM3-B1", "04-YM2-BY5", "08-YM3-BY5")
column_names <- paste0("o25448_1_", sample_names, " [normalized count]")
B1_v_YM5_counts_mat <- B1_v_YM5_df[, column_names]
name_splits <- strsplit(sample_names, "-", fixed = TRUE)
new_sample_names <- vapply(name_splits, function(x) paste0(x[-1], collapse = "_"), "")
new_sample_names <- paste0("File", sapply(name_splits, "[[", 1), "_",
                           "B1-v-YM5_", new_sample_names
                           )
colnames(B1_v_YM5_counts_mat) <- new_sample_names



sample_names <- c("03-YM2-BD5", "07-YM3-BD5", "04-YM2-BY5", "08-YM3-BY5")
column_names <- paste0("o25448_1_", sample_names, " [normalized count]")
D5_v_YM5_counts_mat <- D5_v_YM5_df[, column_names]
name_splits <- strsplit(sample_names, "-", fixed = TRUE)
new_sample_names <- vapply(name_splits, function(x) paste0(x[-1], collapse = "_"), "")
new_sample_names <- paste0("File", sapply(name_splits, "[[", 1), "_",
                           "D5-v-YM5_", new_sample_names
                           )
colnames(D5_v_YM5_counts_mat) <- new_sample_names



differential_df <- data.frame("Sequence" = B1_v_D5_df[, "Gene_Name"],
                              B1_v_D5_df["Gene_symbol"],
                              B1_v_D5_counts_mat,
                              B1_v_YM5_counts_mat,
                              D5_v_YM5_counts_mat,
                              stringsAsFactors = FALSE,
                              check.names = FALSE
                              )




# Create histograms of counts data ----------------------------------------

CountsHistogram(raw_df,
                file_name = "Histograms of original read counts for non-targeting control sgRNAs",
                top_title = "Raw read count distribution of non-targeting sgRNAs (n = 1000)"
                )

CountsHistogram(raw_df,
                file_name = "Histograms of original read counts for all sgRNAs",
                top_title = paste0("Raw read count distribution of all sgRNAs (n = ", nrow(counts_mat), ")"),
                only_NT   = FALSE
                )


CountsHistogram(differential_df,
                file_name = "Histograms of normalized read counts for non-targeting control sgRNAs",
                top_title = "Normalized read count distribution of non-targeting sgRNAs (n = 1000)"
                )

CountsHistogram(differential_df,
                file_name = "Histograms of normalized read counts for all sgRNAs",
                top_title = paste0("Normalized read count distribution of all sgRNAs (n = ", nrow(counts_mat), ")"),
                only_NT   = FALSE
                )


CountsHistogram(norm_by_NT_df,
                file_name = "Histograms of NT-normalized read counts for non-targeting control sgRNAs",
                top_title = "NT-normalized read count distribution of non-targeting sgRNAs (n = 1000)"
                )

CountsHistogram(norm_by_NT_df,
                file_name = "Histograms of NT-normalized read counts for all sgRNAs",
                top_title = paste0("NT-normalized read count distribution of all sgRNAs (n = ", nrow(counts_mat), ")"),
                only_NT   = FALSE
                )




