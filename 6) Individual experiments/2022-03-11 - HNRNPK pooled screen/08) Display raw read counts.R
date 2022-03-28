### 12th March 2022 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("RColorBrewer")
library("beeswarm")
library("edgeR")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-03-11 - HNRNPK pooled screen")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")
plots_directory       <- file.path(file_output_directory, "Plots", "8) Raw counts")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Read in data, and rank hit genes.RData"))
load(file.path(R_objects_directory, "06) Rank hit genes.RData"))
load(file.path(R_objects_directory, "07) Show differential enrichment or depletion.RData"))





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



DataForGene <- function(gene_symbol, counts_df = TMM_df) {
  stopifnot("raw_df" %in% ls(envir = globalenv()))
  are_this_gene <- counts_df[, "GeneSymbol"] == gene_symbol
  if (!(any(are_this_gene))) {
    stop(paste0("The gene symbol '", gene_symbol, "' was not found!"))
  }
  data_columns <- grep("SVS", names(counts_df), fixed = TRUE, value = TRUE)
  data_mat <- as.matrix(counts_df[are_this_gene, data_columns])
  row.names(data_mat) <- NULL
  return(data_mat)
}



CountsForSg <- function(gene_symbol, sg_number) {

  ## Prepare the data
  numeric_vec <- DataForGene(gene_symbol)[sg_number, ]


  ## Define the color scheme
  baseline_color <- brewer.pal(9, "Greys")[[6]]
  control_color  <- brewer.pal(9, "Blues")[[5]]
  HNRNPK_color   <- brewer.pal(9, "Reds")[[7]]

  groups_fac <- factor(rep(1:3, 2))
  group_colors <- c(baseline_color, control_color, HNRNPK_color)
  group_labels <- c("Baseline", "NT", expression(italic("HNRNPK")))
  have_replicates <- TRUE

  group_means <- tapply(numeric_vec, groups_fac, mean)

  ## Prepare the final colors
  if (have_replicates) {
    use_fraction_pale <- 0.4
  } else {
    use_fraction_pale <- 0.4
  }
  bar_colors <- vapply(group_colors, Palify, "", fraction_pale = use_fraction_pale)


  ## Determine group positions
  group_positions <- seq_len(nlevels(groups_fac))

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
           cex = par("cex"),
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



TwoByTwoForGene <- function(gene_symbol) {
  old_par <- par(mfrow = c(2, 2),
                 mar   = c(3, 4, 2.6, 1),
                 oma   = c(0, 1, 3, 1)
                 )
  for (i in 1:4) {
    CountsForSg(gene_symbol, i)
    if (i %in% c(1, 3)) {
      mtext("Read count", side = 2, cex = par("cex"), line = 2.4)
    }
    title(paste0("sg", i), cex.main = 1)
  }
  title(gene_symbol, font.main = 4, outer = TRUE, line = 1)
  par(old_par)
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

  png(filename = file.path(plots_directory,
                           sub_folder,
                           paste0(file_name, ".png")
                           ),
      width = 5, height = 4.2,
      units = "in", res = 600
      )
  TwoByTwoForGene(gene_symbol)
  dev.off()

  return(invisible(NULL))
}


CountsHistogram <- function(use_counts_df,
                            file_name,
                            top_title = "",
                            only_NT = TRUE,
                            set_axis_limits = TRUE,
                            truncate_outliers = TRUE,
                            x_axis_limits = NULL,
                            y_axis_limits = NULL
                            ) {

  assign("delete_use_counts_df", use_counts_df, envir = globalenv())

  ## Prepare the data
  sample_columns <- grep("SVS", names(use_counts_df), value = TRUE, fixed = TRUE)
  use_counts_mat <- as.matrix(use_counts_df[, sample_columns])
  sample_names <- sample_columns
  # sample_names <- sapply(strsplit(sample_columns, "-", fixed = TRUE), "[[", 2)
  if (only_NT) {
    are_NT <- grepl("non-targeting", use_counts_df[, "GeneSymbol"], ignore.case = TRUE)
    use_counts_mat <- use_counts_mat[are_NT, ]
    rownames(use_counts_mat) <- NULL
    if (set_axis_limits) {
      x_axis_limits <- c(0, 6050)
      if (is.null(y_axis_limits)) {
        y_axis_limits <- c(-2, 80)
      }
    }
  } else {
    if (set_axis_limits) {
      x_axis_limits <- c(0, 6050)
      if (is.null(y_axis_limits)) {
        y_axis_limits <- c(-100, 4000)
      }
    }
  }

  if (truncate_outliers) {
    use_counts_mat[use_counts_mat > 6000] <- 6000
  }

  ## Set up the graphics device

  file_path <- file.path(file_output_directory,
                         "Plots", "7) Quality control",
                         paste0(file_name, ".png")
                         )
  png(file_path, height = 13, width = 6, res = 900, units = "in")


  ## Set up the plot layout
  layout_mat <- cbind(1:6, 7:12)
  layout_mat <- rbind(c(1, 1), layout_mat + 1L)
  layout(layout_mat,
         heights = c(0.5, rep(1, 3))
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

  assign("delete_use_counts_mat", use_counts_mat, envir = globalenv())
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
      axes     = FALSE,
      xpd      = NA
    )
    if (is.null(x_axis_limits)) {
      hist_args[["xlim"]] <- NULL
    }
    do.call(hist, hist_args)

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

    if (i %in% c(6, num_samples)) {
      mtext("Read count", side = 1, line = 2, cex = par("cex"), xpd = NA)
    }

  }
  assign("delete_num_samples", num_samples, envir = globalenv())

  for (i in seq_len(12 - num_samples)) {
    MakeEmptyPlot()
  }

  dev.off()

  return(invisible(NULL))
}





# Normalize raw data using non-targeting control sgRNAs -------------------

norm_by_NT_df <- raw_df[, c("Sequence", "GeneSymbol")]
are_NT <- grepl("Non-Targeting", raw_df[, "GeneSymbol"], fixed = TRUE)
for (column_name in grep("SVS", names(raw_df), value = TRUE, fixed = TRUE)) {
  trimmed_mean_NT <- mean(raw_df[are_NT, column_name], trim = 0.05)
  norm_by_NT_df[[column_name]] <- raw_df[, column_name] / trimmed_mean_NT * 1000
}



# Perform TMM normalization of raw data -----------------------------------

counts_object <- DGEList(raw_df[, 6:11])
counts_object <- edgeR::calcNormFactors(counts_object, method = "TMM")
TMM_mat <- edgeR::cpm(counts_object)
TMM_df <- data.frame(raw_df[, 1:5], TMM_mat, stringsAsFactors = FALSE)





# Plot counts data for individual genes -----------------------------------

TwoByTwoForGene("GTF3C6")


for (i in seq_along(top_10_HS)) {
  ExportPNGsForGene(gene_symbol = top_10_HS[[i]], gene_number = i,
                    sub_folder = "a) Top 10 - overall (by hit strength)"
                    )
}


for (i in seq_along(top_10_NE)) {
  ExportPNGsForGene(gene_symbol = top_10_NE[[i]], gene_number = i,
                    sub_folder = "b) Additional 10 (non-essential and outside top 10)"
                    )
}



for (additional_gene in c("PRNP", "CDNF", "HNRNPK")) {
  ExportPNGsForGene(additional_gene)
}



for (i in seq_along(Stefano_top_genes)) {
  ExportPNGsForGene(gene_symbol = Stefano_top_genes[[i]], gene_number = i,
                    sub_folder = "d) Stefano's selection - top genes"
                    )
}


for (i in seq_along(Stefano_non_essential_genes)) {
  ExportPNGsForGene(gene_symbol = Stefano_non_essential_genes[[i]], gene_number = i,
                    sub_folder = "e) Stefano's selection - non-essential"
                    )
}



# Compile normalized data -------------------------------------------------

stopifnot(length(unique(list(base_v_NT_df[, "gene_id"],   base_v_HNRNPK_df[, "gene_id"],   NT_v_HNRNPK_df[, "gene_id"]))) == 1)
stopifnot(length(unique(list(base_v_NT_df[, "Gene_symbol"], base_v_HNRNPK_df[, "Gene_symbol"], NT_v_HNRNPK_df[, "Gene_symbol"]))) == 1)


sample_names <- c("1-SVSI_DAY1", "4-SVSII_DAY1", "2-SVSI_DAY14_NT", "5-SVSII_DAY14_NT")
column_names <- paste0("o25644_1_", sample_names, " [normalized count]")
base_v_NT_counts_mat <- base_v_NT_df[, column_names]


sample_names <- c("1-SVSI_DAY1", "4-SVSII_DAY1", "3-SVSI_DAY14_hnRNPK", "6-SVSII_DAY14_hnRNPK")
column_names <- paste0("o25644_1_", sample_names, " [normalized count]")
base_v_HNRNPK_counts_mat <- base_v_HNRNPK_df[, column_names]


sample_names <- c("2-SVSI_DAY14_NT", "5-SVSII_DAY14_NT", "3-SVSI_DAY14_hnRNPK", "6-SVSII_DAY14_hnRNPK")
column_names <- paste0("o25644_1_", sample_names, " [normalized count]")
NT_v_HNRNPK_counts_mat <- NT_v_HNRNPK_df[, column_names]


differential_df <- data.frame("Sequence" = sapply(strsplit(base_v_NT_df[, "gene_id"], "-", fixed = TRUE), "[[", 2),
                              "GeneSymbol" = base_v_NT_df[, "Gene_symbol"],
                              base_v_NT_counts_mat,
                              base_v_HNRNPK_counts_mat,
                              NT_v_HNRNPK_counts_mat,
                              stringsAsFactors = FALSE,
                              check.names = FALSE
                              )
names(differential_df) <- sub("o25644_1_", "", names(differential_df), fixed = TRUE)
names(differential_df) <- sub(" [normalized count]", "", names(differential_df), fixed = TRUE)






# Create histograms of counts data ----------------------------------------

CountsHistogram(raw_df,
                file_name = "Histograms of original read counts for non-targeting control sgRNAs",
                top_title = "Raw read count distribution of non-targeting sgRNAs (n = 1000)"
                )

CountsHistogram(raw_df,
                file_name = "Histograms of original read counts for all sgRNAs",
                top_title = paste0("Raw read count distribution of all sgRNAs (n = ", nrow(raw_df), ")"),
                only_NT   = FALSE
                )


CountsHistogram(differential_df,
                file_name = "Histograms of normalized read counts for non-targeting control sgRNAs",
                top_title = "Normalized read count distribution of non-targeting sgRNAs (n = 1000)"
                )

CountsHistogram(differential_df,
                file_name = "Histograms of normalized read counts for all sgRNAs",
                top_title = paste0("Normalized read count distribution of all sgRNAs (n = ", nrow(raw_df), ")"),
                only_NT   = FALSE
                )


CountsHistogram(norm_by_NT_df,
                file_name = "Histograms of NT-normalized read counts for non-targeting control sgRNAs",
                top_title = "NT-normalized read count distribution of non-targeting sgRNAs (n = 1000)"
                )

CountsHistogram(norm_by_NT_df,
                file_name = "Histograms of NT-normalized read counts for all sgRNAs",
                top_title = paste0("NT-normalized read count distribution of all sgRNAs (n = ", nrow(raw_df), ")"),
                only_NT   = FALSE,
                y_axis_limits = c(0, 10000)
                )


CountsHistogram(TMM_df,
                file_name = "Histograms of TMM-normalized read counts for non-targeting control sgRNAs",
                top_title = "TMM-normalized read count distribution of non-targeting sgRNAs (n = 1000)",
                x_axis_limits = c(0, 100)
                )

CountsHistogram(TMM_df,
                file_name = "Histograms of TMM-normalized read counts for all sgRNAs",
                top_title = paste0("TMM-normalized read count distribution of all sgRNAs (n = ", nrow(raw_df), ")"),
                only_NT   = FALSE,
                x_axis_limits = c(0, 100)
                )



