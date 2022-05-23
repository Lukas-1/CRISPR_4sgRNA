## 2022-05-23



# Load packages and source code -------------------------------------------

library("RColorBrewer")
library("beeswarm")

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 4sg - first trial")
input_dir   <- file.path(project_dir, "02_input_data")
rdata_dir   <- file.path(project_dir, "03_R_objects")
figures_dir <- file.path(project_dir, "04_output_data", "Figures")
DepMap2020Q2_dir <- file.path(input_dir, "Essential genes", "Used for Nunez et al")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))
load(file.path(rdata_dir, "05_assign_sgRNAs_to_plasmids__counts_df.RData"))




# Read in data ------------------------------------------------------------

essentials_df <- read.table(file.path(DepMap2020Q2_dir, "DepMap_20Q2__common_essentials.csv"),
                            stringsAsFactors = FALSE, sep = "\t", header = TRUE
                            )
non_essentials_df <- read.table(file.path(DepMap2020Q2_dir, "DepMap_20Q2__nonessentials.csv"),
                                stringsAsFactors = FALSE, sep = "\t", header = TRUE
                                )






# Define functions --------------------------------------------------------

GetLog2FC <- function(input_df,
                      allow_switch = TRUE,
                      allow_1MM = TRUE,
                      choose_rep = NULL
                      ) {
  column_prefix <- paste0(if (allow_switch) "MaySwitch" else "NoSwitch", "_",
                          if (allow_1MM) "xMM" else "0MM", "_"
                          )
  column_names <- grep(column_prefix, names(counts_df), fixed = TRUE, value = TRUE)
  counts_mat <- as.matrix(counts_df[, column_names[3:6]])
  if (is.null(choose_rep)) {
    results_vec <- rowMeans(counts_mat[, 3:4]) / rowMeans(counts_mat[, 1:2])
  } else {
    if (choose_rep == 1) {
      use_columns <- c(1, 3)
    } else {
      use_columns <- c(2, 4)
    }
    results_vec <- counts_mat[, use_columns[[2]]] / counts_mat[, use_columns[[1]]]
  }
  results_vec <- log2(results_vec)
  return(results_vec)
}



CountsForPlasmid <- function(plasmid_ID, allow_switch = TRUE, allow_1MM = TRUE) {
  stopifnot("counts_df" %in% ls(envir = globalenv()))
  plasmid_index <- match(plasmid_ID, counts_df[, "Plasmid_ID"])
  if (is.na(plasmid_index)) {
    stop("A plasmid with the ID '", plasmid_ID, "' was not found!")
  }
  column_prefix <- paste0(if (allow_switch) "MaySwitch" else "NoSwitch", "_",
                          if (allow_1MM) "xMM" else "0MM", "_"
                          )
  column_names <- grep(column_prefix, names(counts_df), fixed = TRUE, value = TRUE)
  group_names <- c("Pre-screen", "T0", "T12")
  groups_fac <- factor(rep(group_names, each = 2), levels = group_names)
  counts_mat <- as.matrix(counts_df[, column_names])
  counts_vec <- counts_mat[plasmid_index, ]
  PlotCounts(counts_vec, groups_fac)
  title(plasmid_ID, font.main = 4)
}


PlotCounts <- function(numeric_vec,
                       groups_fac,
                       group_colors = c("Greys", "Blues", "Purples"),
                       group_labels = NULL
                       ) {

  group_means <- tapply(numeric_vec, groups_fac, mean)

  if (is.null(group_labels)) {
    group_labels <- levels(groups_fac)
  }

  ## Determine group positions
  num_groups <- nlevels(groups_fac)
  group_positions <- seq_len(num_groups)
  side_gap <- 0.5
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )
  bar_width <- 0.6 #2 / 3


  ## Prepare the data axis
  use_y_limits <- range(numeric_vec)
  y_space <- (use_y_limits[[2]] - use_y_limits[[1]]) * 0.02
  use_y_limits <- c(use_y_limits[[1]] - y_space, use_y_limits[[2]] + y_space)
  if (use_y_limits[[1]] > 0) {
    use_y_limits[[1]] <- 0
  }

  ## Prepare the colors
  bar_colors   <- vapply(group_colors, function(x) brewer.pal(9, x)[[5]], "")
  point_colors <- vapply(group_colors, function(x) brewer.pal(9, x)[[9]], "")

  ## Set up the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = use_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )

  mtext(group_labels, at = group_positions, side = 1, line = 0.5, cex = par("cex"))


  ## Draw the bars
  rect(xleft   = group_positions - (bar_width / 2),
       xright  = group_positions + (bar_width / 2),
       ybottom = 0,
       ytop    = group_means,
       col     = bar_colors,
       border  = NA,
       xpd     = NA
       )

  ## Draw the points
  beeswarm_df <- beeswarm(numeric_vec ~ groups_fac,
                          at       = group_positions,
                          col      = point_colors,
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

  axis(2,
       mgp      = c(3, 0.54, 0),
       gap.axis = 0,
       tcl      = -0.4,
       las      = 1,
       lwd      = par("lwd"),
       cex      = par("cex"),
       xpd      = NA
       )

  box(bty = "l", lwd = par("lwd"))
}





# Tidy data on essential genes --------------------------------------------

TidyDepMap <- function(input_df) {
  splits_list <- strsplit(input_df[, 1], " (", fixed = TRUE)
  results_df <- data.frame(
    "Gene_symbol" = sapply(splits_list, "[[", 1),
    "Entrez_ID"   = sub(")", "", sapply(splits_list, "[[", 2), fixed = TRUE),
    stringsAsFactors = FALSE
  )
  results_df[, "Entrez_ID"] <- as.integer(results_df[, "Entrez_ID"])
  return(results_df)
}

essentials_df <- TidyDepMap(essentials_df)
non_essentials_df <- TidyDepMap(non_essentials_df)



# Check for plasmids that are not represented in the data -----------------

table(counts_df[, "Sum_MaySwitch_xMM"] == 0)
table(counts_df[, "Sum_MaySwitch_0MM"] == 0)
table(counts_df[, "Sum_NoSwitch_xMM"] == 0)
table(counts_df[, "Sum_NoSwitch_0MM"] == 0)



# Calculate log2FCs -------------------------------------------------------

counts_df[, "Log2FC_MaySwitch_xMM"] <- GetLog2FC(counts_df, allow_switch = TRUE,  allow_1MM = TRUE)
counts_df[, "Log2FC_MaySwitch_0MM"] <- GetLog2FC(counts_df, allow_switch = TRUE,  allow_1MM = FALSE)
counts_df[, "Log2FC_NoSwitch_xMM"]  <- GetLog2FC(counts_df, allow_switch = FALSE, allow_1MM = TRUE)
counts_df[, "Log2FC_NoSwitch_0MM"]  <- GetLog2FC(counts_df, allow_switch = FALSE, allow_1MM = FALSE)

for (ri in 1:2) {
  counts_df[, paste0("Log2FC_MaySwitch_xMM_R", ri)] <- GetLog2FC(counts_df, allow_switch = TRUE,  allow_1MM = TRUE,  choose_rep = ri)
}
for (ri in 1:2) {
  counts_df[, paste0("Log2FC_MaySwitch_0MM_R", ri)] <- GetLog2FC(counts_df, allow_switch = TRUE,  allow_1MM = FALSE, choose_rep = ri)
}
for (ri in 1:2) {
  counts_df[, paste0("Log2FC_NoSwitch_xMM_R", ri)]  <- GetLog2FC(counts_df, allow_switch = FALSE, allow_1MM = TRUE,  choose_rep = ri)
}
for (ri in 1:2) {
  counts_df[, paste0("Log2FC_NoSwitch_0MM_R", ri)]  <- GetLog2FC(counts_df, allow_switch = FALSE, allow_1MM = FALSE, choose_rep = ri)
}




GetAvailableGenes <- function(entrezs_vec, count_column = "NoSwitch_xMM", min_count = 200) {

  stopifnot(all(c("counts_df", "CRISPRoff_df") %in% ls(envir = globalenv())))
  stopifnot(is.numeric(entrezs_vec))
  stopifnot(!(any(duplicated(entrezs_vec))))
  # stopifnot(identical(CRISPRoff_df[, "Entrez_ID"], counts_df[, "Entrez_ID"]))

  matches_vec <- match(entrezs_vec, CRISPRoff_df[, "Entrez_ID"])
  message(paste0("\nOf the ", length(entrezs_vec), " Entrez gene IDs, ",
                 sum(is.na(matches_vec)), " were not available in the library.\n"
                 ))
  matches_vec <- matches_vec[!(is.na(matches_vec))]
  entrezs_vec <- entrezs_vec[!(is.na(matches_vec))]

  have_multiple_plasmids <- CRISPRoff_df[, "Num_plasmids_for_Entrez"][matches_vec]
  message(paste0(sum(is.na(matches_vec)), " genes were represented by multiple ",
                 "plasmids (targeting multiple TSSs).\nOnly the first plasmid ",
                 "was chosen for each gene.\n"
                 ))

  # sample_names <- c("T0_R1", "T0_R2")
  sample_names <- c("Tbefore_R1", "Tbefore_R2")
  column_names <- paste0(count_column, "_", sample_names)
  counts_mat <- as.matrix(counts_df[, column_names])
  total_counts_vec <- rowSums(counts_mat[matches_vec, ])
  too_few_counts <- total_counts_vec < min_count

  message(paste0(sum(too_few_counts), " plasmids were represented by fewer ",
                 "than ", min_count, " reads across both replicates at the ",
                 "T0 timepoint,\nand were excluded from the analysis. ",
                 sum(!(too_few_counts)), " genes remained.\n"
                 ))

  results_vec <- entrezs_vec[!(too_few_counts)]
  return(results_vec)
}


GetAvailableGenes(essentials_df[, "Entrez_ID"])
GetAvailableGenes(non_essentials_df[, "Entrez_ID"])









CountsForPlasmid("PRNP")





