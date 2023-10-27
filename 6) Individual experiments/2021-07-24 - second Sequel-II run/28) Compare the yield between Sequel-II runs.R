### 23rd July 2021 ###



# Import packages and source code -----------------------------------------

library("beeswarm")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions
source(file.path(R_functions_directory, "20) Summarizing data across wells.R"))




# Define folder paths -----------------------------------------------------

s2r1_directory           <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r1_R_objects_directory <- file.path(s2r1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")

file_output_directory    <- file.path(s2r2_directory, "5) Output", "Figures", "Troubleshooting")
PNGs_output_directory    <- file.path(s2r2_directory, "5) Output", "PNGs", "Troubleshooting")




# Load data ---------------------------------------------------------------

load(file.path(s2r1_R_objects_directory, "01) Process and export plate barcodes.RData"))
s2r1_plates_df <- plates_df
load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
s2r2_plates_df <- plates_df
rm(plates_df)

s2r1_plates_df[["Run2_pool"]] <- NA

plates_columns <- intersect(names(s2r2_plates_df), names(s2r1_plates_df))
plates_df <- rbind.data.frame(s2r1_plates_df[, plates_columns],
                              s2r2_plates_df[, plates_columns],
                              stringsAsFactors = FALSE,
                              make.row.names = FALSE
                              )


load(file.path(s2r1_R_objects_directory, "05) Read in PacBio data.RData"))
s2r1_ccs_df <- ccs_df

load(file.path(s2r1_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))
s2r1_ccs3_df_list <- ccs3_df_list
s2r1_ccs5_df_list <- ccs5_df_list
s2r1_ccs7_df_list <- ccs7_df_list
load(file.path(s2r1_R_objects_directory, "13) Process demultiplexed reads - with subsampling.RData"))

s2r1_subsampled_list <- subsampled_list
rm(subsampled_list)


load(file.path(s2r2_R_objects_directory, "05) Read in PacBio data.RData"))
s2r2_ccs_df <- ccs_df
rm(ccs_df)

load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))
s2r2_ccs3_df_list <- ccs3_df_list
s2r2_ccs5_df_list <- ccs5_df_list
s2r2_ccs7_df_list <- ccs7_df_list
rm(ccs3_df_list)
rm(ccs5_df_list)
rm(ccs7_df_list)





# Define functions --------------------------------------------------------

AddMetricColumns <- function(input_ccs_df, use_df_list) {
  input_ccs_df[["Is_hifi_read"]] <- (input_ccs_df[["Read_quality"]] != -1) &
                                    (input_ccs_df[["Num_full_passes"]] >= 3)

  lima_pass_filters <- use_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
  use_zmws <- use_df_list[["individual_reads_df"]][["ZMW"]][lima_pass_filters]
  input_ccs_df[["Passes_barcode_filters"]] <- input_ccs_df[["ZMW"]] %in% use_zmws

  return(input_ccs_df)
}



ReadCountsByRun <- function(numeric_vec, groups_fac, logical_vec) {

  stopifnot(length(numeric_vec) == length(groups_fac))
  stopifnot(length(numeric_vec) == length(logical_vec))


  ## Determine the point colors
  point_blue <- brewer.pal(9, "Blues")[[9]]
  point_red <- brewer.pal(9, "Reds")[[7]]
  colors_vec <- ifelse(logical_vec, point_red, point_blue)


  ## Determine group positions
  num_groups <- length(unique(means_df[["SmrtCell"]]))
  side_gap <- 0.5
  group_positions <- seq_len(num_groups)
  group_limits <- c((min(group_positions) - side_gap) - (num_groups * 0.04),
                     max(group_positions) + side_gap  + (num_groups * 0.04)
                    )


  ## Prepare the data axis
  numeric_axis_pos <- pretty(c(0, max(numeric_vec)), n = 10)
  numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
  numeric_axis_labels <- format(numeric_axis_pos)



  ## Set up the plot canvas and annotations
  use_mar <- c(5.5, 5, 4, 8)
  old_mar <- par("mar"     = use_mar,
                 "lheight" = 1.25
                 )
  plot(1,
       xlim = group_limits,
       ylim = numeric_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )

  DrawGridlines(numeric_limits)

  mtext(text = "Mean number of reads per well", side = 2, line = 3)
  title("Read counts (plate-level summary)")

  axis(2,
       mgp    = c(2.5, 0.65, 0),
       tcl    = -0.45,
       las    = 1,
       at     = numeric_axis_pos,
       labels = numeric_axis_labels
       )
  box(bty = "l")

  mtext(text = c("Run 1",
                 "Run 1\n45% sampled\nSimulation 1",
                 "Run 1\n45% sampled\nSimulation 2",
                 "Run 2\nPool 1",
                 "Run 2\nPool 2"
                 ),
        at = seq_len(num_groups),
        side = 1,
        padj = 1
        )


  ## Draw the superimposed boxplots
  boxplot(numeric_vec ~ groups_fac,
          at         = group_positions,
          boxwex     = 0.8,
          outline    = FALSE,
          names      = rep.int("", length(group_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * 3,
          col        = brewer.pal(9, "Blues")[[2]],
          border     = brewer.pal(9, "Blues")[[7]],
          add        = TRUE,
          axes       = FALSE,
          lwd        = 1
          )


  # Draw the points
  set.seed(1) # For reproducible positions
  beeswarm(numeric_vec ~ groups_fac,
           at       = group_positions,
           pch      = 16,
           cex      = 0.8,
           yaxs     = "i",
           pwcol    = colors_vec,
           bty      = "l",
           xlab     = "",
           ylab     = "Mean number of reads per well",
           labels   = rep("", num_groups),
           ylim     = numeric_limits,
           yaxs     = "i",
           xaxs     = "i",
           add      = TRUE,
           priority = "random"
           )


  ## Prepare for drawing the legend
  y_mid <- 0.5
  large_gap_lines <- 2.5
  large_gap <- diff(grconvertY(c(0, large_gap_lines), from = "lines", to = "npc"))
  small_gap <- large_gap / 2

  labels_list <- list(
    c("At least 75%",
      "of wells are",
      "unproblematic"
    ),
    c("Fewer than",
      "10 reads for",
      "> 25% of wells"
    )
  )

  is_large_gap <- lapply(labels_list, function(x) {
    if (length(x) == 1) {
      FALSE
    } else {
      c(TRUE, rep(FALSE, length(x) - 1))
    }
  })

  gaps_vec <- ifelse(unlist(is_large_gap), large_gap, small_gap)
  gaps_vec[[1]] <- 0
  total_span <- sum(gaps_vec)
  start_y <- y_mid + (total_span / 2)
  y_sequence <- start_y - cumsum(gaps_vec)
  y_pos <- grconvertY(y = y_sequence, from = "npc", to = "user")
  x_start <- 1 + diff(grconvertX(c(0, 0.8), from = "lines", to = "npc"))

  ## Draw the legend
  text(x      = grconvertX(x = x_start, from = "npc", to = "user"),
       y      = y_pos,
       cex    = 1,
       labels = sapply(unlist(labels_list), VerticalAdjust),
       adj    = c(0, 0.5),
       xpd    = NA
       )

  points(x   = rep(grconvertX(x = x_start + 0.003, from = "npc", to = "user"), 2),
         y   = y_pos[c(1, 4)],
         cex = 1,
         pch = 16,
         col = c(point_blue, point_red),
         xpd = NA
         )

  return(invisible(NULL))
}





# Compare the yield between the first and second sequencing runs ----------

table((s2r1_ccs_df[["Read_quality"]] != -1) & (s2r1_ccs_df[["Num_full_passes"]] >= 3))
table((s2r2_ccs_df[["Read_quality"]] != -1) & (s2r2_ccs_df[["Num_full_passes"]] >= 3))

s2r1_ccs_df <- AddMetricColumns(s2r1_ccs_df, s2r1_ccs3_df_list)
s2r2_ccs_df <- AddMetricColumns(s2r2_ccs_df, s2r2_ccs3_df_list)

nrow(s2r1_ccs_df)
nrow(s2r2_ccs_df)

table(s2r1_ccs_df[["Is_hifi_read"]])
table(s2r2_ccs_df[["Is_hifi_read"]])

table(s2r1_ccs_df[["Is_hifi_read"]], s2r1_ccs_df[["Passes_barcode_filters"]])
table(s2r2_ccs_df[["Is_hifi_read"]], s2r2_ccs_df[["Passes_barcode_filters"]])





# Create a combined data frame of read counts -----------------------------

stopifnot(!(anyDuplicated(plates_df[["Plate_number"]])))

matches_vec <- match(s2r2_ccs7_df_list[["filtered_summary_df"]][["Plate_number"]],
                     s2r2_plates_df[["Plate_number"]]
                     )
s2r2_ccs7_df_list[["filtered_summary_df"]][["Pool_number"]] <- s2r2_plates_df[["Run2_pool"]][matches_vec]
are_pool1 <- s2r2_ccs7_df_list[["filtered_summary_df"]][["Pool_number"]] == 1

use_columns <- c("Combined_ID", "Plate_number", "Well_number", "Count_total")

counts_df <- rbind.data.frame(
  data.frame("SmrtCell" = "Run 1",
             s2r1_ccs7_df_list[["filtered_summary_df"]][, use_columns],
             stringsAsFactors = FALSE
             ),
  data.frame("SmrtCell" = "Run 1 (45%, sub-sampled #1)",
             s2r1_subsampled_list[["45% sampled"]][["rep1"]][["ccs7"]][["filtered_summary_df"]][, use_columns],
             stringsAsFactors = FALSE
             ),
  data.frame("SmrtCell" = "Run 1 (45%, sub-sampled #2)",
             s2r1_subsampled_list[["45% sampled"]][["rep2"]][["ccs7"]][["filtered_summary_df"]][, use_columns],
             stringsAsFactors = FALSE
             ),
  data.frame("SmrtCell" = "Run 2, pool 1",
             s2r2_ccs7_df_list[["filtered_summary_df"]][are_pool1, use_columns],
             stringsAsFactors = FALSE
             ),
  data.frame("SmrtCell" = "Run 2, pool 2",
             s2r2_ccs7_df_list[["filtered_summary_df"]][!(are_pool1), use_columns],
             stringsAsFactors = FALSE
             ),
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

names(counts_df)[names(counts_df) == "Count_total"] <- "Num_reads"

matches_vec <- match(counts_df[["Plate_number"]], plates_df[["Plate_number"]])
counts_df[["Colony_picked"]] <- plates_df[["Colony_picked"]][matches_vec]





# Summarize the combined data frame of read counts ------------------------

plates_vec <- paste0(counts_df[["SmrtCell"]], "__", counts_df[["Plate_number"]])

counts_list <- tapply(seq_len(nrow(counts_df)),
                      factor(plates_vec, levels = unique(plates_vec)),
                      function(x) {
                        results_list <- c(
                          as.list(counts_df[x[[1]], c("SmrtCell", "Plate_number", "Colony_picked")]),
                          list("Mean_num_reads"    = mean(counts_df[["Num_reads"]][x]),
                               "Fraction_below_10" = sum(counts_df[["Num_reads"]][x] < 10) / length(x)
                               )
                          )
                      })
means_df <- do.call(rbind.data.frame, c(counts_list, stringsAsFactors = FALSE, make.row.names = FALSE))






# Select plates to show ---------------------------------------------------

are_excluded <- means_df[["Colony_picked"]] |
                (means_df[["Plate_number"]] == 79)

use_df <- means_df[!(are_excluded), ]





# Plot comparisons of read counts per plate -------------------------------

ReadCountsByRun(use_df[["Mean_num_reads"]],
                factor(use_df[["SmrtCell"]], levels = unique(use_df[["SmrtCell"]])),
                use_df[["Fraction_below_10"]] > 0.25
                )


use_width <- 8.9
use_height <- 6


## Set up PDF
pdf(file = file.path(file_output_directory, "Comparison of sequencing yields.pdf"),
    width = use_width, height = use_height
    )
ReadCountsByRun(use_df[["Mean_num_reads"]],
                factor(use_df[["SmrtCell"]], levels = unique(use_df[["SmrtCell"]])),
                use_df[["Fraction_below_10"]] > 0.25
                )
dev.off()



## Set up PNG
png(file = file.path(PNGs_output_directory, "Comparison of sequencing yields.png"),
    width = use_width, height = use_height, units = "in", res = 600
    )
ReadCountsByRun(use_df[["Mean_num_reads"]],
                factor(use_df[["SmrtCell"]], levels = unique(use_df[["SmrtCell"]])),
                use_df[["Fraction_below_10"]] > 0.25
                )
dev.off()







