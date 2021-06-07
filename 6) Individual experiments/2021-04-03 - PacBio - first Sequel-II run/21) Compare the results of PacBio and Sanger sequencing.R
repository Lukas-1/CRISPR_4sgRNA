### 3rd June 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For MakeInvisible and related functions





# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_input_directory     <- file.path(sql2_directory, "2) Input")
sanger_excel_path        <- file.path(file_input_directory,
                                      "Sanger sequencing",
                                      "sanger confirmation of NGS with retransformation and single colony picking.xlsx"
                                      )
file_output_directory    <- file.path(sql2_directory, "5) Output", "Figures", "Comparison of Sanger and PacBio sequencing")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))




# Read in data ------------------------------------------------------------

all_sheets <- excel_sheets(sanger_excel_path)
sanger_list <- lapply(all_sheets, function(x) read_excel(sanger_excel_path, sheet = x))





# Define functions --------------------------------------------------------

AddPercentages <- function(input_df) {
  count_regex <- "^(Count|Num)_"
  count_columns <- setdiff(grep(count_regex, names(input_df), value = TRUE), "Count_total")
  for (column_name in count_columns) {
    perc_column <- sub(count_regex, "Perc_", column_name)
    input_df[[perc_column]] <- input_df[[column_name]] / input_df[["Count_total"]]
  }
  return(input_df)
}


MakeCombinedIDs <- function(plate_numbers, well_numbers) {
  paste0("Plate", formatC(plate_numbers, width = 2, flag = "0"), "_",
         "Well",  formatC(well_numbers,  width = 3, flag = "0")
         )
}





# Process the Excel sheet from Sanger sequencing --------------------------

sanger_list <- lapply(sanger_list, function(x) as.data.frame(x, stringsAsFactors = FALSE))

sanger_df_list <- lapply(sanger_list, function(x) {
  well_name <- names(x)[[1]]
  use_df <- x[, 1:9]
  for (i in 2:9) {
    use_df[[i]] <- as.integer(use_df[[i]])
  }
  use_indices <- seq_len(nrow(use_df) - 3)
  use_df <- data.frame("Original_ID" = well_name,
                       "Read_number" = as.integer(use_df[[1]][use_indices]),
                       use_df[use_indices, 2:9],
                       stringsAsFactors = FALSE
                       )
  return(use_df)
})

sanger_df <- do.call(rbind.data.frame,
                     c(sanger_df_list,
                       list(stringsAsFactors = FALSE, make.row.names = FALSE)
                       )
                     )




# Translate the well names ------------------------------------------------

well_numbers_mat <- matrix(1:384, nrow = 16, ncol = 24, byrow = TRUE)
well_names_vec <- do.call(paste0, expand.grid(1:24, LETTERS[1:16])[, 2:1])
well_names_mat <- matrix(well_names_vec, nrow = 16, ncol = 24, byrow = TRUE)

well_IDs <- unique(sanger_df[["Original_ID"]])
well_ID_splits <- strsplit(well_IDs, " ", fixed = TRUE)
well_names <- sapply(well_ID_splits, "[[", 2)
plate_names <- sapply(well_ID_splits, "[[", 1)

well_numbers <- vapply(well_names, function(x) well_numbers_mat[well_names_mat == x], integer(1))
column_plate_numbers <- plates_df[["Plate_number"]][match(plate_names, plates_df[["Plate_name"]])]
bead_plate_numbers <- plates_df[["Plate_number"]][match(paste0(plate_names, "-beads"), plates_df[["Plate_name"]])]
columns_combined_IDs <- MakeCombinedIDs(column_plate_numbers, well_numbers)
beads_combined_IDs <- ifelse(is.na(bead_plate_numbers),
                             NA,
                             MakeCombinedIDs(bead_plate_numbers, well_numbers)
                             )

matches_vec <- match(sanger_df[["Original_ID"]], well_IDs)

sanger_df <- data.frame("Columns_ID" = columns_combined_IDs[matches_vec],
                        "Beads_ID" = beads_combined_IDs[matches_vec],
                        sanger_df,
                        stringsAsFactors = FALSE
                        )




# Evaluate combined sg_cr regions -----------------------------------------

sg_cr_mat <- do.call(cbind, lapply(1:4, function(x) {
  as.integer(rowSums(as.matrix(sanger_df[, paste0(c("sg", "cr"), x)])) == 2)
}))
colnames(sg_cr_mat) <- paste0("sg", 1:4, "_cr", 1:4)

at_least_mat <- do.call(cbind, lapply(1:3, function(x) {
  as.integer(rowSums(sg_cr_mat) >= x)
}))
colnames(at_least_mat) <- paste0("at_least_", 1:3)

sanger_df <- data.frame(sanger_df,
                        sg_cr_mat,
                        at_least_mat,
                        "all_4" = as.integer(rowSums(sg_cr_mat) == 4),
                        stringsAsFactors = FALSE
                        )

use_columns <- c(
  "sg1_cr1", "sg2_cr2", "sg3_cr3", "sg4_cr4",
  "at_least_1", "at_least_2", "at_least_3", "all_4"
)

summary_list <- by(sanger_df, sanger_df[["Columns_ID"]], function(x) {
  percent_vec <- colSums(as.matrix(x[, use_columns])) / nrow(x)
  names(percent_vec) <- paste0("Perc_", names(percent_vec))
  results_list <- c(
    as.list(unique(x[, c("Columns_ID", "Beads_ID", "Original_ID")])),
    list("Count" = nrow(x)),
    as.list(percent_vec)
  )
}, simplify = FALSE)
sanger_summary_df <- do.call(rbind.data.frame, c(summary_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))





# Estimate the effect of spurious contaminations --------------------------

use_columns <- c(
  "Combined_ID", "Plate_number", "Well_number", "Count_total",
  paste0("Count_", use_columns),
  "Num_contaminated_reads"
)

comparison_plates <- c("HA_11", "HO_1")
columns_plates <- plates_df[["Plate_number"]][plates_df[["Plate_name"]] %in% comparison_plates]
beads_plates   <- plates_df[["Plate_number"]][plates_df[["Plate_name"]] %in% paste0(comparison_plates, "-beads")]

use_df <- ccs7_df_list[["filtered_summary_df"]]

columns_df <- use_df[use_df[["Plate_number"]] %in% columns_plates, use_columns]
beads_df <-  use_df[use_df[["Plate_number"]] %in% beads_plates, use_columns]

columns_df <- AddPercentages(columns_df)
beads_df <- AddPercentages(beads_df)

contam_delta_vec <- columns_df[["Perc_contaminated_reads"]] -
                    beads_df[["Perc_contaminated_reads"]]
are_to_exclude <- (columns_df[["Perc_contaminated_reads"]] >= 0.3) |
                  (beads_df[["Perc_contaminated_reads"]] >= 0.3)

mean_contam_delta <- mean(contam_delta_vec[!(are_to_exclude)])





# Estimate the error rate for single colony-picked controls ---------------

control_plate <- plates_df[["Plate_number"]][plates_df[["Plate_name"]] == "Intctrl"]
controls_df <- use_df[use_df[["Plate_number"]] %in% control_plate, use_columns]

controls_df <- AddPercentages(controls_df)

perc_columns <- grep("^Perc_", names(controls_df), value = TRUE)
mean_percentages <- colMeans(controls_df[, perc_columns])
mean_percentages <- ifelse(perc_columns == "Perc_contaminated_reads",
                           mean_percentages,
                           1 - mean_percentages
                           )
names(mean_percentages) <- perc_columns





# Pull out the PacBio data for the sequenced wells ------------------------

matches_vec <- match(sanger_summary_df[["Columns_ID"]], use_df[["Combined_ID"]])
pacbio_columns_df <- AddPercentages(use_df[matches_vec, use_columns])
row.names(pacbio_columns_df) <- NULL

matches_vec <- match(sanger_summary_df[["Beads_ID"]], use_df[["Combined_ID"]])
pacbio_beads_df <- AddPercentages(use_df[matches_vec, use_columns])
row.names(pacbio_beads_df) <- NULL





# Compare the results of PacBio sequencing with Sanger sequencing ---------

common_columns <- grep("^Perc_", names(sanger_summary_df), value = TRUE)

use_control_percentages <- colMeans(controls_df[, common_columns])
columns_list <- lapply(common_columns, function(x) {
  control_perc_errors <- 1 - use_control_percentages[[x]]
  results_mat <- cbind("Sanger"                    = sanger_summary_df[, x],
                       "PacBio_columns"            = pacbio_columns_df[, x],
                       "PacBio_beads"              = pacbio_beads_df[, x],
                       "Corrected_PacBio_columns"  = pacbio_columns_df[, x] + control_perc_errors + mean_contam_delta,
                       "Corrected_PacBio_beads"    = pacbio_beads_df[, x]   + control_perc_errors,
                       "Corrected_delta_columns"   = sanger_summary_df[, x] - pacbio_columns_df[, x] - control_perc_errors - mean_contam_delta,
                       "Corrected_delta_beads"     = sanger_summary_df[, x] - pacbio_beads_df[, x]   - control_perc_errors,
                       "Uncorrected_delta_columns" = sanger_summary_df[, x] - pacbio_columns_df[, x],
                       "Uncorrected_delta_beads"   = sanger_summary_df[, x] - pacbio_beads_df[, x]
                       )
  colnames(results_mat) <- paste0(colnames(results_mat), "_", tolower(x))
  return(results_mat)
})
combined_df <- data.frame(
  sanger_summary_df[, c("Columns_ID", "Beads_ID", "Original_ID")],
  "Count_Sanger" = sanger_summary_df[["Count"]],
  "Count_PacBio_columns" = pacbio_columns_df[["Count_total"]],
  "Count_PacBio_beads" = pacbio_beads_df[["Count_total"]],
  do.call(cbind, columns_list),
  stringsAsFactors = FALSE
)




# Export tabular data -----------------------------------------------------

export_combined_df <- combined_df

have_NA_columns <- vapply(export_combined_df, anyNA, logical(1))
for (i in which(have_NA_columns)) {
  export_combined_df[[i]] <- ifelse(is.na(export_combined_df[[i]]),
                                    " ",
                                    export_combined_df[[i]]
                                    )
}

write.table(export_combined_df,
            file      = file.path(file_output_directory, "Comparison of Sanger and PacBio sequencing.tsv"),
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE
            )




# Perform statistical tests -----------------------------------------------

test_column <- "Perc_all_4"

uncorrected_t_test_results <- t.test(combined_df[["Sanger_perc_all_4"]],
                                     combined_df[["PacBio_columns_perc_all_4"]],
                                     paired = TRUE
                                     )

corrected_t_test_results   <- t.test(combined_df[["Sanger_perc_all_4"]],
                                     combined_df[["Corrected_PacBio_columns_perc_all_4"]],
                                     paired = TRUE
                                     )


t.test(combined_df[["Sanger_perc_all_4"]],
       ifelse(is.na(combined_df[["Corrected_PacBio_columns_perc_all_4"]]),
              combined_df[["Corrected_PacBio_columns_perc_all_4"]],
              combined_df[["Corrected_PacBio_beads_perc_all_4"]]
              ),
       paired = TRUE
       )





# Plot the comparison of Sanger vs. PacBio sequencing ---------------------

## Prepare data
use_metric <- "Perc_all_4"
use_prefixes <- c("Sanger",
                  "PacBio_columns", "Corrected_PacBio_columns",
                  "PacBio_beads", "Corrected_PacBio_beads"
                  )
use_columns <- paste0(use_prefixes, "_", tolower(use_metric))
numeric_mat <- as.matrix(combined_df[, use_columns])



## Set up PDF
pdf(file = file.path(file_output_directory, "Comparison of Sanger and PacBio sequencing.pdf"),
    width = 8, height = 5
    )

par(mar = c(5, 6, 7, 12))



## Determine group positions
num_groups <- 3
group_positions <- seq_len(num_groups)
width <- 2/3
final_width <- width * ((max(group_positions) - min(group_positions)) / (num_groups - 1))
group_limits <- c((min(group_positions) - 0.2) - (num_groups * 0.04),
                  max(group_positions) + 0.2 + (num_groups * 0.04)
                  )

## Prepare the data axis
use_numeric_limits <- c(0, max(numeric_mat, na.rm = TRUE) * 1.00)
numeric_axis_pos <- pretty(use_numeric_limits)
numeric_limits <- c(numeric_axis_pos[[1]], numeric_axis_pos[[length(numeric_axis_pos)]])
numeric_axis_labels <- paste0(format(numeric_axis_pos * 100), "%")


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

axis(2,
     at       = numeric_axis_pos,
     labels   = numeric_axis_labels,
     mgp      = c(3, 0.38, 0),
     gap.axis = 0,
     tcl      = -0.3,
     las      = 1,
     lwd      = par("lwd")
     )

mtext(text = "Reads with 4 correct sgRNAs + tracRNAs",
      side = 2,
      line = 3.2
      )
box(bty = "l")





## Draw the connecting lines
x_positions <- c(1, 1.9, 2.1, 2.9, 3.1)

lines_alpha <- 0.2
alpha_hex <- substr(rgb(1, 1, 1, lines_alpha), 8, 9)
use_grey <- brewer.pal(9, "Greys")[[7]]
use_color <- paste0(use_grey, alpha_hex)
for (row_index in seq_len(nrow(numeric_mat))) {
  are_NA <- is.na(numeric_mat[row_index, ])
  lines(x   = x_positions[!(are_NA)],
        y   = numeric_mat[row_index, ][!(are_NA)],
        col = use_color,
        lwd = 1,
        xpd = NA
        )
}



## Draw the data points
points_alpha <- 0.8
alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
uncorrected_color <- paste0(brewer.pal(9, "Blues")[[7]], alpha_hex)
corrected_color <- paste0(brewer.pal(9, "Purples")[[7]], alpha_hex)
column_colors <- ifelse(grepl("^Corrected_", colnames(numeric_mat)),
                        corrected_color,
                        uncorrected_color
                        )

for (column_index in seq_len(ncol(numeric_mat))) {
  points(x   = rep(x_positions[[column_index]], nrow(numeric_mat)),
         y   = numeric_mat[, column_index],
         pch = 16,
         cex = 1,
         col = column_colors[[column_index]],
         xpd = NA
         )
}




## Draw the x axis legend
text(x      = 1:3,
     y      = grconvertY(y = -0.05, from = "npc", to = "user"),
     labels = c("Sanger", "PacBio", "PacBio"),
     adj    = c(0.5, 1),
     xpd    = NA
     )

small_gap <- 0.08
text(x      = 1:3,
     y      = grconvertY(y = -0.05 - small_gap, from = "npc", to = "user"),
     labels = c("sequencing", "(columns)", "(beads)"),
     adj    = c(0.5, 1),
     xpd    = NA
     )


use_color_text <- 'color1("original") * color2(" / ") * color3("corrected")'
text_colors <- c(uncorrected_color, "black", corrected_color)
for (i in seq_along(text_colors)) {
  text(x      = 2:3,
       y      = grconvertY(y = -0.24, from = "npc", to = "user"),
       labels = VerticalAdjust(parse(text = MakeInvisible(use_color_text, i))),
       adj    = c(0.5, 1),
       col    = text_colors[[i]],
       cex    = 0.9,
       xpd    = NA
       )
}



## Draw the p values
use_grey <- "gray60"
lwd_grey <- "gray80"

y_start <- 1.16
use_lwd <- 0.75
segments(x0  = 1,
         x1  = 1.9,
         y0  = grconvertY(y = y_start, from = "npc", to = "user"),
         col = lwd_grey,
         lwd = use_lwd,
         xpd = NA
         )

y_tick_length <- 0.025
segments(x0  = c(1, 1.9),
         y0  = grconvertY(y = y_start, from = "npc", to = "user"),
         y1  = grconvertY(y = y_start - y_tick_length, from = "npc", to = "user"),
         col = lwd_grey,
         lwd = use_lwd,
         xpd = NA
         )

text(x      = sum(c(1, 1.9)) / 2,
     y      = grconvertY(y = y_start + 0.035, from = "npc", to = "user"),
     labels = as.expression(bquote(italic("p") * " = " * .(signif(uncorrected_t_test_results[["p.value"]], digits = 1)))),
     cex    = 0.8,
     col    = use_grey,
     lwd    = use_lwd,
     xpd    = NA
     )

y_gap <- 0.13
segments(x0  = c(1, 2.1),
         y0  = grconvertY(y = y_start + y_gap, from = "npc", to = "user"),
         y1  = grconvertY(y = y_start + y_gap - y_tick_length, from = "npc", to = "user"),
         col = lwd_grey,
         lwd = use_lwd,
         xpd = NA
         )
segments(x0  = 1,
         x1  = 2.1,
         y0  = grconvertY(y = y_start + y_gap, from = "npc", to = "user"),
         col = lwd_grey,
         lwd = use_lwd,
         xpd = NA
         )


text(x      = sum(c(1, 2.1)) / 2,
     y      = grconvertY(y = y_start + 0.035 + y_gap, from = "npc", to = "user"),
     labels = as.expression(bquote(italic("p") * " = " * .(signif(corrected_t_test_results[["p.value"]], digits = 1)))),
     cex    = 0.8,
     col    = use_grey,
     lwd    = use_lwd,
     xpd    = NA
     )



## Prepare for drawing the "error rates" legend
y_mid <- 0.5
large_gap <- 0.13

labels_list <- list(
  "Corrections:",
  c("Excess",
    "contaminations",
    "(columns)"
  ),
  c("Error rate in",
    "colony-picked",
    "controls"
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
x_start <- 1.1


## Show the descriptions of the error rates
text(x      = grconvertX(x = x_start,    from = "npc", to = "user"),
     y      = grconvertY(y = y_sequence, from = "npc", to = "user"),
     cex    = 1,
     labels = unlist(labels_list),
     adj    = c(0, 0),
     xpd    = NA
     )



## Show a barplot representation of the error rates
use_control_perc <- (1 - use_control_percentages)[[use_metric]]
rect_width_npc <- 0.005
rect_width <- diff(grconvertX(x = c(0, rect_width_npc), from = "npc", to = "user"))
rect_color <- brewer.pal(9, "Blues")[[9]]
rect_mid <- x_start + 0.33


rect(xleft   = grconvertX(x = rect_mid - rect_width, from = "npc", to = "user") - rect_width,
     xright  = grconvertX(x = rect_mid + rect_width, from = "npc", to = "user") + rect_width,
     ybottom = grconvertY(y = y_sequence[[4]], from = "npc", to = "user"),
     ytop    = grconvertY(y = y_sequence[[4]], from = "npc", to = "user") + mean_contam_delta,
     col     = rect_color,
     border  = NA,
     xpd     = NA
     )

rect(xleft   = grconvertX(x = rect_mid - rect_width, from = "npc", to = "user") - rect_width,
     xright  = grconvertX(x = rect_mid + rect_width, from = "npc", to = "user") + rect_width,
     ybottom = grconvertY(y = y_sequence[[7]], from = "npc", to = "user"),
     ytop    = grconvertY(y = y_sequence[[7]], from = "npc", to = "user") + use_control_perc,
     col     = rect_color,
     border  = NA,
     xpd     = NA
     )


## Label the bars with their numerical values
y_gap <- 0.03
text_cex <- 0.8
text(x      = grconvertX(x = rect_mid, from = "npc", to = "user"),
     y      = grconvertY(y = mean_contam_delta + y_sequence[[4]] + y_gap, from = "npc", to = "user"),
     cex    = text_cex,
     labels = paste0(signif(mean_contam_delta * 100, digits = 2), "%"),
     adj    = c(0.5, 0),
     font   = 2,

     xpd    = NA
     )

text(x      = grconvertX(x = rect_mid, from = "npc", to = "user"),
     y      = grconvertY(y = use_control_perc + y_sequence[[7]] + y_gap, from = "npc", to = "user"),
     cex    = text_cex,
     labels = paste0(signif(use_control_perc * 100, digits = 2), "%"),
     adj    = c(0.5, 0),
     font   = 2,
     xpd    = NA
     )


dev.off()





