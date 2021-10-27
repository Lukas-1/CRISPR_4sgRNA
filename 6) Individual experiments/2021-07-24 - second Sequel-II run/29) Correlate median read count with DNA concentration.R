### 18th August 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")
R_functions_directory      <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))
source(file.path(s2r1_R_functions_directory, "06) Drawing scatter plots of plate-level summary data.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory          <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
concentration_directory <- file.path(s2r2_directory, "2) Input/Metadata/DNA concentrations")
R_objects_directory     <- file.path(s2r2_directory, "3) R objects")
file_output_directory   <- file.path(s2r2_directory, "5) Output", "Figures", "Troubleshooting")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(R_objects_directory, "11) Process demultiplexed PacBio reads - plates_analysis_list.RData"))




# Define functions --------------------------------------------------------

GetValidRows <- function(input_df) {
  are_all_NA <- apply(input_df, 1, function(x) all(is.na(x)))
  first_all_NA <- which(are_all_NA)[[1]]
  use_indices <- seq_len(first_all_NA - 1L)
  results_df <- input_df[use_indices, ]
  return(results_df)
}



RemoveRedundantColumns <- function(input_df) {
  are_all_the_same <- vapply(input_df,
                             function(x) length(unique(x)) == 1,
                             logical(1)
                             )
  input_df <- input_df[, !(are_all_the_same)]
  return(input_df)
}



PoolScatter <- function(input_df, x_column, y_column) {

  ## Define the point colors

  pool1_color  <- brewer.pal(9, "Greens")[[7]]
  pool2_color  <- brewer.pal(9, "Blues")[[7]]
  no_DNA_color <- brewer.pal(9, "Reds")[[3]]
  hi_DNA_color <- brewer.pal(9, "Reds")[[8]]

  colors_vec <- ifelse(input_df[, "Run2_pool"] == 1, pool1_color, pool2_color)
  colors_vec[input_df[, "No_DNA"]] <- no_DNA_color
  colors_vec[input_df[, "Was_corrected"]] <- hi_DNA_color

  input_df[["Point_colors"]] <- colors_vec
  point_cex = 1.1

  use_mai <- c(0.9, 0.9, 0.4, 1.5)


  ## Draw the scatter plot

  RegressionScatter(input_df,
                    x_column,
                    y_column,
                    color_column = "Point_colors",
                    x_label      = axis_labels_list[[x_column]],
                    y_label      = axis_labels_list[[y_column]],
                    use_mai      = use_mai,
                    point_cex    = point_cex
                    )


  ## Draw the legend

  legend_colors <- c(pool1_color, pool2_color, no_DNA_color, hi_DNA_color)
  border_alpha_hex <- substr(rgb(1, 1, 1, 0.75), 8, 9)
  fill_alpha_hex <- substr(rgb(1, 1, 1, 0.5), 8, 9)

  old_mai <- par("mai" = use_mai)

  legend("right",
         inset     = -0.27,
         legend    = c("Pool 1",
                       "Pool 2",
                       "No DNA",
                       "Half amount"
                       ),
         y.intersp = 1.4,
         pch       = 21,
         pt.cex    = point_cex,
         col       = paste0(legend_colors, border_alpha_hex),
         pt.bg     = paste0(legend_colors, fill_alpha_hex),
         bty       = "n",
         xpd       = NA
         )

  par(old_mai)

  return(invisible(NULL))
}



BarcodeScatter <- function(input_df, use_var, same_limits_xy = TRUE) {

  x_column <- paste0("Pool1_", tolower(use_var))
  y_column <- paste0("Pool2_", tolower(use_var))

  use_mai <- c(0.8, 0.9, 0.8, 1)

  RegressionScatter(input_df,
                    x_column,
                    y_column,
                    x_label        = "Pool 1",
                    y_label        = "Pool 2",
                    use_mai        = c(0.8, 0.9, 0.8, 1),
                    label_line     = 2.4,
                    same_limits_xy = same_limits_xy
                    )
  old_mai <- par(mai = use_mai)

  mtext(axis_labels_list[[use_var]], line = 0.9)

  par(old_mai)

  return(invisible(NULL))
}




# Read in data ------------------------------------------------------------

ReadConc <- function(file_name) {
  stopifnot("concentration_directory" %in% ls(envir = globalenv()))
  results_df <- data.frame(read_excel(file.path(concentration_directory,
                                                paste0(file_name, ".xlsx")
                                                ),
                                      skip = 1,
                                      sheet = "Data Summary"
                                      ),
                            stringsAsFactors = FALSE,
                            check.names = FALSE
                            )[, 1:8]
  return(results_df)
}

conc_pool1_df <- ReadConc("ABMM20210714_P1")
conc_pool2_df <- ReadConc("ABMM2021_07_14_P2")




# Merge the concentration measurements ------------------------------------

names(conc_pool1_df)[[6]] <- names(conc_pool2_df)[[6]]
names(conc_pool2_df)[[8]] <- names(conc_pool1_df)[[8]]

conc_df <- rbind.data.frame(GetValidRows(conc_pool1_df),
                            GetValidRows(conc_pool2_df),
                            stringsAsFactors = FALSE,
                            make.row.names = FALSE
                            )

new_order <- order(match(conc_df[["b-fabric number"]], plates_df[["Number_96wp"]]))
conc_df <- conc_df[new_order, ]
row.names(conc_df) <- NULL




# Tidy concentration data -------------------------------------------------

conc_df <- RemoveRedundantColumns(conc_df)
conc_df[[2]] <- as.numeric(conc_df[[2]])

missing_plate_number <- plates_df[["Number_96wp"]][which(plates_df[["Plate_name"]] == "HO_13")]

conc_df[["Was_corrected"]] <- conc_df[[2]] > 25

conc_df[["No_DNA"]] <- conc_df[, "b-fabric number"] == missing_plate_number

conc_df[["Corrected_concentration"]] <- ifelse(conc_df[["Was_corrected"]],
                                               conc_df[[2]] / 2,
                                               conc_df[[2]]
                                               )




# Revert the deconvolution of plates 79 and 82 ----------------------------

new_reads_df <- plates_analysis_list[["individual_reads_df"]]

were_found <- new_reads_df[["ZMW"]] %in% ccs_df[["ZMW"]]
table(new_reads_df[["Plate_number"]][!(were_found)])

new_reads_df <- new_reads_df[were_found, ]
new_reads_df[["Plate_number"]] <- ifelse(new_reads_df[["Plate_number"]] == 79L,
                                         82L,
                                         new_reads_df[["Plate_number"]]
                                         )

new_reads_df[["Combined_ID"]] <- paste0("Plate", formatC(new_reads_df[["Plate_number"]], width = 3, flag = "0"),
                                        "_Well", formatC(new_reads_df[["Well_number"]],  width = 3, flag = "0")
                                        )

new_order <- order(new_reads_df[["Combined_ID"]])
new_reads_df <- new_reads_df[new_order, ]

row.names(new_reads_df) <- NULL




# Create new summary data, after reverting the convolution ----------------

new_analysis_list <- list("individual_reads_df" = new_reads_df)


ccs_df[["Passed_filters"]] <- ccs_df[["Plate_passed_filters"]] &
                              (ccs_df[["Well_passed_filters"]] %in% TRUE)
ccs7_zmws <- GetCCS7_ZMWs(ccs_df)

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)

ccs7_df_list <- SummarizeWells(new_analysis_list,
                               use_zmws           = ccs7_zmws,
                               ID_column          = "Combined_ID",
                               unique_IDs         = sg_sequences_df[["Combined_ID"]],
                               deletions_df       = NULL,
                               aligned_contam_df  = NULL,
                               filter_cross_plate = FALSE
                               )




# Produce plate-level summary data ----------------------------------------

summaries_df <- SummarizePlates(ccs7_df_list[["filtered_summary_df"]])




# Merge the two data frames -----------------------------------------------

matches_vec <- match(summaries_df[["Number_96wp"]], conc_df[["b-fabric number"]])
setdiff(seq_len(96), matches_vec)
matches_vec[duplicated(matches_vec)] <- 96L

use_columns <- c("position on the FGCZ plate",
                 "no of PacBio reads (kbp)",
                 names(conc_df)[[2]],
                 "Corrected_concentration",
                 "No_DNA",
                 "Was_corrected"
                 )
matched_df <- conc_df[matches_vec, use_columns]
names(matched_df)[1:3] <- c("Position_FGCZ_plate", "Thousands_of_reads", "Original_concentration")

merged_df <- data.frame(summaries_df,
                        matched_df,
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                        )



# Look for adapter-specific effects ---------------------------------------

eligible_barcodes <- merged_df[["Barcode_ID"]][!(merged_df[["Colony_picked"]])]
eligible_barcodes <- eligible_barcodes[duplicated(eligible_barcodes)]

pool1_df <- merged_df[merged_df[["Run2_pool"]] %in% 1, ]
pool2_df <- merged_df[merged_df[["Run2_pool"]] %in% 2, ]

use_columns <- c(
  "Number_of_wells", "Median_read_count", "Median_count_adjust_by_no_of_wells",
  "Fraction_wells_with_few_reads", "Fraction_wells_with_no_reads",
  "Sum_counts_in_k", "Thousands_of_reads"
)

pool1_df <- pool1_df[match(eligible_barcodes, pool1_df[["Barcode_ID"]]), use_columns]
pool2_df <- pool2_df[match(eligible_barcodes, pool2_df[["Barcode_ID"]]), use_columns]

names(pool1_df) <- paste0("Pool1_", tolower(names(pool1_df)))
names(pool2_df) <- paste0("Pool2_", tolower(names(pool2_df)))

barcodes_df <- data.frame(
  "Barcode_ID" = eligible_barcodes,
  pool1_df,
  pool2_df,
  stringsAsFactors = FALSE,
  row.names = NULL
)




# Explore results ---------------------------------------------------------

cor.test(merged_df[["Original_concentration"]], merged_df[["Sum_counts_in_k"]])
cor.test(merged_df[["Corrected_concentration"]], merged_df[["Sum_counts_in_k"]])




# Plot the sequencing yield vs. the DNA concentration ---------------------

PoolScatter(merged_df, "Corrected_concentration", "Sum_counts_in_k")

for (make_PDF in c(FALSE, TRUE)) {

  if (make_PDF) {
    pdf(file.path(file_output_directory, "Sequencing yield vs. DNA concentration.pdf"),
        width  = 5 + 2.4,
        height = 5 + 1.3
        )
  }

  for (x_var in c("Corrected_concentration", "Original_concentration")) {
    for (y_var in names(axis_labels_list)[3:(length(axis_labels_list) - 1)]) {
      PoolScatter(merged_df, x_var, y_var)
    }
  }

  if (make_PDF) {
    dev.off()
  }
}




# Correlate the sequencing yields for pairs of barcodes -------------------

BarcodeScatter(barcodes_df, "Sum_counts_in_k")

for (make_PDF in c(FALSE, TRUE)) {

  if (make_PDF) {
    pdf(file.path(file_output_directory, "Sequencing yields for pairs of barcodes.pdf"),
        width  = 5.9,
        height = 5.6
        )
  }

  for (plot_var in names(axis_labels_list)[3:(length(axis_labels_list) - 1)]) {
    print(plot_var)
    BarcodeScatter(barcodes_df, plot_var)
  }

  if (make_PDF) {
    dev.off()
  }
}




# Save data ---------------------------------------------------------------

merged_plates_df <- merged_df

save(list = "merged_plates_df",
     file = file.path(R_objects_directory, "29) Correlate median read count with DNA concentration.RData")
     )




