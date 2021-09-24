### 20th September 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")
R_functions_directory      <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))
source(file.path(s2r1_R_functions_directory, "06) Drawing scatter plots of plate-level summary data.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")
file_output_directory    <- file.path(s2r3_directory, "5) Output", "Figures", "Troubleshooting")




# Load data ---------------------------------------------------------------

load(file.path(s2r3_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r3_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))

load(file.path(s2r2_R_objects_directory, "30) Calculate correction factors to account for mean read counts.RData"))




# Define functions --------------------------------------------------------

RunScatter <- function(input_df, use_var, point_cex = 1.1) {

  x_column <- paste0("Run2_", tolower(use_var))
  y_column <- paste0("Run3_", tolower(use_var))

  use_mai <- c(0.9, 0.9, 0.7, 1.2)


  ## Define the point colors

  three_colors <- c(brewer.pal(9, "Blues")[[8]],
                    brewer.pal(9, "Purples")[[8]],
                    brewer.pal(9, "Reds")[[8]]
                    )
  colors_vec <- ifelse(input_df[, "Correction_factor"] == 4,
                       three_colors[[3]],
                       ifelse(input_df[, "Correction_factor"] == 2,
                              three_colors[[2]],
                              three_colors[[1]]
                              )
                       )
  input_df[, "Point_color"] <- colors_vec


  ## Draw the scatter plot

  RegressionScatter(input_df,
                    x_column,
                    y_column,
                    same_limits_xy = TRUE,
                    color_column   = "Point_color",
                    x_label        = "Run 2",
                    y_label        = "Run 3",
                    use_mai        = use_mai,
                    label_line     = 2.4,
                    point_cex      = point_cex
                    )


  ## Draw the legend

  border_alpha_hex <- substr(rgb(1, 1, 1, 0.75), 8, 9)
  fill_alpha_hex <- substr(rgb(1, 1, 1, 0.5), 8, 9)

  old_mai <- par("mai" = use_mai)

  legend("right",
         inset     = -0.205,
         legend    = c("1\u00d7 DNA",
                       "2\u00d7 DNA"#,
                       #"4\u00d7 DNA"
                       ),
         y.intersp = 1.4,
         pch       = 21,
         pt.cex    = point_cex,
         col       = paste0(three_colors[1:2], border_alpha_hex),
         pt.bg     = paste0(three_colors[1:2], fill_alpha_hex),
         bty       = "n",
         xpd       = NA
         )


  ## Draw the title

  mtext(axis_labels_list[[use_var]], line = 0.9)

  par(old_mai)

  return(invisible(NULL))
}




# Produce plate-level summary data ----------------------------------------

pool3_summaries_df <- SummarizePlates(ccs7_df_list[["filtered_summary_df"]])




# Merge data from run 2 and run 3 -----------------------------------------

use_columns <- c(
  "Number_of_wells", "Median_count", "Median_count_per_well",
  "Fraction_wells_with_few_reads", "Fraction_wells_with_no_reads",
  "Sum_counts_in_k"
)

run3_df <- pool3_summaries_df[, use_columns]
matches_vec <- match(pool3_summaries_df[, "Plate_name"],
                     extended_df[, "Plate_name"]
                     )
stopifnot(!(anyNA(matches_vec)))
run2_df <- extended_df[matches_vec, use_columns]

names(run2_df) <- paste0("Run2_", tolower(names(run2_df)))
names(run3_df) <- paste0("Run3_", tolower(names(run3_df)))

common_columns <- c(
  "Plate_name", "Plate_number", "Barcode_ID",
  "Colony_picked", "Number_of_wells"
)

runs_df <- data.frame(
  pool3_summaries_df[, common_columns],
  extended_df[matches_vec, "Correction_factor", drop = FALSE],
  run2_df,
  run3_df,
  stringsAsFactors = FALSE,
  row.names = NULL
)



# Correlate the sequencing yields between runs 2 and 3 --------------------

RunScatter(runs_df, "Sum_counts_in_k")

for (make_PDF in c(FALSE, TRUE)) {

  if (make_PDF) {
    pdf(file.path(file_output_directory, "Sequencing yields for runs 2 and 3.pdf"),
        width  = 5 + 2.1,
        height = 5 + 1.6
        )
  }

  for (plot_var in names(axis_labels_list)[4:(length(axis_labels_list) - 1)]) {
    print(plot_var)
    RunScatter(runs_df, plot_var)
  }

  if (make_PDF) {
    dev.off()
  }
}



# Save data ---------------------------------------------------------------

save(list = "pool3_summaries_df",
     file = file.path(s2r3_R_objects_directory, "29) Correlate median read count with DNA concentration.RData")
     )




