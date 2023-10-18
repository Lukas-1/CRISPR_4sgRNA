### 28 May 2020 ###




# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR_4sgRNA"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")
file_output_directory            <- file.path(file_directory, "4) Output")






# Load data ---------------------------------------------------------------

load(file.path(intermediate_R_objects_directory, "3) Import data from external tools.RData"))






# Define functions --------------------------------------------------------

MakeEmptyPlot <- function(y_limits = c(0, 1)) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = y_limits, xaxs = "i", yaxs = "i"
       )
}


MakeHist <- function(numeric_vec,
                     use_color,
                     show_axis_labels = TRUE
                     ) {

  histogram_results <- hist(numeric_vec,
                            breaks = 200,
                            plot   = FALSE,
                            )

  MakeEmptyPlot(y_limits = c(0, 30))

  axis_grey <- "black"
  gridlines_grey <- "gray87"

  abline(v = seq(0.1, 0.9, by = 0.1), col = gridlines_grey, lwd = 0.75)
  box(col = axis_grey, lwd = 0.75)

  plot(histogram_results,
       add = TRUE,
       col = use_color,
       border = NA
       )

  tick_locations <- axTicks(1)
  if (show_axis_labels) {
    tick_labels <- paste0(tick_locations * 100, "%")
  } else {
    tick_labels <- rep("", length(tick_locations))
  }
  axis(1,
       labels    = tick_labels,
       at        = tick_locations,
       las       = 1,
       mgp       = c(3, 0.8, 0),
       tcl       = 0,
       cex.axis  = 1.3,
       font.axis = 1,
       lwd       = 0.75,
       col   = axis_grey
       )
  axis(2,
       las = 1,
       mgp = c(3, 0.45, 0),
       tcl = -0.3,
       lwd = 0.75,
       col = axis_grey
       )

}




DrawHistogramPanel <- function(assembly_df) {

  par(mfrow = c(4, 1),
      mar = c(0.5, 2, 1, 2),
      oma = c(4, 4, 2, 3)
      )

  for (i in 1:4) {
    title_text <- paste0("sg", i)
    MakeHist(assembly_df[[paste0("Fraction_correct_sg", i)]],
             use_color = "#003876",
             show_axis_labels = i == 4
             )
    text(x      = -0.067,
         y      = par("usr")[[3]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.5),
         labels = title_text,
         font   = 2,
         cex    = 1.5,
         xpd    = NA,
         adj    = c(1, 0.5)
         )
    if (i == 1) {
      text(x      = -0.07,
           y      = par("usr")[[4]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.11),
           adj    = c(0, 0),
           labels = "Counts",
           xpd    = NA
           )
    }
    if (i == 4) {
      text(x      = 0.5,
           y      = par("usr")[[3]] - ((par("usr")[[4]] - par("usr")[[3]]) * 0.4),
           labels = "Accuracy",
           font   = 1,
           xpd    = NA,
           cex    = 1.5
           )
    }
  }
  return(invisible(NULL))
}








# Draw the histograms -----------------------------------------------------

DrawHistogramPanel(assembly_df)

png(filename = file.path(file_output_directory, "Accuracy histograms.png"),
    res    = 600,
    height = 5.5,
    width  = 6.5,
    units  = "in"
    )
DrawHistogramPanel(assembly_df)
dev.off()



pdf(file.path(file_output_directory, "Accuracy histograms.pdf"),
    height = 5.5,
    width  = 6.5
    )
DrawHistogramPanel(assembly_df)
dev.off()












