### 20th September 2021 ###



# Import packages and source code -----------------------------------------

library("RColorBrewer")




# Define labels -----------------------------------------------------------

axis_labels_list <- list(
  "Original_concentration"        = expression("Measured DNA concentration (Units" * scriptscriptstyle(" ") * "/" * scriptscriptstyle(" ") * mu * "L)"),
  "Corrected_concentration"       = expression("Corrected DNA concentration (Units" * scriptscriptstyle(" ") * "/" * scriptscriptstyle(" ") * mu * "L)"),
  "Thousands_of_reads"            = expression("Read count, demultiplexing by FGCZ (thousands)"),
  "Sum_counts_in_k"               = expression("High-quality read count (thousands)"),
  "Median_count_per_well"         = expression("Median number of high-quality reads per well"),
  "Fraction_wells_with_few_reads" = expression("Fraction of wells with" < "10 reads"),
  "Fraction_wells_with_no_reads"  = expression("Fraction of wells with zero reads")
)



# Define functions --------------------------------------------------------

RegressionScatter <- function(input_df,
                              x_column,
                              y_column,
                              same_limits_xy = FALSE,
                              color_column   = NULL,
                              x_label        = "",
                              y_label        = "",
                              use_mai        = rep(1, 4),
                              point_cex      = 1.1,
                              label_line     = 2.7
                              ) {

  ## Make sure the x values are in order (for the lm model below)

  new_order <- order(input_df[, x_column])
  input_df <- input_df[new_order, ]
  row.names(input_df) <- NULL


  ## Define the point colors

  if (is.null(color_column)) {
    colors_vec <- brewer.pal(9, "Purples")[[8]]
  } else {
    colors_vec <- input_df[, color_column]
  }
  border_alpha_hex <- substr(rgb(1, 1, 1, 0.75), 8, 9)
  fill_alpha_hex <- substr(rgb(1, 1, 1, 0.5), 8, 9)


  ## Perform a linear regression

  model_df <- input_df[, c(x_column, y_column)]
  names(model_df) <- c("x_var", "y_var")
  lm_model <- lm(y_var ~ x_var, data = model_df)
  lm_summary <- summary(lm_model)

  if (same_limits_xy) {
    use_x_max <- max(c(input_df[, x_column], input_df[, y_column]), na.rm = TRUE)
  } else {
    use_x_max <- max(input_df[, x_column], na.rm = TRUE)
  }
  assign("delete_input_df", input_df, envir = globalenv())
  assign("delete_x_column", x_column, envir = globalenv())
  assign("delete_y_column", y_column, envir = globalenv())
  assign("delete_use_x_max", use_x_max, envir = globalenv())
  new_seq <- seq(0, use_x_max, length.out = 200)
  new_df <- data.frame("x_var" = new_seq)
  conf_int_mat <- predict(lm_model,
                          newdata = new_df,
                          interval = "confidence",
                          level = 0.95
                          )


  ## Establish the plot limits
  x_range <- c(0, max(input_df[[x_column]], na.rm = TRUE))
  y_range <- c(0, max(input_df[[y_column]], na.rm = TRUE))
  space_factor <- 0.02
  x_space <- (x_range[[2]] - x_range[[1]]) * space_factor
  y_space <- (y_range[[2]] - y_range[[1]]) * space_factor
  x_limits <- c(x_range[[1]] - x_space, x_range[[2]] + x_space)
  y_limits <- c(y_range[[1]] - y_space, y_range[[2]] + y_space)

  if (same_limits_xy) {
    common_limits <- c(min(x_limits[[1]], y_limits[[1]]),
                       max(x_limits[[2]], y_limits[[2]])
                       )
    x_limits <- common_limits
    y_limits <- common_limits
  }

  ## Prepare the plot region
  old_mai <- par("mai" = use_mai)
  plot(1,
       type = "n",
       las  = 1,
       mgp  = c(label_line, 0.55, 0),
       tcl  = -0.35,
       xaxs = "i",
       yaxs = "i",
       xlim = x_limits,
       ylim = y_limits,
       xlab = x_label,
       ylab = y_label
       )


  ## Add the R^2
  r2_text <- bquote(italic("R") * ""^2  ~ "=" ~
                    .(round(lm_summary[["r.squared"]], digits = 2))
                    )

  text(x      = par("usr")[[2]] + grconvertX(0.05, from = "npc", to = "user"),
       y      = par("usr")[[4]] - grconvertY(0.05, from = "npc", to = "user"),
       labels = r2_text,
       xpd    = NA,
       adj    = 0
       )


  ## Draw the regression line
  polygon(c(new_df[, 1], rev(new_df[, 1])),
          c(conf_int_mat[, 2], rev(conf_int_mat[, 3])),
          col = "gray92", border = NA
          )
  lines(new_df[, 1], conf_int_mat[, 1], col = "gray70", lwd = 2)


  ## Draw the points
  points(input_df[[x_column]],
         input_df[[y_column]],
         cex  = point_cex,
         pch  = 21,
         col  = paste0(colors_vec, border_alpha_hex),
         bg   = paste0(colors_vec, fill_alpha_hex)
         )


  ## Final steps
  box()
  par(old_mai)
  return(invisible(NULL))
}




SummarizePlates <- function(summary_df) {

  plates_fac <- factor(summary_df[["Plate_number"]])

  median_counts <- tapply(summary_df[["Count_total"]], plates_fac, median)
  sum_counts <- tapply(summary_df[["Count_total"]], plates_fac, sum)

  num_wells_with_few_reads <- tapply(summary_df[["Count_total"]],
                                     plates_fac,
                                     function(x) sum(x < 10)
                                     )
  num_wells_with_no_reads  <- tapply(summary_df[["Count_total"]],
                                     plates_fac,
                                     function(x) sum(x == 0)
                                     )

  matches_vec <- match(levels(plates_fac), plates_df[["Plate_number"]])
  number_of_wells <- tabulate(plates_fac)

  results_df <- data.frame(
    plates_df[matches_vec, names(plates_df) != "Barcode_sequence"],
    "Number_of_wells"               = number_of_wells,
    "Median_count"                  = median_counts,
    "Median_count_per_well"         = median_counts / 384 * number_of_wells,
    "Fraction_wells_with_few_reads" = num_wells_with_few_reads / number_of_wells,
    "Fraction_wells_with_no_reads"  = num_wells_with_no_reads / number_of_wells,
    "Sum_counts"                    = sum_counts,
    "Sum_counts_in_k"               = sum_counts / 1000,
    stringsAsFactors                = FALSE,
    row.names                       = NULL
  )
  return(results_df)
}




