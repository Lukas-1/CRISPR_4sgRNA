## 2022-06-04


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define functions --------------------------------------------------------

MakeROCDf <- function(input_df, numeric_column, logical_column = "Is_essential") {
  input_df[, "Specificity"] <- vapply(input_df[, numeric_column], function(x) {
    pass_threshold <- input_df[, numeric_column] <= x
    count_TN <- sum((!(input_df[, logical_column])) & (!(pass_threshold)))
    count_FP <- sum((!(input_df[, logical_column])) & pass_threshold)
    count_TN / (count_TN + count_FP)
  }, numeric(1))
  input_df[, "Sensitivity"] <- vapply(input_df[, numeric_column], function(x) {
    pass_threshold <- input_df[, numeric_column] <= x
    count_TP <- sum(input_df[, logical_column] & pass_threshold)
    count_FN <- sum(input_df[, logical_column] & (!(pass_threshold)))
    count_TP / (count_TP + count_FN)
  }, numeric(1))
  return(input_df)
}



GetROCMat <- function(input_df) {
  sens_vec <- input_df[, "Sensitivity"]
  spec_vec <- input_df[, "Specificity"]
  stopifnot(identical(sens_vec, sort(sens_vec)))
  stopifnot(identical(spec_vec, sort(spec_vec, decreasing = TRUE)))
  if (!(any((sens_vec == 0) & (spec_vec == 1)))) {
    sens_vec <- c(0, sens_vec)
    spec_vec <- c(1, spec_vec)
  }
  return(cbind("Sensitivity" = sens_vec, "Specificity" = spec_vec))
}



GetAUC <- function(use_ROC_df) {
  use_mat <- GetROCMat(use_ROC_df)
  x_vec <- 1 - use_mat[, "Specificity"]
  y_vec <- use_mat[, "Sensitivity"]
  indices_vec <- seq_len(length(x_vec) - 1)
  x_starts <- x_vec[indices_vec]
  x_ends <- x_vec[indices_vec + 1]
  rect_areas <- (x_ends - x_starts) * y_vec[indices_vec]
  auc_result <- sum(rect_areas)
  return(auc_result)
}



PlotROCDf <- function(use_ROC_df,
                      flip       = FALSE,
                      add        = FALSE,
                      line_color = RColorBrewer::brewer.pal(9, "Blues")[[7]],
                      show_AUC   = !(add),
                      use_title  = "Gene essentiality with CRISPRoff"
                      ) {

  ROC_mat <- GetROCMat(use_ROC_df)

  if (flip) {
    sens_vec <- ROC_mat[, "Sensitivity"]
    spec_vec <- ROC_mat[, "Specificity"]
    ROC_mat[, "Sensitivity"] <- spec_vec
    ROC_mat[, "Specificity"] <- sens_vec
  }

  AUC_value <- GetAUC(use_ROC_df)

  if (!(add)) {
    axis_limits <- c(0, 1)
    use_tcl <- -0.36

    plot(1, type = "n", ann = FALSE, axes = FALSE,
         xlim = axis_limits, ylim = axis_limits, xaxs = "i", yaxs = "i"
         )
    axis(1, mgp = c(3, 0.4, 0), tcl = use_tcl)
    mtext("False positive rate", side = 1, line = 1.8, cex = par("cex"))
    axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1)
    mtext("True positive rate", side = 2, line = 2.3, cex = par("cex"))
    abline(a = 0, b = 1, col = "gray78", lty = "dashed")
    box()
  }

  lines(1 - ROC_mat[, "Specificity"], ROC_mat[, "Sensitivity"],
        lwd = 2, col = line_color, xpd = NA
        )

  if (show_AUC) {
    space_in_lines <- 1.1
    text(x      = par("usr")[[2]] - diff(grconvertY(c(0, space_in_lines), from = "lines", to = "user")),
         y      = par("usr")[[3]] + diff(grconvertY(c(0, space_in_lines), from = "lines", to = "user")),
         labels = paste0("AUC = ", round(AUC_value, digits = 2)),
         adj    = c(1, 0),
         xpd    = NA
         )
  }
  return(invisible(AUC_value))
}



