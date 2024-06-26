## 2022-06-04


# Load packages and source code -------------------------------------------

library("RColorBrewer")
library("PRROC") # for AUCForROCdf



# Define functions --------------------------------------------------------

MakeROCDf <- function(input_df, numeric_column, logical_column = "Is_essential") {
  original_accuracy_mat <- BasicROCMat(input_df[, numeric_column], input_df[, logical_column])
  flipped_accuracy_mat  <- BasicROCMat(input_df[, numeric_column], !(input_df[, logical_column]), reverse_threshold = FALSE)
  colnames(flipped_accuracy_mat) <- paste0(colnames(flipped_accuracy_mat), "_flipped")
  results_df <- data.frame(input_df, original_accuracy_mat,
                           flipped_accuracy_mat[, "Precision_flipped", drop = FALSE], row.names = NULL
                           )
  return(results_df)
}


BasicROCMat <- function(numeric_vec, logical_vec, reverse_threshold = TRUE) {
  stopifnot(length(numeric_vec) == length(logical_vec))
  accuracy_list <- lapply(numeric_vec, function(x) {
    if (reverse_threshold) {
      pass_threshold <- numeric_vec <= x
    } else {
      pass_threshold <- numeric_vec > x
    }
    count_TN <- sum((!(logical_vec)) & (!(pass_threshold)))
    count_TP <- sum(logical_vec & pass_threshold)
    count_FN <- sum(logical_vec & (!(pass_threshold)))
    count_FP <- sum((!(logical_vec)) & pass_threshold)
    c("Specificity" = count_TN / (count_TN + count_FP),
      "Sensitivity" = count_TP / (count_TP + count_FN),
      "Precision"   = count_TP / (count_TP + count_FP)
      )
  })
  accuracy_mat <- do.call(rbind, accuracy_list)
  return(accuracy_mat)
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
                      xlab_line  = 1.8,
                      ylab_line  = 2.3,
                      ROC_lwd    = 2
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
    axis(1, mgp = c(3, 0.4, 0), tcl = use_tcl, lwd = par("lwd"))
    mtext("False positive rate", side = 1, line = xlab_line, cex = par("cex"))
    axis(2, mgp = c(3, 0.5, 0), tcl = use_tcl, las = 1, lwd = par("lwd"))
    mtext("True positive rate", side = 2, line = ylab_line, cex = par("cex"))
    abline(a = 0, b = 1, col = "gray78", lty = "dashed")
    box()
  }

  lines(1 - ROC_mat[, "Specificity"], ROC_mat[, "Sensitivity"],
        lwd = ROC_lwd, col = line_color, xpd = NA
        )

  if (show_AUC) {
    space_in_lines <- 1.1
    text(x      = par("usr")[[2]] - diff(grconvertY(c(0, space_in_lines), from = "lines", to = "user")),
         y      = par("usr")[[3]] + diff(grconvertY(c(0, space_in_lines), from = "lines", to = "user")),
         labels = paste0("AUC = ", format(round(AUC_value, digits = 2), nsmall = 2)),
         adj    = c(1, 0),
         xpd    = NA
         )
  }
  return(invisible(AUC_value))
}



RemoveInfinite <- function(numeric_vec) {
  finite_vec <- numeric_vec[is.finite(numeric_vec)]
  min_value <- min(finite_vec)
  max_value <- max(finite_vec)
  data_span <- max_value - min_value
  to_add <- data_span * 0.01
  numeric_vec[numeric_vec == Inf] <- max_value + to_add
  numeric_vec[numeric_vec == -Inf] <- min_value - to_add
  return(numeric_vec)
}



AUCForROCdf <- function(numeric_vec, logical_vec, flip = TRUE, precision_recall = TRUE, return_curve = FALSE) {
  numeric_vec <- RemoveInfinite(numeric_vec)
  if (flip) {
    scores_class_A <- numeric_vec[!(logical_vec)]
    scores_class_B <- numeric_vec[logical_vec]
  } else {
    numeric_vec <- -(numeric_vec) # Because the PRROC::pr.curve function assumes that controls < cases, the values have to be inverted
    scores_class_A <- numeric_vec[logical_vec]
    scores_class_B <- numeric_vec[!(logical_vec)]
  }
  if (precision_recall) {
    PRROC_object <- PRROC::pr.curve(scores_class_A, scores_class_B, curve = return_curve)
    AUC_value <- PRROC_object[["auc.integral"]]
  } else {
    PRROC_object <- PRROC::roc.curve(scores_class_A, scores_class_B, curve = return_curve)
    AUC_value <- PRROC_object[["auc"]]
  }
  if (return_curve) {
    return(PRROC_object[["curve"]])
  } else {
    return(AUC_value)
  }
}

