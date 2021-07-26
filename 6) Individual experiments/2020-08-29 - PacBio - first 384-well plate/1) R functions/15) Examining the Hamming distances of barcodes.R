### 4th April 2021 ###




# Import packages and source code -----------------------------------------

library("beeswarm")




# Define functions --------------------------------------------------------

GetDistances <- function(sequence_vec) {

  common_length <- unique(nchar(sequence_vec))
  stopifnot(length(common_length) == 1)

  results_df <- data.frame(t(combn(seq_along(sequence_vec), 2)))
  names(results_df) <- paste0("Index_", 1:2)

  results_df[["Sequence_1"]] <- sequence_vec[results_df[["Index_1"]]]
  results_df[["Sequence_2"]] <- sequence_vec[results_df[["Index_2"]]]

  results_df[["Physical_distance"]] <- results_df[[2]] - results_df[[1]]

  seq_chars <- strsplit(sequence_vec, "")
  seq_common <- mapply(function(x, y) sum(seq_chars[[x]] == seq_chars[[y]]),
                       results_df[["Index_1"]],
                       results_df[["Index_2"]]
                       )

  results_df[["Sequence_distance"]] <- common_length - seq_common

  return(results_df)
}



PlotDistances <- function(distances_df, title_unit = "", point_cex = 0.45, use_spacing = 1) {
  beeswarm(split(distances_df[, "Physical_distance"],
                 distances_df[, "Sequence_distance"]
                 ),
           las     = 1,
           cex     = point_cex,
           spacing = use_spacing,
           pch     = 16,
           mgp     = c(2.4, 0.6, 0),
           tcl     = -0.4,
           ylab    = "Physical distance",
           xlab    = "Hamming distance of the sequences (number of differing bases)",
           ylim    = c(0, max(distances_df[["Physical_distance"]]) + 1),
           yaxs    = "i",
           method  = "center",
           col     = toupper("#3056db")
           )
  cor_results <- suppressWarnings(cor.test(distances_df[["Physical_distance"]],
                                           distances_df[["Sequence_distance"]],
                                           method = "spearman"
                                           ))

  spearmans_rho <- format(signif(cor_results[["estimate"]][[1]], 1),
                          scientific = 10
                          )
  p_value <- format(signif(cor_results[["p.value"]], 1),
                    scientific = 10
                    )

  title_string <- bquote(plain(.(as.character(nrow(distances_df))) *
                         .(title_unit) * "(Spearman's " *
                         italic(rho) * " = " * .(spearmans_rho) *
                         ", " * italic("p") * " = " * .(p_value)) * ")"
                         )
  title(title_string, cex.main = 1.1)
}


