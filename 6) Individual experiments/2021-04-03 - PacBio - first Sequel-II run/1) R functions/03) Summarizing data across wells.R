### 13th June 2021 ###





# Functions for calculating additional summary statistics -----------------

AddContaminationSummary <- function(ccs_df_list, contamin_df, summary_df_name = "filtered_summary_df") {

  pass_filters <- ccs_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
  use_zmws <- ccs_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]
  are_to_use <- contamin_df[["ZMW"]] %in% use_zmws

  cross_well_contam_df  <- contamin_df[are_to_use & !(contamin_df[["Are_cross_plate"]]), ]
  cross_plate_contam_df <- contamin_df[are_to_use & contamin_df[["Are_cross_plate"]], ]


  summary_df <- ccs_df_list[[summary_df_name]]

  cross_well_contam_df[["Reference_ID"]]  <- factor(cross_well_contam_df[["Reference_ID"]],
                                                    levels = summary_df[["Combined_ID"]]
                                                    )
  cross_plate_contam_df[["Reference_ID"]] <- factor(cross_plate_contam_df[["Reference_ID"]],
                                                    levels = summary_df[["Combined_ID"]]
                                                    )

  well_contam_splits  <- split(cross_well_contam_df[["ZMW"]], cross_well_contam_df[["Reference_ID"]])
  plate_contam_splits <- split(cross_plate_contam_df[["ZMW"]], cross_plate_contam_df[["Reference_ID"]])

  well_contam_splits <- lapply(well_contam_splits, unique)
  plate_contam_splits <- lapply(plate_contam_splits, unique)

  summary_df[["Num_contaminated_reads_aligned"]] <- lengths(well_contam_splits)
  summary_df[["Num_cross_plate_contaminated"]] <- lengths(plate_contam_splits)


  for (column_name in c("Num_contaminated_reads_aligned", "Num_cross_plate_contaminated")) {
    summary_df[[column_name]] <- ifelse(summary_df[["Count_total"]] == 0,
                                        NA,
                                        summary_df[[column_name]]
                                        )
  }

  ccs_df_list[[summary_df_name]] <- summary_df
  return(ccs_df_list)
}




AddDeletionSummary <- function(ccs_df_list, del_df, summary_df_name = "filtered_summary_df") {

  pass_filters <- ccs_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
  use_zmws <- ccs_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]
  are_to_use <- del_df[["ZMW"]] %in% use_zmws

  del_df  <- del_df[are_to_use, ]

  summary_df <- ccs_df_list[[summary_df_name]]

  del_df[["Combined_ID"]] <- factor(del_df[["Combined_ID"]],
                                    levels = summary_df[["Combined_ID"]]
                                    )

  del_splits <- split(del_df[["ZMW"]], del_df[["Combined_ID"]])
  del_splits <- lapply(del_splits, unique)
  summary_df[["Num_reads_with_deletions_exceeding_20bp"]] <- lengths(del_splits)

  span_columns <- c("Span_tracrRNAs", "Span_promoters", "Span_sg_cr", "Span_sgRNAs")
  span_list <- lapply(span_columns, function(x) {
    are_spanning <- del_df[[x]]
    span_splits <- split(del_df[["ZMW"]][are_spanning],
                         del_df[["Combined_ID"]][are_spanning]
                         )
    span_splits <- lapply(span_splits, unique)
    return(lengths(span_splits))
  })
  span_mat <- do.call(cbind, span_list)
  colnames(span_mat) <- sub("Span_", "Num_reads_with_deletions_spanning_", span_columns, fixed = TRUE)
  summary_df <- data.frame(summary_df, span_mat, stringsAsFactors = FALSE)
  ccs_df_list[[summary_df_name]] <- summary_df
  return(ccs_df_list)
}



AddNoContamCounts <- function(ccs_df_list, contamin_df, summary_df_name = "filtered_summary_df") {

  indiv_df <- ccs_df_list[["individual_reads_df"]]
  summary_df <- ccs_df_list[[summary_df_name]]

  are_contaminated <- (indiv_df[["ZMW"]] %in% contamin_df[["ZMW"]]) |
                      ((indiv_df[["Num_contaminating_guides"]] >= 1) %in% TRUE)

  pass_filters <- indiv_df[["Passes_filters"]] == 1
  no_contam_df <- indiv_df[(!(are_contaminated)) & pass_filters, ]

  binary_columns <- c(paste0("sg", 1:4, "_cr", 1:4),
                      paste0("at_least_", 1:3),
                      "all_4", "all_4_promoters", "whole_plasmid"
                      )
  no_contam_mat <- as.matrix(no_contam_df[, binary_columns])

  wells_fac <- factor(no_contam_df[["Combined_ID"]],
                      levels = summary_df[["Combined_ID"]]
                      )

  summary_vec_list <- tapply(seq_len(nrow(no_contam_mat)),
                             wells_fac,
                             function(x) colSums(no_contam_mat[x, , drop = FALSE])
                             )
  are_empty <- vapply(summary_vec_list, is.null, logical(1))
  if (any(are_empty)) {
    summary_vec_list[are_empty] <- list(integer(ncol(no_contam_mat)))
  }

  summary_counts_mat <- do.call(rbind, summary_vec_list)
  mode(summary_counts_mat) <- "integer"

  colnames(summary_counts_mat) <- paste0("Count_no_contam_", colnames(summary_counts_mat))
  summary_df <- data.frame(summary_df,
                           "Count_total_no_contam" = tabulate(wells_fac),
                           summary_counts_mat,
                           stringsAsFactors = FALSE
                           )

  ccs_df_list[[summary_df_name]] <- summary_df
  return(ccs_df_list)
}







# Helper functions for creating plots -------------------------------------

DrawGridlines <- function(y_limits, extra_grid_lines = TRUE) {
  y_range <- y_limits[[2]] - y_limits[[1]]
  divide_by <- 20
  grid_seq <- seq(from = y_limits[[1]], to = y_limits[[2]], by = y_range / divide_by)
  are_main <- rep(c(TRUE, FALSE), divide_by)[seq_along(grid_seq)]
  segments(x0   = par("usr")[[1]],
           x1   = par("usr")[[2]],
           y0   = grid_seq[are_main],
           col  = "gray88",
           lend = "butt",
           xpd  = NA
           )
  if (extra_grid_lines) {
    segments(x0   = par("usr")[[1]],
             x1   = par("usr")[[2]],
             y0   = grid_seq[!(are_main)],
             col  = "gray95",
             lend = "butt",
             xpd  = NA
             )
  }
  return(invisible(NULL))
}



