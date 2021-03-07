### 27th January 2021 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

file_directory        <- "~/CRISPR/6) Individual experiments/2021-01-27 - calculate gene-wise p values for a pooled CRISPR screen"
file_input_directory  <- file.path(file_directory, "1) Input")
file_output_directory <- file.path(file_directory, "2) Output")
RData_directory       <- file.path(file_directory, "3) R objects")




# Read in data ------------------------------------------------------------

file_names_1 <- list.files(file.path(file_input_directory, "Sent by Tingting"))

df_list_1 <- lapply(file_names_1, function(x) {
  read.csv(file.path(file_input_directory,"Sent by Tingting", x),
           check.names = FALSE,
           stringsAsFactors = FALSE
           )
})
names(df_list_1) <- sub(".csv", "", file_names_1, fixed = TRUE)


file_names_2 <- list.files(file.path(file_input_directory, "Sent by Davide", "2) After TMM"))
file_names_2_splits <- sapply(strsplit(file_names_2, "_without", fixed = TRUE), "[", 2)
are_selected <- grepl("NBH", file_names_2_splits, fixed = TRUE) &
                grepl("RML6", file_names_2_splits, fixed = TRUE)
file_names_2 <- file_names_2[are_selected]

df_list_2 <- lapply(file_names_2, function(x) {
  read.csv(file.path(file_input_directory,"Sent by Davide", "2) After TMM", x),
           check.names = FALSE,
           stringsAsFactors = FALSE
           )
})



# Define functions --------------------------------------------------------

SummarizeSgRNAs <- function(sub_df, control_p_values, num_top_guides = 2L) {
  assign("delete_sub_df", sub_df, envir = globalenv())
  unique_transcripts <- unique(sub_df[["transcript_id"]])
  num_transcripts <- length(unique_transcripts)
  num_guides <- nrow(sub_df)
  gene_symbol <- unique(sub_df[["gene_name"]])
  stopifnot(length(gene_symbol) == 1)
  gene_ID <- unique(sub_df[["gene_id"]])
  stopifnot(length(gene_ID) == 1)
  new_order <- order(sub_df[["pValue"]], decreasing = FALSE)
  use_num_guides <- min(num_top_guides, num_guides)
  logfc <- mean(sub_df[["log2_dcare"]][new_order][seq_len(use_num_guides)])
  wilcox_results <- wilcox.test(sub_df[["pValue"]][new_order][seq_len(use_num_guides)],
                                control_p_values,
                                alternative = "less"
                                )
  p_value <- wilcox_results[["p.value"]]
  results_list <- list(
    "Ensembl_gene_ID"       = gene_ID,
    "Gene_symbol"           = gene_symbol,
    "Transcript_IDs"        = paste0(unique_transcripts, collapse = ", "),
    "Number_of_guides"      = num_guides,
    "Number_of_transcripts" = num_transcripts,
    "Log_fold_change"       = logfc,
    "Unadjusted_p_value"    = p_value
  )
  return(results_list)
}



SummarizeDf <- function(input_df,
                        control_p_vec,
                        use_num_guides = 2L,
                        split_by_column = "transcript_id"
                        ) {

  include_columns <- c("gene_name", "gene_id", "transcript_id",
                       "log2_dcare", "pValue"
                       )
  split_df_list <- split(input_df[, include_columns],
                         input_df[, split_by_column]
                         )

  summarized_list <- lapply(split_df_list,
                            SummarizeSgRNAs,
                            control_p_values = control_p_vec,
                            num_top_guides = use_num_guides
                            )

  results_df <- do.call(rbind.data.frame,
                        c(summarized_list,
                          list(stringsAsFactors = FALSE,
                               make.row.names = FALSE
                          )
                        ))

  results_df[["Hit_strength"]] <- results_df[["Log_fold_change"]] *
                                  (-log10(results_df[["Unadjusted_p_value"]]))

  return(results_df)
}




SimulateAndSummarize <- function(targeting_df,
                                 control_df,
                                 use_num_guides = 3L,
                                 use_num_repetitions = 2L,
                                 split_by_column = "transcript_id",
                                 use_prefix = NULL,
                                 calculate_q_values = TRUE,
                                 hit_strength_q_values = TRUE
                                 ) {
  results_df <- SummarizeDf(targeting_df,
                            control_df[["pValue"]],
                            use_num_guides = use_num_guides,
                            split_by_column = split_by_column
                            )
  if (calculate_q_values) {
    sim_hit_strengths <- GetSimulatedHitStrengths(results_df[, "Number_of_guides"],
                                                  control_df,
                                                  use_num_guides = use_num_guides,
                                                  num_repetitions = use_num_repetitions
                                                  )

    HS_q_values_df <- GetHitStrengthQValues(results_df[, "Hit_strength"],
                                            sim_hit_strengths
                                            )

    sim_p_values <- GetSimulatedPValues(results_df[, "Number_of_guides"],
                                        control_df[["pValue"]],
                                        use_num_guides = use_num_guides,
                                        num_repetitions = use_num_repetitions
                                        )

    p_q_values_df <- GetQValues(results_df[, "Unadjusted_p_value"],
                                sim_p_values
                                )

    results_df <- data.frame(results_df, HS_q_values_df, p_q_values_df,
                             stringsAsFactors = FALSE
                             )
  }

  if (split_by_column == "transcript_id") {
    results_df <- results_df[, names(results_df) != "Number_of_transcripts"]
  }
  if (!(is.null(use_prefix))) {
    prefix_columns <- c("Log_fold_change", "Unadjusted_p_value", "Hit_strength",
                        "Uncorrected_HS_q_value", "Hit_strength_q_value",
                        "Uncorrected_q_value", "Q_value"
                        )
    prefix_columns <- intersect(prefix_columns, names(results_df))
    new_prefix_columns <- paste0(use_prefix, tolower(prefix_columns))
    names(results_df)[match(prefix_columns, names(results_df))] <- new_prefix_columns
  }
  results_list <- list("summary_df" = results_df,
                       "simulated_p_values" = sim_p_values,
                       "simulated_hit_strengths" = sim_hit_strengths
                       )
  return(results_list)
}




GetSimulatedPValues <- function(num_guides_vec,
                                control_p_vec,
                                use_num_guides = 2L,
                                num_repetitions = 5
                                ) {
  set.seed(1)
  simulated_p_values <- vapply(rep(num_guides_vec, num_repetitions),
                               function(x) {
                                 random_sample <- sample(control_p_vec, use_num_guides)
                                 wilcox_results <- wilcox.test(random_sample,
                                                               control_p_vec,
                                                               alternative = "less"
                                                               )
                                 return(wilcox_results[["p.value"]])
                               },
                               numeric(1)
                               )
  simulated_p_values <- sort(simulated_p_values)
  return(simulated_p_values)
}




GetSimulatedHitStrengths <- function(num_guides_vec,
                                     control_df,
                                     use_num_guides = 2L,
                                     num_repetitions = 5
                                     ) {
  set.seed(1)
  index_vec <- seq_len(nrow(control_df))
  simulated_hit_strengths <- vapply(rep(num_guides_vec, num_repetitions),
                                    function(x) {
                                      random_indices <- sample(index_vec, use_num_guides)
                                      log_fc <- mean(control_df[["log2_dcare"]][random_indices])
                                      wilcox_results <- wilcox.test(control_df[["pValue"]][random_indices],
                                                                    control_df[["pValue"]],
                                                                    alternative = "less"
                                                                    )
                                      log_p <- -(log10(wilcox_results[["p.value"]]))
                                      return(log_fc * log_p)
                                    },
                                    numeric(1)
                                    )
  simulated_hit_strengths <- sort(simulated_hit_strengths)
  return(simulated_hit_strengths)
}



GetHitStrengthQValues <- function(observed_HS, sim_HS) {
  num_ratio <- length(observed_HS) / length(sim_HS)
  observed_HS <- abs(observed_HS)
  sim_HS <- abs(sim_HS)
  q_values_vec <- vapply(observed_HS,
                         function(x) {
                           num_observed <- sum(observed_HS >= x)
                           num_sim <- sum(sim_HS >= x) * num_ratio
                           ratio <- num_sim / (num_sim + num_observed)
                           return(ratio)
                         },
                         numeric(1)
                         )
  corrected_values <- vapply(seq_along(observed_HS), function(x) {
    are_lower_ranked <- observed_HS <= observed_HS[[x]]
    min(q_values_vec[are_lower_ranked])
  }, numeric(1))
  results_df <- data.frame(
    "Uncorrected_HS_q_value" = q_values_vec,
    "Hit_strength_q_value" = corrected_values
  )
  return(results_df)
}




GetQValues <- function(observed_p, sim_p) {
  num_ratio <- length(observed_p) / length(sim_p)
  q_values_vec <- vapply(observed_p,
                         function(x) {
                           num_observed <- sum(observed_p <= x)
                           num_sim <- sum(sim_p <= x) * num_ratio
                           ratio <- num_sim / (num_sim + num_observed)
                           return(ratio)
                         },
                         numeric(1)
                         )

  corrected_values <- vapply(seq_along(observed_p), function(x) {
    are_lower_ranked <- observed_p >= observed_p[[x]]
    min(q_values_vec[are_lower_ranked])
  }, numeric(1))
  results_df <- data.frame(
    "Uncorrected_q_value" = q_values_vec,
    "Q_value" = corrected_values
  )
  return(results_df)
}



CompleteAnalysis <- function(use_df,
                             split_by_column = "transcript_id",
                             num_repetitions = 20L,
                             calculate_q_values = FALSE,
                             order_by_top2 = FALSE
                             ) {

  are_present <- use_df[["isPresent"]] == "TRUE"
  are_controls <- use_df[["transcript_id"]] == "nontargetting.control"

  gene_IDs_df <- unique(use_df[!(are_controls), c("gene_id", "gene_name")])
  have_multiple_symbols <- (table(gene_IDs_df[["gene_id"]])[gene_IDs_df[["gene_id"]]] > 1) %in% TRUE
  duplicates_df <- gene_IDs_df[have_multiple_symbols, ]
  duplicates_df <- duplicates_df[order(duplicates_df[["gene_id"]]), ]
  stopifnot(!(any(duplicated(duplicates_df[["gene_name"]]))))
  num_ambiguous_IDs <- length(unique(duplicates_df[["gene_id"]]))
  message(paste0(num_ambiguous_IDs,
                 " Ensembl gene IDs were mapped to multiple gene symbols!"
                 )
          )


  have_no_ID <- !(are_controls) & is.na(use_df[["gene_name"]])
  message(paste0(length(unique(use_df[["transcript_id"]][have_no_ID])),
                 " unique transcripts were not associated with a gene symbol!"
                 )
          )

  control_df <- use_df[are_controls & are_present, c("pValue", "log2_dcare")]
  targeting_df <- use_df[!(are_controls) & are_present, ]

  top2_list <- SimulateAndSummarize(targeting_df,
                                    control_df,
                                    use_num_guides = 2L,
                                    use_num_repetitions = num_repetitions,
                                    split_by_column = split_by_column,
                                    use_prefix = "Top2_",
                                    calculate_q_values = calculate_q_values
                                    )

  top3_list <- SimulateAndSummarize(targeting_df,
                                    control_df,
                                    use_num_guides = 3L,
                                    use_num_repetitions = num_repetitions,
                                    split_by_column = split_by_column,
                                    use_prefix = "Top3_",
                                    calculate_q_values = calculate_q_values
                                    )

  results_df <- data.frame(
    top2_list[["summary_df"]],
    top3_list[["summary_df"]][, !(names(top3_list[["summary_df"]]) %in% names(top2_list[["summary_df"]]))],
    stringsAsFactors = FALSE
  )
  if (order_by_top2) {
    order_columns <- c("Top2_hit_strength", "Top2_unadjusted_p_value")
  } else {
    order_columns <- c("Top3_hit_strength", "Top3_unadjusted_p_value")
  }
  new_order <- order(-(abs(results_df[, order_columns[[1]]])),
                     results_df[, order_columns[[2]]]
                     )
  results_df <- results_df[new_order, ]
  row.names(results_df) <- NULL

  results_list <- list(
    "summary_df"                   = results_df,
    "top2_simulated_p_values"      = top2_list[["simulated_p_values"]],
    "top2_simulated_hit_strengths" = top2_list[["simulated_hit_strengths"]],
    "top3_simulated_p_values"      = top3_list[["simulated_p_values"]],
    "top3_simulated_hit_strengths" = top3_list[["simulated_hit_strengths"]]
  )

  return(results_list)
}




WriteTable <- function(use_df, file_path) {
  write.table(use_df,
              file = file_path,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE
              )
}


ReorderByTop3 <- function(top_df) {
  new_order <- order(top_df[, "Top3_hit_strength_q_value"],
                     top_df[, "Top3_unadjusted_p_value"]
                     )
  top_df <- top_df[new_order, ]
  row.names(top_df) <- NULL
  return(top_df)
}








PlotPValueHistograms <- function(use_list,
                                 num_guides = 2L,
                                 num_repetitions = 20L,
                                 main_title = "",
                                 passage_y_position = 0.7
                                 ) {

  simulated_hist <- lapply(use_list, function(x) {
    hist_results <- hist(x[[paste0("top", num_guides, "_simulated_p_values")]],
                         breaks = 50,
                         plot = FALSE
                         )
    hist_results[["counts"]] <- hist_results[["counts"]] / num_repetitions
    return(hist_results)
  })
  real_hist <- lapply(use_list, function(x) {
    hist_results <- hist(x[["summary_df"]][[paste0("Top", num_guides, "_unadjusted_p_value")]],
                         breaks = 50,
                         plot = FALSE
                         )
    return(hist_results)
  })

  combined_hist <- c(real_hist, simulated_hist)

  max_count <- max(vapply(combined_hist, function(x) max(x[["counts"]]), numeric(1)))
  # use_limits <- GetFlexibleAxisLimits(c(0, max_count))
  use_limits <- range(0, max_count)
  use_limits[[1]] <- -(use_limits[[2]] * 0.01)
  use_color <- brewer.pal(9, "Blues")[[7]]

  use_hist <- combined_hist[[1]]


  layout_mat <- matrix(c(1:2, rep(3:7, each = 2)), ncol = 2, byrow = TRUE)
  layout_mat <- rbind(layout_mat[1:2, ],
                      max(layout_mat) + c(1, 5),
                      layout_mat[3, ],
                      max(layout_mat) + c(2, 6),
                      layout_mat[4, ],
                      max(layout_mat) + c(3, 7),
                      layout_mat[5, ],
                      max(layout_mat) + c(4, 8),
                      layout_mat[6, ]
                      )

  layout_mat <- cbind(max(layout_mat) + 1, layout_mat)

  heights_vec <- rep(1, nrow(layout_mat))
  heights_vec[[1]] <- 1.5
  heights_vec[[10]] <- 0.8
  heights_vec[c(3, 5, 7, 9)] <- 2.5


  layout(layout_mat,
         heights = heights_vec,
         widths = c(0.2, 1, 1)
         )
  par(mar = c(0, 2, 0, 2))

  for (use_text in c("Real", "Simulated")) {
    MakeEmptyPlot()
    text(as.expression(bquote(bold(.(use_text) ~ bolditalic("p") ~ "values"))),
         y = -0.3, x = 0.5, cex = 1.6, xpd = NA
         )
  }

  MakeEmptyPlot()
  text(main_title, x = 0.5, y = 1.75, xpd = NA, cex = 1.7, font = 1)

  for (i in 1:4) {
    MakeEmptyPlot()
    # text(paste0("Passage ", c(0, 2, 4, 8)[[i]]),
    #      y = 0.1, x = 0.5, cex = 1.4, font = 2,
    #      col = "gray40", xpd = NA
    #      )
  }

  for (i in seq_along(combined_hist)) {
    plot(combined_hist[[i]],
         ylim   = use_limits,
         las    = 1,
         col    = use_color,
         border = NA,
         tcl    = -0.3,
         xaxs   = "i",
         xlim   = c(-0.02, 1.02),
         ylab   = "",
         xlab   = "",
         main   = "",
         xpd    = NA,
         axes   = FALSE
         )

    if (i %in% 1:4) {
      text(paste0("Passage ", c(0, 2, 4, 8))[[i]],
           x   = par("usr")[[1]] - ((par("usr")[[2]] - par("usr")[[1]]) * 0.2),
           y   = par("usr")[[3]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.5),
           cex = 1.5,
           font = 2,
           col = "gray20",
           xpd = NA,
           srt = 90,
           )
    }
    axis_lwd <- 0.75

    axis(1, mgp = c(2, 0.45, 0), tcl = -0.35, lwd = axis_lwd)
    axis(2, mgp = c(2, 0.45, 0), las = 1, tcl = -0.35, lwd = axis_lwd,
         col.axis = if (i %in% 1:4) "black" else "black"
         )
    segments(x0  = par("usr")[[1]],
             x1  = par("usr")[[1]],
             y0  = 0,
             y1  = par("usr")[[3]],
             xpd = NA,
             lwd = axis_lwd
             )
    segments(x0  = par("usr")[[1]],
             x1  = 0,
             y0  = par("usr")[[3]],
             y1  = par("usr")[[3]],
             xpd = NA,
             lwd = axis_lwd
             )
  }

  par(mar = rep(0, 4))
  MakeEmptyPlot()
}









# Perform calculations ----------------------------------------------------

num_simulations <- 100L

by_transcript_list_1 <- lapply(df_list_1,
                               function(x) CompleteAnalysis(x,
                                                            split_by_column = "transcript_id",
                                                            num_repetitions = num_simulations,
                                                            calculate_q_values = TRUE,
                                                            order_by_top2 = TRUE
                                                            )
                               )

by_symbol_list_1 <- lapply(df_list_2,
                           function(x) CompleteAnalysis(x,
                                                        split_by_column = "gene_name",
                                                        num_repetitions = num_simulations,
                                                        calculate_q_values = TRUE,
                                                        order_by_top2 = TRUE
                                                        )
                           )


by_transcript_list_2 <- lapply(df_list_2,
                               function(x) CompleteAnalysis(x,
                                                            split_by_column = "transcript_id",
                                                            num_repetitions = num_simulations,
                                                            calculate_q_values = TRUE,
                                                            order_by_top2 = TRUE
                                                            )
                               )




# Create plots ------------------------------------------------------------

for (make_PDF in c(FALSE, TRUE)) {

  if (make_PDF) {
    pdf(file = file.path(file_output_directory, "Sample selection 1", "P value distribution.pdf"),
        height = 8.1, width = 7.5
        )
  }

  PlotPValueHistograms(by_transcript_list_1,
                       num_guides = 2,
                       main_title = "Top 2 sgRNAs per transcript",
                       num_repetitions = num_simulations
                       )

  PlotPValueHistograms(by_transcript_list_1,
                       num_guides = 3,
                       main_title = "Top 3 sgRNAs per transcript",
                       passage_y_position = 0.92,
                       num_repetitions = num_simulations
                       )

  PlotPValueHistograms(by_symbol_list_1,
                       num_guides = 2,
                       main_title = "Top 2 sgRNAs per gene symbol",
                       num_repetitions = num_simulations
                       )

  PlotPValueHistograms(by_symbol_list_2,
                       num_guides = 3,
                       main_title = "Top 3 sgRNAs per gene symbol",
                       passage_y_position = 0.92,
                       num_repetitions = num_simulations
                       )

  if (make_PDF) {
    dev.off()
  }

  if (make_PDF) {
    pdf(file = file.path(file_output_directory, "Sample selection 2", "P value distribution.pdf"),
        height = 8.1, width = 7.5
        )
  }

  PlotPValueHistograms(by_transcript_list_2,
                       num_guides = 2,
                       main_title = "Top 2 sgRNAs per transcript",
                       num_repetitions = num_simulations
                       )

  PlotPValueHistograms(by_transcript_list_2,
                       num_guides = 3,
                       main_title = "Top 3 sgRNAs per transcript",
                       passage_y_position = 0.92,
                       num_repetitions = num_simulations
                       )


  if (make_PDF) {
    dev.off()
  }
}




# Export data -------------------------------------------------------------

for (i in seq_along(df_list_1)) {
  WriteTable(by_transcript_list_1[[i]][["summary_df"]],
             file = file.path(file_output_directory,
                              "Sample selection 1",
                              "Ordered by top 2 sgRNAs",
                              "By transcript",
                              paste0(names(df_list_1)[[i]], "_top2_by_transcript.tsv")
                              )
             )
}


for (i in seq_along(df_list_1)) {
  WriteTable(by_symbol_list_1[[i]][["summary_df"]],
             file = file.path(file_output_directory,
                              "Sample selection 1",
                              "Ordered by top 2 sgRNAs",
                              "By gene symbol",
                              paste0(names(df_list_1)[[i]], "_top2_by_gene_symbol.tsv")
                              )
             )
}


for (i in seq_along(df_list_1)) {
  WriteTable(ReorderByTop3(by_transcript_list_1[[i]][["summary_df"]]),
             file = file.path(file_output_directory,
                              "Sample selection 1",
                              "Ordered by top 3 sgRNAs",
                              "By transcript",
                              paste0(names(df_list_1)[[i]], "_top3_by_transcript.tsv")
                              )
             )
}


for (i in seq_along(df_list_1)) {
  WriteTable(ReorderByTop3(by_symbol_list_1[[i]][["summary_df"]]),
             file = file.path(file_output_directory,
                              "Sample selection 1",
                              "Ordered by top 3 sgRNAs",
                              "By gene symbol",
                              paste0(names(df_list_1)[[i]], "_top3_by_gene_symbol.tsv")
                              )
             )
}




for (i in seq_along(df_list_2)) {
  WriteTable(by_transcript_list_2[[i]][["summary_df"]],
             file = file.path(file_output_directory,
                              "Sample selection 1",
                              "Ordered by top 2 sgRNAs",
                              "By transcript",
                              paste0(names(df_list_2)[[i]], "_top2_by_transcript.tsv")
                              )
             )
}






# Save data ---------------------------------------------------------------

save(list = c("by_transcript_list_1", "by_symbol_list_1",
              "by_transcript_list_2", "by_symbol_list_2"
              ),
     file = file.path(RData_directory, "1) Calculate gene-wise p values for a pooled CRISPR screen.RData")
     )







#
#
#
# # Try stuff ---------------------------------------------------------------
#
#
# use_df <- df_list[[1]]
# split_by_column = "transcript_id"
# num_repetitions = 20L
# calculate_q_values = FALSE
# order_by_top2 = FALSE
#
# are_present <- use_df[["isPresent"]] == "TRUE"
# are_controls <- use_df[["transcript_id"]] == "nontargetting.control"
#
# gene_IDs_df <- unique(use_df[!(are_controls), c("gene_id", "gene_name")])
# have_multiple_symbols <- (table(gene_IDs_df[["gene_id"]])[gene_IDs_df[["gene_id"]]] > 1) %in% TRUE
# duplicates_df <- gene_IDs_df[have_multiple_symbols, ]
# duplicates_df <- duplicates_df[order(duplicates_df[["gene_id"]]), ]
# stopifnot(!(any(duplicated(duplicates_df[["gene_name"]]))))
# num_ambiguous_IDs <- length(unique(duplicates_df[["gene_id"]]))
# message(paste0(num_ambiguous_IDs,
#                " Ensembl gene IDs were mapped to multiple gene symbols!"
#                ))
#
#
# have_no_ID <- !(are_controls) & is.na(use_df[["gene_name"]])
# message(paste0(length(unique(use_df[["transcript_id"]][have_no_ID])),
#                " unique transcripts were not associated with a gene symbol!"
#                ))
#
#
# control_df <- use_df[are_controls & are_present, c("pValue", "log2_dcare")]
# targeting_df <- use_df[!(are_controls) & are_present, ]
#
#
#
# use_num_guides = 2L
# use_num_repetitions = 20L
# split_by_column = "transcript_id"
# use_prefix = NULL
# calculate_q_values = TRUE
#
#
# results_df <- SummarizeDf(targeting_df,
#                           control_df[["pValue"]],
#                           use_num_guides = use_num_guides,
#                           split_by_column = split_by_column
#                           )
# if (calculate_q_values) {
#   sim_hit_strengths <- GetSimulatedHitStrengths(results_df[, "Number_of_guides"],
#                                                 control_df,
#                                                 use_num_guides = use_num_guides,
#                                                 num_repetitions = use_num_repetitions
#                                                 )
#
#   HS_q_values_df <- GetHitStrengthQValues(results_df[, "Hit_strength"],
#                                           sim_hit_strengths
#                                           )
#
#   sim_p_values <- GetSimulatedPValues(results_df[, "Number_of_guides"],
#                                       control_df[["pValue"]],
#                                       use_num_guides = use_num_guides,
#                                       num_repetitions = use_num_repetitions
#                                       )
#
#   p_q_values_df <- GetQValues(results_df[, "Unadjusted_p_value"],
#                               sim_p_values
#                               )
#
#   results_df <- data.frame(results_df, HS_q_values_df, p_q_values_df,
#                            stringsAsFactors = FALSE
#                            )
# }
#
#
#
#
# if (split_by_column == "transcript_id") {
#   results_df <- results_df[, names(results_df) != "Number_of_transcripts"]
# }
# if (!(is.null(use_prefix))) {
#   prefix_columns <- intersect(c("Log_fold_change", "Unadjusted_p_value", "Q_value"), names(results_df))
#   new_prefix_columns <- paste0(use_prefix, tolower(prefix_columns))
#   names(results_df)[match(prefix_columns, names(results_df))] <- new_prefix_columns
# }
#
#




