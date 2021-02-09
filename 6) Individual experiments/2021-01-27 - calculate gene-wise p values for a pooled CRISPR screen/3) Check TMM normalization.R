### 6th February 2021 ###





# Import packages and source code -----------------------------------------

library("RColorBrewer")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))





# Define folder paths -----------------------------------------------------

file_directory        <- "~/CRISPR/6) Individual experiments/2021-01-27 - calculate gene-wise p values for a pooled CRISPR screen"
file_input_directory  <- file.path(file_directory, "1) Input", "Sent by Davide")
file_output_directory <- file.path(file_directory, "2) Output")
RData_directory       <- file.path(file_directory, "3) R objects")




# Read in data ------------------------------------------------------------

tsv_file_names <- list.files(file.path(file_input_directory, "1) Before TMM"))
before_df_list <- lapply(tsv_file_names, function(x) {
  read.table(file.path(file_input_directory, "1) Before TMM", x),
             header = TRUE,
             check.names = FALSE,
             stringsAsFactors = FALSE,
             sep = "\t"
             )
})
names(before_df_list) <- sub(".txt", "", tsv_file_names, fixed = TRUE)


csv_file_names <- list.files(file.path(file_input_directory, "2) After TMM"))
after_df_list <- lapply(csv_file_names, function(x) {
  read.csv(file.path(file_input_directory, "2) After TMM", x),
           check.names = FALSE,
           stringsAsFactors = FALSE
           )
})
names(after_df_list) <- sub(".csv", "", csv_file_names, fixed = TRUE)






# Define functions --------------------------------------------------------

CountBarPlots <- function(before_after_list) {

  combined_list <- unlist(before_after_list, recursive = FALSE)
  combined_list <- c(combined_list[grepl("before", names(combined_list), fixed = TRUE)],
                     combined_list[grepl("after", names(combined_list), fixed = TRUE)]
                     )

  use_limits <- c(0, 50)
  use_color <- brewer.pal(9, "Blues")[[7]]

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
         widths = c(0.15, 1, 1)
         )
  par(mar = c(0, 3, 0, 2))

  for (use_text in c("Before normalization", "After normalization")) {
    MakeEmptyPlot()
    text(labels = use_text, font = 2,
         y = -0.3, x = 0.5, cex = 1.6, xpd = NA
         )
  }

  MakeEmptyPlot()
  text("Sum of all sgRNA counts",
       x = 0.5, y = 1.75, xpd = NA, cex = 1.8, font = 1
       )

  for (i in 1:4) {
    MakeEmptyPlot()
  }

  for (i in seq_along(combined_list)) {
    pos <- barplot(combined_list[[i]] / 10^6,
                   ylim   = use_limits,
                   las    = 1,
                   col    = ifelse(grepl("NBH", names(combined_list[[i]]), fixed = TRUE),
                                   brewer.pal(9, "Blues")[[7]],
                                   brewer.pal(9, "Purples")[[7]]
                   ),
                   border = NA,
                   tcl    = -0.3,
                   # xaxs   = "i",
                   # xlim   = c(-0.02, 1.02),
                   ylab   = "",
                   xlab   = "",
                   main   = "",
                   xpd    = NA,
                   axes   = FALSE,
                   space  = 0.5,
                   names.arg = rep("", length(combined_list[[i]]))
                   )
    if (i %in% 1:4) {
      text(paste0("Passage ", c(0, 2, 4, 8))[[i]],
           x   = par("usr")[[1]] - ((par("usr")[[2]] - par("usr")[[1]]) * 0.215),
           y   = par("usr")[[3]] + ((par("usr")[[4]] - par("usr")[[3]]) * 0.5),
           cex = 1.5,
           font = 2,
           col = "gray20",
           xpd = NA,
           srt = 90,
           )
    }
    axis_positions <- axTicks(2)
    mtext(sub("_rep", "\u2013", names(combined_list[[i]])),
          line = 0.5, side = 1, at = pos[, 1], xpd = NA, cex = par("cex"),
          )

    axis(2, mgp = c(2, 0.45, 0), las = 1, tcl = -0.35,
         col = "gray50",
         col.axis = if (i %in% 1:4) "black" else "black",
         at = axis_positions, labels = paste0(axis_positions, "M"),
         lwd = 0.75
         )
    box(bty = "l", col = "gray50", lwd = 0.75)
  }
  par(mar = rep(0, 4))
  MakeEmptyPlot()

}



# Draw scatter plots ------------------------------------------------------

plot(after_df_list[["P0_RML6NHB_comparison"]][["2-P0-NBH_GT_[normalized count]"]],
     before_df_list[["2-P0-NBH_GT-result"]][["Count"]]
     )

plot(after_df_list[["P2_RML6NHB_comparison"]][["1-P2-RML6_GT_[normalized count]"]],
     before_df_list[["1-P2-RML6_GT-result"]][["Count"]]
     )

plot(after_df_list[["P4_RML6NHB_comparison"]][["1-P4-RML6_GT_[normalized count]"]],
     before_df_list[["1-P4-RML6_GT-result"]][["Count"]]
     )

plot(after_df_list[["P8_RML6NHB_comparison_without1NBH"]][["1-P8-RML6_GT_[normalized count]"]],
     before_df_list[["1-P8-RML6_GT-result"]][["Count"]]
     )





# Prepare bar plots -------------------------------------------------------

sum_of_counts_list <- list(
  "P0" = list("before" = c("NBH_rep1"  = NA,
                           "RML6_rep1" = NA,
                           "NBH_rep2"  = sum(before_df_list[["2-P0-NBH_GT-result"]][["Count"]]),
                           "RML6_rep2" = sum(before_df_list[["2-P0-RML6_GT-result"]][["Count"]]),
                           "NBH_rep3"  = sum(before_df_list[["3-P0-NBH_GT-result"]][["Count"]]),
                           "RML6_rep3" = sum(before_df_list[["3-P0-RML6_GT-result"]][["Count"]])
                           ),
              "after" =  c("NBH_rep1"  = NA,
                           "RML6_rep1" = NA,
                           "NBH_rep2"  = sum(after_df_list[["P0_RML6NHB_comparison"]][["2-P0-NBH_GT_[normalized count]"]]),
                           "RML6_rep2" = sum(after_df_list[["P0_RML6NHB_comparison"]][["2-P0-RML6_GT_[normalized count]"]]),
                           "NBH_rep3"  = sum(after_df_list[["P0_RML6NHB_comparison"]][["3-P0-NBH_GT_[normalized count]"]]),
                           "RML6_rep3" = sum(after_df_list[["P0_RML6NHB_comparison"]][["3-P0-RML6_GT_[normalized count]"]])
                           )
              ),
  "P2" = list("before" = c("NBH_rep1"  = NA,
                           "RML6_rep1" = NA,
                           "NBH_rep2"  = sum(before_df_list[["2-P2-NBH_GT-result"]][["Count"]]),
                           "RML6_rep2" = sum(before_df_list[["2-P2-RML6_GT-result"]][["Count"]]),
                           "NBH_rep3"  = sum(before_df_list[["3-P2-NBH_GT-result"]][["Count"]]),
                           "RML6_rep3" = sum(before_df_list[["3-P2-RML6_GT-result"]][["Count"]])
                           ),
              "after" =  c("NBH_rep1"  = NA,
                           "RML6_rep1" = NA,
                           "NBH_rep2"  = sum(after_df_list[["P2_RML6NHB_comparison_without1NBH1RML6"]][["2-P2-NBH_GT_[normalized count]"]]),
                           "RML6_rep2" = sum(after_df_list[["P2_RML6NHB_comparison_without1NBH1RML6"]][["2-P2-RML6_GT_[normalized count]"]]),
                           "NBH_rep3"  = sum(after_df_list[["P2_RML6NHB_comparison_without1NBH1RML6"]][["3-P2-NBH_GT_[normalized count]"]]),
                           "RML6_rep3" = sum(after_df_list[["P2_RML6NHB_comparison_without1NBH1RML6"]][["3-P2-RML6_GT_[normalized count]"]])
                           )
              ),
  "P4" = list("before" = c("NBH_rep1"  = NA,
                           "RML6_rep1" = sum(before_df_list[["1-P4-RML6_GT-result"]][["Count"]]),
                           "NBH_rep2"  = sum(before_df_list[["2-P4-NBH_GT-result"]][["Count"]]),
                           "RML6_rep2" = sum(before_df_list[["2-P4-RML6_GT-result"]][["Count"]]),
                           "NBH_rep3"  = sum(before_df_list[["3-P4-NBH_GT-result"]][["Count"]]),
                           "RML6_rep3" = sum(before_df_list[["3-P4-RML6_GT-result"]][["Count"]])
                           ),
              "after" =  c("NBH_rep1"  = NA,
                           "RML6_rep1" = sum(after_df_list[["P4_RML6NHB_comparison_without1NBH"]][["1-P4-RML6_GT_[normalized count]"]]),
                           "NBH_rep2"  = sum(after_df_list[["P4_RML6NHB_comparison_without1NBH"]][["2-P4-NBH_GT_[normalized count]"]]),
                           "RML6_rep2" = sum(after_df_list[["P4_RML6NHB_comparison_without1NBH"]][["2-P4-RML6_GT_[normalized count]"]]),
                           "NBH_rep3"  = sum(after_df_list[["P4_RML6NHB_comparison_without1NBH"]][["3-P4-NBH_GT_[normalized count]"]]),
                           "RML6_rep3" = sum(after_df_list[["P4_RML6NHB_comparison_without1NBH"]][["3-P4-RML6_GT_[normalized count]"]])
                           )
              ),
  "P8" = list("before" = c("NBH_rep1"  = NA,
                           "RML6_rep1" = sum(before_df_list[["1-P8-RML6_GT-result"]][["Count"]]),
                           "NBH_rep2"  = sum(before_df_list[["2-P8-NBH_GT-result"]][["Count"]]),
                           "RML6_rep2" = sum(before_df_list[["2-P8-RML6_GT-result"]][["Count"]]),
                           "NBH_rep3"  = sum(before_df_list[["3-P8-NBH_GT-result"]][["Count"]]),
                           "RML6_rep3" = sum(before_df_list[["3-P8-RML6_GT-result"]][["Count"]])
                           ),
              "after" =  c("NBH_rep1"  = NA,
                           "RML6_rep1" = sum(after_df_list[["P8_RML6NHB_comparison_without1NBH"]][["1-P8-RML6_GT_[normalized count]"]]),
                           "NBH_rep2"  = sum(after_df_list[["P8_RML6NHB_comparison_without1NBH"]][["2-P8-NBH_GT_[normalized count]"]]),
                           "RML6_rep2" = sum(after_df_list[["P8_RML6NHB_comparison_without1NBH"]][["2-P8-RML6_GT_[normalized count]"]]),
                           "NBH_rep3"  = sum(after_df_list[["P8_RML6NHB_comparison_without1NBH"]][["3-P8-NBH_GT_[normalized count]"]]),
                           "RML6_rep3" = sum(after_df_list[["P8_RML6NHB_comparison_without1NBH"]][["3-P8-RML6_GT_[normalized count]"]])
                           )
              )
)




CountBarPlots(sum_of_counts_list)

pdf(file = file.path(file_output_directory, "Sum of sgRNA counts.pdf"),
    height = 8.1, width = 7.5
    )
CountBarPlots(sum_of_counts_list)
dev.off()







