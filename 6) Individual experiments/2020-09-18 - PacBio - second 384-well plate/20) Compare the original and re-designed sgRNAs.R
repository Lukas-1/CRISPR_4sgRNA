### 28th October 2020 ###





# Import packages and source code -----------------------------------------

library("beeswarm")
library("RColorBrewer")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")
file_output_directory  <- file.path(plate2_directory, "3) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
p1_sg_sequences_df <- sg_sequences_df
load(file.path(p1_R_objects_directory, "11) Process demultiplexed PacBio reads.RData"))
p1_sl7_ccs5_df_list <- sl7_ccs5_df_list

load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
p2_sg_sequences_df <- sg_sequences_df
load(file.path(p2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))
p2_sl7_ccs5_df_list <- sl7_ccs5_df_list

rm(sg_sequences_df)
rm(sl7_ccs5_df_list)






# Define functions --------------------------------------------------------

PlotSelectionBias <- function(p1_df_list, p2_df_list, summary_df_name, show_column) {

  p1_vec <- p1_df_list[[summary_df_name]][[show_column]]
  p2_vec <- p2_df_list[[summary_df_name]][[show_column]]

  are_eligible <- (p2_sg_sequences_df[["Longest_subsequence"]] <= 7) &
                  (p2_sg_sequences_df[["Block"]] == 2)

  well_matches <- match(p2_sg_sequences_df[["Target_gene"]][are_eligible],
                        p1_sg_sequences_df[["Target_gene"]]
                        )

  were_picked <- p1_sg_sequences_df[["Target_gene"]] %in% p2_sg_sequences_df[["Target_gene"]]


  swarm_list <- lapply(1:20, function(x) {
    have_this_homology <- p1_sg_sequences_df[["Longest_subsequence"]] == x
    list(p1_vec[have_this_homology & !(were_picked)],
         p1_vec[have_this_homology & were_picked]
         )
  })


  plot(1, type = "n", axes = FALSE, ann = FALSE,
       xlim = c(3, 20), ylim = c(0, 100), xaxs = "i", yaxs = "i"
       )

  tick_locations <- axTicks(2)
  tick_labels <- paste0(tick_locations, "%")

  abline(h = seq(10, 90, 20), col = "gray98", lwd = 0.75)
  abline(h = tick_locations, col = "gray94", lwd = 0.75)

  box(bty = "l", lwd = 0.75)

  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.35,
       lwd      = 0.75,
       cex.axis = 0.9
       )

  for (i in 1:20) {
    for (j in 1:2) {
      sub_vec <- swarm_list[[i]][[j]]
      if (length(sub_vec > 1)) {
        x_position <- i
        x_dodge <- 1 / 6
        if (j == 1) {
          x_position <- x_position - x_dodge
          color_scheme <- "Blues"
        } else {
          x_position <- x_position + x_dodge
          color_scheme <- "Reds"
        }
        boxplot(x         = sub_vec,
                at        = x_position,
                boxwex    = 0.5,
                outline   = FALSE,
                whisklty  = "blank",
                staplewex = 0,
                axes      = FALSE,
                whisklwd  = 0,
                staplelty = 0,
                col       = brewer.pal(9, color_scheme)[[2]],
                boxlwd    = 0.75,
                medlwd    = par("lwd") * 2,
                add       = TRUE
        )

        point_cex <- 0.4
        beeswarm_df <- beeswarm(sub_vec,
                                at       = x_position,
                                priority = "random",
                                spacing  = 0.7,
                                cex      = point_cex,
                                do.plot  = FALSE
                                )
        points(beeswarm_df[["x"]],
               beeswarm_df[["y"]],
               pch = 16,
               cex = point_cex,
               col = brewer.pal(9, color_scheme)[[7]],
               xpd = NA
               )

      }
    }
  }

  hom_vec <- c(p1_sg_sequences_df[["Longest_subsequence"]],
               p2_sg_sequences_df[["Longest_subsequence"]]
               )

  present_lengths <- seq(from = min(hom_vec), to = max(hom_vec))

  text(x      = present_lengths,
       y      = par("usr")[[3]] - ((par("usr")[[4]] - par("usr")[[3]]) * 0.04),
       labels = present_lengths,
       font   = 1,
       cex    = 0.8,
       xpd    = NA
       )

  title("Accuracy vs. number of homologies \u2013 selected for the second plate?", cex.main = 0.9, font.main = 1)

  legend("bottom",
         fill = c(brewer.pal(9, "Blues")[[3]], brewer.pal(9, "Reds")[[3]]),
         legend = c("Not selected", "Selected"),
         horiz = TRUE,
         inset = -0.21,
         cex = 0.9,
         bty = "n",
         xpd = NA
         )
  return(invisible(NULL))
}




# Plot selection bias -----------------------------------------------------

PlotSelectionBias(p1_sl7_ccs5_df_list,
                  p2_sl7_ccs5_df_list,
                  "filtered_summary_df",
                  "Perc_all_4"
                  )


pdf_height <- 5.4
pdf_width <- 7.2
use_mar <- c(4.5, 4, 4, 2.5)



pdf(file = file.path(plots_output_directory,
                     "Bias in selection for re-sequencing.pdf"
                     ),
    width = pdf_width,
    height = pdf_height
    )
par("mar" = use_mar)
PlotSelectionBias(p1_sl7_ccs5_df_list,
                  p2_sl7_ccs5_df_list,
                  "filtered_summary_df",
                  "Perc_all_4"
                  )
dev.off()





