### 21st October 2020 ###




# Define plot dimensions --------------------------------------------------

use_height <- 6.5
use_width <- 7




# Define functions --------------------------------------------------------

PlotBarplotMat <- function(barplot_mat,
                           colors_vec,
                           positions_vec = seq_len(ncol(barplot_mat))
                           ) {

  num_categories <- nrow(barplot_mat)
  stopifnot(length(colors_vec) == num_categories)

  lower_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    if (x == 1) {
      rep(0, ncol(barplot_mat))
    } else {
      colSums(barplot_mat[1:(x - 1), , drop = FALSE])
    }
  })

  upper_borders_vec_list <- lapply(seq_len(num_categories), function(x) {
    colSums(barplot_mat[seq_len(x), , drop = FALSE])
  })

  for (i in seq_len(ncol(barplot_mat))) {
    for (j in seq_len(num_categories)) {
      rect(xleft   = (positions_vec[[i]] - 1) / max(positions_vec),
           xright  = (positions_vec[[i]] / max(positions_vec)),
           ybottom = lower_borders_vec_list[[j]][[i]],
           ytop    = upper_borders_vec_list[[j]][[i]],
           col     = colors_vec[[j]],
           border  = NA,
           xpd     = NA
           )
    }
  }
  return(invisible(NULL))
}



SideTextAndAxes <- function(side_text) {

  tick_locations <- axTicks(2)
  tick_labels <- paste0(tick_locations * 100, "%")
  axis(2,
       labels   = tick_labels,
       at       = tick_locations,
       las      = 1,
       mgp      = c(3, 0.45, 0),
       tcl      = -0.3,
       lwd      = 0.75,
       cex.axis = 0.8
       )

  text(x      = -0.06,
       y      = 0.5,
       adj    = c(1, 0.5),
       labels = side_text,
       xpd    = NA,
       cex    = 1.3
       )

  return(invisible(NULL))

}




DrawAlterationBarplot <- function(summary_df,
                                  main_title    = NULL,
                                  reorder_wells = FALSE,
                                  gap_weight    = 2L
                                  ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  percent_columns <- paste0("Perc_sg", 1:4, "_cr", 1:4)
  if (reorder_wells) {
    mean_accuracies <- rowMeans(as.matrix(summary_df[, percent_columns]))
    new_order <- order(summary_df[["Perc_all_4"]],
                       summary_df[["Perc_at_least_3"]],
                       summary_df[["Perc_at_least_2"]],
                       summary_df[["Perc_at_least_1"]],
                       mean_accuracies
                       )
  } else {
    new_order <- seq_len(nrow(summary_df))
  }

  if ("Empty_well" %in% names(sg_sequences_df)) {
    are_to_include <- !(sg_sequences_df[["Empty_well"]])
  } else {
    are_to_include <- rep(TRUE, nrow(sg_sequences_df))
  }

  use_indices <- new_order[are_to_include[new_order]]
  summary_df <- summary_df[use_indices, ]
  rownames(summary_df) <- NULL

  num_wells <- length(use_indices)


  ## Set up the layout

  barplot_height <- 2.8

  space_height <- 1

  layout(cbind(rep(1, 9),
               3:(9 + 3 - 1),
               rep(2, 9)
               ),
         widths  = c(0.13, 0.8, 0.07),
         heights = c(3,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height,
                     barplot_height,
                     space_height
                     )
         )

  par(mar = rep(0, 4))

  for (i in 1:3) {
    MakeEmptyPlot()
  }
  if (!(is.null(main_title))) {
    text(x = 0.5,
         y = 0.65,
         labels = main_title,
         cex = 1.1,
         xpd = NA
         )

  }

  four_colors <- c("#F9F4EC", "#EE442F", "#63ACBE", "#601A4A")
  four_alterations <- c("Correct", "Mutation", "Deletion", "Contamination")

  four_colors <- rev(four_colors)
  four_alterations <- rev(four_alterations)

  color_text_vec <- c('bold(color1("% plasmids with") * ',
                      'color1(" ") * color2("mutations,") * ',
                      'color1(" ") * color3("deletions,") * ',
                      'color1(" or ") * color4("contaminations"))'
                      )

  use_colors <- c("#000000",
                  rev(four_colors)[2:4]
                  )

  color_indices <- seq_along(use_colors)
  text_indices <- seq_along(color_text_vec)
  color_text <- paste0(color_text_vec[text_indices], collapse = "")


  add_gap <- (!(reorder_wells)) && ("Block" %in% names(sg_sequences_df))
  if (add_gap) {
    block_vec <- sg_sequences_df[["Block"]][are_to_include]
    wells_seq <- rep(NA, num_wells)
    current_index <- 0L
    current_block <- block_vec[[1]]
    for (i in seq_len(num_wells)) {
      if (current_block != block_vec[[i]]) {
        current_block <- block_vec[[i]]
        current_index <- current_index + gap_weight
      }
      current_index <- current_index + 1L
      wells_seq[[i]] <- current_index
    }
  } else {
    wells_seq <- seq_len(num_wells)
  }



  ## Draw a vertical barplot for the mutation categories

  column_name_list <- lapply(1:4,
                             function(x) {
                               paste0(four_alterations, "_sg", x, "_cr", x)
                             })

  column_mat_list <- lapply(column_name_list, function(x) {
    use_mat <- as.matrix(summary_df[, x])
    for (i in seq_len(ncol(use_mat))) {
      use_mat[, i] <- use_mat[, i] / summary_df[["Count_total"]]
    }
    return(t(use_mat))
  })


  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[1]], four_colors, wells_seq)
  SideTextAndAxes("sg1")


  for (j in color_indices) {
    text(x      = 0.5,
         y      = 1.25,
         labels = VerticalAdjust(parse(text = MakeInvisible(color_text, j))),
         adj    = c(0.5, 0),
         xpd    = NA,
         col    = use_colors[[j]],
         cex    = 1.2
         )
  }

  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[2]], four_colors, wells_seq)
  SideTextAndAxes("sg2")
  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[3]], four_colors, wells_seq)
  SideTextAndAxes("sg3")
  MakeEmptyPlot()

  MakeEmptyPlot()
  PlotBarplotMat(column_mat_list[[4]], four_colors, wells_seq)
  SideTextAndAxes("sg4")
  MakeEmptyPlot()


  return(invisible(NULL))
}
