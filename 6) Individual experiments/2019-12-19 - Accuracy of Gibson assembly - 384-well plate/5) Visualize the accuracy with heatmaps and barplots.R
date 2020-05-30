### 28 May 2020 ###




# Import packages and source code -----------------------------------------

library("RColorBrewer")
library("squash")
library("scales")
library("viridis")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")
file_output_directory            <- file.path(file_directory, "4) Output")




# Load data ---------------------------------------------------------------

load(file.path(intermediate_R_objects_directory, "3) Import data from external tools.RData"))





# Re-order assembly_df ----------------------------------------------------

new_order <- order(assembly_df[["Fraction_correct_4sg"]])







# Define functions --------------------------------------------------------

MakeEmptyPlot <- function(y_limits = c(0, 1)) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = y_limits, xaxs = "i", yaxs = "i"
       )
}


parula_colors <- c("#352A87", "#0F5CDD", "#1481D6", "#06A4CA", "#2EB7A4",
                   "#87BF77", "#D1BB59", "#FEC832", "#F9FB0E"
                   )


magma_colors <- c("#000004FF", "#02020DFF", "#06051BFF", "#0D0A28FF", "#140E36FF",
                  "#1C1044FF", "#241254FF", "#2E1162FF", "#390F6DFF", "#440F76FF",
                  "#4E117BFF", "#58157EFF", "#611980FF", "#6B1D81FF", "#752181FF",
                  "#7F2482FF", "#882781FF", "#932B80FF", "#9C2E7FFF", "#A6317DFF",
                  "#B1357BFF", "#BB3978FF", "#C53C74FF", "#CF4070FF", "#D8456CFF",
                  "#E14D66FF", "#E95462FF", "#EF5D5EFF", "#F4675CFF", "#F7735CFF",
                  "#FA7E5EFF", "#FC8A62FF", "#FD9668FF", "#FEA16EFF", "#FEAD77FF",
                  "#FEB87FFF", "#FEC488FF", "#FECF92FF", "#FDDB9DFF", "#FDE6A8FF",
                  "#FCF2B4FF", "#FCFDBFFF"
                  )

ParulaFunction <- colorRampPalette(rev(parula_colors))

MagmaFunction <- colorRampPalette(rev(magma_colors))

CividisFunction <- function(x) viridis::cividis(x, direction = -1)



original_cividis_100 <- CividisFunction(100)

percentage_orange_vec <- c(seq(50, 0, by = -(50 / 30)),
                           rep(0, 69)
                           )


MixInColor <- function(color_1, color_2, fraction_color_2 = 0.5) {
  color_index <- max(1, round(1000 * fraction_color_2))
  colorRampPalette(c(color_1, color_2))(1000)[[color_index]]
}


use_orange <- toupper("#e6b135")

orange_cividis_100 <- vapply(1:100,
                             function(x) MixInColor(original_cividis_100[[x]],
                                                    use_orange,
                                                    percentage_orange_vec[[x]] / 100
                                                    ),
                             ""
                             )



OrangeCividisFunction <- colorRampPalette(orange_cividis_100)






DrawColorTrapezoid <- function(start_x     = 0,
                               start_y     = 0,
                               end_x       = 1,
                               end_y       = 1,
                               trapezoid_y = NULL,
                               my_colors   = OrangeCividisFunction(100)
                               ) {

  my_x_positions <- seq(start_x, end_x, by = (end_x - start_x) / length(my_colors))
  num_positions <- length(my_x_positions)

  bottom_left_corners_x  <- my_x_positions[2:(num_positions - 1L)]
  top_left_corners_x     <- bottom_left_corners_x
  top_right_corners_x    <- my_x_positions[3:num_positions]
  bottom_right_corners_x <- top_right_corners_x

  draw_triangle <- !(!(is.null(trapezoid_y)) && (trapezoid_y == end_y))

  if (draw_triangle) {
    if (is.null(trapezoid_y)) {
      triangle_start_y <- start_y
    } else {
      triangle_start_y <- trapezoid_y
    }

    my_y_positions <- seq(triangle_start_y,
                          end_y,
                          by = (end_y - triangle_start_y) / length(my_colors)
                          )


    bottom_left_corners_y  <- rep_len(start_y, num_positions - 1L)
    top_left_corners_y     <- my_y_positions[2:(num_positions - 1L)]
    top_right_corners_y    <- my_y_positions[3:num_positions]
    bottom_right_corners_y <- rep_len(start_y, num_positions - 1L)
  }




  if (!(is.null(trapezoid_y))) {
    rect(xleft   = my_x_positions[[1]],
         xright  = my_x_positions[[2]],
         ybottom = start_y,
         ytop    = trapezoid_y,
         col     = my_colors[[1]],
         border  = NA,
         xpd     = NA
         )

    for (i in seq_len(length(my_colors) - 1L)) {
      my_color <- my_colors[[i + 1L]]
      rect(xleft   = bottom_left_corners_x[[i]],
           xright  = top_right_corners_x[[i]],
           ybottom = start_y,
           ytop    = trapezoid_y,
           col     = my_color,
           border  = NA,
           xpd     = NA
           )
    }


  }

  if (!(!(is.null(trapezoid_y)) && (trapezoid_y == end_y))) {
    polygon(x      = c(my_x_positions[[1]], my_x_positions[[2]], my_x_positions[[2]]),
            y      = c(my_y_positions[[1]], my_y_positions[[2]], my_y_positions[[1]]),
            col    = my_colors[[1]],
            border = NA,
            xpd    = NA
            )

    for (i in seq_len(length(my_colors) - 1L)) {
      my_color <- my_colors[[i + 1L]]
      polygon(x      = c(bottom_left_corners_x[[i]], top_left_corners_x[[i]], top_right_corners_x[[i]], bottom_right_corners_x[[i]]),
              y      = c(bottom_left_corners_y[[i]], top_left_corners_y[[i]], top_right_corners_y[[i]], bottom_right_corners_y[[i]]),
              col    = my_color,
              border = NA,
              xpd    = NA
              )
    }
  }
}





# Necessary for vertical bars and input
Transposify <- function(numeric_input) {
  t(apply(as.matrix(numeric_input), 2, rev))
}


Do_cimage <- function(color_input) {
# If a vector is provided, by default it will be turned into a HORIZONTAL matrix
  if (!(is.matrix(color_input))) {
    input_matrix <- as.matrix(color_input)
  } else if (nrow(color_input) == 1) {
    input_matrix <- t(color_input)
  } else {
    input_matrix <- Transposify(color_input)
  }
  cimage(zcol = input_matrix, axes = FALSE, xlab = "", ylab = "")
}











# Set up the layout -------------------------------------------------------

png(filename = file.path(file_output_directory, "Accuracy heatmap.png"),
    res    = 600,
    height = 4.5,
    width  = 6.5,
    units  = "in"
    )

space_height <- 1

golden_ratio <- (1 + sqrt(5)) / 2

layout(cbind(rep(1, 9),
             3:(9 + 3 - 1),
             rep(2, 9)
             ),
       widths  = c(0.1, 0.8, 0.1),
       heights = c(space_height * 1.2,
                   3.5,
                   space_height * 0.13,
                   space_height * 0.5,
                   space_height * 1.2,
                   3 * golden_ratio,
                   space_height * 1.5,
                   space_height * 0.5,
                   space_height * 0.6
                   )
       )

par(mar = rep(0, 4))





# Draw a vertical barplot -------------------------------------------------

for (i in 1:4) {
  MakeEmptyPlot()
}

# text(x = 0.5,
#      y = 0.5,
#      labels = "Sequencing of Gibson assembly \u2013 384 genes, re-ordered by accuracy",
#      font = 2,
#      xpd = NA
#      )


fraction_correct_vec <- assembly_df[["Fraction_correct_4sg"]][new_order]
are_80perc_or_more <- fraction_correct_vec > 0.8
num_above_80 <- sum(are_80perc_or_more)
over_80_fraction <- num_above_80 / 384


darker_colors <- toupper("#15396d")
lighter_color <- "#ECEBE9"



rect(xleft = 0,
     xright = 1,
     ybottom = 0,
     ytop = 1,
     col = lighter_color,
     border = NA
     )


# grid_lines_grey <- "gray75"
#
# dark_line_positions <- seq(0.2, 1, by = 0.2)
# light_line_positions <- seq(0.1, 0.9, by = 0.2)
# segments(x0 = 0, x1 = 1, y0 = dark_line_positions, y1 = dark_line_positions,
#          col = grid_lines_grey, lend = "butt"
#          )
# segments(x0 = 0, x1 = 1, y0 = light_line_positions, y1 = light_line_positions,
#          col = grid_lines_grey, lend = "butt"
#          )


for (i in 1:384) {
  rect(xleft   = (i - 1) / 384,
       xright  = (i / 384),
       ybottom = 0,
       ytop    = fraction_correct_vec[[i]],
       col     = darker_colors[[1]],
       border  = NA
       )
}

# segments(x0 = 0, x1 = 1, y0 = 0, y1 = 0,
#          col = grid_lines_grey, lend = "butt"
#          )



# abline(h = 0.8, col = "gray60")
# abline(v = over_80_fraction, col = "white")

tick_locations <- axTicks(2)
tick_labels <- paste0(tick_locations * 100, "%")
axis(2,
     labels   = tick_labels,
     at       = tick_locations,
     las      = 1,
     mgp      = c(3, 0.45, 0),
     tcl      = -0.3
     )

text(x      = -0.062,
     y      = 1.13,
     adj    = c(0, 0),
     labels = "Mean accuracy",
     xpd    = NA,
     font   = 2
     )

# box()




# Draw the % accuracy horizontal barplot ----------------------------------

for (i in 1:2) {
  MakeEmptyPlot()
}




two_grey_colors <- c("gray77", "gray40")

rect(xleft   = c(0, 1 - over_80_fraction),
     xright  = c(1 - over_80_fraction, 1),
     ybottom = 0,
     ytop    = 1,
     border  = NA, #"gray20",
     lwd     = 0.5,
     col     = two_grey_colors
     )

x_mid <- (num_above_80 / 2) / 384

mid_x <- par("usr")[[1]] + ((par("usr")[[2]] - (par("usr")[[1]])) * x_mid)
mid_y <- 0.5

text(x      = 0.5,
     y      = -0.65,
     adj    = c(0.5, 0.5),
     labels = bquote(bold(.(as.character(num_above_80)) *  " / 384 genes " *
                          "were " >= "80% accurate"
                          )),
     font   = 2,
     col    = "black",
     xpd    = NA
     )


MakeEmptyPlot()





# Draw the heatmap --------------------------------------------------------

numeric_mat <- t(as.matrix(assembly_df[, paste0("Fraction_correct_sg", 1:4)]))

my_breaks <- seq(0, 1, by = 0.01)
my_breaks[c(1, length(my_breaks))] <- c(0, ceiling(my_breaks[length(my_breaks)]))

my_cmap <- makecmap(numeric_mat, colFn = OrangeCividisFunction, breaks = my_breaks)

my_color_mat <- cmap(numeric_mat, my_cmap)

Do_cimage(my_color_mat)

x_range <- par("usr")[[2]] - par("usr")[[1]]

text(x      = par("usr")[[1]] - (x_range * 0.018),
     y      = 4:1,
     labels = paste0("sg", 1:4),
     xpd    = NA,
     adj    = c(1, 0.5)
     )
# box()

MakeEmptyPlot()





# Draw the color indicator (trapezoid) ------------------------------------

trapezoid_start_x <- 0.85
trapezoid_start_y <- 0.65
trapezoid_end_y   <- 0.9

DrawColorTrapezoid(start_x     = trapezoid_start_x,
                   end_x       = 1,
                   start_y     = trapezoid_start_y,
                   end_y       = trapezoid_end_y,
                   trapezoid_y = 0.72
                   )

text(x      = trapezoid_start_x - 0.009,
     y      = trapezoid_start_y + ((trapezoid_end_y - trapezoid_start_y) * 0.3),
     adj    = c(1, 0.5),
     labels = "Accuracy",
     xpd    = NA,
     font   = 2,
     cex    = 0.9
     )


trapezoid_seq <- seq(trapezoid_start_x, 1, by = ((1 - trapezoid_start_x) / 5))

text(x      = trapezoid_seq,
     y      = trapezoid_start_y - 0.16,
     labels = paste0(seq(0, 100, by = 20), "%"),
     font   = 2,
     cex    = 0.6,
     xpd    = NA
     )
segments(x0   = trapezoid_seq,
         x1   = trapezoid_seq,
         y0   = trapezoid_start_y - 0.088,
         y1   = trapezoid_start_y - 0.02,
         xpd  = NA,
         lwd  = 0.5,
         lend = "butt"
         )

segments(x0   = trapezoid_start_x,
         x1   = 1,
         y0   = trapezoid_start_y - 0.02,
         y1   = trapezoid_start_y - 0.02,
         lwd  = 0.5,
         lend = "butt"
         )

MakeEmptyPlot()




# Draw the strip indicating genes with homologies >8 bp -------------------

longest_subsequence_vec <- assembly_df[["Longest_subsequence"]][new_order]
have_homology <- longest_subsequence_vec >= 8

# homology_colors <- brewer.pal(8, "Paired")[c(3, 4)]
homology_colors <- two_grey_colors

for (i in 1:384) {
  rect(xleft   = (i - 1) / 384,
       xright  = (i / 384),
       ybottom = 0,
       ytop    = 1,
       col     = homology_colors[[as.integer(have_homology[[i]]) + 1]],
       border  = NA
       )
}


text(x      = 0.5,
     y      = -0.65,
     adj    = c(0.5, 0.5),
     labels = bquote(bold(.(as.character(sum(have_homology))) *  " / 384 genes " *
                          "had " >= "8bp homology"
                          )),
     font   = 2,
     col    = "black",
     xpd    = NA
     )

MakeEmptyPlot()
dev.off()






# Compute some statistics -------------------------------------------------

table(are_80perc_or_more, have_homology)
fisher.test(are_80perc_or_more, have_homology)






