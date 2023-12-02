### 2023-10-08


# Define functions --------------------------------------------------------

AnnotatePrepoolPostpoolPlot <- function(top_title = NULL, top_labels = FALSE) {

  bar_positions <- PrepoolPostpoolBarPositions()

  if (top_labels) {
    if (!(is.null(top_title))) {
      mtext(top_title, line = 1.6, padj = 0, cex = par("cex"))
    }
    segments(x0  = bar_positions[c(1, 5)] - 0.25,
             x1  = bar_positions[c(4, 8)] + 0.25,
             y0  = par("usr")[[4]] + diff(grconvertY(c(0, 0.35), from = "lines", to = "user")),
             col = "gray50",
             xpd = NA
             )
    mtext(c("prepool", "postpool"),
          at = c(mean(bar_positions[1:4]), mean(bar_positions[5:8])),
          line = 0.525, padj = 0, cex = par("cex")
          )
  }

  ## Draw the bottom labels
  mtext(text = sapply(as.character(rep(1:2, times = 4)), VerticalAdjust),
        at = bar_positions, side = 1, line = 0.925, cex = par("cex")
        )
  are_rep1 <- rep(c(TRUE, FALSE), times = 4)
  segments(x0  = bar_positions[are_rep1] - 0.1,
           x1  = bar_positions[!(are_rep1)] + 0.1,
           y0  = par("usr")[[3]] - diff(grconvertY(c(0, 1.725), from = "lines", to = "user")),
           col = "gray50",
           xpd = NA
           )
  mtext(text = sapply(c("T0", "end"), VerticalAdjust),
        at   = tapply(bar_positions, rep(1:4, each = 2), mean),
        side = 1,
        line = 2.125,
        cex  = par("cex")
        )

  if (!(top_labels)) {
    segments(x0  = bar_positions[c(1, 5)] - 0.1,
             x1  = bar_positions[c(4, 8)] + 0.1,
             y0  = par("usr")[[3]] - diff(grconvertY(c(0, 2.925), from = "lines", to = "user")),
             col = "gray50",
             xpd = NA
             )
    mtext(text = sapply(c("prepool", "postpool"), VerticalAdjust),
          at   = tapply(bar_positions, rep(1:2, each = 4), mean),
          side = 1,
          line = 3.325,
          cex  = par("cex")
          )

    mtext(VerticalAdjust("lentivirus:"), side = 1, line = 3.325, cex = par("cex"),
        at = par("usr")[[1]] + diff(grconvertX(c(0, 1), from = "lines", to = "user")),
        adj = 1
        )
  }

  mtext(VerticalAdjust("replicate:"), side = 1, line = 0.925, cex = par("cex"),
        at = par("usr")[[1]] + diff(grconvertX(c(0, 1), from = "lines", to = "user")),
        adj = 1
        )
  mtext(VerticalAdjust("timepoint:"), side = 1, line = 2.125, cex = par("cex"),
        at = par("usr")[[1]] + diff(grconvertX(c(0, 1), from = "lines", to = "user")),
        adj = 1
        )

  return(invisible(NULL))
}



PrepoolPostpoolBarPositions <- function() {
  bar_positions <- RepositionByGroups(c(1, 1, 2, 2), gap_ratio = 1.25)
  bar_positions <- c(scales::rescale(bar_positions, to = c(1, 3.7)),
                     scales::rescale(bar_positions, to = c(5.3, 8))
                     )
  return(bar_positions)
}


