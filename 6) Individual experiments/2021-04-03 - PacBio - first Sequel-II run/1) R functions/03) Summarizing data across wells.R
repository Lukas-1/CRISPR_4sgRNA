### 13th June 2021 ###





# Define functions --------------------------------------------------------

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
