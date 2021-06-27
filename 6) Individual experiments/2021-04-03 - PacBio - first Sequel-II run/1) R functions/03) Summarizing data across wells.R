### 13th June 2021 ###






# Helper functions for creating plots -------------------------------------

GetPlateSelection <- function(plate_names) {
  require_objects <- c("plate_selection_list", "plate_selection_titles_list")
  stopifnot(all(require_objects %in% ls(envir = globalenv())))
  if (is.null(plate_names)) {
    plate_names <- "All plates"
  }
  if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
    plate_names <- plate_selection_list[[plate_names]]
  } else if ((length(plate_names) == 1) && (plate_names %in% names(plate_selection_titles_list))) {
    main_title <- plate_selection_titles_list[[plate_names]]
  } else {
    main_title <- "PacBio sequencing"
  }
  results_list <- list("plate_names" = plate_names,
                       "title" = main_title
                       )
  return(results_list)
}



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



