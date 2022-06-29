### 20th November 2021 ###



# Define functions --------------------------------------------------------

SetUpMappings <- function() {

  plate_names <- c("A1", "A2", "B1", "B2")

  scheme_2rows <- rbind(rep(plate_names[1:2], times = 12), rep(plate_names[3:4], times = 12))
  scheme_mat <- do.call(rbind, lapply(1:8, function(x) scheme_2rows))

  indices_384_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)

  indices_mat_list <- sapply(plate_names, function(x) {
    are_this_scheme <- scheme_mat == x
    rows_list <- lapply(1:16, function(y) {
      indices_384_mat[y, ][are_this_scheme[y, ]]
    })
    do.call(rbind, rows_list)
  }, simplify = FALSE)


  coords_96wp_mat <- do.call(rbind, lapply(LETTERS[1:8], function(x) {
    paste0(x, 1:12)
  }))

  long_df_list <- lapply(1:4, function(x) {
    data.frame("Well_number_384"  = as.vector(t(indices_mat_list[[x]])),
               "Plate_96"         = names(indices_mat_list)[[x]],
               "Well_coords_96"   = as.vector(t(coords_96wp_mat)),
               stringsAsFactors = FALSE
               )
    })
  long_df <- do.call(rbind.data.frame,
                     c(long_df_list,
                       stringsAsFactors = FALSE,
                       make.row.names = FALSE
                     ))

  long_df[, "String_96wp"] <- paste0("Plate ", long_df[, "Plate_96"],
                                     " - well ", long_df[, "Well_coords_96"]
                                     )
  return(long_df)
}




ConvertWellNumbers <- function(well_numbers) {

  are_valid <- well_numbers %in% seq_len(384)

  if (!(all(are_valid))) {
    stop("The input must consist of whole numbers from 1 to 384!")
  }
  conversion_df <- SetUpMappings()

  matches_vec <- match(well_numbers, conversion_df[, "Well_number_384"])
  results_vec <- conversion_df[matches_vec, "String_96wp"]
  return(results_vec)
}



