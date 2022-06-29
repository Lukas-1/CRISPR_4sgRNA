### 16th March 2022 ###



# Define functions --------------------------------------------------------

MakeDistanceList <- function(manhattan_distance = FALSE, matrix_format = FALSE, pad_sides = NULL) {

  num_rows <- 16
  num_columns <- 24
  if (!(is.null(pad_sides))) {
    num_rows <- num_rows + (pad_sides * 2)
    num_columns <- num_columns + (pad_sides * 2)
  }
  num_cells <- num_rows * num_columns

  plate_mat <- matrix(seq_len(num_cells), nrow = num_rows, ncol = num_columns, byrow = TRUE)

  current_well <- 1
  distance_list <- lapply(seq_len(num_cells), function(current_well) {
    distance_vec <- vapply(seq_len(num_cells), function(target_well) {

      are_current_well <- plate_mat == current_well
      are_target_well <- plate_mat == target_well

      current_row_index <- which(rowSums(are_current_well) == 1)
      target_row_index  <- which(rowSums(are_target_well) == 1)
      current_col_index <- which(colSums(are_current_well) == 1)
      target_col_index  <- which(colSums(are_target_well) == 1)

      col_distance <- abs(current_col_index - target_col_index)
      row_distance <- abs(current_row_index - target_row_index)

      if (manhattan_distance) {
        distance <- col_distance + row_distance
      } else {
        distance <- sqrt(col_distance^2 + row_distance^2) # Compute the Euclidian distance
      }
      return(distance)
    }, if (manhattan_distance) integer(1) else numeric(1))
    if (matrix_format) {
      return(matrix(distance_vec, nrow = num_rows, ncol = num_columns, byrow = TRUE))
    } else {
      return(distance_vec)
    }
  })
  return(distance_list)
}







ReplaceByNearest <- function(numeric_mat, distance_list, expand_by = 3L) {

  stopifnot((nrow(numeric_mat) * ncol(numeric_mat)) == length(distance_list))

  if (!(is.null(expand_by))) {
    left_right_mat <- matrix(nrow = nrow(numeric_mat) + (expand_by * 2),
                             ncol = expand_by
                             )
    top_bottom_mat <- matrix(nrow = expand_by, ncol = ncol(numeric_mat))
    left_right_mat[] <- NA
    top_bottom_mat[] <- NA
    numeric_mat <- rbind(top_bottom_mat, numeric_mat, top_bottom_mat)
    numeric_mat <- cbind(left_right_mat, numeric_mat, left_right_mat)
  }

  numeric_vec <- as.vector(t(numeric_mat))
  are_eligible <- !(is.na(numeric_vec))
  all_seq <- seq_along(numeric_vec)
  substituted_vec <- vapply(all_seq, function(x) {
    distance_vec <- distance_list[[x]]
    are_min_distance <- are_eligible & (distance_vec == min(distance_vec[are_eligible]))
    return(mean(numeric_vec[are_min_distance]))
  }, numeric(1))





}

