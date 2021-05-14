### 15th October 2020




# Define functions --------------------------------------------------------

GetMeanQuality <- function(qualities, rescale = TRUE) {
  if (!(identical("PhredQuality", as.character(class(qualities))))) {
    qualities <- PhredQuality(qualities)
  }
  mean_qualities <- ShortRead::alphabetScore(qualities) / width(qualities)
  if (rescale) {
    mean_qualities <- mean_qualities / 93 * 100
  }
  return(mean_qualities)
}



PassBasicFilters <- function(input_ccs_df, wells_vec) {
  are_valid <- input_ccs_df[, "Passed_filters"] &
               (input_ccs_df[, "Well_number"] %in% wells_vec)
  if ("Well_exists" %in% names(input_ccs_df)) {
    are_valid <- are_valid & (input_ccs_df[["Well_exists"]] %in% TRUE)
  }
  return(are_valid)
}



GetCCS3_ZMWs <- function(input_ccs_df, wells_vec = seq_len(384)) {
  are_CCS3 <- PassBasicFilters(input_ccs_df, wells_vec) &
              (input_ccs_df[, "Read_quality"] >= 0.99) &
              (input_ccs_df[, "Num_full_passes"] >= 3)
  return(input_ccs_df[["ZMW"]][are_CCS3])
}


GetCCS5_ZMWs <- function(input_ccs_df, wells_vec = seq_len(384)) {
  are_eligible <- PassBasicFilters(input_ccs_df, wells_vec)
  if ("Pass_CCS5" %in% names(input_ccs_df)) {
    are_CCS5 <- are_eligible & input_ccs_df[["Pass_CCS5"]]
  } else {
    are_CCS5 <- are_eligible &
                (input_ccs_df[, "Read_quality"] >= 0.999) &
                (input_ccs_df[, "Num_full_passes"] >= 5)
  }
  return(input_ccs_df[["ZMW"]][are_CCS5])
}


GetCCS7_ZMWs <- function(input_ccs_df, wells_vec = seq_len(384)) {
  are_CCS7 <- PassBasicFilters(input_ccs_df, wells_vec) &
              (input_ccs_df[, "Read_quality"] >= 0.9999) &
              (input_ccs_df[, "Num_full_passes"] >= 7)
  return(input_ccs_df[["ZMW"]][are_CCS7])
}



