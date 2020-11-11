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


GetCCS5ZMWs <- function(ccs_df, wells_vec = seq_len(384)) {
  are_CCS5 <- ccs_df[["Passed_filters"]] &
              ccs_df[["Pass_CCS5"]] &
              (ccs_df[["Well_number"]] %in% wells_vec)
  return(ccs_df[["ZMW"]][are_CCS5])
}





