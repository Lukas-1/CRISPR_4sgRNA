### 15th October 2020




# Define functions --------------------------------------------------------

GetWellNumbers <- function(lima_report_df) {
  combo_IDs_vec <- paste0(lima_report_df[["IdxLowestNamed"]], "--",
                          lima_report_df[["IdxHighestNamed"]]
                          )
  barcodes_to_wells_map[combo_IDs_vec]
}


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

