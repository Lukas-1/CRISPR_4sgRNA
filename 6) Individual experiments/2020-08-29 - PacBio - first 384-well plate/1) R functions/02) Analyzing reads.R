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



GetWellNumbers <- function(lima_report_df) {
  combo_IDs_vec <- paste0(lima_report_df[["IdxLowestNamed"]], "--",
                          lima_report_df[["IdxHighestNamed"]]
                          )
  barcodes_to_wells_map[combo_IDs_vec]
}



# Functions for reading in PacBio output ----------------------------------

ExtractZMWs <- function(report_df, lima_list, wells_vec = seq_len(384)) {
  ccs_well_numbers <- GetWellNumbers(report_df)
  are_valid_wells <- ccs_well_numbers %in% wells_vec

  report_zmws <- as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]])))
  lima_zmws <- as.integer(substr(lima_list[["qname"]], 22, nchar(lima_list[["qname"]]) - 4))

  result_zmws <- report_zmws[are_valid_wells & (report_zmws %in% lima_zmws)]
  return(result_zmws)
}


ReadLimaReport <- function(file_path) {
  read.table(file_path,
             header = TRUE, sep = "\t",
             quote = "", stringsAsFactors = FALSE
             )
}


