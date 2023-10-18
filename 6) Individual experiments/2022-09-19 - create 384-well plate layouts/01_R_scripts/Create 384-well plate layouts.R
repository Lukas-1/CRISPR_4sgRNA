## 2022-09-19



# Load packages and source code -------------------------------------------

library("writexl")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir           <- file.path(experiments_directory, "2022-09-19 - create 384-well plate layouts")
input_dir             <- file.path(project_dir, "02_input_data")
output_dir            <- file.path(project_dir, "03_output_data")



# Read in data ------------------------------------------------------------

CRISPRa_df  <- read.delim(file.path(input_dir, "CRISPRa_4sg_full_ordered_by_well.tsv"),
                          stringsAsFactors = FALSE
                          )
CRISPRko_df <- read.delim(file.path(input_dir, "CRISPRko_4sg_full_ordered_by_well.tsv"),
                          stringsAsFactors = FALSE
                          )



# Define functions --------------------------------------------------------

ConvertPlateNumbers <- function(plate_strings) {
  splits_list <- strsplit(plate_strings, "_", fixed = TRUE)
  plate_numbers_vec <- sapply(splits_list, "[[", 2)
  plate_numbers_vec <- sub("tf", "", plate_numbers_vec, fixed = TRUE)
  plate_numbers_vec <- sub("+", "_plus", plate_numbers_vec, fixed = TRUE)
  return(plate_numbers_vec)
}


SplitIntoPlates <- function(input_df) {
  plasmids_list <- split(input_df[, "Plasmid_name"],
                         factor(input_df[, "Plate_ID"],
                                levels = unique(input_df[, "Plate_ID"])
                                )
                         )
  mat_list <- lapply(plasmids_list, function(x) {
    if (length(x) != 384) {
      x <- c(x, rep("", 384 - length(x)))
    }
    results_mat <- matrix(x, nrow = 16, ncol = 24, byrow = TRUE)
    results_mat <- cbind(paste0("Row_", LETTERS[1:16]), results_mat)
    colnames(results_mat) <- c(" ", paste0("Col_", 1:24))
    return(results_mat)
  })
  df_list <- lapply(mat_list, function(x) {
    x <- data.frame(x, stringsAsFactors = FALSE)
    names(x)[[1]] <- ""
    return(x)
  })
  return(df_list)
}



# Prepare data ------------------------------------------------------------

## CRISPRa
CRISPRa_df[, "Plasmid_name"] <- paste0(CRISPRa_df[, "Gene_symbol"],
                                       ifelse(CRISPRa_df[, "TSS_ID"] == "", "", "_"),
                                       CRISPRa_df[, "TSS_ID"]
                                       )

CRISPRa_df[, "Plate_ID"] <- paste0("HA_", ConvertPlateNumbers(CRISPRa_df[, "Plate_string"]))
CRISPRa_df[, "Plasmid_ID"] <- paste0(CRISPRa_df[, "Plate_ID"], "_", CRISPRa_df[, "Plasmid_name"])
use_columns <- c("Plate_ID", "Plasmid_name", "Well_number")
stopifnot(all(table(CRISPRa_df[, "Plasmid_ID"]) == 4))
CRISPRa_df <- CRISPRa_df[!(duplicated(CRISPRa_df[, "Plasmid_ID"])), use_columns]
row.names(CRISPRa_df) <- NULL

## CRISPRko
CRISPRko_df[, "Plasmid_name"] <- CRISPRko_df[, "Gene_symbol"]
CRISPRko_df[, "Plate_ID"] <- paste0("HO_", ConvertPlateNumbers(CRISPRko_df[, "Plate_string"]))
CRISPRko_df[, "Plasmid_ID"] <- paste0(CRISPRko_df[, "Plate_ID"], "_", CRISPRko_df[, "Plasmid_name"])
stopifnot(all(table(CRISPRko_df[, "Plasmid_ID"]) == 4))
CRISPRko_df <- CRISPRko_df[!(duplicated(CRISPRko_df[, "Plasmid_ID"])), use_columns]
row.names(CRISPRko_df) <- NULL



# Split into plates -------------------------------------------------------

CRISPRa_df_list <- SplitIntoPlates(CRISPRa_df)
CRISPRko_df_list <- SplitIntoPlates(CRISPRko_df)



# Export Excel sheets -----------------------------------------------------

write_xlsx(CRISPRa_df_list,
           path = file.path(output_dir, "CRISPRa_sheets.xlsx"),
           format_headers = FALSE
           )
write_xlsx(CRISPRko_df_list,
           path = file.path(output_dir, "CRISPRko_sheets.xlsx"),
           format_headers = FALSE
           )




