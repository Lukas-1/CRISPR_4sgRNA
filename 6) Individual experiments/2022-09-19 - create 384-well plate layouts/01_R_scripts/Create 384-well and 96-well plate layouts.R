## 2022-09-19



# Load packages and source code -------------------------------------------

library("writexl")

CRISPR_root_dir   <- "~/CRISPR_4sgRNA"
experiments_dir   <- file.path(CRISPR_root_dir, "6) Individual experiments")
pick_plasmids_dir <- file.path(experiments_dir, "2022-01-04 - pick genes from the libraries")
source(file.path(pick_plasmids_dir, "1) R functions", "01) Converting plate layouts.R"))


# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_dir, "2022-09-19 - create 384-well plate layouts")
input_dir   <- file.path(project_dir, "02_input_data")
output_dir  <- file.path(project_dir, "03_output_data")



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


SplitInto384WellPlates <- function(input_df) {
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



SplitInto96WellPlates <- function(input_df) {

  plates_vec <- input_df[, "Plate_ID"]
  plates_fac <- factor(plates_vec, levels = unique(plates_vec))
  plates_df_list <- split(input_df[, c("Plasmid_name", "Well_number")], plates_fac)

  coords_vec <- ConvertWellNumbers(1:384)
  coords_splits <- strsplit(coords_vec, ", ", fixed = TRUE)

  four_96wp_vec <- sub("Plate ", "", sapply(coords_splits, "[[", 1), fixed = TRUE)
  well_96wp_vec <- sub("well ", "", sapply(coords_splits, "[[", 2), fixed = TRUE)
  row_96wp_vec <- substr(well_96wp_vec, 1, 1)
  column_96wp_vec <- as.integer(substr(well_96wp_vec, 2, nchar(well_96wp_vec)))
  map_df <- data.frame(
    "Well_number_384" = 1:384,
    "Sub_plate_96"    = four_96wp_vec,
    "Row_96"          = row_96wp_vec,
    "Column_96"       = column_96wp_vec
  )

  mat_96wp <- matrix("", nrow = 8, ncol = 12)
  rownames(mat_96wp) <- LETTERS[1:8]
  mat_list <- rep(list(mat_96wp), 4)
  names(mat_list) <- c("A1", "A2", "B1", "B2")

  plates_list_list <- lapply(plates_df_list, function(x) {
    result_list <- mat_list
    matches_vec <- match(x[, "Well_number"], map_df[, "Well_number_384"])
    for (i in seq_along(matches_vec)) {
      use_index <- matches_vec[[i]]
      sub_plate <- map_df[, "Sub_plate_96"][[i]]
      row_name <- map_df[, "Row_96"][[i]]
      column_index <- map_df[, "Column_96"][[i]]
      result_list[[sub_plate]][row_name, column_index] <- x[, "Plasmid_name"][[i]]
    }
    result_list
  })

  mod_plates_list_list <- lapply(plates_list_list, function(x) {
    lapply(1:4, function(y) {
      results_mat <- rbind(rep("", 12), 1:12, x[[y]])
      results_mat <- cbind(c(paste0("Plate ", names(x)[[y]]), "", LETTERS[1:8]), results_mat)
      dimnames(results_mat) <- NULL
      return(results_mat)
    })
  })

  mat_list <- lapply(mod_plates_list_list, function(x) {
    upper_mat <- cbind(x[[1]], matrix("", nrow = nrow(x[[1]]), ncol = 2), x[[2]])
    lower_mat <- cbind(x[[3]], matrix("", nrow = nrow(x[[3]]), ncol = 2), x[[4]])
    rbind(upper_mat,
          matrix("", nrow = 3, ncol = ncol(upper_mat)),
          lower_mat
          )
  })

  df_list <- lapply(mat_list, function(x) {
    results_df <- data.frame(x)
    names(results_df) <- rep("", ncol(results_df))
    return(results_df)
  })
  return(df_list)
}



# Prepare data ------------------------------------------------------------

## CRISPRa
CRISPRa_df[, "Plasmid_name"] <- paste0(CRISPRa_df[, "Gene_symbol"],
                                       ifelse(nzchar(CRISPRa_df[, "TSS_ID"]), "_", ""),
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

CRISPRa_384_df_list  <- SplitInto384WellPlates(CRISPRa_df)
CRISPRko_384_df_list <- SplitInto384WellPlates(CRISPRko_df)

CRISPRa_96_df_list   <- SplitInto96WellPlates(CRISPRa_df)
CRISPRko_96_df_list  <- SplitInto96WellPlates(CRISPRko_df)




# Export Excel sheets -----------------------------------------------------

write_xlsx(CRISPRa_384_df_list,
           path = file.path(output_dir, "CRISPRa_384wp_sheets.xlsx"),
           format_headers = FALSE
           )
write_xlsx(CRISPRko_384_df_list,
           path = file.path(output_dir, "CRISPRko_384wp_sheets.xlsx"),
           format_headers = FALSE
           )

write_xlsx(CRISPRa_96_df_list,
           path = file.path(output_dir, "CRISPRa_96wp_sheets.xlsx"),
           format_headers = FALSE
           )
write_xlsx(CRISPRko_96_df_list,
           path = file.path(output_dir, "CRISPRko_96wp_sheets.xlsx"),
           format_headers = FALSE
           )



