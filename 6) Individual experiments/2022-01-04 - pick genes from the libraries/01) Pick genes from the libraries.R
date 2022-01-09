### 5th January 2022 ##


# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_directory     <- file.path(experiments_directory, "2022-01-04 - pick genes from the libraries")
R_functions_directory <- file.path(project_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Converting plate layouts.R"))



# Define folder paths -----------------------------------------------------

RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")

project_output_directory <- file.path(project_directory, "4) Output")




# Load data ---------------------------------------------------------------

CRISPRa_objects <- load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRa_sgRNA_df <- full_4sg_by_well_df

CRISPRko_objects <- load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_sgRNA_df <- full_4sg_by_well_df

rm(list = union(CRISPRa_objects, CRISPRko_objects))





# Pick genes --------------------------------------------------------------

CRISPRa_columns <- c("Plate_number", "Well_number", "Coords_96wp",
                     "Entrez_ID", "Gene_symbol",
                     "Is_control", "Sublibrary_4sg",
                     "Is_main_TSS", "TSS_number", "TSS_ID"
                     )

TidySgDf <- function(sg_df) {
  sg_df[, "Gene_symbol"] <- ifelse(is.na(sg_df[, "Gene_symbol"]),
                                   sg_df[, "Combined_ID"],
                                   sg_df[, "Gene_symbol"]
                                   )
  plate_string_splits <- strsplit(sg_df[, "Plate_string"], "_", fixed = TRUE)
  sg_df[, "Coords_96wp"] <- ConvertWellNumbers(sg_df[, "Well_number"])
  sg_df[, "Plate_number"] <- sapply(plate_string_splits, "[[", 2)
  use_columns <- intersect(CRISPRa_columns, names(sg_df))
  sg_df <- sg_df[sg_df[, "Rank"] == 1, use_columns]
  row.names(sg_df) <- NULL
  return(sg_df)
}


ExportTable <- function(sg_df, file_name) {
  write.table(sg_df, row.names = FALSE,
              file = file.path(project_output_directory, paste0(file_name, ".tsv")),
              sep = "\t"
              )
}


CRISPRa_df <- TidySgDf(CRISPRa_sgRNA_df)
CRISPRko_df <- TidySgDf(CRISPRko_sgRNA_df)

CRISPRko_df[CRISPRko_df[, "Is_control"] == "Yes", ]

use_genes <- c("PRNP", "IKBKG", "Control_33", "Control_64")

ExportTable(CRISPRko_df[CRISPRko_df[, "Gene_symbol"] %in% use_genes, ], "Four_genes_Lukas")













