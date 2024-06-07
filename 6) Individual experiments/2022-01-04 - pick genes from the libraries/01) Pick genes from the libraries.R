### 5th January 2022 ##


# Import packages and source code -----------------------------------------

library("readxl")
library("Biostrings")

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_directory     <- file.path(experiments_directory, "2022-01-04 - pick genes from the libraries")
R_functions_directory <- file.path(project_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Converting plate layouts.R"))



# Define folder paths -----------------------------------------------------

RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")

project_input_directory  <- file.path(project_directory, "2) Input")
project_output_directory <- file.path(project_directory, "4) Output")



# Load data ---------------------------------------------------------------

CRISPRa_objects <- load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRa_sgRNA_df <- full_4sg_by_well_df

CRISPRko_objects <- load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_sgRNA_df <- full_4sg_by_well_df

rm(list = union(CRISPRa_objects, CRISPRko_objects))




# Read in data ------------------------------------------------------------

FLAER_df <- data.frame(read_excel(file.path(project_input_directory, "25 genes for ko efficiency 1sg vs 4sg FLAER assay.xlsx")),
                       stringsAsFactors = FALSE
                       )
CRISPRa_genes_df <- data.frame(read_excel(file.path(project_input_directory, "12 genes activation 1sg vs 4sg qPCR.xlsx")),
                               stringsAsFactors = FALSE
                               )



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
  sg_df[, "Coords_96wp"] <- ConvertWellNumbers(sg_df[, "Well_number"])
  plate_string_splits <- strsplit(sg_df[, "Plate_string"], "_", fixed = TRUE)
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





# Pick GPI genes for CRISPRko (for the revision) --------------------------

GPI_genes <- FLAER_df[, 1]
GPI_genes[GPI_genes == "GAA1"] <- "GPAA1"

CRISPRko_sg_mat <- do.call(cbind, lapply(1:4, function(x) {
  CRISPRko_sgRNA_df[CRISPRko_sgRNA_df[, "Rank"] == x, "sgRNA_sequence"]
}))
colnames(CRISPRko_sg_mat) <- paste0("Sequence_sg", 1:4)


fwd_primer_mat <- CRISPRko_sg_mat
fwd_primer_mat <- paste0("accg", fwd_primer_mat)
dim(fwd_primer_mat) <- dim(CRISPRko_sg_mat)
colnames(fwd_primer_mat) <- paste0("Fwd_primer_sg", 1:4)

rev_primer_mat <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(CRISPRko_sg_mat)))
rev_primer_mat <- paste0("AAAC", rev_primer_mat)
dim(rev_primer_mat) <- dim(CRISPRko_sg_mat)
colnames(rev_primer_mat) <- paste0("Rev_primer_sg", 1:4)

use_df <- data.frame("Library" = "T.spiezzo",
                     CRISPRko_df, CRISPRko_sg_mat,
                     fwd_primer_mat, rev_primer_mat,
                     stringsAsFactors = FALSE
                     )


ExportTable(use_df[use_df[, "Gene_symbol"] %in% GPI_genes, names(use_df) != "Is_control"], "Revisions/GPI_genes")




# Pick CRISPRa genes (for the revision) -----------------------------------

CRISPRa_genes <- CRISPRa_genes_df[, "Gene"]

CRISPRa_sg_mat <- do.call(cbind, lapply(1:4, function(x) {
  CRISPRa_sgRNA_df[CRISPRa_sgRNA_df[, "Rank"] == x, "sgRNA_sequence"]
}))
colnames(CRISPRa_sg_mat) <- paste0("Sequence_sg", 1:4)


fwd_primer_mat <- CRISPRa_sg_mat
fwd_primer_mat <- paste0("accg", fwd_primer_mat)
dim(fwd_primer_mat) <- dim(CRISPRa_sg_mat)
colnames(fwd_primer_mat) <- paste0("Fwd_primer_sg", 1:4)

rev_primer_mat <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(CRISPRa_sg_mat)))
rev_primer_mat <- paste0("AAAC", rev_primer_mat)
dim(rev_primer_mat) <- dim(CRISPRa_sg_mat)
colnames(rev_primer_mat) <- paste0("Rev_primer_sg", 1:4)


use_df <- data.frame("Library" = "T.gonfio",
                     CRISPRa_df, CRISPRa_sg_mat,
                     fwd_primer_mat, rev_primer_mat,
                     stringsAsFactors = FALSE
                     )
use_df[, "TSS_ID"] <- ifelse(is.na(use_df[, "TSS_ID"]), " ", use_df[, "TSS_ID"])

use_order <- order(match(use_df[, "Gene_symbol"], use_df[, "Gene_symbol"]),
                   match(use_df[, "Is_main_TSS"], "Main")
                   )
use_df <- use_df[use_order, ]
row.names(use_df) <- NULL

export_df <- use_df[use_df[, "Gene_symbol"] %in% CRISPRa_genes, !(names(use_df) %in% c("Is_control", "TSS_number"))]

ExportTable(export_df, "Revisions/CRISPRa_genes")







