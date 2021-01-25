### 25th January 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")

file_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-01-25 - provide coordinates of obsolete TF wells")
file_output_directory    <- file.path(file_directory, "2) Output")





# Load data ---------------------------------------------------------------

object_names <- load(file.path(CRISPRa_RData_directory,  "28) Distribute sgRNAs for the whole genome onto plates.RData"))
by_well_CRISPRa_df <- full_4sg_by_well_df
load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
by_well_CRISPRko_df <- full_4sg_by_well_df
rm(list = object_names)





# Filter only obsolete genes ----------------------------------------------

FilterCRISPRDf <- function(CRISPR_df) {
  are_selected <- (CRISPR_df[["Is_obsolete"]] %in% "Yes") &
                  (CRISPR_df[["Rank"]] %in% 1)
  CRISPR_df <- CRISPR_df[are_selected, ]
  row.names(CRISPR_df) <- NULL
  CRISPR_df[["Plate_string"]] <- sub("_sg1", "", CRISPR_df[["Plate_string"]], fixed = TRUE)
  return(CRISPR_df)
}

obsolete_CRISPRa_df <-FilterCRISPRDf(by_well_CRISPRa_df)
obsolete_CRISPRko_df <- FilterCRISPRDf(by_well_CRISPRko_df)





# Export obsolete genes ---------------------------------------------------

ExportPlates(obsolete_CRISPRa_df,  "CRISPRa_obsolete_TF_genes", sub_folder = "")
ExportPlates(obsolete_CRISPRko_df, "CRISPRo_obsolete_TF_genes", sub_folder = "")







