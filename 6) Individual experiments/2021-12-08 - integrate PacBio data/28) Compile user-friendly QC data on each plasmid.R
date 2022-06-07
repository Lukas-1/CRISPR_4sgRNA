### 17th April 2022



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")
source(file.path(p1_R_functions_directory, "08) Processing demultiplexed PacBio reads.R")) # For ExportTable
source(file.path(s2r1_R_functions_directory, "03) Tidying CRISPR gRNA data frames.R"))



# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")
file_output_directory    <- file.path(s2rI_directory, "5) Output")
tables_output_directory  <- file.path(file_output_directory, "Tables")
pretty_tables_directory  <- file.path(tables_output_directory, "Summary tables", "User-friendly tables")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(library_RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(library_RData_directory, "3) CRISPRko")



# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Select the quality control data used ------------------------------------

use_summary_df <- ccs3_df_list[["original_summary_df"]]
stopifnot(identical(use_summary_df[, "Combined_ID"], library_df[, "Combined_ID"]))
use_reads_df <- ccs3_df_list[["individual_reads_df"]]



# Summarize the major contaminating wells ---------------------------------

wells_fac <- factor(use_reads_df[, "Combined_ID"], levels = unique(use_reads_df[, "Combined_ID"]))
contam_wells_list <- tapply(use_reads_df[, "Contaminating_well"],
                            wells_fac,
                            function(x) unlist(strsplit(x, ", ", fixed = TRUE)),
                            simplify = FALSE
                            )
contam_wells_list <- lapply(contam_wells_list, function(x) x[!(is.na(x))])
num_contam_list <- lapply(contam_wells_list, function(x) if (length(x) == 0) integer(0) else table(x, dnn = NULL))
num_contam_list <- lapply(num_contam_list, sort, decreasing = TRUE)
counts_vec <- tabulate(wells_fac)
perc_contam_list <- mapply("/", num_contam_list, counts_vec)
perc_contam_list <- lapply(perc_contam_list, function(x) x[x >= 0.05])
contaminating_wells_vec <- vapply(perc_contam_list, function(x) {
  if (length(x) == 0) {
    NA_character_
  } else {
    paste0(names(x), collapse = ", ")
  }
}, "")

stopifnot(all(use_reads_df[, "Combined_ID"] %in% use_summary_df[, "Combined_ID"]))
matches_vec <- match(use_summary_df[, "Combined_ID"], levels(wells_fac))
aligned_contam_vec <- contaminating_wells_vec[matches_vec]




# Integrate data ----------------------------------------------------------

library_columns <- c(
  "Entrez_ID", "Gene_symbol", "TSS_ID", "Modality",
  "Plate_string_4sg", "Well_4sg"
)

count_columns <- c(
  "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
  "Num_reads_with_sgRNA_deletion", "Num_contaminated_reads"
)
perc_mat <- apply(use_summary_df[, count_columns], 2, function(x) x / use_summary_df[, "Count_total"])
colnames(perc_mat) <- sub("^(Num|Count)", "Percent", colnames(perc_mat))

combined_df <- data.frame(
  library_df[, library_columns],
  use_summary_df["Count_total"],
  perc_mat,
  "Contaminating_wells" = aligned_contam_vec,
  stringsAsFactors = FALSE,
  row.names = NULL
)

combined_df[, "Gene_symbol"] <- ifelse(is.na(combined_df[, "Gene_symbol"]),
                                       library_df[, "Combined_ID"],
                                       combined_df[, "Gene_symbol"]
                                       )

rename_columns <- c(
  "Count_total"      = "Number_of_reads",
  "Plate_string_4sg" = "Library_plate",
  "Well_4sg"         = "Library_well",
  "TSS_ID"           = "TSS"
)

for (column_name in names(rename_columns)) {
  names(combined_df)[names(combined_df) == column_name] <- rename_columns[[column_name]]
}



# Re-order data -----------------------------------------------------------

new_order <- order(combined_df[, "Modality"],
                   as.integer(sub("H[AO]_", "", sub("+", "", combined_df[, "Library_plate"], fixed = TRUE))),
                   combined_df[, "Library_plate"]
                   )
combined_df <- combined_df[new_order, ]
row.names(combined_df) <- NULL



# Export data -------------------------------------------------------------

CRISPRa_df  <- combined_df[combined_df[, "Modality"] %in% "CRISPRa", ]
CRISPRko_df <- combined_df[combined_df[, "Modality"] %in% "CRISPRko", ]

row.names(CRISPRa_df) <- NULL
row.names(CRISPRko_df) <- NULL

CRISPRa_df <- CRISPRa_df[, names(CRISPRa_df) != "Modality"]
CRISPRko_df <- CRISPRko_df[, names(CRISPRko_df) != "Modality"]

ExportTable(CRISPRa_df,
            file_name_only = "T.gonfio (CRISPRa) QC data",
            file_directory = pretty_tables_directory
            )

ExportTable(CRISPRko_df,
            file_name_only = "T.spiezzo (CRISPRko) QC data",
            file_directory = pretty_tables_directory
            )



