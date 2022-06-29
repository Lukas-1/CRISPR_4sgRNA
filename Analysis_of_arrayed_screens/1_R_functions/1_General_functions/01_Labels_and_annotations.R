# 2021-12-29


# Load packages and source code -------------------------------------------

library("readxl")



# Define folder path ------------------------------------------------------

project_dir <- "~/R_projects/CRISPRa_TF"
input_dir   <- file.path(project_dir, "2_input")



# Read in data ------------------------------------------------------------

columns_df <- data.frame(read_excel(file.path(input_dir, "All_metrics.xlsx")),
                         check.names = FALSE, stringsAsFactors = FALSE
                         )



# Define functions --------------------------------------------------------

AdjustLabels <- function(original_string = "GBA activity",
                         replacement_string = "PrPc levels"
                         ) {


  label_names <- c("column_file_names", "long_column_labels", "short_column_labels")
  for (label_name in label_names) {
    new_vec <- sub(original_string, replacement_string, get(label_name), fixed = TRUE)
    assign(label_name, new_vec, envir = globalenv())
  }
  return(invisible(NULL))
}



# Re-order the columns ----------------------------------------------------

## Order by the type of metric

metrics_in_order <- unique(columns_df[, "Metric"])
are_after <- seq_along(metrics_in_order) > which(metrics_in_order == "Fold-NT")
metrics_in_order <- unique(c(metrics_in_order[!(are_after)], "Log2FC", metrics_in_order[are_after]))


## Place the Glo-normalized version of each metric immediately after its unnormalized equivalent

column_names_stripped <- sub("_Glo", "", columns_df[, "Column name"], fixed = FALSE)
new_order <- order(match(columns_df[, "Metric"], metrics_in_order),
                   match(column_names_stripped, column_names_stripped)
                   )
column_names <- columns_df[, "Column name"][new_order]


## Place the CellTiter-Glo viability values immediately before the SSMD columns

Glo_columns <- c("CellTiterGlo_raw", "CellTiterGlo_foldNT")

are_before <- seq_len(nrow(columns_df)) < which(column_names == "SSMD_deltaNT")
column_names <- unique(c(column_names[are_before], Glo_columns,
                         column_names[!(are_before)])
                         )
columns_df <- columns_df[order(match(columns_df[, "Column name"], column_names)), ]


## Exclude columns that are not needed

are_to_exclude <- ((grepl("_foldNT", columns_df[, "Column name"], fixed = TRUE) &
                   (columns_df[, "Column name"] != "CellTiterGlo_foldNT"))) |  # p values calculated using foldNT values are not plausible
                  (columns_df[, "Metric"] == "t value")                        # t values show a simple linear relationship with SSMD values, and are thus redundant
columns_df <- columns_df[!(are_to_exclude), ]
row.names(columns_df) <- NULL



# Produce named vectors ---------------------------------------------------

column_names <- ifelse(columns_df[, "Replicates"] == "Yes",
                       paste0(columns_df[, "Column name"], "_rep1"),
                       columns_df[, "Column name"]
                       )

column_file_names <- columns_df[, "File name"]
names(column_file_names) <- column_names

long_column_labels <- columns_df[, "Long label"]
names(long_column_labels) <- column_names

short_column_labels <- columns_df[, "Short label"]
names(short_column_labels) <- column_names


