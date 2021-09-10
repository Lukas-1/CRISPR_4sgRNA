### 9 September 2021 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2021-08-18 - find non-essential hit genes")
file_input_directory  <- file.path(file_directory, "1) Input")
screen_data_directory <- file.path(file_input_directory, "Pooled screen", "Data")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Find non-essential hit genes.RData"))





# Read in data ------------------------------------------------------------

baseline_df <- read.csv(file.path(screen_data_directory,
                                  "DMSO5--over--B1_DOWN_GeneName.csv"
                                  ),
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                        )




# Tidy data ---------------------------------------------------------------

are_empty_columns <- vapply(baseline_df, function(x) {
  length(unique(x)) == 1
}, logical(1))

baseline_df <- baseline_df[, !(are_empty_columns)]





# Set a threshold ---------------------------------------------------------

baseline_df[["Pass_threshold"]] <- (baseline_df[, "log2 Ratio"] < -1) &
                                   (baseline_df[, "fdr"] < 0.01)

genes_fac <- factor(baseline_df[["Gene_symbol"]])

baseline_genes_df <- data.frame(
  "Gene_symbol"           = levels(genes_fac),
  "Mean_log2FC_essential" = tapply(baseline_df[["log2 Ratio"]], genes_fac, mean),
  "Num_essential_gRNAs"   = tapply(baseline_df[["Pass_threshold"]], genes_fac, sum),
  "Min_log2FC_essential"  = tapply(baseline_df[["log2 Ratio"]], genes_fac, min),
  "Min_FDR_essential"     = tapply(baseline_df[["fdr"]], genes_fac, min),
  stringsAsFactors        = FALSE,
  row.names               = NULL
)

matches_vec <- match(combined_df[, "Mouse_symbol"],
                     baseline_genes_df[, "Gene_symbol"]
                     )

combined_df <- data.frame(combined_df,
                          baseline_genes_df[matches_vec, 2:5],
                          stringsAsFactors = FALSE,
                          row.names = NULL
                          )

column_indices <- c(1:7,
                    (ncol(combined_df) - 4):ncol(combined_df),
                    8:(ncol(combined_df) - 5)
                    )
combined_df <- combined_df[, column_indices]




# Prepare for export ------------------------------------------------------

export_df <- combined_df
NA_empty_columns <- c("CRISPR_mean_probability", "CRISPR_num_essential",
                      "CRISPR_num_cell_lines", "CRISPR_mean_effect",
                      "Achilles_mean_probability", "Achilles_num_essential", "Achilles_num_cell_lines",
                      "DEMETER2_mean_probability", "DEMETER2_num_essential", "DEMETER2_num_cell_lines",
                      "Human_Entrez_ID", "Human_symbol",
                      "Mouse_symbol", "Mouse_Entrez_ID", "Mapping",
                      "Expressed_in_GT17",
                      "Num_essential_gRNAs", "Mean_log2FC_essential",
                      "Min_log2FC_essential", "Min_FDR_essential"
                      )
for (column_name in NA_empty_columns) {
  export_df[, column_name] <- ifelse(is.na(export_df[, column_name]),
                                     "",
                                     export_df[, column_name]
                                     )
}

omit_columns <- c("Achilles_mean_probability", "Achilles_num_essential",
                  "Achilles_num_cell_lines",
                  "Combined_category"
                  )

export_df <- export_df[, !(names(export_df) %in% omit_columns)]





# Export data -------------------------------------------------------------

write.csv(export_df,
          file = file.path(file_output_directory, "Up_genes_essential_B1_v_DMSO5.csv"),
          row.names = FALSE,
          na = "N/A"
          )



