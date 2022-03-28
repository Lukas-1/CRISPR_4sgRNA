### 12th March 2022 ###


# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-03-11 - HNRNPK pooled screen")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Read in data, and rank hit genes.RData"))
load(file.path(R_objects_directory, "03) Find non-essential hit genes.RData"))



# Define functions --------------------------------------------------------

MakeSummaryDf <- function(input_df, use_suffix = NULL) {
  genes_fac <- factor(input_df[, "Gene_symbol"])
  genes_df <- data.frame(
    "Gene_symbol"           = levels(genes_fac),
    "Category"              = NA,
    "Mean_log2FC_essential" = tapply(input_df[["log2 Ratio"]], genes_fac, mean),
    "Num_essential_gRNAs"   = tapply(input_df[["Is_fitness_relevant"]], genes_fac, sum),
    "Min_log2FC_essential"  = tapply(input_df[["log2 Ratio"]], genes_fac, min),
    "Min_FDR_essential"     = tapply(input_df[["fdr"]], genes_fac, min),
    stringsAsFactors        = FALSE,
    row.names               = NULL
  )
  genes_df[["Category"]] <- ifelse(genes_df[["Min_log2FC_essential"]] < log2(0.5),
                                   "Essential",
                                   ifelse(genes_df[["Min_log2FC_essential"]] < log2(0.8),
                                          "Intermediate",
                                          "Non-essential"
                                          )
                                   )
  if (!(is.null(use_suffix))) {
    names(genes_df)[2:6] <- paste0(names(genes_df)[2:6], use_suffix)
  }
  return(genes_df)
}




# Add data on "hit" status ------------------------------------------------

hit_gRNAs <- NT_v_HNRNPK_df[["gene_id"]][NT_v_HNRNPK_df[["Is_hit"]]]
base_v_NT_df[["Is_hit"]] <- base_v_NT_df[["gene_id"]] %in% hit_gRNAs




# Create gene-level summary data frames -----------------------------------

base_v_NT_df[["Is_fitness_relevant"]] <- (base_v_NT_df[, "log2 Ratio"] < -1) &
                                         (base_v_NT_df[, "fdr"] < 0.01)

hits_base_v_NT_df <- base_v_NT_df[base_v_NT_df[["Is_hit"]], ]

all_df <- MakeSummaryDf(base_v_NT_df, use_suffix = "_all_gRNAs")
hits_df <- MakeSummaryDf(hits_base_v_NT_df, use_suffix = "_hit_gRNAs")



# Combine the data --------------------------------------------------------

all_matches_vec <- match(combined_df[, "Gene_symbol"], all_df[, "Gene_symbol"])
hits_matches_vec <- match(combined_df[, "Gene_symbol"], hits_df[, "Gene_symbol"])

combined_df <- data.frame(combined_df[, 1:4],
                          all_df[all_matches_vec, 2:6],
                          hits_df[hits_matches_vec, 2:6],
                          combined_df[, 5:ncol(combined_df)],
                          stringsAsFactors = FALSE,
                          row.names = NULL
                          )



# Prepare for export ------------------------------------------------------

export_df <- combined_df
export_df[["Category"]] <- as.character(export_df[["Category"]])
keep_NA_columns <- c("CRISPR_common", "Achilles_common", "BlomenHart_intersect",
                     "BlomenHart_intersect_DepMap",  "Hart_3_or_more_lines",
                     "Hart_HeLa", "Blomen_HAP1_KBM7_intersect",
                     "Category", "Category_all_gRNAs", "Category_hit_gRNAs"
                      )
for (column_name in setdiff(names(export_df), keep_NA_columns)) {
  export_df[, column_name] <- ifelse(is.na(export_df[, column_name]),
                                     "",
                                     export_df[, column_name]
                                     )
}

omit_columns <- c("Achilles_mean_probability", "Achilles_num_essential",
                  "Achilles_num_cell_lines",
                  "Combined_category", "Four_categories"
                  )

export_df <- export_df[, !(names(export_df) %in% omit_columns)]




# Export data -------------------------------------------------------------

write.csv(export_df,
          file = file.path(file_output_directory, "Up_genes_essential_base_v_NT.csv"),
          row.names = FALSE,
          na = "N/A"
          )



# Save data ---------------------------------------------------------------

save(list = "base_v_NT_df",
     file = file.path(R_objects_directory, "04) Add data on the fitness effect with control treatment.RData")
     )



