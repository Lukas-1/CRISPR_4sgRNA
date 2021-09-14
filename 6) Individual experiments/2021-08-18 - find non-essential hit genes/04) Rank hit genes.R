### 12 September 2021 ###


# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2021-08-18 - find non-essential hit genes")
file_input_directory  <- file.path(file_directory, "1) Input")
screen_data_directory <- file.path(file_input_directory, "Pooled screen", "Data")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Identify and rank hit genes.RData"))
load(file.path(R_objects_directory, "03) Add data on the fitness effect with control treatment.RData"))




# Integrate 'baseline vs. DMSO5' data with 'YM5 vs. DMSO5' data -----------

matches_vec <- match(YM5_v_DMSO5_df[, "Gene_Name"], baseline_df[, "Gene_Name"])

stopifnot(!(anyNA(matches_vec)))

stopifnot(identical(baseline_df[matches_vec, "Gene_symbol"],
                    YM5_v_DMSO5_df[, "Gene_symbol"]
                    ))

integrated_df <- data.frame(
  "sgRNA_sequence" = YM5_v_DMSO5_df[["Gene_Name"]],
  YM5_v_DMSO5_df["Gene_symbol"],
  "B1vD5_log2FC"        = baseline_df[matches_vec, "log2 Ratio"],
  "B1vD5_p_value"       = baseline_df[matches_vec, "pValue"],
  "B1vD5_FDR"           = baseline_df[matches_vec, "fdr"],

  "YM5vD5_log2FC"       = YM5_v_DMSO5_df[, "log2 Ratio"],
  "YM5vD5_hit_strength" = YM5_v_DMSO5_df[["log2 Ratio"]] *
                          (-log10(YM5_v_DMSO5_df[["pValue"]])),
  "YM5vD5_p_value"      = YM5_v_DMSO5_df[, "pValue"],
  "YM5vD5_FDR"          = YM5_v_DMSO5_df[, "fdr"],

  "Is_intersect_hit"    = YM5_v_DMSO5_df[["Is_hit"]],
  stringsAsFactors      = FALSE,
  row.names             = NULL
)



# Re-order hit genes (by hit strength, +- essentiality) -------------------

inters_ES_df <- integrated_df[integrated_df[["Is_intersect_hit"]], ]

inters_ES_df <- inters_ES_df[, names(inters_ES_df) != "Is_intersect_hit"]

fitness_log2FC_cutoff <- log2(0.5)

inters_ES_df[, "Is_not_essential"] <- inters_ES_df[, "B1vD5_log2FC"] >
                                      fitness_log2FC_cutoff

new_order <- order(inters_ES_df[, "Is_not_essential"],
                   decreasing = TRUE
                   )
inters_ES_df <- inters_ES_df[new_order, ]
row.names(inters_ES_df) <- NULL

new_order <- order(inters_ES_df[, "Is_not_essential"],
                   inters_ES_df[, "YM5vD5_hit_strength"],
                   decreasing = TRUE
                   )
inters_NE_df <- inters_ES_df[new_order, ]
row.names(inters_NE_df) <- NULL




# Identify top 10 overall genes and top 10 non-essential genes ------------

number_of_genes <- 10L

top_10_ES <- unique(inters_ES_df[, "Gene_symbol"])[seq_len(number_of_genes)]

intersect(top_10_ES, inters_NE_df[seq_len(number_of_genes), "Gene_symbol"])


unique(inters_NE_df[["Gene_symbol"]][inters_NE_df[["Is_not_essential"]]])



top_10_NE <- setdiff(unique(inters_NE_df[["Gene_symbol"]][inters_NE_df[["Is_not_essential"]]]),
                     top_10_ES
                     )[seq_len(number_of_genes)]



# Export data -------------------------------------------------------------

write.csv(inters_ES_df,
          file = file.path(file_output_directory, "Top genes",
                           "Top sgRNAs - 1) ranked by effect strength.csv"
                           ),
          row.names = FALSE,
          na = ""
          )

write.csv(inters_NE_df,
          file = file.path(file_output_directory, "Top genes",
                           "Top sgRNAs - 2) non-essential preferred.csv"
                           ),
          row.names = FALSE,
          na = ""
          )


# Save data ---------------------------------------------------------------

save(list = c("top_10_ES", "top_10_NE",
              "integrated_df",
              "inters_ES_df", "inters_NE_df"
              ),
     file = file.path(R_objects_directory, "04) Rank hit genes.RData")
     )



