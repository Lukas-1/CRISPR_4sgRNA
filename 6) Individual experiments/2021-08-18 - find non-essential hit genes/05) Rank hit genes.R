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
load(file.path(R_objects_directory, "04) Merge mouse gene essentiality data with 4sg library coordinates.RData"))




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

table(integrated_df[["Is_intersect_hit"]])




# Add data on library coordinates -----------------------------------------

use_columns <- c("Mouse_symbol",
                 "Human_symbol", "Human_Entrez_ID",
                 "Sublibrary_4sg", "Plate_ID", "Well_number"
                 )

df_list <- lapply(integrated_df[, "Gene_symbol"], function(x) {
  if (is.na(x)) {
    return(mouse_expanded_df[NA_integer_, ])
  }
  are_this_gene <- mouse_expanded_df[["Mouse_symbol"]] %in% x
  if (any(are_this_gene)) {
    return(mouse_expanded_df[are_this_gene, use_columns])
  } else {
    return(mouse_expanded_df[NA_integer_, ])
  }
})

bound_df <- do.call(rbind.data.frame,
                    c(df_list,
                      stringsAsFactors = FALSE,
                      make.row.names = FALSE
                      ))

expanded_indices <- rep(seq_len(nrow(integrated_df)),
                        vapply(df_list, nrow, integer(1))
                        )

integrated_df <- data.frame(
  integrated_df[expanded_indices, 1:2],
  bound_df[, !(names(bound_df) %in% "Mouse_symbol")],
  integrated_df[expanded_indices, 3:ncol(integrated_df)],
  stringsAsFactors = FALSE,
  row.names = NULL
)



# Re-order hit genes (by hit strength, +- essentiality) -------------------

inters_HS_df <- integrated_df[integrated_df[["Is_intersect_hit"]], ]

exclude_columns <- c("TSS_number", "Is_intersect_hit")
inters_HS_df <- inters_HS_df[, !(names(inters_HS_df) %in% exclude_columns)]

fitness_log2FC_cutoff <- log2(0.5)

inters_HS_df[, "Is_not_essential"] <- inters_HS_df[, "B1vD5_log2FC"] >
                                      fitness_log2FC_cutoff

new_order <- order(inters_HS_df[, "YM5vD5_hit_strength"],
                   decreasing = TRUE
                   )
inters_HS_df <- inters_HS_df[new_order, ]
row.names(inters_HS_df) <- NULL

new_order <- order(inters_HS_df[, "Is_not_essential"],
                   inters_HS_df[, "YM5vD5_hit_strength"],
                   decreasing = TRUE
                   )
inters_NE_df <- inters_HS_df[new_order, ]
row.names(inters_NE_df) <- NULL



# Identify top 10 overall genes and top 10 non-essential genes ------------

number_of_genes <- 10L

top_10_HS <- unique(inters_HS_df[, "Gene_symbol"])[seq_len(number_of_genes)]

intersect(top_10_HS, inters_NE_df[seq_len(number_of_genes), "Gene_symbol"])

unique(inters_NE_df[["Gene_symbol"]][inters_NE_df[["Is_not_essential"]]])

top_10_NE <- setdiff(unique(inters_NE_df[["Gene_symbol"]][inters_NE_df[["Is_not_essential"]]]),
                     top_10_HS
                     )[seq_len(number_of_genes)]



# Export data -------------------------------------------------------------

write.csv(inters_HS_df,
          file = file.path(file_output_directory, "Top genes",
                           "Top sgRNAs - 1) ranked by hit strength.csv"
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

save(list = c("top_10_HS", "top_10_NE",
              "integrated_df",
              "inters_HS_df", "inters_NE_df"
              ),
     file = file.path(R_objects_directory, "04) Rank hit genes.RData")
     )



