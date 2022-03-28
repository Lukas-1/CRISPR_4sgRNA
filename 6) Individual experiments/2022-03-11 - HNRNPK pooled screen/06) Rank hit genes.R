### 12th March 2022 ###


# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-03-11 - HNRNPK pooled screen")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Read in data, and rank hit genes.RData"))
load(file.path(R_objects_directory, "04) Add data on the fitness effect with control treatment.RData"))
load(file.path(R_objects_directory, "05) Merge gene essentiality data with 4sg library coordinates.RData"))




# Integrate 'baseline vs. NT' data with 'NT vs. HNRNPK' data --------------

matches_vec <- match(NT_v_HNRNPK_df[, "gene_id"], base_v_HNRNPK_df[, "gene_id"])

stopifnot(!(anyNA(matches_vec)))

stopifnot(identical(base_v_HNRNPK_df[matches_vec, "Gene_symbol"],
                    NT_v_HNRNPK_df[, "Gene_symbol"]
                    ))

integrated_df <- data.frame(
  "sgRNA_sequence" = sapply(strsplit(NT_v_HNRNPK_df[["gene_id"]], "-", fixed = TRUE), "[[", 2),
  NT_v_HNRNPK_df["Gene_symbol"],
  "BasevNT_log2FC"       = base_v_HNRNPK_df[matches_vec, "log2 Ratio"],
  "BasevNT_hit_strength" = base_v_HNRNPK_df[matches_vec, "log2 Ratio"] *
                           (-log10(base_v_HNRNPK_df[matches_vec, "pValue"])),
  "BasevNT_p_value"      = base_v_HNRNPK_df[matches_vec, "pValue"],
  "BasevNT_FDR"          = base_v_HNRNPK_df[matches_vec, "fdr"],

  "NTvHN_log2FC"         = NT_v_HNRNPK_df[, "log2 Ratio"],
  "NTvHN_hit_strength"   = NT_v_HNRNPK_df[["log2 Ratio"]] *
                           (-log10(NT_v_HNRNPK_df[["pValue"]])),
  "NTvHN_p_value"        = NT_v_HNRNPK_df[, "pValue"],
  "NTvHN_FDR"            = NT_v_HNRNPK_df[, "fdr"],

  "Is_intersect_hit"     = NT_v_HNRNPK_df[["Is_hit"]],
  stringsAsFactors       = FALSE,
  row.names              = NULL
)

table(integrated_df[["Is_intersect_hit"]])




# Add data on library coordinates -----------------------------------------

use_columns <- c("Gene_symbol", "Entrez_ID",
                 "Sublibrary_4sg", "Plate_ID", "Well_number"
                 )

df_list <- lapply(integrated_df[, "Gene_symbol"], function(x) {
  if (is.na(x)) {
    return(expanded_df[NA_integer_, use_columns])
  }
  are_this_gene <- expanded_df[["Gene_symbol"]] %in% x
  if (any(are_this_gene)) {
    return(expanded_df[are_this_gene, use_columns])
  } else {
    return(expanded_df[NA_integer_, use_columns])
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
  bound_df[, !(names(bound_df) %in% "Gene_symbol")],
  integrated_df[expanded_indices, 3:ncol(integrated_df)],
  stringsAsFactors = FALSE,
  row.names = NULL
)



# Re-order hit genes (by hit strength, +- essentiality) -------------------

inters_HS_df <- integrated_df[integrated_df[["Is_intersect_hit"]], ]

exclude_columns <- c("TSS_number", "Is_intersect_hit")
inters_HS_df <- inters_HS_df[, !(names(inters_HS_df) %in% exclude_columns)]

fitness_log2FC_cutoff <- log2(0.5)

inters_HS_df[, "Is_not_essential"] <- inters_HS_df[, "BasevNT_log2FC"] >
                                      fitness_log2FC_cutoff

new_order <- order(inters_HS_df[, "NTvHN_hit_strength"],
                   decreasing = TRUE
                   )
inters_HS_df <- inters_HS_df[new_order, ]
row.names(inters_HS_df) <- NULL

new_order <- order(inters_HS_df[, "Is_not_essential"],
                   inters_HS_df[, "NTvHN_hit_strength"],
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
     file = file.path(R_objects_directory, "06) Rank hit genes.RData")
     )



