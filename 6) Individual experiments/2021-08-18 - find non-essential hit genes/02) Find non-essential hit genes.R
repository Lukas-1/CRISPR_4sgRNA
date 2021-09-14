### 18 August 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

essential_genes_path  <- file.path(experiments_directory,
                                   "2021-08-16 - annotate the Brie library with data on essential genes",
                                   "2) R objects",
                                   "Annotate the Brie library with data on essential genes.RData"
                                   )

file_directory        <- file.path(experiments_directory, "2021-08-18 - find non-essential hit genes")
file_input_directory  <- file.path(file_directory, "1) Input")
top_genes_directory   <- file.path(file_input_directory, "Pooled screen", "Top genes")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")




# Load data ---------------------------------------------------------------

load(essential_genes_path)

load(file.path(R_objects_directory, "01) Identify and rank hit genes.RData"))




# Read in data ------------------------------------------------------------

gt_df <- read.delim(file.path(file_input_directory, "Other", "scGT vs NGT.txt"),
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                    )




# Explore data ------------------------------------------------------------

table(mouse_essential_df[["CRISPR_common"]], useNA = "ifany")





# Create categories -------------------------------------------------------

are_common_essential <- mouse_essential_df[["CRISPR_common"]] == "Essential"
mouse_essential_df[["Category"]] <- ifelse(are_common_essential,
                                           "Essential",
                                           ifelse(mouse_essential_df[["CRISPR_num_essential"]] == 0,
                                                  "Non-essential",
                                                  "Intermediate"
                                                  )
                                           )

other_columns <- c("Achilles_common", "BlomenHart_intersect",
                   "BlomenHart_intersect_DepMap",  "Hart_3_or_more_lines",
                   "Hart_HeLa", "Blomen_HAP1_KBM7_intersect"
                   )

other_mat <- as.matrix(mouse_essential_df[, other_columns])

consensus_vec <- apply(other_mat, 1, function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else if (all(c("Essential", "Non-essential") %in% x)) {
    return("Both")
  } else {
    results_vec <- unique(x)
    results_vec <- results_vec[!(is.na(results_vec))]
    return(results_vec)
  }
})

alternative_categories <- vapply(seq_along(consensus_vec), function(x) {
  if (is.na(are_common_essential[[x]])) {
    if (consensus_vec[[x]] %in% "Both") {
      "Intermediate"
    } else {
      consensus_vec[[x]]
    }
  } else if (!(are_common_essential[[x]])) {
    if (consensus_vec[[x]] %in% "Non-essential") {
      "Non-essential"
    } else {
      "Intermediate"
    }
  } else {
    if (consensus_vec[[x]] %in% "Essential") {
      "Essential"
    } else {
      "Intermediate"
    }
  }
}, "")

alternative_categories <- ifelse((mouse_essential_df[["DEMETER2_mean_probability"]] > 0.5) %in% TRUE,
                                 "Essential",
                                 alternative_categories
                                 )

alternative_categories <- ifelse((alternative_categories == "Non-essential") &
                                 ((mouse_essential_df[["DEMETER2_num_essential"]] > 0) %in% TRUE),
                                 "Intermediate",
                                 alternative_categories
                                 )


are_NA <- is.na(mouse_essential_df[["Category"]])

mouse_essential_df[["Category"]][are_NA] <- alternative_categories[are_NA]

mouse_essential_df[["Category"]] <- factor(mouse_essential_df[["Category"]],
                                           levels = c("Non-essential",
                                                      "Intermediate",
                                                      "Essential"
                                                      )
                                           )




# Combine the data frame --------------------------------------------------

df_list <- lapply(up_genes_df[, "Gene_symbol"], function(x) {
  are_this_gene <- mouse_essential_df[["Mouse_symbol"]] %in% x
  if (!(any(are_this_gene))) {
    sub_df <- mouse_essential_df[NA_integer_, ]
  } else {
    sub_df <- mouse_essential_df[are_this_gene, ]
  }
  categories_fac <- sub_df[["Category"]]
  are_NA <- is.na(categories_fac)
  if (all(are_NA)) {
    category <- NA_character_
  } else {
    categories_fac <- unique(categories_fac[!(are_NA)])
    category <- paste0(as.character(categories_fac), collapse = "/")
  }
  sub_df[["Combined_category"]] <- category
  row.names(sub_df) <- NULL
  return(sub_df)
})

combined_df <- do.call(rbind.data.frame,
                       c(df_list,
                         stringsAsFactors = FALSE,
                         make.row.names = FALSE
                         )
                       )

indices_vec <- rep(seq_len(nrow(up_genes_df)),
                   vapply(df_list, nrow, integer(1))
                   )
combined_df <- data.frame(
  up_genes_df[indices_vec, ],
  combined_df[, names(combined_df) != "Mouse_symbol"],
  stringsAsFactors = FALSE,
  row.names = NULL
)
names(combined_df)[[1]] <- "Mouse_symbol"

column_indices <- c(1, 3, 2, 4, ncol(combined_df) - 1,
                    5:(ncol(combined_df) - 2), ncol(combined_df)
                    )

combined_df <- combined_df[, column_indices]




# Re-order the data frame -------------------------------------------------

combined_categories <- c("Non-essential",
                         "Intermediate",
                         "Essential/Non-essential",
                         "Essential"
                         )

stopifnot(all(combined_df[["Combined_category"]] %in% c(NA, combined_categories)))

new_order <- order(match(combined_df[["Combined_category"]], combined_categories))
combined_df <- combined_df[new_order, ]
row.names(combined_df) <- NULL

table(combined_df[["Category"]], useNA = "ifany")




# Add data on expression in the GT1-7 cell line ---------------------------

matches_vec <- match(combined_df[, "Mouse_symbol"], gt_df[, "gene_name"])
combined_df[["Expressed_in_GT17"]] <- gt_df[matches_vec, "isPresent"]

column_indices <- c(1:7, ncol(combined_df), 8:(ncol(combined_df) - 1))
combined_df <- combined_df[, column_indices]

combined_df[["Expressed_in_GT17"]] <- ifelse(combined_df[["Expressed_in_GT17"]],
                                             "Yes", "No"
                                             )



# Prepare for export ------------------------------------------------------

export_df <- combined_df
export_df[["Category"]] <- as.character(export_df[["Category"]])
keep_NA_columns <- c("CRISPR_common", "Achilles_common", "BlomenHart_intersect",
                     "BlomenHart_intersect_DepMap",  "Hart_3_or_more_lines",
                     "Hart_HeLa", "Blomen_HAP1_KBM7_intersect",
                     "Category"
                      )
for (column_name in setdiff(names(export_df), keep_NA_columns)) {
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
          file = file.path(file_output_directory, "Up_genes_essentiality.csv"),
          row.names = FALSE,
          na = "N/A"
          )



# Save data ---------------------------------------------------------------

save(combined_df,
     file = file.path(R_objects_directory, "02) Find non-essential hit genes.RData")
     )



