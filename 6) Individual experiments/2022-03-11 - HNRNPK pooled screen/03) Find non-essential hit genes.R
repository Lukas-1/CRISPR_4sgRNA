### 12th March 2022 ###



# Import packages and source code -----------------------------------------

library("readxl")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-03-11 - HNRNPK pooled screen")
R_objects_directory   <- file.path(file_directory, "2) R objects")
file_output_directory <- file.path(file_directory, "3) Output")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Read in data, and rank hit genes.RData"))
load(file.path(R_objects_directory, "02) Compile data on essential genes - essential_df.RData"))





# Create categories -------------------------------------------------------

are_common_essential <- essential_df[["CRISPR_common"]] == "Essential"

essential_df[is.na(are_common_essential), ]

essential_df[["Category"]] <- ifelse(are_common_essential,
                                     "Essential",
                                     ifelse(essential_df[["CRISPR_num_essential"]] == 0,
                                            "Non-essential",
                                            "Intermediate"
                                            )
                                     )

other_columns <- c("Achilles_common", "BlomenHart_intersect",
                   "BlomenHart_intersect_DepMap",  "Hart_3_or_more_lines",
                   "Hart_HeLa", "Blomen_HAP1_KBM7_intersect"
                   )

other_mat <- as.matrix(essential_df[, other_columns])

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

alternative_categories <- ifelse((essential_df[["DEMETER2_mean_probability"]] > 0.5) %in% TRUE,
                                 "Essential",
                                 alternative_categories
                                 )

alternative_categories <- ifelse((alternative_categories == "Non-essential") &
                                 ((essential_df[["DEMETER2_num_essential"]] > 0) %in% TRUE),
                                 "Intermediate",
                                 alternative_categories
                                 )

are_NA <- is.na(essential_df[["Category"]])

essential_df[["Category"]][are_NA] <- alternative_categories[are_NA]

essential_df[["Category"]] <- factor(essential_df[["Category"]],
                                     levels = c("Non-essential",
                                                "Intermediate",
                                                "Essential"
                                                )
                                     )




# Combine the data frame --------------------------------------------------

df_list <- lapply(up_genes_df[, "Gene_symbol"], function(x) {
  are_this_gene <- essential_df[["Gene_symbol"]] %in% x
  if (!(any(are_this_gene))) {
    sub_df <- essential_df[NA_integer_, ]
  } else {
    sub_df <- essential_df[are_this_gene, ]
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
  combined_df[, names(combined_df) != "Gene_symbol"],
  stringsAsFactors = FALSE,
  row.names = NULL
)

column_indices <- c(1, 3, 2, ncol(combined_df) - 1,
                    4:(ncol(combined_df) - 2), ncol(combined_df)
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
                  "Combined_category", "Four_categories"
                  )

export_df <- export_df[, !(names(export_df) %in% omit_columns)]



# Export data -------------------------------------------------------------

write.csv(export_df,
          file = file.path(file_output_directory, "Up_genes_essentiality.csv"),
          row.names = FALSE,
          na = "N/A"
          )



# Save data ---------------------------------------------------------------

save(list = c("combined_df", "essential_df"),
     file = file.path(R_objects_directory, "03) Find non-essential hit genes.RData")
     )



