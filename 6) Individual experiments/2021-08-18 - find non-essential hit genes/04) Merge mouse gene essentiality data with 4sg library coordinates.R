### 14th September 2021 ###


# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
experiments_directory   <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory          <- file.path(experiments_directory, "2021-08-18 - find non-essential hit genes")
R_objects_directory     <- file.path(file_directory, "2) R objects")

library_RData_directory <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(library_RData_directory, "2) CRISPRa")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "02) Find non-essential hit genes.RData"))
load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))




# Try stuff ---------------------------------------------------------------

plate_IDs <- sapply(strsplit(full_4sg_by_well_df[["Plate_string"]], "_", fixed = TRUE), "[[", 2)

full_4sg_by_well_df[["Plate_ID"]] <- plate_IDs

are_obsolete <- full_4sg_by_well_df[["Is_obsolete"]] %in% "Yes"
full_4sg_by_well_df <- full_4sg_by_well_df[!(are_obsolete), ]
row.names(full_4sg_by_well_df) <- NULL

use_columns <- c("Sublibrary_4sg", "Plate_ID", "Well_number", "Gene_symbol",
                 "Entrez_ID", "TSS_number", "TSS_ID", "Is_main_TSS"
                 )
coordinates_df <- unique(full_4sg_by_well_df[, use_columns])

num_TSSs <- as.integer(table(coordinates_df[["Entrez_ID"]])[coordinates_df[["Entrez_ID"]]])
coordinates_df[["Is_main_TSS"]] <- ifelse(num_TSSs == 1,
                                          "Single",
                                          coordinates_df[["Is_main_TSS"]]
                                          )

for (column_name in c("TSS_ID", "TSS_number")) {
  coordinates_df[[column_name]] <- ifelse(num_TSSs == 1,
                                          NA,
                                          coordinates_df[[column_name]]
                                          )
}


df_list <- lapply(mouse_essential_df[, "Human_Entrez_ID"], function(x) {
  if (is.na(x)) {
    return(coordinates_df[NA_integer_, ])
  }
  are_this_gene <- coordinates_df[["Entrez_ID"]] %in% x
  if (any(are_this_gene)) {
    return(coordinates_df[are_this_gene, ])
  } else {
    return(coordinates_df[NA_integer_, ])
  }
})

bound_df <- do.call(rbind.data.frame,
                    c(df_list,
                      stringsAsFactors = FALSE,
                      make.row.names = FALSE
                      ))

expanded_indices <- rep(seq_len(nrow(mouse_essential_df)),
                        vapply(df_list, nrow, integer(1))
                        )

mouse_expanded_df <- data.frame(
  mouse_essential_df[expanded_indices, 1:5],
  bound_df[, !(names(bound_df) %in% c("Entrez_ID", "Gene_symbol"))],
  mouse_essential_df[expanded_indices, 6:ncol(mouse_essential_df)],
  stringsAsFactors = FALSE,
  row.names = NULL
)




# Save data ---------------------------------------------------------------

save(list = "mouse_expanded_df",
     file = file.path(R_objects_directory, "04) Merge mouse gene essentiality data with 4sg library coordinates.RData")
     )




