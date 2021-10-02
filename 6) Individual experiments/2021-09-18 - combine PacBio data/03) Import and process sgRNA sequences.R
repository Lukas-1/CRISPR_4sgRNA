### 21st September 2021 ###


# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "19) Annotating duplicated gRNAs.R"))



# Define folder paths -----------------------------------------------------

s2r1_directory           <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")

s2r1_R_objects_directory <- file.path(s2r1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(s2r1_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
run1_library_df <- library_df

load(file.path(s2r2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
run2_library_df <- library_df



# Filter library_df -------------------------------------------------------

non_library_plates <- c("Vac-1", "PD_A", "PD_O")
beads_plates <- unique(grep("beads", run1_library_df[["Plate_name"]], value = TRUE, fixed = TRUE))
exclude_plates <- c(non_library_plates, sub("-beads", "", beads_plates, fixed = TRUE))

are_included <- !(run1_library_df[["Plate_name"]] %in% exclude_plates)
run1_library_df <- run1_library_df[are_included, ]

run1_library_df[["Plate_name"]] <- sub("-beads", "", run1_library_df[["Plate_name"]], fixed = TRUE)

use_columns <- grep("Num_occurrences", names(run2_library_df), value = TRUE, invert = TRUE, fixed = TRUE)

library_df <- rbind.data.frame(
  run1_library_df[, use_columns],
  run2_library_df[, use_columns],
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

library_df <- AddNumOccurrences(library_df)




# Standardize the "Combined_ID" column ------------------------------------

IDs_vec <- paste0("Plate",
                  formatC(library_df[["Plate_number"]], width = 3, flag = "0"),
                  "_Well",
                  formatC(library_df[["Well_number"]], width = 3, flag = "0")
                  )
library_df[["Combined_ID"]] <- IDs_vec




# Save data ---------------------------------------------------------------

save(list = "library_df",
     file = file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData")
     )



