### 20th September 2021 ###


# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2r1_directory           <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")

s2r1_R_objects_directory <- file.path(s2r1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(s2r1_R_objects_directory, "01) Process and export plate barcodes.RData"))
run1_plates_df <- plates_df

load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
run2_plates_df <- plates_df

load(file.path(s2r2_R_objects_directory, "30) Calculate correction factors to account for mean read counts.RData"))

rm(plates_df)



# Create the new plates_df ------------------------------------------------

run1_plates_df[["DNA_isolation"]] <- ifelse(grepl("beads", run1_plates_df[["Plate_name"]], fixed = TRUE),
                                            "Beads", "Columns"
                                            )

run1_plates_df[["Run2_pool"]] <- 0L
run1_plates_df[["Number_96wp"]] <- NA

run2_plates_df[["DNA_isolation"]] <- "Beads"

common_columns <- intersect(names(run2_plates_df), names(run1_plates_df))

plates_df <- rbind.data.frame(
  run1_plates_df[, common_columns],
  run2_plates_df[, common_columns],
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)



# Save data ---------------------------------------------------------------

save(list = "plates_df",
     file = file.path(s2rC_R_objects_directory, "01) Process and export plate barcodes.RData")
     )



