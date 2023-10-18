### 19th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2rC_directory        <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
R_functions_directory <- file.path(s2rC_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Reassigning read unique IDs.R"))




# Define folder paths -----------------------------------------------------

s2r1_directory           <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")
s2r4_directory           <- file.path(experiments_directory, "2021-12-07 - fourth Sequel-II run")
s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")

s2r1_R_objects_directory <- file.path(s2r1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")
s2r4_R_objects_directory <- file.path(s2r4_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "30) Calculate plate weightings for the next round of sequencing.RData"))

load(file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData"))

load(file.path(s2r1_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
run1_alignments_df <- alignments_df

load(file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - alignments_df.RData"))
run2_alignments_df <- alignments_df

load(file.path(s2r3_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
run3_alignments_df <- alignments_df

load(file.path(s2r4_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
run4_alignments_df <- alignments_df

rm(alignments_df)

load(file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - ccs_df.RData"))
run2_ccs_df <- ccs_df

load(file.path(s2rI_R_objects_directory, "05) Read in PacBio data.RData"))





# Standardize the "Combined_ID" column ------------------------------------

run1_alignments_df[["Combined_ID"]] <- sub("Plate", "Plate0", run1_alignments_df[["Combined_ID"]], fixed = TRUE)

are_included <- run1_alignments_df[, "Combined_ID"] %in% library_df[, "Combined_ID"]
run1_alignments_df <- run1_alignments_df[are_included, ]
row.names(run1_alignments_df) <- NULL




# Assign original ZMWs ----------------------------------------------------

run1_alignments_df[["Original_ZMW"]] <- run1_alignments_df[["ZMW"]]
run3_alignments_df[["Original_ZMW"]] <- run3_alignments_df[["ZMW"]]

matches_vec <- match(run2_alignments_df[["ZMW"]], run2_ccs_df[["ZMW"]])
run2_alignments_df[["Original_ZMW"]] <- run2_ccs_df[["Original_ZMW"]][matches_vec]




# Split up the reads from run 2 into the two pools ------------------------

run2_alignments_df[["Pool"]] <- match(run2_ccs_df[matches_vec, "SmrtCell"],
                                      c("Sequel2_run2_pool1", "Sequel2_run2_pool2")
                                      )

R2P1_align_df <- run2_alignments_df[run2_alignments_df[["Pool"]] %in% 1, ]
R2P2_align_df <- run2_alignments_df[run2_alignments_df[["Pool"]] %in% 2, ]
row.names(R2P2_align_df) <- NULL
rm(run2_alignments_df)




# Assign the new ZMWs to the alignment data frames ------------------------

run1_new_align_df <- AddNewZMWs(ccs_df,
                                run1_alignments_df,
                                pool_number = 0L
                                )

R2P1_new_align_df <- AddNewZMWs(ccs_df,
                                R2P1_align_df,
                                pool_number = 1L
                                )

R2P2_new_align_df <- AddNewZMWs(ccs_df,
                                R2P2_align_df,
                                pool_number = 2L
                                )

run3_new_align_df <- AddNewZMWs(ccs_df,
                                run3_alignments_df,
                                pool_number = 3L
                                )




# Combine the alignment data frames from runs 1 to 3 ----------------------

runs1to3_alignments_df <- rbind.data.frame(
  run1_new_align_df,
  R2P1_new_align_df,
  R2P2_new_align_df,
  run3_new_align_df,
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)



# Filter out unneeded reads from runs 1 to 3 ------------------------------

replace_plates <- run4_df[run4_df[, "Religation"], "Plate_number"]

ID_splits <- strsplit(runs1to3_alignments_df[, "Combined_ID"], "_", fixed = TRUE)
plate_numbers <- as.integer(sub("Plate", "", sapply(ID_splits, "[[", 1)))
are_right_plates <- !(plate_numbers %in% replace_plates)
runs1to3_alignments_df <- runs1to3_alignments_df[are_right_plates, ]
row.names(runs1to3_alignments_df) <- NULL





# Filter out unneeded reads from run 4 ------------------------------------

ID_splits <- strsplit(run4_alignments_df[, "Combined_ID"], "_", fixed = TRUE)
plate_numbers <- as.integer(sub("Plate", "", sapply(ID_splits, "[[", 1)))
are_eligible <- (!(plate_numbers %in% c(46, 47, 71))) # These are the three plates that were wrongly re-sequenced
run4_alignments_df <- run4_alignments_df[are_eligible, ]
row.names(run4_alignments_df) <- NULL



# Add additional columns to run4_alignments_df ----------------------------

stopifnot(!(any(duplicated(ccs_df[, "ZMW"]))))
matches_vec <- match(run4_alignments_df[, "ZMW"], ccs_df[, "ZMW"])
for (column_name in c("Run", "Pool", "Original_ZMW")) {
  run4_alignments_df[[column_name]] <- ccs_df[[column_name]][matches_vec]
}




# Combine all alignment data frames ---------------------------------------

alignments_df <- rbind.data.frame(
  runs1to3_alignments_df,
  run4_alignments_df,
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)


stopifnot(!(any(duplicated(alignments_df[, "ZMW"]))))




# Filter out unneeded reads -----------------------------------------------

are_present <- alignments_df[, "ZMW"] %in% ccs_df[, "ZMW"]
alignments_df <- alignments_df[are_present, ]
row.names(alignments_df) <- NULL

stopifnot(all(ccs_df[, "ZMW"] %in% alignments_df[, "ZMW"]))




# Save data ---------------------------------------------------------------

save(list = "alignments_df",
     file = file.path(s2rI_R_objects_directory,
                      "06) Perform pairwise alignments with the reference sequence.RData"
                      )
     )




