### 24th September 2021 ###


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

s2r1_R_objects_directory <- file.path(s2r1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))

load(file.path(s2r1_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
run1_alignments_df <- alignments_df

load(file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - alignments_df.RData"))
run2_alignments_df <- alignments_df

load(file.path(s2r3_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
run3_alignments_df <- alignments_df

rm(alignments_df)

load(file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - ccs_df.RData"))
run2_ccs_df <- ccs_df

load(file.path(s2rC_R_objects_directory, "05) Read in PacBio data.RData"))




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



# Combine the alignment data frames ---------------------------------------

alignments_df <- rbind.data.frame(
  run1_new_align_df,
  R2P1_new_align_df,
  R2P2_new_align_df,
  run3_new_align_df,
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

stopifnot(!(anyDuplicated(alignments_df[, "ZMW"])))




# Save data ---------------------------------------------------------------

save(list = "alignments_df",
     file = file.path(s2rC_R_objects_directory,
                      "06) Perform pairwise alignments with the reference sequence.RData"
                      )
     )




