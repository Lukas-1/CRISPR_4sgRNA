### 21st September 2021 ###



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

load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))

load(file.path(s2r1_R_objects_directory, "05) Read in PacBio data.RData"))
run1_ccs_df <- ccs_df

load(file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - ccs_df.RData"))
run2_ccs_df <- ccs_df

load(file.path(s2r3_R_objects_directory, "05) Read in PacBio data.RData"))
run3_ccs_df <- ccs_df

rm(ccs_df)



# Assign run and pool numbers ---------------------------------------------

run1_ccs_df[["Run"]]  <- 1L
run1_ccs_df[["Pool"]] <- 0L

run2_ccs_df[["Run"]]  <- 2L
run2_ccs_df[["Pool"]] <- match(run2_ccs_df[, "SmrtCell"],
                               c("Sequel2_run2_pool1", "Sequel2_run2_pool2")
                               )

run3_ccs_df[["Run"]]  <- 3L
run3_ccs_df[["Pool"]] <- 3L



# Count the number of ZMWs per pool ---------------------------------------

run2_pools <- table(run2_ccs_df[["Pool"]])

num_reads_vec <- c(
  "Pool0" = nrow(run1_ccs_df),
  "Pool1" = run2_pools[[1]],
  "Pool2" = run2_pools[[2]],
  "Pool3" = nrow(run3_ccs_df)
)



# Assign new read identifiers ---------------------------------------------

run1_ccs_df[["Original_ZMW"]] <- run1_ccs_df[["ZMW"]]
run3_ccs_df[["Original_ZMW"]] <- run3_ccs_df[["ZMW"]]

num_chars_zmw <- nchar(as.character(max(num_reads_vec)))

add_to_zmws <- as.integer(paste0(c(9, 1:3),
                                 paste0(rep("0", num_chars_zmw), collapse = "")
                                 )
                          )

run1_ccs_df[["ZMW"]] <- add_to_zmws[[1]] + seq_len(num_reads_vec[[1]])
run3_ccs_df[["ZMW"]] <- add_to_zmws[[4]] + seq_len(num_reads_vec[[4]])
run2_ccs_df[["ZMW"]] <- ifelse(run2_ccs_df[["Pool"]] == 2,
                               run2_ccs_df[["ZMW"]] + add_to_zmws[[2]],
                               run2_ccs_df[["ZMW"]]
                               )



# Combine the data frames -------------------------------------------------

use_columns <- union(c("Run", "Pool", "ZMW", "Original_ZMW"),
                     names(run1_ccs_df)
                     )

ccs_df <- rbind.data.frame(
  run1_ccs_df[, use_columns],
  run2_ccs_df[, use_columns],
  run3_ccs_df[, use_columns],
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)



# Filter out unneeded reads -----------------------------------------------

table(ccs_df[["Well_exists"]], useNA = "ifany")
table(ccs_df[["Read_quality"]] > 0)

are_eligible <- (ccs_df[, "Well_exists"] %in% TRUE) &
                (ccs_df[, "Read_quality"] > 0) &
                (ccs_df[, "Plate_number"] %in% library_df[, "Plate_number"]) &
                (ccs_df[["Read_quality"]] > 0.999) &
                (ccs_df[["Num_full_passes"]] >= 5)

ccs_df <- ccs_df[are_eligible, ]
row.names(ccs_df) <- NULL




# Standardize the "Combined_ID" column ------------------------------------

IDs_vec <- paste0("Plate",
                  formatC(ccs_df[["Plate_number"]], width = 3, flag = "0"),
                  "_Well",
                  formatC(ccs_df[["Well_number"]], width = 3, flag = "0")
                  )
ccs_df[["Combined_ID"]] <- IDs_vec




# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(s2rC_R_objects_directory, "05) Read in PacBio data.RData")
     )



