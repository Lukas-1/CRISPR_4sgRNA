### 27th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rI_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(s2rI_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
load(file.path(s2rI_R_objects_directory, "07) Extract barcode sequences and quality scores.RData"))




# Create read filters -----------------------------------------------------

reads_df <- ccs_df
names(reads_df)[names(reads_df) == "Well_clip_start"] <- "Clip_start"
names(reads_df)[names(reads_df) == "Well_clip_end"] <- "Clip_end"
reads_df <- PassQualityFilters(reads_df, barcodes_df)

ccs_df[["Passed_filters"]] <- ccs_df[["Plate_passed_filters"]] &
                              (ccs_df[["Well_passed_filters"]] %in% TRUE)
ccs7_zmws <- GetCCS7_ZMWs(ccs_df)

reads_df[["Passes_CCS7"]] <- reads_df[, "ZMW"] %in% ccs7_zmws
ccs_df[["Passes_all_filters"]] <- reads_df[, "Passes_filters"] & reads_df[["Passes_CCS7"]]




# Filter out reads that do not meet barcode criteria ----------------------

originally_empty_wells <- table(library_df[, "Combined_ID"] %in% ccs_df[, "Combined_ID"])

ccs_df <- ccs_df[reads_df[, "Passes_barcode_filters"] == 1, ]
new_order <- order(match(ccs_df[, "Combined_ID"], library_df[, "Combined_ID"]))
ccs_df <- ccs_df[new_order, ]




# Select CCS7 reads for wells where hi-fi reads are available -------------

IDs_fac <- factor(ccs_df[, "Combined_ID"], levels = unique(ccs_df[, "Combined_ID"]))
have_hifi_splits <- split(ccs_df[, "Passes_all_filters"], IDs_fac)
have_hifi_read <- vapply(have_hifi_splits, any, logical(1))
have_hifi_long <- rep(have_hifi_read, lengths(have_hifi_splits))
are_selected <- ifelse(have_hifi_long, ccs_df[["Passes_all_filters"]], TRUE)
ccs_df <- ccs_df[are_selected, ]
row.names(ccs_df) <- NULL

now_empty_wells <- table(library_df[, "Combined_ID"] %in% ccs_df[, "Combined_ID"])




# Filter unused reads from alignments_df ----------------------------------

are_used <- alignments_df[, "ZMW"] %in% ccs_df[, "ZMW"]
alignments_df <- alignments_df[are_used, ]
new_order <- order(match(alignments_df[, "Combined_ID"], library_df[, "Combined_ID"]))
alignments_df <- alignments_df[new_order, ]
row.names(alignments_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(s2rI_R_objects_directory,
                      "07.5) Pre-filter reads - ccs_df.RData"
                      )
     )

save(list = "alignments_df",
     file = file.path(s2rI_R_objects_directory,
                      "07.5) Pre-filter reads - alignments_df.RData"
                      )
     )





