### 27th October 2020 ###




# Import packages and source code -----------------------------------------

library("Rsamtools")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory  <- file.path(file_directory, "2) Input")
raw_data_directory    <- file.path(file_input_directory, "Raw data", "384-well-1_2_B01")
R_objects_directory   <- file.path(file_directory, "3) R objects")

run_name <- "m54073_190912_185203"
subreads_bam_file     <- file.path(raw_data_directory, paste0(run_name, ".subreads.bam"))
subreads_stats_file   <- file.path(raw_data_directory, paste0(run_name, ".subreads.stats.csv"))




# Read in data ------------------------------------------------------------

subreads_stats_df <- read.csv(subreads_stats_file, stringsAsFactors = FALSE)

subreads_bam <- scanBam(subreads_bam_file,
                        param = ScanBamParam(what = c("qname", "seq"))
                        )[[1]]



# Explore data ------------------------------------------------------------

stopifnot(identical(subreads_stats_df[["id"]],
                    subreads_bam[["qname"]]
                    ))

table(subreads_stats_df[["passes"]])
table(subreads_stats_df[["qual"]])

exclude_columns <- c("passes", "qual", "mean.qual")
for (column_name in exclude_columns) {
  stopifnot(length(unique(subreads_stats_df[[column_name]])) == 1)
}

id_splits <- strsplit(subreads_stats_df[["id"]], "/", fixed = TRUE)

stopifnot(identical(as.integer(sapply(id_splits, "[[", 2)), subreads_stats_df[["hole"]]))
stopifnot(length(unique(sapply(id_splits, "[[", 1))) == 1)
stopifnot(identical(sapply(id_splits, "[[", 3),
                    paste0(subreads_stats_df[["start"]], "_", subreads_stats_df[["end"]])
                    )
          )
stopifnot(identical(subreads_stats_df[["length"]], subreads_stats_df[["end"]] - subreads_stats_df[["start"]]))




# Drop unnecessary columns ------------------------------------------------

subreads_stats_df <- subreads_stats_df[, !(names(subreads_stats_df) %in% c(exclude_columns, "id"))]

# table(subreads_stats_df[["hole"]] %in% sl7_ccs_df[["ZMW"]])




# Save data ---------------------------------------------------------------

save(list = c("subreads_bam", "subreads_stats_df"),
     file = file.path(R_objects_directory, "19) Import raw subreads.RData")
     )




