### 11th April 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "05) Performing pairwise alignments with the reference sequence.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData"))
load(file.path(sql2_R_objects_directory, "05) Read in PacBio data.RData"))





# Filter reads ------------------------------------------------------------

are_eligible <- (ccs_df[["Well_exists"]] %in% TRUE) &
                (ccs_df[["Read_quality"]] > 0)

stopifnot(all(library_df[["Combined_ID"]] %in% ccs_df[["Combined_ID"]][are_eligible]))

# are_ccs3 <- are_eligible & (ccs_df[["Read_quality"]] > 0.99) & (ccs_df[["Num_full_passes"]] > 3)




# Extract barcodes --------------------------------------------------------

plate_numbers <- unique(library_df[["Plate_number"]])
chunk_size <- 3
num_chunks <- ceiling(length(plate_numbers) / chunk_size)
chunk_vec <- rep(seq_len(num_chunks), each = chunk_size)[seq_len(length(plate_numbers))]
chunk_list <- split(plate_numbers, chunk_vec)

RData_prefix <- "06) Perform pairwise alignments with the reference sequence"

previous_time <- Sys.time()

for (i in seq_len(num_chunks)) {

  message("Processing chunk #", i, " of ", num_chunks, "...")

  are_these_plates <- ccs_df[["Plate_number"]] %in% chunk_list[[i]]
  are_selected <- are_eligible & are_these_plates
  selected_IDs <- library_df[["Combined_ID"]][library_df[["Plate_number"]] %in% chunk_list[[i]]]

  chunk_list[[i]] <- ExtractAlignedSequences(ccs_df[are_selected, ],
                                             ID_column = "Combined_ID",
                                             unique_IDs = selected_IDs,
                                             parallel_mode = TRUE
                                             )
  chunk_number <- formatC(i, width = 2, flag = "0")
  file_name <- paste0(RData_prefix, " - chunk ", chunk_number, ".RData")
  object_name <- paste0("align_chunk", chunk_number, chunk_list[[i]])
  assign(object_name, chunk_list[[i]])
  message("Saving data...")
  save(list = object_name,
       file = file.path(R_objects_directory, file_name)
       )
  rm(object_name)
  gc()
  message(format(previous_time - Sys.time(), digits = 3), " have elapsed.")
  message("")
  previous_time <- Sys.time()
}

alignments_df <- data.table::rbindlist(chunk_list)
data.table::setDF(alignments_df)




# Save data ---------------------------------------------------------------

save(list = "alignments_df",
     file = file.path(R_objects_directory, paste0(file_name, ".RData"))
     )




