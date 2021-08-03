### 3rd August 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R")) # For the "features_list" object
source(file.path(R_functions_directory, "18) Characterizing deletions.R"))
source(file.path(R_functions_directory, "21) Simulating random deletions (for comparison with observed deletions).R"))



# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p1_R_objects_directory, "04) Create reference sequences for each well - constant sequences.RData"))
load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "10) Identify and characterize deletions.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Prepare plasmid coordinates ---------------------------------------------

plate_bc <- unique(nchar(plates_df[["Barcode_sequence"]])) + 1L
well_bc <- unique(nchar(c(row_barcodes, column_barcodes)))
whole_plasmid <- nchar(plasmid_string)
barcode_length <- plate_bc + well_bc

deletion_size <- 20L




# Select only one example read per well (avoid PCR duplicates) ------------

all_reads_del_df <- FilterDeletionsDf(deletions_df, ccs7_df_list)
one_read_del_df <- all_reads_del_df[all_reads_del_df[["Well_random_rank"]] %in% 1, ]
row.names(one_read_del_df) <- NULL




# Run simulations ---------------------------------------------------------

one_read_simul_mat  <- GetSimulatedMat(one_read_del_df[["Deletion_size"]],
                                       whole_plasmid, barcode_length,
                                       num_simulations = 10000
                                       )
all_reads_simul_mat <- GetSimulatedMat(all_reads_del_df[["Deletion_size"]],
                                       whole_plasmid, barcode_length,
                                       num_simulations = 10000
                                       )



# Save data ---------------------------------------------------------------

save(list = c("one_read_simul_mat", "all_reads_simul_mat",
              "one_read_del_df", "all_reads_del_df"
              ),
     file = file.path(s2r2_R_objects_directory, "21) Simulate random deletions (to compare with observed deletions).RData")
     )




