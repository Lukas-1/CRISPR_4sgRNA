### 6th June 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(R_functions_directory, "17) Characterizing deletions.R"))





# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")





# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p1_R_objects_directory, "04) Create reference sequences for each well - constant sequences.RData"))
load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))
load(file.path(sql2_R_objects_directory, "19) Identify and characterize large deletions.RData"))





# Define functions --------------------------------------------------------

MakeRandomDeletionDf <- function(range_spans, deletion_sizes, barcode_length = 27L) {

  stopifnot(identical(length(range_spans), length(deletion_sizes)))

  random_ranges <- vapply(range_spans, function(x) sample(seq_len(x), 1), integer(1))
  random_starts <- random_ranges + barcode_length

  random_df <- data.frame(
    "Deletion_start" = random_starts,
    "Deletion_end"   = random_starts + deletion_sizes
  )

  random_df <- AnnotateDeletions(random_df)
  return(random_df)
}



GetSimulatedMat <- function(deletion_sizes,
                            plasmid_length = 2225L,
                            barcode_length = 27L,
                            set_seed = TRUE,
                            num_simulations = 1000L
                            ) {

  range_spans <- plasmid_length - deletion_sizes

  if (set_seed) {
    set.seed(1)
  }

  random_del_df_list <- lapply(seq_len(num_simulations), function(x) {
    if (x %in% c(1, seq(20, num_simulations, by = 20))) {
      message(paste0("Working on iteration #", x, "..."))
    }
    MakeRandomDeletionDf(range_spans, deletion_sizes, barcode_length = barcode_length)
  })

  use_columns <- c("Span_tracrRNAs", "Span_sgRNAs", "Span_sg_cr", "Span_promoters")
  num_spanning_mat <- do.call(rbind, lapply(random_del_df_list, function(x) {
    colSums(as.matrix(x[, use_columns]))
  }))
  percent_spanning_mat <- num_spanning_mat / length(deletion_sizes)
  return(percent_spanning_mat)
}



# Prepare plasmid coordinates ---------------------------------------------

plate_bc <- unique(nchar(plates_df[["Barcode_sequence"]])) + 1L
well_bc <- unique(nchar(c(row_barcodes, column_barcodes)))
whole_plasmid <- nchar(plasmid_string)
barcode_length <- plate_bc + well_bc

deletion_size <- 20L




# Filter reads ------------------------------------------------------------

pass_filters <- ccs7_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
use_zmws <- ccs7_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]
are_selected <- deletions_df[["ZMW"]] %in% use_zmws

all_reads_del_df <- deletions_df[are_selected, ]
row.names(all_reads_del_df) <- NULL




# Select only one example read per well (avoid PCR duplicates) ------------

all_reads_del_df <- PrioritizeDeletions(all_reads_del_df)
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
     file = file.path(sql2_R_objects_directory, "20) Simulate random deletions (to compare with observed deletions).RData")
     )




