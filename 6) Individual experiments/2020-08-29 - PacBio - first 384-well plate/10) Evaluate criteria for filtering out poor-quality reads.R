### 15th October 2020 ###



# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory   <- file.path(file_directory, "2) Input")
file_output_directory  <- file.path(file_directory, "5) Output")
R_objects_directory    <- file.path(file_directory, "3) R objects")

plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "9) Process demultiplexed PacBio reads.RData"))







# Check various cutoffs ---------------------------------------------------

use_df <- sl7_ccs3_df_list[["individual_reads_df"]]

combined_score_cutoff <- 80
score_lead_cutoff <- 40

mean_quality_cutoff <- 85


table(use_df[["BC_combined_score"]] < 80, useNA = "ifany")
table(use_df[["BC_score_lead"]] < 40, useNA = "ifany")
table(use_df[["Mean_quality"]] < 85, useNA = "ifany")


are_poor_barcodes <- (use_df[["BC_combined_score"]] < combined_score_cutoff) |
                     (use_df[["BC_score_lead"]] < score_lead_cutoff)

are_poor_quality <- use_df[["Mean_quality"]] < mean_quality_cutoff






# Try something -----------------------------------------------------------

ccs3_df <- sl7_ccs3_df_list[["individual_reads_df"]]
ccs5_df <- sl7_ccs5_df_list[["individual_reads_df"]]

matches_vec <- match(ccs5_df[["ZMW"]], ccs3_df[["ZMW"]])

replica_ccs5_df <- ccs3_df[matches_vec, ]
row.names(replica_ccs5_df) <- NULL

identical(replica_ccs5_df[, names(replica_ccs5_df) != "Random_distance"], ccs5_df[, names(replica_ccs5_df) != "Random_distance"])

are_identical <- vapply(seq_len(nrow(replica_ccs5_df)),
                        function(x) identical(replica_ccs5_df[x, ], ccs5_df[x, ]),
                        logical(1)
                        )

replica_ccs5_df[!(are_identical), ]




# Check for associations with contaminations ------------------------------

use_reads_df <- sl7_ccs3_df_list[["individual_reads_df"]]

bc_comb_cutoff <- 80
bc_lead_cutoff <- 40
close_well_range <- 3L
are_poor_barcodes <- (use_reads_df[["BC_combined_score"]] < bc_comb_cutoff) &
                     (use_reads_df[["BC_score_lead"]] < bc_lead_cutoff)

are_close_contams <- use_reads_df[["Random_distance"]] <= close_well_range

are_standard_lengths <- use_reads_df[["Length"]] %in% 2223:2229

are_contaminated <- use_reads_df[["Contam_guides"]] >= 1


fisher.test(table("are_contam" = are_contaminated,
                  "bad_bc"     = are_poor_barcodes
                  ))

fisher.test(table("are_close" = are_close_contams,
                  "bad_bc"    = are_poor_barcodes
                  ))
wilcox.test(use_reads_df[["Random_distance"]] ~
              are_poor_barcodes
            )



fisher.test(table("are_contam"      = are_contaminated,
                  "standard_length" = are_standard_lengths
                  ))

fisher.test(table("standard_length" = are_standard_lengths,
                  "bad_bc"          = are_poor_barcodes
                  ))













