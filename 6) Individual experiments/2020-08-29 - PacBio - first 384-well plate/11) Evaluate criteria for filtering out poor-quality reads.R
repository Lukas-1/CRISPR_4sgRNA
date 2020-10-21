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

load(file.path(R_objects_directory, "08) Extract barcode sequences and quality scores.RData"))
load(file.path(R_objects_directory, "10) Process demultiplexed PacBio reads.RData"))






# Define the "individual reads" data frame to be used ---------------------

use_df <- sl7_ccs5_df_list[["individual_reads_df"]]
use_df <- use_df[order(use_df[["Well_number"]]), ]
row.names(use_df) <- NULL





# Check various cutoffs ---------------------------------------------------

table(use_df[["BC_combined_score"]])
table(use_df[["BC_score_lead"]])


combined_score_cutoff <- 80
score_lead_cutoff <- 40
mean_quality_cutoff <- 85
close_well_range <- 3L

combined_score_low_cutoff <- 60
score_lead_low_cutoff <- 30

table(use_df[["BC_combined_score"]] < 100, useNA = "ifany")

table(use_df[["BC_combined_score"]] < 80, useNA = "ifany")
table(use_df[["BC_score_lead"]] < 40, useNA = "ifany")
table(use_df[["Mean_quality"]] < 85, useNA = "ifany")

are_poor_barcodes <- (use_df[["BC_combined_score"]] < combined_score_cutoff) |
                     (use_df[["BC_score_lead"]] < score_lead_cutoff)

are_very_poor_barcodes <- (use_df[["BC_combined_score"]] < combined_score_low_cutoff) |
                          (use_df[["BC_score_lead"]] < score_lead_low_cutoff)

are_poor_quality <- use_df[["Mean_quality"]] < mean_quality_cutoff




# Check for associations with contaminations ------------------------------

are_close_contams <- use_df[["Random_distance"]] <= close_well_range

are_standard_lengths <- use_df[["Length"]] %in% 2223:2229

are_contaminated <- use_df[["Contam_guides"]] >= 1


fisher.test(table("are_contam" = are_contaminated,
                  "bad_bc"     = are_poor_barcodes
                  ))

fisher.test(table("are_close" = are_close_contams,
                  "bad_bc"    = are_poor_barcodes
                  ))
wilcox.test(use_df[["Random_distance"]] ~ are_poor_barcodes)



fisher.test(table("are_contam"      = are_contaminated,
                  "standard_length" = are_standard_lengths
                  ))

fisher.test(table("standard_length" = are_standard_lengths,
                  "bad_bc"          = are_poor_barcodes
                  ))





# Evaluate my own barcode criteria ----------------------------------------

table(use_df[["Correct_barcodes"]])
table(use_df[["Starts_with_row_barcode"]] & use_df[["Ends_with_column_barcode"]])

table(use_df[["Starts_with_row_barcode"]])
table(use_df[["Ends_with_column_barcode"]])

table(use_df[["Correct_row_flank"]])
table(use_df[["Correct_column_flank"]])

table(use_df[["Row_bc_length"]])
table(use_df[["Column_bc_length"]])

table(use_df[["Row_bc_length"]] >= 8)
table(use_df[["Column_bc_length"]] >= 8)

have_correct_row <- use_df[["Starts_with_row_barcode"]] |
                    ((use_df[["Row_bc_length"]] >= 8) &
                      use_df[["Correct_row_flank"]]
                    )
have_correct_column <- use_df[["Ends_with_column_barcode"]] |
                       ((use_df[["Column_bc_length"]] >= 8) &
                         use_df[["Correct_column_flank"]]
                        )


are_correct <- have_correct_row & have_correct_column

table(are_correct)
table(are_correct & !(are_very_poor_barcodes))
table(are_correct & !(are_poor_barcodes))

only_8_bp <- (use_df[["Column_bc_length"]] <= 8) | (use_df[["Row_bc_length"]] <= 8)

use_df[are_correct & are_poor_barcodes, ]
use_df[are_correct & are_poor_barcodes & only_8_bp, ]


fisher.test(table("are_contam" = are_contaminated[are_correct],
                  "bad_bc"     = are_poor_barcodes[are_correct]
                  ))

fisher.test(table("are_contam" = are_contaminated[are_correct & !(are_very_poor_barcodes)],
                  "bad_bc"     = are_poor_barcodes[are_correct & !(are_very_poor_barcodes)]
                  ))









# Analyze short barcodes --------------------------------------------------

table(use_df[["Row_bc_length"]])
table(use_df[["Column_bc_length"]])

short_barcode_cutoff <- 7

are_short_bc <- (use_df[["Row_bc_length"]] <= short_barcode_cutoff) |
                (use_df[["Column_bc_length"]] <= short_barcode_cutoff)


table(use_df[["BC_combined_score"]][are_short_bc])


use_df[are_short_bc & !(are_poor_barcodes), ]
















