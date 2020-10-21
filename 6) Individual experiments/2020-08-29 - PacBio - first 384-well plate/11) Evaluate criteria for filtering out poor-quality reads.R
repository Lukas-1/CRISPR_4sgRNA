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






# Analyze just the barcodes data frame ------------------------------------

barcodes_df <- sl7_barcodes_df

table(sl7_barcodes_df[["Starts_with_column_barcode"]] &
      !(sl7_barcodes_df[["Ends_with_column_barcode"]])
      )
table(sl7_barcodes_df[["Ends_with_column_barcode"]] &
      !(sl7_barcodes_df[["Starts_with_column_barcode"]])
      )
table(sl7_barcodes_df[["Starts_with_row_barcode"]] &
      !(sl7_barcodes_df[["Ends_with_row_barcode"]])
      )
table(sl7_barcodes_df[["Ends_with_row_barcode"]] &
      !(sl7_barcodes_df[["Starts_with_row_barcode"]])
      )

table(sl7_barcodes_df[["Contains_column_barcode"]] &
      !(sl7_barcodes_df[["Ends_with_column_barcode"]])
      )
table(sl7_barcodes_df[["Contains_row_barcode"]] &
      !(sl7_barcodes_df[["Starts_with_row_barcode"]])
      )



# Define the "individual reads" data frame to be used ---------------------

use_df <- sl7_ccs5_df_list[["individual_reads_df"]]
use_df <- use_df[order(use_df[["Well_number"]]), ]
row.names(use_df) <- NULL





# Check various cutoffs ---------------------------------------------------


combined_score_cutoff <- 80
score_lead_cutoff <- 40

mean_quality_cutoff <- 85


table(use_df[["BC_combined_score"]] < 100, useNA = "ifany")

table(use_df[["BC_combined_score"]] < 80, useNA = "ifany")
table(use_df[["BC_score_lead"]] < 40, useNA = "ifany")
table(use_df[["Mean_quality"]] < 85, useNA = "ifany")


are_poor_barcodes <- (use_df[["BC_combined_score"]] < combined_score_cutoff) |
                     (use_df[["BC_score_lead"]] < score_lead_cutoff)

are_poor_quality <- use_df[["Mean_quality"]] < mean_quality_cutoff






# Analyze short barcodes --------------------------------------------------

table(use_df[["Row_bc_length"]])
table(use_df[["Column_bc_length"]])

short_barcode_cutoff <- 7

are_short_bc <- (use_df[["Row_bc_length"]] <= short_barcode_cutoff) |
                (use_df[["Column_bc_length"]] <= short_barcode_cutoff)


table(use_df[["BC_combined_score"]][are_short_bc])


use_df[are_short_bc & !(are_poor_barcodes), ]

goo



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













