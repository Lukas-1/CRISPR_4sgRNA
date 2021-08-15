### 9th August 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(R_functions_directory, "17) Characterizing contaminations (using aligned reads).R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r2_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(s2r2_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))
load(file.path(s2r2_R_objects_directory, "07) Extract barcode sequences and quality scores.RData"))
load(file.path(s2r2_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(s2r2_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData"))




# Define functions --------------------------------------------------------

ReassignZMWs79to82 <- function(input_df, assignments_df, new_zmws_mat) {

  has_plate_number <- "Plate_number" %in% names(input_df)
  if (has_plate_number) {
    plate_numbers <- input_df[["Plate_number"]]
  } else {
    plate_numbers <- as.integer(substr(input_df[["Combined_ID"]], 6, 8))
  }

  are_preceding <- (plate_numbers < 79) %in% TRUE
  are_intervening <- ((plate_numbers > 79) &
                      (plate_numbers < 82)) %in% TRUE
  are_following <- (plate_numbers > 82) %in% TRUE
  are_plate79or82 <- plate_numbers %in% 82

  plate79_only_zmws <- assignments_df[assignments_df[, "Assigned_plate"] == "79", "ZMW"]
  both_zmws <- assignments_df[assignments_df[, "Assigned_plate"] == "both", "ZMW"]
  plate82_zmws <- assignments_df[assignments_df[, "Assigned_plate"] %in% c("82", "both"), "ZMW"]

  are_plate79_only <- input_df[, "ZMW"] %in% plate79_only_zmws
  are_both <- input_df[, "ZMW"] %in% both_zmws
  are_plate82 <- input_df[, "ZMW"] %in% plate82_zmws

  plate79_only_df <- input_df[are_plate79_only, ]

  plate79_both_df <- input_df[are_both, ]
  new_zmw_matches <- match(plate79_both_df[, "ZMW"], new_zmws_mat[, "Plate82_ZMW"])
  stopifnot(!(any(duplicated(new_zmw_matches))))
  plate79_both_df[, "ZMW"] <- new_zmws_mat[new_zmw_matches, "Plate79_ZMW"]

  if ("Combined_ID" %in% names(input_df)) {
    plate79_only_df[, "Combined_ID"] <- sub("Plate082", "Plate079", plate79_only_df[, "Combined_ID"], fixed = TRUE)
    plate79_both_df[, "Combined_ID"] <- sub("Plate082", "Plate079", plate79_both_df[, "Combined_ID"], fixed = TRUE)
  }
  if (has_plate_number) {
    plate79_only_df[, "Plate_number"] <- 79L
    plate79_both_df[, "Plate_number"] <- 79L
  }

  results_df <- rbind.data.frame(
    input_df[are_preceding, ],
    plate79_only_df,
    plate79_both_df,
    input_df[are_intervening, ],
    input_df[are_plate82, ],
    input_df[are_following, ],
    make.row.names = FALSE,
    stringsAsFactors = FALSE
  )

  stopifnot(!(any(duplicated(results_df[, "ZMW"]))))
  return(results_df)
}





# Re-assign the reads mapped to plate 82 ----------------------------------


## Select only reads from plate 82
are_82 <- ccs_df[["Plate_number"]] %in% 82

keep_columns <- c("ZMW", "Original_ZMW", "Combined_ID", "Plate_number",
                  "Well_number", "Well_exists", "Num_full_passes"
                  )

plate82_df <- ccs_df[are_82, keep_columns]
row.names(plate82_df) <- NULL



## Identify reads from plate 82 with correct guide RNAs

are_sg <- extracted_df[["Feature"]] %in% paste0("sg", 1:4)
are_correct_82 <- (extracted_df[["Plate_number"]] %in% 82) &
                   are_sg & extracted_df[["Is_correct"]]

correct_82_zmws <- unique(extracted_df[["ZMW"]][are_correct_82])

plate82_df[["Has_correct_guide_plate82"]] <- plate82_df[["ZMW"]] %in% extracted_df[["ZMW"]][are_correct_82]



## Identify reads from plate 79 with correct guide RNAs

contam_82_df <- contam_df[contam_df[["Reference_plate_number"]] %in% 82, ]
row.names(contam_82_df) <- NULL
contam_82_df[["Is_same_well"]] <- contam_82_df[["Contaminating_well_number"]] == contam_82_df[["Reference_well_number"]]

are_79 <- contam_82_df[["Contaminating_plate_number"]] == 79
are_correct_79 <- are_79 & contam_82_df[["Is_same_well"]]
correct_79_zmws <- contam_82_df[["ZMW"]][are_correct_79]
plate82_df[["Has_correct_guide_plate79"]] <- plate82_df[["ZMW"]] %in% correct_79_zmws



## Identify contaminated reads from plates 79 and 82

are_82 <- contam_82_df[["Contaminating_plate_number"]] == 82
contaminated_82_zmws <- contam_82_df[["ZMW"]][are_82]
plate82_df[["Is_contamination_from_plate82"]] <- plate82_df[["ZMW"]] %in% contaminated_82_zmws

contaminated_79_zmws <- contam_82_df[["ZMW"]][are_79 & !(contam_82_df[["Is_same_well"]])]
plate82_df[["Is_contamination_from_plate79"]] <- plate82_df[["ZMW"]] %in% contaminated_79_zmws



## Check for reads that fall outside the expected categories

stopifnot(!(any(plate82_df[["Has_correct_guide_plate79"]] & plate82_df[["Has_correct_guide_plate82"]])))
stopifnot(!(any(plate82_df[["Has_correct_guide_plate82"]] & plate82_df[["Is_contamination_from_plate79"]])))
stopifnot(!(any(plate82_df[["Has_correct_guide_plate79"]] & plate82_df[["Is_contamination_from_plate82"]])))



## Examine reads that seem to be a mix of guides from plates 79 and 82

table(plate82_df[["Is_contamination_from_plate79"]] & plate82_df[["Is_contamination_from_plate82"]])
mixed_zmws <- plate82_df[["ZMW"]][plate82_df[["Is_contamination_from_plate79"]] & plate82_df[["Is_contamination_from_plate82"]]]
contam_82_df[contam_82_df[["ZMW"]] %in% mixed_zmws, ]



## Assign reads to either plate 79, plate 82, or both

are_79 <- plate82_df[["Has_correct_guide_plate79"]] |
          (plate82_df[["Is_contamination_from_plate79"]] &
           !(plate82_df[["Is_contamination_from_plate82"]]))

are_82 <-  plate82_df[["Has_correct_guide_plate82"]] |
          (plate82_df[["Is_contamination_from_plate82"]] &
           !(plate82_df[["Is_contamination_from_plate79"]]))

are_both <- !(plate82_df[["Has_correct_guide_plate79"]] |
              plate82_df[["Has_correct_guide_plate82"]] |
              plate82_df[["Is_contamination_from_plate79"]] |
              plate82_df[["Is_contamination_from_plate82"]]
              ) | (
                plate82_df[["Is_contamination_from_plate79"]] &
                plate82_df[["Is_contamination_from_plate82"]]
              )
plate82_df[["Assigned_plate"]] <- ifelse(are_both, "both",
                                         ifelse(are_79,
                                                "79",
                                                ifelse(are_82, "82", "error")
                                                )
                                         )




# Assign new ZMWs to the reads assigned to both plates --------------------

new_zmws_mat <- cbind("Plate82_ZMW" = plate82_df[["ZMW"]][plate82_df[["Assigned_plate"]] == "both"])

max_zmw <- max(ccs_df[["ZMW"]])

updated_zmws_mat <- cbind(new_zmws_mat,
                          "Plate79_ZMW" = seq_len(nrow(new_zmws_mat)) + max_zmw
                          )




# Re-structure existing data to fit the new plate assignments -------------

new_alignments_df <- ReassignZMWs79to82(alignments_df,
                                        plate82_df,
                                        updated_zmws_mat
                                        )

new_ccs_df <- ReassignZMWs79to82(ccs_df,
                                 plate82_df,
                                 updated_zmws_mat
                                 )

new_barcodes_df <- ReassignZMWs79to82(barcodes_df,
                                      plate82_df,
                                      updated_zmws_mat
                                      )




# Repeat analyses using the new plate assignments -------------------------

new_extracted_df <- ExtractAlignedSequences(new_ccs_df,
                                            new_alignments_df,
                                            ID_column = "Combined_ID",
                                            unique_IDs = sg_sequences_df[["Combined_ID"]]
                                            )


new_contam_df <- CharacterizeContaminations(new_extracted_df, sg_sequences_df)




# Update object names -----------------------------------------------------

ccs_df        <- new_ccs_df
alignments_df <- new_alignments_df
barcodes_df   <- new_barcodes_df
extracted_df  <- new_extracted_df
contam_df     <- new_contam_df




# Save data ---------------------------------------------------------------

save(list = "ccs_df",
     file = file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - ccs_df.RData")
     )

save(list = "alignments_df",
     file = file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - alignments_df.RData")
     )

save(list = "barcodes_df",
     file = file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - barcodes_df.RData")
     )

save(list = "extracted_df",
     file = file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - extracted_df.RData")
     )

save(list = "contam_df",
     file = file.path(s2r2_R_objects_directory, "09.5) Deconvolve the plates with a barcoding error - contam_df.RData")
     )


