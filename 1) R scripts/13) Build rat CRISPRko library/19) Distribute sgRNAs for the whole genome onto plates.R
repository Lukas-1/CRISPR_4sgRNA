### 14th April 2020 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "10) Rat - General")
CRISPRko_RData_directory <- file.path(RData_directory, "12) Rat - CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Rat - CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Read in gene lists.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Create an empty SNP column ----------------------------------------------

merged_CRISPRko_df[[preferred_AF_max_column]] <- NA_real_






# Assign pharmacology sgRNAs to plates ------------------------------------

are_pharma <- (merged_CRISPRko_df[["Combined_ID"]] %in% pharma_df[["Entrez_ID"]]) |
              ((merged_CRISPRko_df[["Is_control"]] %in% "Yes") & (merged_CRISPRko_df[["Num_0MM"]] == 0) & (merged_CRISPRko_df[["Num_1MM"]] == 0))
pharma_CRISPRa_df <- merged_CRISPRko_df[are_pharma, ]





# Add the control guides to the pharma plates -----------------------------

pharma_4sg_df <- AllocateAllGuidesToPlates(pharma_CRISPRa_df,
                                           list("Pharma" = pharma_CRISPRa_df[["Entrez_ID"]]),
                                           num_control_wells = 6,
                                           reorder_df = FALSE
                                           )


# Re-order the data frame according to the plate layout -------------------

pharma_4sg_df[["Old_order"]] <- seq_len(nrow(pharma_4sg_df))
pharma_4sg_reordered_df <- ReorderPlates(pharma_4sg_df)




# Move the control guides to the last plate -------------------------------

are_controls <- pharma_4sg_reordered_df[["Is_control"]] == "Yes"
last_plate_number <- max(pharma_4sg_reordered_df[["Plate_number"]])
max_well_last_plate <- max(pharma_4sg_reordered_df[["Well_number"]][(pharma_4sg_reordered_df[["Plate_number"]] == last_plate_number) & !(are_controls)])
control_wells_seq <- rep(seq(from = max_well_last_plate + 1L,
                             to = max_well_last_plate + sum(are_controls) / 4,
                             ),
                         each = 4
                         )
pharma_4sg_reordered_df[["Well_number"]][are_controls] <- control_wells_seq
pharma_4sg_reordered_df[["Plate_number"]][are_controls] <- last_plate_number
pharma_4sg_reordered_df <- AssignPlateStrings(pharma_4sg_reordered_df, use_prefix = "ph_")





# Restore the original plate order ----------------------------------------

pharma_4sg_df <- pharma_4sg_reordered_df[order(pharma_4sg_reordered_df[["Old_order"]]), ]
row.names(pharma_4sg_df) <- NULL
pharma_4sg_reordered_df <- pharma_4sg_reordered_df[, names(pharma_4sg_reordered_df) != "Old_order"]
pharma_4sg_df <- pharma_4sg_df[, names(pharma_4sg_df) != "Old_order"]





# Export the pharmacology gene plate layouts ------------------------------

pharma_4sg_df <- pharma_4sg_df[, names(pharma_4sg_df) != "Sublibrary_4sg"]
pharma_4sg_reordered_df <- pharma_4sg_reordered_df[, names(pharma_4sg_reordered_df) != "Sublibrary_4sg"]


omit_columns <- c("GuideScan_Num_2MM", "GuideScan_Num_3MM",
                  "GuideScan_specificity", "GuideScan_efficiency",
                  "GuideScan_offtarget_category",
                  preferred_AF_max_column,
                  "Original_symbol"
                  )

export_columns <- setdiff(export_columns, omit_columns)

ExportPlates(pharma_4sg_df, "Pharmacology_4sg_original_order", sub_folder = "Plate layout - pharma")
ExportPlates(pharma_4sg_reordered_df, "Pharmacology_4sg_reordered", sub_folder = "Plate layout - pharma")

for (i in 1:4) {
  use_df <- pharma_4sg_reordered_df[pharma_4sg_reordered_df[["Rank"]] %in% i, ]
  ExportPlates(use_df,
               paste0("Pharmacology_4sg_reordered_sg", i),
               sub_folder = "Plate layout - pharma",
               add_padding_between_plates = TRUE
               )
}




# Remove the empty SNP column ---------------------------------------------

pharma_4sg_df <- pharma_4sg_df[, names(pharma_4sg_df) != preferred_AF_max_column]
pharma_4sg_reordered_df <- pharma_4sg_reordered_df[, names(pharma_4sg_reordered_df) != preferred_AF_max_column]




# Save data ---------------------------------------------------------------

save(list = c("pharma_4sg_df", "pharma_4sg_df"),
     file = file.path(CRISPRko_RData_directory, "19) Distribute sgRNAs for the whole genome onto plates - PD genes.RData")
     )








