### 4th April 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")
CRISPR_root_directory <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts", "1) R functions")
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))         # For exporting the library as a reference
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))      # For ExportPlates



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
file_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
file_input_directory     <- file.path(file_directory, "2) Input")
file_output_directory    <- file.path(file_directory, "5) Output", "Library reference")
R_objects_directory      <- file.path(file_directory, "3) R objects")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(library_RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(library_RData_directory, "3) CRISPRko")
CRISPRi_RData_directory  <- file.path(library_RData_directory, "4) CRISPRi")

controls_file <- file.path(file_input_directory, "Metadata", "Internal Control NGS plate.xlsx")



# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
rm(vac_by_gene_df)

load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates - PD genes.RData"))

CRISPRa_df <- full_4sg_by_well_df
PD_CRISPRa_df <- PD_4sg_by_well_df

CRISPRko_objects    <- load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_PD_objects <- load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates - PD genes.RData"))

CRISPRko_df <- full_4sg_by_well_df
PD_CRISPRko_df <- PD_4sg_by_well_df

rm(list = c(CRISPRko_objects, CRISPRko_PD_objects))
rm(list = c("CRISPRko_objects", "CRISPRko_PD_objects"))

load(file.path(R_objects_directory, "01) Process and export plate barcodes.RData"))




# Read in data ------------------------------------------------------------

controls_df <- as.data.frame(read_excel(controls_file, col_names = FALSE),
                             stringsAsFactors = FALSE, check.names = FALSE
                             )




# Assemble a data frame of sgRNA sequences --------------------------------

TidySequencesDf <- function(CRISPR_df) {
  results_df <- CRISPR_df[, c("Plate_string", "Well_number", "Rank", "sgRNA_sequence")]
  results_df[["Plate_string"]] <- sub("_tf", "_", results_df[["Plate_string"]], fixed = TRUE)
  results_df[["Plate_string"]] <- sub("_sg[1-4]$", "", results_df[["Plate_string"]])
  results_df[["Plate_string"]] <- toupper(results_df[["Plate_string"]])
  results_df[["Plate_string"]] <- sub("+", "plus", results_df[["Plate_string"]], fixed = TRUE)
  names(results_df) <- c("Plate_name", "Well_number", "Sg_number", "Sequence")
  return(results_df)
}

CRISPRa_seq_df     <- TidySequencesDf(CRISPRa_df)
CRISPRko_seq_df    <- TidySequencesDf(CRISPRko_df)
PD_CRISPRa_seq_df  <- TidySequencesDf(PD_CRISPRa_df)
PD_CRISPRko_seq_df <- TidySequencesDf(PD_CRISPRko_df)
CRISPRi_seq_df     <- TidySequencesDf(vac_by_well_df)
CRISPRi_seq_df[["Plate_name"]] <- sub("VACI_", "Vac-", CRISPRi_seq_df[["Plate_name"]], fixed = TRUE)

combined_df <- rbind.data.frame(CRISPRi_seq_df,
                                PD_CRISPRa_seq_df,
                                PD_CRISPRko_seq_df,
                                CRISPRa_seq_df,
                                CRISPRko_seq_df,
                                stringsAsFactors = FALSE,
                                make.row.names = FALSE
                                )

combined_df <- combined_df[combined_df[["Plate_name"]] %in% plates_df[["Plate_name"]], ]
combined_df[["Sequence"]] <- toupper(combined_df[["Sequence"]])




# Process colony-picked controls ------------------------------------------

controls_mat <- as.matrix(controls_df[15:(15 + 16 - 1), 3:25])
controls_mat <- cbind(controls_mat, rep(NA, 16))
colnames(controls_mat) <- NULL
wells_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
are_present <- !(is.na(controls_mat))
control_wells <- wells_mat[are_present]

are_controls <- (CRISPRko_seq_df[["Plate_name"]] %in% "HO_38") &
                (CRISPRko_seq_df[["Well_number"]] %in% control_wells)
control_plate_df <- CRISPRko_seq_df[are_controls, ]
control_plate_df[["Plate_name"]] <- "Intctrl"




# Add the bead-purified plates --------------------------------------------

beads_df <- combined_df[combined_df[["Plate_name"]] %in% c("HA_11", "HO_1"), ]
beads_df[["Plate_name"]] <- paste0(beads_df[["Plate_name"]], "-beads")

combined_df <- rbind.data.frame(combined_df,
                                control_plate_df,
                                beads_df,
                                stringsAsFactors = FALSE
                                )



# Re-format the data frame ------------------------------------------------

sequences_list <- lapply(1:4, function(x) combined_df[["Sequence"]][combined_df[["Sg_number"]] %in% x])
sequences_mat <- do.call(cbind, sequences_list)
colnames(sequences_mat) <- paste0("Sequence_sg", 1:4)

use_columns <- setdiff(names(combined_df), c("Sequence", "Sg_number"))

library_df <- data.frame("Combined_ID" = NA,
                         "Plate_number" = NA,
                         combined_df[combined_df[["Sg_number"]] %in% 1, use_columns],
                         sequences_mat,
                         stringsAsFactors = FALSE,
                         row.names = NULL
                         )

plate_matches <- match(library_df[["Plate_name"]], plates_df[["Plate_name"]])
plate_numbers <- plates_df[["Plate_number"]][plate_matches]
well_names <- paste0("Plate", formatC(plate_numbers, width = 2, flag = "0"), "_",
                     "Well", formatC(library_df[["Well_number"]], width = 3, flag = "0")
                     )
library_df[["Combined_ID"]] <- well_names
library_df[["Plate_number"]] <- plate_numbers





# Correct for the mix-up on plate 13 --------------------------------------

are_plate_13 <- library_df[["Plate_number"]] == 13

scheme_mat <- rbind(rep(c("A", "B"), times = 12), rep(c("C", "D"), times = 12))
full_scheme_mat <- do.call(rbind, lapply(1:8, function(x) scheme_mat))

scheme_vec <- unlist(lapply(1:16, function(x) full_scheme_mat[x, ]))

original_13_mat <- as.matrix(library_df[are_plate_13, paste0("Sequence_sg", 1:4)])
rownames(original_13_mat) <- NULL

new_13_mat <- matrix(nrow = 384, ncol = 4)
for (i in 1:4) {
  new_13_mat[scheme_vec == "B", ] <- original_13_mat[scheme_vec == "D", ]
  new_13_mat[scheme_vec == "D", ] <- original_13_mat[scheme_vec == "B", ]
  new_13_mat[scheme_vec == "A", ] <- original_13_mat[scheme_vec == "A", ]
  new_13_mat[scheme_vec == "C", ] <- original_13_mat[scheme_vec == "C", ]
}

for (i in 1:4) {
  library_df[[paste0("Sequence_sg", i)]][are_plate_13] <- new_13_mat[, i]
}





# Annotate duplicated sgRNAs ----------------------------------------------

all_sequences <- unlist(library_df[, paste0("Sequence_sg", 1:4)])
num_occurrences_mat <- do.call(cbind,
                               lapply(library_df[, paste0("Sequence_sg", 1:4)],
                                      function(x) vapply(x, function(y) sum(y == all_sequences), integer(1), USE.NAMES = FALSE)
                               ))
colnames(num_occurrences_mat) <- paste0("Num_occurrences_sg", 1:4)

library_df <- data.frame(library_df, num_occurrences_mat, stringsAsFactors = FALSE)





# Add data on (known) empty wells -----------------------------------------

empty_wells <- c("Plate08_Well024", "Plate10_Well024")

library_df[["Empty_well"]] <- library_df[["Combined_ID"]] %in% empty_wells






# Create a large combined data frame, for reference purposes --------------

df_names <- c("vac_by_well_df", "PD_CRISPRa_df", "PD_CRISPRko_df",
              "CRISPRa_df", "CRISPRko_df"
              )

CRISPRa_df[["Modality"]]     <- "CRISPRa"
CRISPRko_df[["Modality"]]    <- "CRISPRko"
PD_CRISPRa_df[["Modality"]]  <- "CRISPRa"
PD_CRISPRko_df[["Modality"]] <- "CRISPRko"
vac_by_well_df[["Modality"]] <- "CRISPRi"

PD_CRISPRa_df[["Sublibrary_4sg"]]  <- "PD-a"
PD_CRISPRko_df[["Sublibrary_4sg"]] <- "PD-o"
vac_by_well_df[["Sublibrary_4sg"]] <- "Vacuolation"

common_columns <- Reduce(intersect, lapply(df_names, function(x) names(get(x))))

df_list <- lapply(df_names, function(x) get(x)[, common_columns])

all_guides_df <- do.call(rbind.data.frame, c(df_list, stringsAsFactors = FALSE, make.row.names = FALSE))

export_columns <- c("Modality", export_columns)

ExportPlates(all_guides_df, sub_folder = "", file_name = "All_guides",
             no_modality = TRUE, add_primers = FALSE
             )




# Save data ---------------------------------------------------------------

save(list = "library_df",
     file = file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData")
     )




