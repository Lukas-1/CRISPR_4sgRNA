### 4th April 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")
CRISPR_root_directory <- "~/CRISPR"

general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts", "1) R functions")
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))         # For exporting the library as a reference
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))      # For ExportPlates

experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "19) Annotating duplicated gRNAs.R"))
source(file.path(s2r1_R_functions_directory, "03) Tidying CRISPR gRNA data frames.R"))




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

CRISPRa_df[["Modality"]]     <- "CRISPRa"
CRISPRko_df[["Modality"]]    <- "CRISPRko"
PD_CRISPRa_df[["Modality"]]  <- "CRISPRa"
PD_CRISPRko_df[["Modality"]] <- "CRISPRko"
vac_by_well_df[["Modality"]] <- "CRISPRi"


df_names <- c(
  "vac_by_well_df" = "CRISPRi_seq_df",
  "PD_CRISPRa_df"  = "PD_CRISPRa_seq_df",
  "PD_CRISPRko_df" = "PD_CRISPRko_seq_df",
  "CRISPRa_df"     = "CRISPRa_seq_df",
  "CRISPRko_df"    = "CRISPRko_seq_df"
)

for (df_name in names(df_names)) {
  assign(df_names[[df_name]], TidySequencesDf(get(df_name)))
}

CRISPRi_seq_df[["Plate_name"]] <- sub("VACI_", "Vac-", CRISPRi_seq_df[["Plate_name"]], fixed = TRUE)

combined_df <- do.call(rbind.data.frame,
                       c(lapply(df_names, get),
                         list(stringsAsFactors = FALSE, make.row.names = FALSE)
                         )
                       )

combined_df <- combined_df[combined_df[["Plate_name"]] %in% plates_df[["Plate_name"]], ]




# Process colony-picked controls ------------------------------------------

control_plate_df <- CreateControlsDf(controls_df, CRISPRko_seq_df)




# Add the bead-purified plates --------------------------------------------

beads_df <- combined_df[combined_df[["Plate_name"]] %in% c("HA_11", "HO_1"), ]
beads_df[["Plate_name"]] <- paste0(beads_df[["Plate_name"]], "-beads")

combined_df <- rbind.data.frame(combined_df,
                                control_plate_df,
                                beads_df,
                                stringsAsFactors = FALSE
                                )



# Re-format the data frame ------------------------------------------------

library_df <- MakeLibraryDf(plates_df, combined_df)





# Correct for the mix-up on plate 13 --------------------------------------

are_plate_13 <- library_df[["Plate_number"]] == 13

scheme_mat <- rbind(rep(c("A", "B"), times = 12), rep(c("C", "D"), times = 12))
full_scheme_mat <- do.call(rbind, lapply(1:8, function(x) scheme_mat))

scheme_vec <- unlist(lapply(1:16, function(x) full_scheme_mat[x, ]))

original_13_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)

new_13_mat <- matrix(nrow = 16, ncol = 24)
for (i in 1:4) {
  new_13_mat[full_scheme_mat == "B"] <- original_13_mat[full_scheme_mat == "D"]
  new_13_mat[full_scheme_mat == "D"] <- original_13_mat[full_scheme_mat == "B"]
  new_13_mat[full_scheme_mat == "A"] <- original_13_mat[full_scheme_mat == "A"]
  new_13_mat[full_scheme_mat == "C"] <- original_13_mat[full_scheme_mat == "C"]
}


are_preceding <- library_df[["Plate_number"]] < 13
are_following <- library_df[["Plate_number"]] > 13

new_13_df <- library_df[are_plate_13, ]
new_13_indices <- unlist(lapply(1:16, function(x) new_13_mat[x, ]))
new_13_df <- new_13_df[new_13_indices, ]

new_13_df[["Well_number"]] <- seq_len(384)
new_13_df[["Combined_ID"]] <- paste0("Plate13_Well",
                                     formatC(seq_len(384), width = 3, flag = "0")
                                     )

library_df <- rbind.data.frame(
  library_df[are_preceding, ],
  new_13_df,
  library_df[are_following, ],
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)



# Annotate duplicated sgRNAs ----------------------------------------------

library_df <- AddNumOccurrences(library_df)





# Add data on (known) empty wells -----------------------------------------

empty_wells <- c("Plate08_Well024", "Plate10_Well024")

library_df[["Empty_well"]] <- library_df[["Combined_ID"]] %in% empty_wells





# Create a large combined data frame, for reference purposes --------------

PD_CRISPRa_df[["Sublibrary_4sg"]]  <- "PD-a"
PD_CRISPRko_df[["Sublibrary_4sg"]] <- "PD-o"
vac_by_well_df[["Sublibrary_4sg"]] <- "Vacuolation"

common_columns <- Reduce(intersect, lapply(names(df_names), function(x) names(get(x))))

df_list <- lapply(names(df_names), function(x) get(x)[, common_columns])

all_guides_df <- do.call(rbind.data.frame, c(df_list, stringsAsFactors = FALSE, make.row.names = FALSE))




# Export the combined data frame ------------------------------------------

export_columns <- c("Modality", export_columns)

ExportPlates(all_guides_df, sub_folder = "", file_name = "All_guides",
             no_modality = TRUE, add_primers = FALSE
             )




# Save data ---------------------------------------------------------------

save(list = "library_df",
     file = file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData")
     )



