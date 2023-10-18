### 26th July 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")
CRISPR_root_directory <- "~/CRISPR_4sgRNA"

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

s2r2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-07-24 - second Sequel-II run")
file_output_directory    <- file.path(s2r2_directory, "5) Output", "Library reference")
R_objects_directory      <- file.path(s2r2_directory, "3) R objects")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(library_RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(library_RData_directory, "3) CRISPRko")

controls_file <- file.path(s2r1_directory, "2) Input", "Metadata", "Internal Control NGS plate.xlsx")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))

CRISPRa_df <- full_4sg_by_well_df
CRISPRko_objects <- load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_df <- full_4sg_by_well_df

rm(list = CRISPRko_objects)
rm(list = "CRISPRko_objects")

load(file.path(R_objects_directory, "01) Process and export plate barcodes.RData"))





# Read in data ------------------------------------------------------------

controls_df <- as.data.frame(read_excel(controls_file, col_names = FALSE),
                             stringsAsFactors = FALSE, check.names = FALSE
                             )




# Assemble a data frame of sgRNA sequences --------------------------------

CRISPRa_df[["Modality"]]  <- "CRISPRa"
CRISPRko_df[["Modality"]] <- "CRISPRko"

CRISPRa_seq_df <- TidySequencesDf(CRISPRa_df)
CRISPRko_seq_df <- TidySequencesDf(CRISPRko_df)

combined_df <- rbind.data.frame(CRISPRa_seq_df,
                                CRISPRko_seq_df,
                                stringsAsFactors = FALSE,
                                make.row.names = FALSE
                                )

all_plate_names <- unlist(strsplit(plates_df[["Plate_name"]], " and ", fixed = TRUE), use.names = FALSE)

combined_df <- combined_df[combined_df[["Plate_name"]] %in% all_plate_names, ]




# Process plate HO-53 -----------------------------------------------------

plate384_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
plate384_rot180_mat <- abs(plate384_mat - 385L)

df_are_HO_53 <- combined_df[["Plate_name"]] %in% "HO_53"
HO_53_wells <- combined_df[["Well_number"]][df_are_HO_53]

are_HO_53_mat <- apply(plate384_rot180_mat, 2, function(x) x %in% HO_53_wells)

HO_53_map_mat <- cbind("old" = seq_len(sum(are_HO_53_mat)),
                       "new" = rev(t(plate384_mat)[t(are_HO_53_mat)])
                       )
HO_53_matches <- match(combined_df[["Well_number"]][df_are_HO_53],
                       HO_53_map_mat[, "old"]
                       )

combined_df[["Well_number"]][df_are_HO_53] <- HO_53_map_mat[HO_53_matches, "new"]

combined_df[["Plate_name"]][combined_df[["Plate_name"]] %in% c("HO_5", "HO_53")] <- "HO_5 and HO_53"

plate_matches <- match(combined_df[["Plate_name"]], plates_df[["Plate_name"]])
stopifnot(!(anyNA(plate_matches)))
new_order <- order(plate_matches, combined_df[["Well_number"]])
combined_df <- combined_df[new_order, ]
row.names(combined_df) <- NULL





# Process colony-picked controls ------------------------------------------

control_plate_df <- CreateControlsDf(controls_df, CRISPRko_seq_df)
pool1_controls_df <- control_plate_df
pool2_controls_df <- control_plate_df

pool1_controls_df[["Plate_name"]] <- "IntCtl_pool1"
pool2_controls_df[["Plate_name"]] <- "IntCtl_pool2"





# Add the colony-picked controls ------------------------------------------

combined_df <- rbind.data.frame(combined_df,
                                pool1_controls_df,
                                pool2_controls_df,
                                stringsAsFactors = FALSE
                                )



# Re-format the data frame ------------------------------------------------

library_df <- MakeLibraryDf(plates_df, combined_df)




# Annotate duplicated sgRNAs ----------------------------------------------

library_df <- AddNumOccurrences(library_df)





# Add data on (known) empty wells -----------------------------------------





# Create a large combined data frame, for reference purposes --------------

df_names <- c("CRISPRa_df", "CRISPRko_df")
common_columns <- Reduce(intersect, lapply(df_names, function(x) names(get(x))))
df_list <- lapply(df_names, function(x) get(x)[, common_columns])
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



