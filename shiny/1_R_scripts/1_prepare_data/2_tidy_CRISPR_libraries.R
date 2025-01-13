### 2024-09-22


# Import packages and source code -----------------------------------------

root_dir <- "~/CRISPR_4sgRNA"
general_functions_directory <- file.path(root_dir, "1) R scripts/1) R functions")
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # for LocationStringToDf
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R")) # for GetCutLocations
source(file.path(root_dir, "6) Individual experiments",
                 "2022-01-04 - pick genes from the libraries",
                 "1) R functions", "01) Converting plate layouts.R"
                 ))



# Define folder paths -----------------------------------------------------

project_dir   <- file.path("~", "CRISPR_4sgRNA", "shiny")
libraries_dir <- file.path(project_dir, "2_input", "our_CRISPR_libraries")
rdata_dir     <- file.path(project_dir, "3_RData")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "1_read_in_CRISPR_libraries.RData"))



# Define functions --------------------------------------------------------

AddPlasmidIDs <- function(input_df) {
  input_df[, "Plate_ID"] <- sapply(strsplit(input_df[, "Plate_string"], "_"), "[[", 2)
  input_df[, "Plate_ID"] <- sub("tf", "", input_df[, "Plate_ID"])
  input_df[, "Plate_ID"] <- ifelse(input_df[, "Plate_ID"] == "5+", "05+",
                                   formatC(as.integer(sub("+", "", input_df[, "Plate_ID"], fixed = TRUE)), width = 2, flag = "0")
                                   )
  plasmids_vec <- paste0(ifelse(is.na(input_df[, "Entrez_ID"]), input_df[, "Gene_symbol"], input_df[, "Entrez_ID"]),
                         "__", input_df[, "Plate_ID"],
                         "__", input_df[, "Well_number"]
                         )
  input_df[, "Plasmid_ID"] <- plasmids_vec
  stopifnot(all(table(input_df[, "Plasmid_ID"]) == 4))
  return(input_df)
}



ReformatLocations <- function(char_vec) {
  are_truncated <- grepl("Truncated", char_vec, fixed = TRUE)
  truncated_list <- strsplit(char_vec[are_truncated], "... The first 3 are: ", fixed = TRUE)
  num_entries <- sapply(truncated_list, "[[", 1)
  have_location <- !(is.na(char_vec) | (char_vec %in% c("", " ")))
  locations_vec <- char_vec
  locations_vec[are_truncated] <- sapply(truncated_list, "[[", 2)
  locations_list <- strsplit(locations_vec[have_location], "; ", fixed = TRUE)
  locations_long_vec <- unlist(locations_list)
  locations_df <- LocationStringToDf(locations_long_vec)
  locations_df[, "Cut_location"] <- GetCutLocations(locations_df)
  locations_long_vec <- paste0(locations_df[, "Chromosome"], ":",
                               locations_df[, "Cut_location"], ":",
                               locations_df[, "Strand"]
                               )
  locations_new_list <- split(locations_long_vec, rep(seq_along(locations_list), lengths(locations_list)))
  stopifnot(identical(unname(lengths(locations_list)), unname(lengths(locations_new_list))))
  locations_new_vec <- vapply(locations_new_list, paste0, collapse = "; ", "")
  results_vec <- char_vec
  results_vec[have_location] <- locations_new_vec
  results_vec[are_truncated] <- paste0(num_entries, "... The first 3 are: ", results_vec[are_truncated])
  return(results_vec)
}



MakeAdjustments <- function(use_df) {

  use_df[, "Plasmid_ID"] <- as.integer(factor(use_df[, "Plasmid_ID"], levels = unique(use_df[, "Plasmid_ID"])))

  are_controls <- use_df[, "Sublibrary_4sg"] %in% c("Controls", "Misc / controls")
  use_df[are_controls, "Rank"] <- rep(1:4, times = sum(are_controls) / 4)

  if ("TSS_ID" %in% names(use_df)) {
    use_df[, "Source"][use_df[, "Source"] == "GPP, Cal, hC-v2"] <- "GPP, Cal, hC"
    use_df[, "Plate_ID"] <- paste0("ha", use_df[, "Plate_ID"])
  } else {
    use_df[, "Source"][use_df[, "Source"] == "GPP, Bru, tk3"] <- "GPP, Bru, TK3"
    use_df[, "Plate_ID"] <- paste0("ho", use_df[, "Plate_ID"])
  }

  use_df[, "all22_SNP_IDs_vcf"][use_df[, "all22_SNP_IDs_vcf"] == " "] <- NA
  use_df[, "Locations_0MM"][are_controls & (use_df[, "Locations_0MM"] == " ")] <- NA
  use_df[, "Locations_1MM"][use_df[, "Locations_1MM"] == " "] <- NA

  use_df[, "Locations_0MM"] <- ReformatLocations(use_df[, "Locations_0MM"])
  use_df[, "Locations_1MM"] <- ReformatLocations(use_df[, "Locations_1MM"])

  use_vec <- use_df[, "Locations_0MM"]
  use_vec[use_vec == " "] <- "?"
  have_location <- grepl(":", use_vec, fixed = TRUE)
  have_multiple <- grepl(";", use_vec, fixed = TRUE)
  use_vec[have_location & !(have_multiple)] <- NA
  use_df[, "Locations_0MM"] <- use_vec

  use_df[, "Coords_96wp"] <- ConvertWellNumbers(use_df[, "Well_number"])
  use_df[, "Coords_96wp"] <- sub("Plate", "plate", use_df[, "Coords_96wp"], fixed = TRUE)

  first_columns <- c("Sublibrary_4sg", "Plate_ID", "Well_number", "Coords_96wp")
  use_df <- use_df[, c(first_columns, setdiff(names(use_df), first_columns))]


  remove_columns <- c("Sequences_1MM", "GuideScan_offtarget_category",
                      "Transcript_ID", "Genomic_sequence_ID",
                      "Plate_string"
                      )
  use_df <- use_df[, !(names(use_df) %in% remove_columns)]


  return(use_df)
}



Summarize4Guides <- function(input_df) {

  ### Create a new data frame with one row per CRISPR library plasmid

  plasmids_fac <- factor(input_df[, "Plasmid_ID"], levels = unique(input_df[, "Plasmid_ID"]))
  are_plasmid_specific <- vapply(names(input_df), function(x) {
    message("Checking column '", x, "' whether it relates to the plasmid (rather than to sgRNAs)... ")
    all(tapply(input_df[, x], plasmids_fac, function(y) length(unique(y)) == 1))
  }, logical(1))


  use_columns <- names(input_df)[are_plasmid_specific]
  plasmids_df_list <- split(input_df[, use_columns], plasmids_fac)

  plasmids_df_list <- lapply(plasmids_df_list, function(x) {
    results_list <- as.list(x[1, ])
    return(results_list)
  })

  sg_sequences_df <- do.call(rbind.data.frame,
                             c(plasmids_df_list,
                               stringsAsFactors = FALSE,
                               make.row.names = FALSE
                             ))

  return(sg_sequences_df)
}



# Add plate and plasmid IDs -----------------------------------------------

CRISPRa_sgRNA_df <- MakeAdjustments(AddPlasmidIDs(CRISPRa_sgRNA_df))
CRISPRko_sgRNA_df <- MakeAdjustments(AddPlasmidIDs(CRISPRko_sgRNA_df))



# # Create data frames with one row per plasmid -----------------------------
#
# CRISPRa_df <- Summarize4Guides(CRISPRa_sgRNA_df)
# CRISPRko_df <- Summarize4Guides(CRISPRko_sgRNA_df)
#
#
# # Try stuff ---------------------------------------------------------------
#
# library("org.Hs.eg.db")
# BimapToList <- function(Bimap_object) {
#   as.list(Bimap_object[mappedkeys(Bimap_object)])
# }
# entrez_to_symbol_vec <- unlist(BimapToList(org.Hs.egSYMBOL))
#
# new_symbols_vec <- entrez_to_symbol_vec[as.character(CRISPRa_df[, "Entrez_ID"])]
# are_different <- !(mapply(identical, new_symbols_vec, CRISPRa_df[, "Gene_symbol"]))
# data.frame(
#   old = CRISPRa_df[, "Gene_symbol"],
#   new = new_symbols_vec
# )[are_different, ]



# Save data ---------------------------------------------------------------

save(list = c("CRISPRa_sgRNA_df", "CRISPRko_sgRNA_df"),
     file = file.path(rdata_dir, "2_tidy_CRISPR_libraries.RData")
     )



