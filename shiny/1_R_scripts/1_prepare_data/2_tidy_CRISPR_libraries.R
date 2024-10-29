### 2024-09-22



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
  plasmids_vec <- paste0(ifelse(is.na(input_df[, "Entrez_ID"]), input_df[, "Gene_symbol"], input_df[, "Entrez_ID"]),
                         "__", input_df[, "Plate_ID"],
                         "__", input_df[, "Well_number"]
                         )
  input_df[, "Plasmid_ID"] <- plasmids_vec
  return(input_df)
}



# Add plate and plasmid IDs -----------------------------------------------

CRISPRa_sgRNA_df <- AddPlasmidIDs(CRISPRa_sgRNA_df)
CRISPRko_sgRNA_df <- AddPlasmidIDs(CRISPRko_sgRNA_df)



# Save data ---------------------------------------------------------------

save(list = c("CRISPRa_sgRNA_df", "CRISPRko_sgRNA_df"),
     file = file.path(rdata_dir, "2_tidy_CRISPR_libraries.RData")
     )




