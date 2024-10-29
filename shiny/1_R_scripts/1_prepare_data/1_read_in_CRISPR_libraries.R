### 2024-09-22



# Define folder paths -----------------------------------------------------

project_dir   <- file.path("~", "CRISPR_4sgRNA", "shiny")
libraries_dir <- file.path(project_dir, "2_input", "our_CRISPR_libraries")
rdata_dir     <- file.path(project_dir, "3_RData")



# Read in data ------------------------------------------------------------

CRISPRa_sgRNA_df  <- read.delim(file.path(libraries_dir, "CRISPRa_all_sublibraries_ordered_by_well.tsv"),
                                quote = "", stringsAsFactors = FALSE, na.strings = c("NA", "")
                                )
CRISPRko_sgRNA_df <- read.delim(file.path(libraries_dir, "CRISPRko_all_sublibraries_ordered_by_well.tsv"),
                                quote = "", stringsAsFactors = FALSE, na.strings = c("NA", "")
                                )


# Save data ---------------------------------------------------------------

save(list = c("CRISPRa_sgRNA_df", "CRISPRko_sgRNA_df"),
     file = file.path(rdata_dir, "1_read_in_CRISPR_libraries.RData")
     )




