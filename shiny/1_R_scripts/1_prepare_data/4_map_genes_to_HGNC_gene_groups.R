### 2024-09-22


# Define folder paths -----------------------------------------------------

project_dir <- file.path("~", "CRISPR_4sgRNA", "shiny")
anno_dir    <- file.path(project_dir, "2_input", "gene_annotation")
rdata_dir   <- file.path(project_dir, "3_RData")



# Read in data ------------------------------------------------------------

hgnc_df <- read.delim(file.path(anno_dir, "HGNC", "hgnc_complete_set_2024-08-23"),
                      quote = "", stringsAsFactors = FALSE
                      )


# Filter data -------------------------------------------------------------


