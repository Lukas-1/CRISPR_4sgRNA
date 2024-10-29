### 2024-09-22


# Define folder paths -----------------------------------------------------

project_dir <- file.path("~", "CRISPR_4sgRNA", "shiny")
anno_dir    <- file.path(project_dir, "2_input", "gene_annotation")
rdata_dir   <- file.path(project_dir, "3_RData")



# Read in data ------------------------------------------------------------

GO_df <- read.delim(file.path(anno_dir, "gene2go_2024-09-22.gz"),
                    quote = "", stringsAsFactors = FALSE, check.names = FALSE
                    )

reactome_df <- read.delim(file.path(anno_dir, "Reactome", "NCBI2Reactome_All_Levels_2024-09-22.txt"),
                          quote = "", stringsAsFactors = FALSE, header = FALSE
                          )



# Filter data -------------------------------------------------------------

GO_df <- GO_df[GO_df[, "#tax_id"] %in% "9606", names(GO_df) != "#tax_id"]
row.names(GO_df) <- NULL

names(reactome_df) <- c("NCBI_gene_ID", "Reactome_pathway_ID", "URL",
                        "Pathway_name", "Evidence_code", "Species"
                        )
reactome_df <- reactome_df[reactome_df[, "Species"] %in% "Homo sapiens", names(reactome_df) != "Species"]
row.names(reactome_df) <- NULL



# Save data ---------------------------------------------------------------

save(list = c("GO_df", "reactome_df"),
     file = file.path(rdata_dir, "3_filter_annotation_data.RData")
     )

