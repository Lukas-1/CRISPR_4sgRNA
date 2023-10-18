### 11th March 2019 ###



# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")
HGNC_directory          <- file.path(CRISPR_input_directory, "Human genome", "HGNC")





# Define functions --------------------------------------------------------

Read_HGNC_Group <- function(file_name) {
  read.table(file.path(HGNC_directory, file_name),
             sep = "\t", quote = "", header = TRUE, row.names = NULL,
             check.names = FALSE, stringsAsFactors = FALSE
             )
}


Tidy_HGNC_df <- function(HGNC_df) {
  results_df <- data.frame(
    "Gene_type"       = HGNC_df[["Locus type"]],
    "Entrez_ID"       = HGNC_df[["NCBI Gene ID"]],
    "Gene_symbol"     = HGNC_df[["Approved symbol"]],
    "Ensembl_gene_ID" = HGNC_df[["Ensembl gene ID"]],
    "HGNC_ID"         = HGNC_df[["HGNC ID"]],
    "HGNC_gene_group" = HGNC_df[["Group name"]],
    stringsAsFactors  = FALSE
  )
  results_df <- results_df[results_df[["Gene_type"]] != "pseudogene", ]
  rownames(results_df) <- NULL
  return(results_df)
}






# Read in data ------------------------------------------------------------

GPCR_output_df <- Read_HGNC_Group("GPCR_group__2020_03_11.tsv")
ion_channel_output_df <- Read_HGNC_Group("Ion_channels_group__2020_03_11.tsv")





# Examine GPCR columns ----------------------------------------------------

table(GPCR_output_df[["Status"]],     useNA = "ifany")
table(GPCR_output_df[["Locus type"]], useNA = "ifany")
table(GPCR_output_df[["Group name"]], useNA = "ifany")

anyNA(GPCR_output_df[["NCBI Gene ID"]])





# Curate a GPCR data frame ------------------------------------------------

GPCR_df <- Tidy_HGNC_df(GPCR_output_df)
ion_channel_df <- Tidy_HGNC_df(ion_channel_output_df)





# Save data ---------------------------------------------------------------

save(list = c("GPCR_df", "ion_channel_df"),
     file = file.path(general_RData_directory, "11) Compile gene groups from HGNC.RData")
     )

















