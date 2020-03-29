### 25th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R")) # For the 'human_genes_GRanges' object
source(file.path(general_functions_directory, "08) Checking for overlaps with genetic polymorphisms.R"))
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "05) Compile data on genetic polymorphisms.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))
strict_CRISPRko_df <- merged_CRISPRko_df
load(file.path(CRISPRko_RData_directory, "10) Find overlaps between CRISPRko sgRNA sequences and genetic polymorphisms.RData"))





# Assign locations for guides with multiple 0MM hits ----------------------

lax_CRISPRko_df <- FindBest0MMLocations(merged_CRISPRko_df)





# Identify genes affected by relaxed locations ----------------------------

are_identical <- mapply(identical, lax_CRISPRko_df[["Start"]], lax_CRISPRko_df[["Start_lax"]])
changed_IDs <- unique(lax_CRISPRko_df[["Combined_ID"]][!(are_identical)])





# Re-order and rename columns ---------------------------------------------

lax_CRISPRko_df <- ReorderLaxColumns(lax_CRISPRko_df)





# Re-assign cut locations -------------------------------------------------

lax_CRISPRko_df[["Cut_location"]] <- GetCutLocations(lax_CRISPRko_df)





# Find SNPs that overlap with sgRNAs, based on relaxed locations ----------

all_columns <- colnames(lax_CRISPRko_df)
lax_CRISPRko_df <- SNPDataForCRISPRdf(lax_CRISPRko_df[, grep("SNP", colnames(lax_CRISPRko_df), fixed = TRUE, invert = TRUE)])[, all_columns]





# Rank sgRNAs -------------------------------------------------------------

lax_CRISPRko_df <- RankCRISPRDf(lax_CRISPRko_df, ID_column = "Combined_ID")






# Find non-overlapping guides, based on relaxed location criteria ---------

lax_CRISPRko_changed_df <- PrioritizeNonOverlapping(lax_CRISPRko_df[lax_CRISPRko_df[["Combined_ID"]] %in% changed_IDs, ],
                                                    ID_column = "Combined_ID",
                                                    parallel_mode = TRUE,
                                                    tolerate_divergent_chromosomes = TRUE
                                                    )





# Merge the changed and unchanged sgRNAs ----------------------------------

position_columns <- c("Chromosome", "Strand", "Start", "End")
for (column_name in position_columns) {
  strict_CRISPRko_df[[paste0(column_name, "_strict")]] <- strict_CRISPRko_df[[column_name]]
}

lax_CRISPRko_df <- rbind.data.frame(
  lax_CRISPRko_changed_df,
  strict_CRISPRko_df[!(strict_CRISPRko_df[["Combined_ID"]] %in% changed_IDs), colnames(lax_CRISPRko_changed_df)],
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

lax_CRISPRko_df <- lax_CRISPRko_df[order(match(lax_CRISPRko_df[["Combined_ID"]], strict_CRISPRko_df[["Combined_ID"]])), ]
rownames(lax_CRISPRko_df) <- NULL







# Save data ---------------------------------------------------------------

save(list = "lax_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory,
                      "12) Pick 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"
                      )
     )







