### 25th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "08) Checking for overlaps with genetic polymorphisms.R"))
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "05) Compile data on genetic polymorphisms.RData"))
load(file.path(CRISPRa_RData_directory, "18) Pick 4 guides per TSS.RData"))
strict_CRISPRa_df <- merged_replaced_CRISPRa_df
load(file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData"))





# Assign locations for guides with multiple 0MM hits ----------------------

lax_CRISPRa_df <- FindBest0MMLocations(merged_replaced_CRISPRa_df)





# Identify genes affected by relaxed locations ----------------------------

are_identical <- mapply(identical, lax_CRISPRa_df[["Start"]], lax_CRISPRa_df[["Start_lax"]])
changed_IDs <- lax_CRISPRa_df[["Combined_ID"]][!(are_identical)]





# Re-order and rename columns ---------------------------------------------

lax_CRISPRa_df <- ReorderLaxColumns(lax_CRISPRa_df)

drop_columns <- c("TSS_number", "Allocated_TSS", "Num_TSSs", "TSS_ID", "AltTSS_ID")
lax_CRISPRa_df <- lax_CRISPRa_df[, !(colnames(lax_CRISPRa_df) %in% drop_columns)]





# Re-assign cut locations -------------------------------------------------

lax_CRISPRa_df[["Cut_location"]] <- GetCutLocations(lax_CRISPRa_df)





# Find SNPs that overlap with sgRNAs, based on relaxed locations ----------

all_columns <- colnames(lax_CRISPRa_df)
lax_CRISPRa_df <- SNPDataForCRISPRdf(lax_CRISPRa_df[, grep("SNP", colnames(lax_CRISPRa_df), fixed = TRUE, invert = TRUE)])[, all_columns]





# Re-assign TSS based on locations assigned using relaxed criteria --------

lax_CRISPRa_df[["Distance_from_TSS"]] <- GetDistanceFromTSS(lax_CRISPRa_df)
lax_CRISPRa_df <- AllocateTSSforAllGenes(lax_CRISPRa_df, parallel_mode = TRUE, tolerate_divergent_chromosomes = TRUE)






# Rank sgRNAs -------------------------------------------------------------

lax_CRISPRa_df <- RankCRISPRDf(lax_CRISPRa_df, ID_column = "AltTSS_ID")






# Find non-overlapping guides, based on relaxed location criteria ---------

lax_CRISPRa_changed_df <- PrioritizeNonOverlapping(lax_CRISPRa_df[lax_CRISPRa_df[["Combined_ID"]] %in% changed_IDs, ],
                                                   ID_column = "AltTSS_ID", parallel_mode = TRUE,
                                                   tolerate_divergent_chromosomes = TRUE
                                                   )





# Merge the changed and unchanged sgRNAs ----------------------------------

position_columns <- c("Chromosome", "Strand", "Start", "End")
for (column_name in position_columns) {
  strict_CRISPRa_df[[paste0(column_name, "_strict")]] <- strict_CRISPRa_df[[column_name]]
}

lax_CRISPRa_df <- rbind.data.frame(
  lax_CRISPRa_changed_df,
  strict_CRISPRa_df[!(strict_CRISPRa_df[["Combined_ID"]] %in% changed_IDs), colnames(lax_CRISPRa_changed_df)],
  stringsAsFactors = FALSE,
  make.row.names = FALSE
)

lax_CRISPRa_df <- lax_CRISPRa_df[order(match(lax_CRISPRa_df[["Combined_ID"]], strict_CRISPRa_df[["Combined_ID"]])), ]
rownames(lax_CRISPRa_df) <- NULL






# Save data ---------------------------------------------------------------

save(list = "lax_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "19) Pick 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"
                      )
     )







