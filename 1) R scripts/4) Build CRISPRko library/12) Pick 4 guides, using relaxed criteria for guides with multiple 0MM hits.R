### 25th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R")) # For the 'human_genes_GRanges' object
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "10) Find overlaps between CRISPRko sgRNA sequences and genetic polymorphisms.RData"))





# Assign locations for guides with multiple 0MM hits ----------------------

lax_CRISPRko_df <- FindBest0MMLocations(merged_CRISPRko_df)





# Re-order and rename columns ---------------------------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End")
columns_vec <- colnames(lax_CRISPRko_df)

are_strict_columns <- columns_vec %in% location_columns
are_relaxed_columns <- columns_vec %in% paste0(location_columns, "_lax")

columns_vec[are_strict_columns] <- paste0(location_columns, "_lax")
columns_vec[are_relaxed_columns] <- location_columns

lax_CRISPRko_df <- lax_CRISPRko_df[, columns_vec]

are_strict_columns <- columns_vec %in% location_columns
are_relaxed_columns <- columns_vec %in% paste0(location_columns, "_lax")

colnames(lax_CRISPRko_df)[are_strict_columns] <- paste0(location_columns, "_strict")
colnames(lax_CRISPRko_df)[are_relaxed_columns] <- location_columns

colnames(lax_CRISPRko_df) <- sub("_lax$", "", colnames(lax_CRISPRko_df))






# Re-assign cut locations -------------------------------------------------

lax_CRISPRko_df[["Cut_location"]] <- GetCutLocations(lax_CRISPRko_df)






# Rank sgRNAs -------------------------------------------------------------

lax_CRISPRko_df <- RankCRISPRDf(lax_CRISPRko_df, ID_column = "Combined_ID")






# Find non-overlapping guides, based on relaxed location criteria ---------

lax_CRISPRko_df <- PrioritizeNonOverlapping(lax_CRISPRko_df,
                                            ID_column = "Combined_ID",
                                            parallel_mode = TRUE,
                                            tolerate_divergent_chromosomes = TRUE
                                            )






# Save data ---------------------------------------------------------------

save(list = "lax_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory,
                      "12) Pick 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"
                      )
     )







