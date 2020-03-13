### 25th February 2020 ###







# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "10) Ranking sgRNAs.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))
source(file.path(general_functions_directory, "13) Separating sgRNAs that target distinct TSSs.R"))
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "18) Re-order the library to prioritize non-overlapping sgRNAs.RData"))





# Prepare for isolating problematic genes ---------------------------------

are_controls <- !(merged_replaced_CRISPRa_df[["Is_control"]] == "No")
combined_IDs_vec <- merged_replaced_CRISPRa_df[["Combined_ID"]][!(are_controls)]
combined_IDs_fac <- factor(combined_IDs_vec, levels = unique(combined_IDs_vec))

stopifnot(identical(length(unique(combined_IDs_fac)),
                    length(rle(as.integer(combined_IDs_fac))[["lengths"]])
                    )
          )




# Identify genes that are to be replaced ----------------------------------

are_spaced <- tapply(merged_replaced_CRISPRa_df[["Spacing"]][!(are_controls)],
                     combined_IDs_fac,
                     function(x) any(x %in% c(12L, 50L))
                     )

unspaced_gene_IDs <- names(which(!(are_spaced)))





# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "17) Integrate the output from CRISPOR.RData"))






# Assign locations for guides with multiple 0MM hits ----------------------

lax_CRISPRa_df <- merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Combined_ID"]] %in% unspaced_gene_IDs, ]

lax_CRISPRa_df <- FindBest0MMLocations(lax_CRISPRa_df)




# Re-order and rename columns ---------------------------------------------

location_columns <- c("Chromosome", "Strand", "Start", "End")
columns_vec <- colnames(lax_CRISPRa_df)

are_strict_columns <- columns_vec %in% location_columns
are_relaxed_columns <- columns_vec %in% paste0(location_columns, "_lax")

columns_vec[are_strict_columns] <- paste0(location_columns, "_lax")
columns_vec[are_relaxed_columns] <- location_columns

lax_CRISPRa_df <- lax_CRISPRa_df[, columns_vec]

are_strict_columns <- columns_vec %in% location_columns
are_relaxed_columns <- columns_vec %in% paste0(location_columns, "_lax")

colnames(lax_CRISPRa_df)[are_strict_columns] <- paste0(location_columns, "_strict")
colnames(lax_CRISPRa_df)[are_relaxed_columns] <- location_columns

drop_columns <- c("TSS_number", "Allocated_TSS", "Num_TSSs", "TSS_ID", "AltTSS_ID")

lax_CRISPRa_df <- lax_CRISPRa_df[, !(colnames(lax_CRISPRa_df) %in% drop_columns)]
colnames(lax_CRISPRa_df) <- sub("_lax$", "", colnames(lax_CRISPRa_df))






# Re-assign TSS based on locations assigned using relaxed criteria --------

lax_CRISPRa_df[["Cut_location"]] <- GetCutLocations(lax_CRISPRa_df)

lax_CRISPRa_df <- AllocateTSSforAllGenes(lax_CRISPRa_df, parallel_mode = TRUE, tolerate_divergent_chromosomes = TRUE)






# Rank sgRNAs -------------------------------------------------------------

lax_CRISPRa_df <- RankCRISPRDf(lax_CRISPRa_df, ID_column = "AltTSS_ID")






# Find non-overlapping guides, based on relaxed location criteria ---------

lax_CRISPRa_df <- PrioritizeNonOverlapping(lax_CRISPRa_df, ID_column = "AltTSS_ID", parallel_mode = FALSE,
                                           tolerate_divergent_chromosomes = TRUE
                                           )







# Save data ---------------------------------------------------------------

save(list = "lax_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory,
                      "18.5.2) Pick the top 4 guides, using relaxed criteria for guides with multiple 0MM hits.RData"
                      )
     )







