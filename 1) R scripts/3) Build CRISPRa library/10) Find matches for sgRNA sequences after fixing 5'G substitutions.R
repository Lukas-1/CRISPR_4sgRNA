### 9th September 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "06) Find matches for sgRNA sequences in the human genome - all_sequences_df.RData"))
load(file.path(CRISPRa_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData"))





# Search the human genome for matches to sgRNAs ---------------------------

replaced_unique_sequences <- unique(toupper(replaced_merged_CRISPRa_df[["sgRNA_sequence"]]))
new_or_not_found_sequences <- setdiff(replaced_unique_sequences, all_sequences_df[["Reference"]])

rm(replaced_merged_CRISPRa_df)

# All replaced sequences should be from the hCRISPRa-v2 library and should be 20-mers.
# However, in the future, some of the "not found" sequences may _NOT_ be 20-mers, in which case, the code would have to be modified.
stopifnot(all(nchar(new_or_not_found_sequences) == 20))

new_sequences_df <- FindSequences(new_or_not_found_sequences)





# Extend sequence matches with additional data ----------------------------

# new_sequences_df <- FindNearestGenes(new_sequences_df)
new_sequences_df[["PAM"]] <- GetNGGPAM(new_sequences_df)






# Merge the results from the new and old previous genome searches ---------

are_still_present <- all_sequences_df[["Reference"]] %in% replaced_unique_sequences

replaced_all_sequences_df <- rbind.data.frame(new_sequences_df,
                                              all_sequences_df[are_still_present, ],
                                              stringsAsFactors = FALSE,
                                              make.row.names = FALSE
                                              )





# Summarize data on all matches for a given sequence ----------------------

rm(all_sequences_df)
rm(new_sequences_df)

replaced_genome_search_df <- SummarizeFoundSequencesDf(replaced_all_sequences_df, all_sequences = replaced_unique_sequences)





# Save data ---------------------------------------------------------------

save(list = "replaced_all_sequences_df",
     file = file.path(CRISPRa_RData_directory, "10) Find matches for sgRNA sequences after fixing 5'G substitutions - replaced_all_sequences_df.RData")
     )
save(list = "replaced_genome_search_df",
     file = file.path(CRISPRa_RData_directory, "10) Find matches for sgRNA sequences after fixing 5'G substitutions - replaced_genome_search_df.RData")
     )







