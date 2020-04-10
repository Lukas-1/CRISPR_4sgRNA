### 9th April 2020 ###




# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRi_RData_directory <- file.path(RData_directory, "4) CRISPRi")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRi_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRi-v2 - CRISPRi_df.RData"))
load(file.path(CRISPRi_RData_directory, "07) Assign genomic locations to sgRNA sequences.RData"))





# Replace artificial 5' G nucleotides from the hCRISPRi-v2 database -------

replaced_merged_CRISPRi_df <- Exchange5PrimeG(merged_CRISPRi_df)




# Replace 5' G nucleotides in control sgRNAs from hCRISPRi-v2 -------------

are_Dolcetto <- grepl("Dolcetto", merged_CRISPRi_df[["Source"]], fixed = TRUE)
are_hCRISPRi_v2 <- grepl("hCRISPRi-v2", merged_CRISPRi_df[["Source"]], fixed = TRUE)
are_controls <- merged_CRISPRi_df[["Is_control"]] == "Yes"
are_to_replace  <- are_controls & are_hCRISPRi_v2

nucleotide_bag <- substr(merged_CRISPRi_df[["sgRNA_sequence"]][are_Dolcetto & !(are_controls)], 1, 1)

set.seed(1)
new_5p_nucleotides <- sample(nucleotide_bag, sum(are_to_replace))

new_5p_nucleotides <- ifelse(new_5p_nucleotides == "G", "G", tolower(new_5p_nucleotides)) # Indicate the replacement

replaced_merged_CRISPRi_df[["sgRNA_sequence"]][are_to_replace] <- paste0(new_5p_nucleotides,
                                                                         substr(merged_CRISPRi_df[["sgRNA_sequence"]][are_to_replace], 2, 20)
                                                                         )




# Remove new duplicates that result from replacing 5'G nucleotides --------

replaced_merged_CRISPRi_df <- ResolveDuplicates(replaced_merged_CRISPRi_df,
                                                concatenate_columns = c("Sublibrary", "hCRISPRi_v2_ID", "hCRISPRi_TSS_source")
                                                )





# Save data ---------------------------------------------------------------

save(list = "replaced_merged_CRISPRi_df",
     file = file.path(CRISPRi_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData")
     )







