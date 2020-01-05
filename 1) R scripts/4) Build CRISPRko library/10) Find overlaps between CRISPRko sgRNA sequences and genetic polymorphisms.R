### 29th October 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "08) Checking for overlaps with genetic polymorphisms.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory    <- file.path(RData_directory, "3) CRISPRko")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "05) Compile data on genetic polymorphisms.RData"))
load(file.path(CRISPRko_RData_directory, "09) Integrate the output from CRISPOR.RData"))






# Search the human genome for matches to sgRNAs ---------------------------

are_mapped <- !(is.na(merged_CRISPRko_df[, "Start"]))

mapped_indices <- rep.int(NA_integer_, length(are_mapped))
mapped_indices[are_mapped] <- seq_len(sum(are_mapped))

sgRNA_polymorphisms_df <- AllPolymorphisms(merged_CRISPRko_df[are_mapped, ])

merged_CRISPRko_df <- data.frame(merged_CRISPRko_df,
                                 sgRNA_polymorphisms_df[mapped_indices, ],
                                 stringsAsFactors = FALSE,
                                 row.names = NULL
                                 )




# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "10) Find overlaps between CRISPRko sgRNA sequences and genetic polymorphisms.RData")
     )


















