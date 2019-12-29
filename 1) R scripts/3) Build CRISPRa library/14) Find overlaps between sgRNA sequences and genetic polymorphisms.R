### 26 July 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "08) Checking for overlaps with genetic polymorphisms.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "05) Compile data on genetic polymorphisms.RData"))
load(file.path(CRISPRa_RData_directory, "13) Integrate the output from GuideScan for individual sgRNA locations.RData"))





# Search the human genome for matches to sgRNAs ---------------------------

are_mapped <- !(is.na(merged_replaced_CRISPRa_df[, "Start"]))

mapped_indices <- rep(NA_integer_, length(are_mapped))
mapped_indices[are_mapped] <- seq_len(sum(are_mapped))

sgRNA_polymorphisms_df <- AllPolymorphisms(merged_replaced_CRISPRa_df[are_mapped, ])

merged_replaced_CRISPRa_df <- data.frame(merged_replaced_CRISPRa_df,
                                         sgRNA_polymorphisms_df[mapped_indices, ],
                                         stringsAsFactors = FALSE,
                                         row.names = NULL
                                         )



# Save data ---------------------------------------------------------------

save(list = "merged_replaced_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "14) Find overlaps between sgRNA sequences and genetic polymorphisms.RData")
     )
























