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

merged_CRISPRko_df <- SNPDataForCRISPRdf(merged_CRISPRko_df)







# Save data ---------------------------------------------------------------

save(list = "merged_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "10) Find overlaps between CRISPRko sgRNA sequences and genetic polymorphisms.RData")
     )


















