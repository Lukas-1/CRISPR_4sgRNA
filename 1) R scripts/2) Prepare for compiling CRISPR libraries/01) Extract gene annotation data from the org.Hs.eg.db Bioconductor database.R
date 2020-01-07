### 22nd July 2019 ###


# Import packages and source code -----------------------------------------

library("org.Hs.eg.db")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")




# Define functions --------------------------------------------------------

BimapToList <- function(Bimap_object) {
  as.list(Bimap_object[mappedkeys(Bimap_object)])
}




# Define gene annotation lookup lists -------------------------------------

entrez_to_symbol_list     <- BimapToList(org.Hs.egSYMBOL)
entrez_to_genename_list   <- BimapToList(org.Hs.egGENENAME)
entrez_to_chromosome_list <- BimapToList(org.Hs.egCHR)

symbol_to_entrez_list     <- BimapToList(org.Hs.egSYMBOL2EG)
alias_to_entrez_list      <- BimapToList(org.Hs.egALIAS2EG)

ensembl_to_entrez_list    <- BimapToList(org.Hs.egENSEMBL2EG)
entrez_to_ensembl_list    <- BimapToList(org.Hs.egENSEMBL)




# Save data ---------------------------------------------------------------

save(list = grep("^(entrez_to_.+_list|.+_to_entrez_list)$", ls(), value = TRUE),
     file = file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData")
     )





