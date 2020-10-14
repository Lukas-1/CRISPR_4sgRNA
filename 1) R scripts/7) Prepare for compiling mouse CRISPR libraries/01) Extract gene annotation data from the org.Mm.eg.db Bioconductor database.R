### 6 October 2020 ###


# Import packages and source code -----------------------------------------

library("org.Mm.eg.db")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "6) Mouse - General")




# Define functions --------------------------------------------------------

BimapToList <- function(Bimap_object) {
  as.list(Bimap_object[mappedkeys(Bimap_object)])
}




# Define gene annotation lookup lists -------------------------------------

entrez_to_symbol_list     <- BimapToList(org.Mm.egSYMBOL)
entrez_to_genename_list   <- BimapToList(org.Mm.egGENENAME)
entrez_to_chromosome_list <- BimapToList(org.Mm.egCHR)

symbol_to_entrez_list     <- BimapToList(org.Mm.egSYMBOL2EG)
alias_to_entrez_list      <- BimapToList(org.Mm.egALIAS2EG)

ensembl_to_entrez_list    <- BimapToList(org.Mm.egENSEMBL2EG)
entrez_to_ensembl_list    <- BimapToList(org.Mm.egENSEMBL)

mgi_to_entrez_list        <- BimapToList(org.Mm.egMGI2EG)
names(mgi_to_entrez_list) <- sub("^MGI:MGI:", "MGI:", names(mgi_to_entrez_list))




# Save data ---------------------------------------------------------------

save(list = grep("^(entrez_to_.+_list|.+_to_entrez_list)$", ls(), value = TRUE),
     file = file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData")
     )





