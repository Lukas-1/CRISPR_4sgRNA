### 22nd July 2019 ###


# Import packages and source code -----------------------------------------

library("org.Hs.eg.db")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Define gene annotation lookup lists -------------------------------------

my_x <- org.Hs.egSYMBOL
entrez_to_symbol_list <- as.list(my_x[mappedkeys(my_x)])

my_x <- org.Hs.egSYMBOL2EG
symbol_to_entrez_list <- as.list(my_x[mappedkeys(my_x)])

my_x <- org.Hs.egGENENAME
entrez_to_genename_list <- as.list(my_x[mappedkeys(my_x)])

my_x <- org.Hs.egCHR
entrez_to_chromosome_list <- as.list(my_x[mappedkeys(my_x)])

my_x <- org.Hs.egALIAS2EG
alias_to_entrez_list <- as.list(my_x[mappedkeys(my_x)])

my_x <- org.Hs.egENSEMBL2EG
ensembl_to_entrez_list <- as.list(my_x[mappedkeys(my_x)])

my_x <- org.Hs.egENSEMBL
entrez_to_ensembl_list <- as.list(my_x[mappedkeys(my_x)])

rm(my_x)





# Save data ---------------------------------------------------------------

save(list = grep("^(entrez_to_.+_list|.+_to_entrez_list)$", ls(), value = TRUE),
     file = file.path(general_RData_directory, "01) Extract gene annotation data from the org.Hs.eg.db Bioconductor database.RData")
     )


