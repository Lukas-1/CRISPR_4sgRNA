### 22nd December 2019 ###



# Import packages and source code -----------------------------------------

library("LncFinder")
library("seqinr")



# Define folder paths -----------------------------------------------------

file_directory       <- "~/Documents/ViennaRNA"
file_input_directory <- file.path(file_directory, "1) Input")
R_objects_directory  <- file.path(file_directory, "2) R objects")




# Read in data ------------------------------------------------------------

sgRNAs_vec <- read.table(file.path(file_input_directory, "Single sequences for self-annealing prediction.txt"),
                         stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                         )[, 1]



# Calculate folding energies ----------------------------------------------

folding_list <- vector(mode = "list", length = length(sgRNAs_vec))

for (i in seq_along(sgRNAs_vec)) {
  folding_list[[i]] <- LncFinder::run_RNAfold(seqinr::as.SeqFastadna(sgRNAs_vec[[i]]), parallel.cores = 4)
}

minimum_free_energies_vec <- vapply(folding_list, function(x) as.numeric(x[3, 1]), numeric(1))




# Save data ---------------------------------------------------------------

save(list = "minimum_free_energies_vec",
     file = file.path(R_objects_directory, "ViennaRNA_minimum_free_energies.RData")
     )




















