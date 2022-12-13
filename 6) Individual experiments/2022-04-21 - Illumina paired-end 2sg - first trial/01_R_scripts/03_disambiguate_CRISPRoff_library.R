## 2022-05-03


# Load packages and source code -------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts", "1) R functions")
experiments_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments")
project_dir                 <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir             <- file.path(project_dir, "01_R_scripts", "R_functions")
source(file.path(R_functions_dir, "03_disambiguating_and_annotating_guides.R"))
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))



# Define folder paths -----------------------------------------------------

rdata_dir <- file.path(project_dir, "03_R_objects")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_annotate_CRISPRoff_library.RData"))
load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))



# Search the human genome for matches to sgRNA sequences ------------------

unique_sequences <- unique(substr(toupper(c(CRISPRoff_df[["protospacer_A"]], CRISPRoff_df[["protospacer_B"]])), 2, 20))
stopifnot(all(nchar(unique_sequences) == 19))

sequences_df <- FindSequences(unique_sequences, max.mismatch = 0)



# Summarize sgRNA matches -------------------------------------------------

sequences_df[["PAM"]] <- GetNGGPAM(sequences_df)

genome_search_df <- SummarizeFoundSequencesDf(sequences_df, all_sequences = unique_sequences)

matches_A <- match(toupper(substr(CRISPRoff_df[, "protospacer_A"], 2, 20)), genome_search_df[, "Sequence"])
matches_B <- match(toupper(substr(CRISPRoff_df[, "protospacer_B"], 2, 20)), genome_search_df[, "Sequence"])

use_columns <- c("Num_0MM", "Locations_0MM")
loci_A_df <- genome_search_df[matches_A, use_columns]
row.names(loci_A_df) <- NULL
names(loci_A_df) <- paste0(names(loci_A_df), "_sg1")

loci_B_df <- genome_search_df[matches_B, use_columns]
row.names(loci_B_df) <- NULL
names(loci_B_df) <- paste0(names(loci_B_df), "_sg2")

CRISPRoff_df <- data.frame(CRISPRoff_df, loci_A_df, loci_B_df)



# Find all TSSs targeted by each sgRNA ------------------------------------

affected_A_df <- FindAffectedGenes(CRISPRoff_df, 1)
affected_B_df <- FindAffectedGenes(CRISPRoff_df, 2)

CRISPRoff_df <- data.frame(CRISPRoff_df, affected_A_df, affected_B_df)


are_ambiguous <- grepl(",", CRISPRoff_df[, "Gene_symbols"], fixed = TRUE)
CRISPRoff_df[, "Gene_symbol"] <- CRISPRoff_df[, "Gene_symbols"]
new_symbols <- Disambiguate_IDs(CRISPRoff_df[are_ambiguous, ],
                                "Gene_symbols", "Affected_gene_symbols_sg1",
                                "Affected_gene_symbols_sg2"
                                )
CRISPRoff_df[, "Gene_symbol"][are_ambiguous] <- new_symbols


are_ambiguous <- grepl(",", CRISPRoff_df[, "Entrez_IDs"], fixed = TRUE)
CRISPRoff_df[, "Entrez_ID"] <- CRISPRoff_df[, "Entrez_IDs"]
new_entrezs <- Disambiguate_IDs(CRISPRoff_df[are_ambiguous, ],
                                "Entrez_IDs", "Affected_Entrez_IDs_sg1",
                                "Affected_Entrez_IDs_sg2"
                                )
CRISPRoff_df[, "Entrez_ID"][are_ambiguous] <- new_entrezs



# Annotate with GC content ------------------------------------------------

CRISPRoff_df[, "Num_GC_sg1"] <- GetGCcontent(CRISPRoff_df[, "protospacer_A"])
CRISPRoff_df[, "Num_GC_sg2"] <- GetGCcontent(CRISPRoff_df[, "protospacer_B"])



# Annotate plasmids that share sgRNAs with other plasmids -----------------

sg1_vec <- toupper(CRISPRoff_df[, "protospacer_A"])
sg2_vec <- toupper(CRISPRoff_df[, "protospacer_B"])
num_occurrences_sg1 <- table(sg1_vec)[toupper(sg1_vec)]
num_occurrences_sg2 <- table(sg2_vec)[toupper(sg2_vec)]
CRISPRoff_df[, "Has_shared_sgRNA"] <- (num_occurrences_sg1 > 1) | (num_occurrences_sg2 > 1)



# Choose between multiple plasmids for some Entrez IDs --------------------

CRISPRoff_df <- ChooseBetweenMultiplePlasmids(CRISPRoff_df)



# Save data ---------------------------------------------------------------

save(list = "CRISPRoff_df",
     file = file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData")
     )



