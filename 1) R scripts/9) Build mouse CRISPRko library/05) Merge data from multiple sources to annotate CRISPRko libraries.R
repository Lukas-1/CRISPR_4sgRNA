### 21st February 2021 ###




# Import packages and source code -----------------------------------------

library("TxDb.Mmusculus.UCSC.mm10.knownGene")

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "6) Mouse - General")
CRISPRko_RData_directory <- file.path(RData_directory, "8) Mouse - CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "01) Extract gene annotation data from the org.Mm.eg.db Bioconductor database.RData"))
load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData"))
load(file.path(CRISPRko_RData_directory, "02) Map CRISPRko sgRNA sequences to the mouse genome - all_sequences_df.RData"))
load(file.path(CRISPRko_RData_directory, "02) Map CRISPRko sgRNA sequences to the mouse genome - genome_search_df.RData"))




# Load global variables ---------------------------------------------------

mouse_genes_mm10_GRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
rm(human_genes_GRanges)




# Include the results of a genome search for the sgRNA sequences ----------

extended_CRISPRko_df <- ExtendWithGenomeSearch(CRISPRko_df, genome_search_df)





# Rename some of the columns ----------------------------------------------

rename_columns_vec <- c(
  "Hits_chromosome" = "Chromosome",
  "Hits_strand"     = "Strand",
  "Hits_start"      = "Start",
  "Hits_end"        = "End"
)

for (column_name in names(rename_columns_vec)) {
  names(extended_CRISPRko_df)[names(extended_CRISPRko_df) == column_name] <- rename_columns_vec[[column_name]]
}




# Assign a single location to sgRNAs targeting duplicated genes -----------

extended_CRISPRko_df[["Entrez_chromosome"]] <- EntrezIDsToChromosomes(extended_CRISPRko_df[["Entrez_ID"]])

extended_CRISPRko_df <- FindBest0MMLocations(extended_CRISPRko_df, mouse_genes_mm10_GRanges)





# Examine ambiguous Entrez IDs --------------------------------------------

are_ambiguous <- grepl(",", extended_CRISPRko_df[["Entrez_ID"]], fixed = TRUE)

if (any(are_ambiguous)) {
  display_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "Source",
                       "sgRNA_sequence", "Chromosome", "Entrez_overlapping_0MM", "Symbol_overlapping_0MM"
                       )
  ambiguous_df <- extended_CRISPRko_df[are_ambiguous, ]
  ambiguous_df[, display_columns]
}




# Assign missing Entrez IDs using the genomic locations of sgRNAs ---------

are_NA <- is.na(extended_CRISPRko_df[["Entrez_ID"]])
are_mapped <- !(is.na(extended_CRISPRko_df[["Start"]]))
are_to_replace <- are_mapped & are_NA & (extended_CRISPRko_df[["Is_control"]] == "No")

if (any(are_to_replace)) {
  reassigned_df <- extended_CRISPRko_df[are_to_replace, ]

  replaced_entrezs <- rep(NA_character_, nrow(reassigned_df))
  overlapping_entrezs_list <- strsplit(extended_CRISPRko_df[["Entrez_overlapping_0MM"]][are_to_replace], ", ", fixed = TRUE)
  are_replaceable <- !(is.na(overlapping_entrezs_list)) & (lengths(overlapping_entrezs_list) == 1)

  for (combined_ID in unique(extended_CRISPRko_df[["Combined_ID"]][are_to_replace][are_replaceable])) {
    are_this_ID <- extended_CRISPRko_df[["Combined_ID"]][are_to_replace] == combined_ID
    this_list <- unique(overlapping_entrezs_list[are_this_ID])
    this_list <- this_list[!(is.na(this_list))]
    if (length(this_list) == 1) {
      replaced_entrezs[are_this_ID] <- this_list[[1]]
    }
  }

  were_replaced <- !(is.na(replaced_entrezs))

  reassigned_df[["Entrez_ID"]][were_replaced] <- replaced_entrezs[were_replaced]
  reassigned_df[["Combined_ID"]][were_replaced] <- ifelse(is.na(reassigned_df[["Entrez_ID"]][were_replaced]),
                                                          toupper(reassigned_df[["Original_symbol"]][were_replaced]),
                                                          reassigned_df[["Entrez_ID"]][were_replaced]
                                                          )
  reassigned_df[["Gene_symbol"]][were_replaced] <- MapToEntrezs(entrez_IDs_vec = reassigned_df[["Entrez_ID"]][were_replaced],
                                                                is_mouse = TRUE
                                                                )[["Gene_symbol"]]

  for (column_name in c("Combined_ID", "Entrez_ID", "Gene_symbol")) {
    extended_CRISPRko_df[[column_name]][are_to_replace] <- reassigned_df[[column_name]]
  }
}




# Remove duplicate sgRNAs -------------------------------------------------

extended_CRISPRko_df <- ResolveDuplicates(extended_CRISPRko_df, concatenate_columns = c("Exon_number_GPP"))





# Remove locations for sgRNAs with discordant chromosome mappings ---------

extended_CRISPRko_df[["Entrez_chromosome"]] <- EntrezIDsToChromosomes(extended_CRISPRko_df[["Entrez_ID"]])

chromosome_vec <- extended_CRISPRko_df[["Entrez_chromosome"]]
chromosome_vec <- ifelse(chromosome_vec == "chrX, chrY", "chrX", chromosome_vec)

are_discordant <- !(is.na(extended_CRISPRko_df[["Chromosome"]])) &
                  !(is.na(chromosome_vec)) &
                  (extended_CRISPRko_df[["Chromosome"]] != chromosome_vec)
extended_CRISPRko_df[are_discordant, ]

for (column in c(location_columns, paste0(location_columns, "_strict"))) {
  extended_CRISPRko_df[[column]][are_discordant] <- NA
}




# Add the cut location ----------------------------------------------------

extended_CRISPRko_df[["Cut_location"]] <- GetCutLocations(extended_CRISPRko_df)




# Eliminate the "5 prime G nucleotide" column -----------------------------

extended_CRISPRko_df[["Num_1MM"]] <- rowSums(extended_CRISPRko_df[, c("Num_5G_MM", "Num_1MM")])
extended_CRISPRko_df <- extended_CRISPRko_df[, names(extended_CRISPRko_df) != "Num_5G_MM"]




# Check for remaining inconsistent chromosome mappings --------------------

inconsistent_IDs <- CheckForInconsistentChromosomes(extended_CRISPRko_df)
if (length(inconsistent_IDs) > 1) {
  message(paste0("The following combined IDs had sgRNAs that mapped to more than one chromosome: ",
                 paste0(inconsistent_IDs, collapse = ", "), "!"
                 )
          )
  selected_gene_columns <- c("Source", "Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "Chromosome") #"Symbol_overlapping_0MM")
  unique(extended_CRISPRko_df[extended_CRISPRko_df[["Combined_ID"]] %in% inconsistent_IDs, selected_gene_columns])
} else {
  message("No inconsistent mappings (sgRNAs from the same gene, but located on different chromosomes) were found!")
}




# Save data ---------------------------------------------------------------

save(list = "extended_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData")
     )



