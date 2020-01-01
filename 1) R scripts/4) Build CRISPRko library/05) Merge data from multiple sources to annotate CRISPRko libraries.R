### 30th October 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "11) Merging data from multiple sources to annotate CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
RData_directory             <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory     <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory    <- file.path(RData_directory, "3) CRISPRko")
file_output_directory       <- file.path(CRISPR_root_directory, "5) Output")
CRISPRko_output_directory   <- file.path(file_output_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRko_RData_directory, "03) Disambiguate gene IDs by mapping hg19 genomic coordinates to genes - CRISPRko_df.RData"))
load(file.path(CRISPRko_RData_directory, "04) Annotate mapped CRISPRko sequences with additional information.RData"))





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
  colnames(extended_CRISPRko_df)[colnames(extended_CRISPRko_df) == column_name] <- rename_columns_vec[[column_name]]
}




# Assign a single location to sgRNAs targeting duplicated genes -----------

extended_CRISPRko_df[, "Entrez_chromosome"] <- EntrezIDsToChromosomes(extended_CRISPRko_df[, "Entrez_ID"])

extended_CRISPRko_df <- ReplaceDuplicatedGenes(extended_CRISPRko_df)







# Examine ambiguous Entrez IDs --------------------------------------------

are_ambiguous <- grepl(",", extended_CRISPRko_df[, "Entrez_ID"], fixed = TRUE)

display_columns <- c("Entrez_ID", "Gene_symbol", "Original_symbol", "Source", "sgRNA_sequence", "Chromosome", "Entrez_overlapping_0MM", "Symbol_overlapping_0MM")

ambiguous_df <- extended_CRISPRko_df[are_ambiguous, ]
ambiguous_df[, display_columns] # There seems like there is nothing that can be done about these ambiguous Entrez IDs...






# Assign missing Entrez IDs using the genomic locations of sgRNAs ---------

are_NA         <- is.na(extended_CRISPRko_df[, "Entrez_ID"])
are_mapped     <- !(is.na(extended_CRISPRko_df[, "Start"]))
are_to_replace <- are_mapped & are_NA & (extended_CRISPRko_df[, "Is_control"] == "No")

reassigned_df <- extended_CRISPRko_df[are_to_replace, ]

replaced_entrezs <- rep(NA_character_, nrow(reassigned_df))
overlapping_entrezs_list <- strsplit(extended_CRISPRko_df[are_to_replace, "Entrez_overlapping_0MM"], ", ", fixed = TRUE)
are_replaceable <- !(is.na(overlapping_entrezs_list)) & (lengths(overlapping_entrezs_list) == 1)

for (combined_ID in unique(extended_CRISPRko_df[are_to_replace, "Combined_ID"][are_replaceable])) {
  are_this_ID <- extended_CRISPRko_df[are_to_replace, "Combined_ID"] == combined_ID
  this_list <- unique(overlapping_entrezs_list[are_this_ID])
  this_list <- this_list[!(is.na(this_list))]
  if (length(this_list) == 1) {
    replaced_entrezs[are_this_ID] <- this_list[[1]]
  }
}

were_replaced <- !(is.na(replaced_entrezs))

reassigned_df[were_replaced, "Entrez_ID"] <- replaced_entrezs[were_replaced]
reassigned_df[were_replaced, "Combined_ID"] <- ifelse(is.na(reassigned_df[were_replaced, "Entrez_ID"]),
                                                      toupper(reassigned_df[were_replaced, "Original_symbol"]),
                                                      reassigned_df[were_replaced, "Entrez_ID"]
                                                      )
reassigned_df[were_replaced, "Gene_symbol"] <- MapToEntrezs(entrez_IDs_vec = reassigned_df[were_replaced, "Entrez_ID"])[, "Gene_symbol"]


for (column_name in c("Combined_ID", "Entrez_ID", "Gene_symbol")) {
  extended_CRISPRko_df[are_to_replace, column_name] <- reassigned_df[, column_name]
}






# Remove duplicate sgRNAs -------------------------------------------------

CRISPRko_df <- ResolveDuplicates(CRISPRko_df, concatenate_columns = "TKOv3_ID")







# Add the cut location ----------------------------------------------------

extended_CRISPRko_df[, "Cut_location"] <- GetCutLocations(extended_CRISPRko_df)





# Remove locations for sgRNAs with discordant chromosome mappings ---------

extended_CRISPRko_df[, "Entrez_chromosome"] <- EntrezIDsToChromosomes(extended_CRISPRko_df[, "Entrez_ID"])

are_discordant <- !(is.na(extended_CRISPRko_df[, "Chromosome"])) &
                  !(is.na(extended_CRISPRko_df[, "Entrez_chromosome"])) &
                  (extended_CRISPRko_df[, "Chromosome"] != extended_CRISPRko_df[, "Entrez_chromosome"])
extended_CRISPRko_df[are_discordant, ]

location_columns <- c("Chromosome", "Strand", "Start", "End")
for (column in c(location_columns, "PAM")) {
  extended_CRISPRko_df[are_discordant, column] <- NA
}






# Eliminate the "5 prime G nucleotide" column -----------------------------

extended_CRISPRko_df[, "Num_1MM"] <- rowSums(extended_CRISPRko_df[, c("Num_5G_MM", "Num_1MM")])
extended_CRISPRko_df <- extended_CRISPRko_df[, colnames(extended_CRISPRko_df) != "Num_5G_MM"]







# Check for remaining inconsistent chromosome mappings --------------------

inconsistent_IDs <- CheckForInconsistentChromosomes(extended_CRISPRko_df)
if (length(inconsistent_IDs) > 1) {
  message(paste0("The following combined IDs had sgRNAs that mapped to more than one chromosome: ",
                 paste0(inconsistent_IDs, collapse = ", "), "!"
                 )
          )
  selected_gene_columns <- c("Source", "Combined_ID", "Entrez_ID", "Gene_symbol", "Original_symbol", "Chromosome", "Symbol_overlapping_0MM")
  unique(extended_CRISPRko_df[extended_CRISPRko_df[, "Combined_ID"] %in% inconsistent_IDs, selected_gene_columns])
} else {
  message("No inconsistent mappings (sgRNAs from the same gene, but located on different chromosomes) were found!")
}





# Save data ---------------------------------------------------------------

save(list = "extended_CRISPRko_df",
     file = file.path(CRISPRko_RData_directory, "05) Merge data from multiple sources to annotate CRISPRko libraries.RData")
     )







