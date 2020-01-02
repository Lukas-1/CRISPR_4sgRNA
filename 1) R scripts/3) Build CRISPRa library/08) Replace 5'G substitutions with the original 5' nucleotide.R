### 9th September 2019 ###




# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
RData_directory         <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory <- file.path(RData_directory, "2) CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "02) Extract the original sequences for sgRNAs from hCRISPRa-v2 - CRISPRa_df.RData"))
load(file.path(CRISPRa_RData_directory, "07) Assign genomic locations to sgRNA sequences.RData"))






# Define functions --------------------------------------------------------

Exchange5PrimeG <- function(CRISPR_df) {
  are_5prime_G <- !(is.na(CRISPR_df[, "Start"])) &
                  (CRISPR_df[, "Num_0MM"] == 0) & (CRISPR_df[, "Num_5G_MM"] == 1) &
                  (grepl("hCRISPRa-v2", CRISPR_df[, "Source"], fixed = TRUE)) &
                  (CRISPR_df[, "Is_control"] != "Yes") &
                  (CRISPR_df[, "Exchanged_5pG"] %in% "No")
  GRanges_object <- GRanges(
    seqnames = CRISPR_df[are_5prime_G, "Chromosome"],
    ranges   = IRanges(start = CRISPR_df[are_5prime_G, "Start"], end = CRISPR_df[are_5prime_G, "End"]),
    strand   = CRISPR_df[are_5prime_G, "Strand"]
  )
  nucleotide_5p_vec <- substr(as.character(motifRG:::getSequence(GRanges_object, BSgenome.Hsapiens.UCSC.hg38)), 1, 1)
  CRISPR_df[, "Exchanged_5pG"] <- ifelse(are_5prime_G, "Yes", CRISPR_df[, "Exchanged_5pG"])
  CRISPR_df[are_5prime_G, "sgRNA_sequence"] <- paste0(nucleotide_5p_vec,
                                                      substr(CRISPR_df[are_5prime_G, "sgRNA_sequence"], 2, nchar(CRISPR_df[are_5prime_G, "sgRNA_sequence"]))
                                                      )
  return(CRISPR_df)
}








# Replace artificial 5' G nucleotides from the hCRISPRa-v2 database -------

replaced_merged_CRISPRa_df <- Exchange5PrimeG(merged_CRISPRa_df)





# Replace 5' G nucleotides in control sgRNAs from hCRISPRa-v2 -------------

are_Calabrese   <- grepl("Calabrese", merged_CRISPRa_df[, "Source"], fixed = TRUE)
are_hCRISPRa_v2 <- grepl("hCRISPRa-v2", merged_CRISPRa_df[, "Source"], fixed = TRUE)
are_controls    <- merged_CRISPRa_df[, "Is_control"] == "Yes"
are_to_replace  <- are_controls & are_hCRISPRa_v2

nucleotide_bag <- substr(merged_CRISPRa_df[are_Calabrese & !(are_controls), "sgRNA_sequence"], 1, 1)

set.seed(1)
new_5p_nucleotides <- sample(nucleotide_bag, sum(are_to_replace))

new_5p_nucleotides <- ifelse(new_5p_nucleotides == "G", "G", tolower(new_5p_nucleotides)) # Indicate the replacement

replaced_merged_CRISPRa_df[are_to_replace, "sgRNA_sequence"] <- paste0(new_5p_nucleotides,
                                                                       substr(merged_CRISPRa_df[are_to_replace, "sgRNA_sequence"], 2, 20)
                                                                       )





# Remove new duplicates that resulted from replacing 5' nucleotides -------

replaced_merged_CRISPRa_df <- ResolveDuplicates(replaced_merged_CRISPRa_df[, c(colnames(CRISPRa_df), "Exchanged_5pG")])






# Save data ---------------------------------------------------------------

save(list = "replaced_merged_CRISPRa_df",
     file = file.path(CRISPRa_RData_directory, "08) Replace 5'G substitutions with the original 5' nucleotide.RData")
     )







