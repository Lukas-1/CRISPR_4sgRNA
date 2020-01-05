### 12 December 2019 ###





# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("liftOver")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "03) Compiling CRISPR libraries.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
liftOver_input_directory <- file.path(CRISPR_input_directory, "Human genome", "liftOver")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "02) Map gene symbols to Entrez IDs.RData"))
load(file.path(CRISPRko_RData_directory, "01) Compile predefined CRISPRko libraries.RData"))

hg19tohg38_chain <- import.chain(file.path(liftOver_input_directory, "hg19ToHg38.over.chain"))




# Loading global variables ------------------------------------------------

human_genes_hg19_GRanges <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)





# Process the genomic locations of sgRNAs in the TKOv3 library ------------

have_locations <- !(is.na(CRISPRko_df[, "TKOv3_ID"])) & (CRISPRko_df[, "Is_control"] == "No")

TKOv3_IDs_vec <- CRISPRko_df[have_locations, "TKOv3_ID"]

TKOv3_IDs_splits        <- strsplit(TKOv3_IDs_vec, "_", fixed = TRUE)
TKOv3_chromosome_splits <- strsplit(sapply(TKOv3_IDs_splits, "[[", 1), ":", fixed = TRUE)
TKOv3_range_splits      <- strsplit(sapply(TKOv3_chromosome_splits, "[[", 2), "-", fixed = TRUE)

TKOv3_locations_df <- data.frame(
  "Chromosome"     = sapply(TKOv3_chromosome_splits, "[[", 1),
  "Strand"         = sapply(TKOv3_IDs_splits, "[[", 3),
  "Start"          = as.integer(sapply(TKOv3_range_splits, "[[", 1)),
  "End"            = as.integer(sapply(TKOv3_range_splits, "[[", 2)),
  stringsAsFactors = FALSE,
  row.names        = NULL
)




# Annotate these locations (using a liftOver from hg19 to hg38) -----------

location_annotations_df <- LiftOverAndAnnotate(TKOv3_locations_df)
gene_overlaps_df <- FindOverlappingGenes(TKOv3_locations_df, gene_models_GRanges = human_genes_hg19_GRanges)

expanded_indices <- rep(NA_integer_, nrow(CRISPRko_df))
expanded_indices[have_locations] <- seq_len(sum(have_locations))

lifted_CRISPRko_df <- data.frame(
  CRISPRko_df,
  location_annotations_df[expanded_indices, ],
  gene_overlaps_df[expanded_indices, ],
  stringsAsFactors = FALSE,
  row.names = NULL
)



# Perform checks on the resulting data frame ------------------------------

selected_columns <- c("Source",
                      "Gene_symbol", "Overlapping_symbols", "Original_symbol",
                      "Entrez_ID", "Overlapping_Entrez_IDs",
                      "TKOv3_ID",
                      "sgRNA_sequence", "Sequence_hg19", "Sequence_liftOver",
                      "Num_overlapping_genes"
                      )
lifted_CRISPRko_df[, selected_columns]

num_sequences_hg19_hg38 <- vapply(seq_len(nrow(lifted_CRISPRko_df)), function(x) {
  sequences <- unique(as.character(lifted_CRISPRko_df[x, c("sgRNA_sequence", "Sequence_hg19", "Sequence_liftOver")]))
  sequences <- sequences[!(is.na(sequences))]
  num_sequences <- length(sequences)
  return(num_sequences)
}, integer(1))

num_sequences_hg19 <- vapply(seq_len(nrow(lifted_CRISPRko_df)), function(x) {
  sequences <- unique(as.character(lifted_CRISPRko_df[x, c("sgRNA_sequence", "Sequence_hg19")]))
  sequences <- sequences[!(is.na(sequences))]
  num_sequences <- length(sequences)
  return(num_sequences)
}, integer(1))

stopifnot(all(num_sequences_hg19 == 1))





# Disambiguate Entrez IDs (based on the location of sgRNAs) ---------------

# Perform checks
are_ambiguous <- grepl(", ", lifted_CRISPRko_df[, "Entrez_ID"], fixed = TRUE)
stopifnot(!(any(are_ambiguous & (lifted_CRISPRko_df[, "Source"] == "Brunello"))))
stopifnot(all(lifted_CRISPRko_df[have_locations, "sgRNA_sequence"] == lifted_CRISPRko_df[have_locations, "Sequence_hg19"]))

# Disambiguate Entrez IDs
disambiguated_df <- ReassignEntrezsByLocations(lifted_CRISPRko_df[have_locations, ])

# Check the results of disambiguation
table(disambiguated_df[, "Entrez_assignment"])
per_gene_columns <- c("Combined_ID", "Old_Entrez", "Entrez_ID", "Old_symbol", "Gene_symbol", "Original_symbol",
                      "Were_replaced", "Exon_number_Brunello", "Exon_number_TKOv3", "Chromosome"
                      )
unique(disambiguated_df[disambiguated_df[, "Entrez_assignment"] %in% "Unambiguous chromosome", per_gene_columns])
unique(disambiguated_df[disambiguated_df[, "Entrez_assignment"] %in% "Overlaps with gene", per_gene_columns])
unique(disambiguated_df[disambiguated_df[, "Entrez_assignment"] %in% "Ambiguous chromosome, and ambiguous overlaps", per_gene_columns])


for (column_name in c("Entrez_ID", "Combined_ID", "Gene_symbol")) {
  CRISPRko_df[have_locations, column_name] <- disambiguated_df[, column_name]
}



# Add the genomic locations based on liftOver to hCRISPRko_df -------------

are_identical_sequences <- (lifted_CRISPRko_df[have_locations, "sgRNA_sequence"] == lifted_CRISPRko_df[have_locations, "Sequence_liftOver"]) %in% TRUE

liftOver_columns <- c("Chromosome_liftOver", "Strand_liftOver", "Start_liftOver", "End_liftOver")
locations_liftOver_vec <- MakeLocationStrings(setNames(lifted_CRISPRko_df[have_locations, liftOver_columns],
                                                       c("Chromosome", "Strand", "Start", "End")
                                                       )
                                              )

locations_liftOver_vec[!(are_identical_sequences)] <- NA_character_
CRISPRko_df[, "Location_liftOver"] <- NA_character_
CRISPRko_df[have_locations, "Location_liftOver"] <- locations_liftOver_vec

CRISPRko_df[have_locations, "Original_PAM"] <- ifelse(are_identical_sequences & is.na(lifted_CRISPRko_df[have_locations, "Original_PAM"]),
                                                      lifted_CRISPRko_df[have_locations, "PAM_liftOver"],
                                                      lifted_CRISPRko_df[have_locations, "Original_PAM"]
                                                      )




# Remove new duplicates after correcting ambiguous Entrez IDs -------------

CRISPRko_df <- ResolveDuplicates(CRISPRko_df, concatenate_columns = "TKOv3_ID")




# Save data ---------------------------------------------------------------

save(list = "CRISPRko_df", file = file.path(CRISPRko_RData_directory, "03) Disambiguate gene IDs by mapping hg19 genomic coordinates to genes - CRISPRko_df.RData"))























