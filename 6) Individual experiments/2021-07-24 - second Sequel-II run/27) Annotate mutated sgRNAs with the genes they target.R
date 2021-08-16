### 15th August 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts/1) R functions")
experiments_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory            <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory       <- file.path(plate1_directory, "1) R functions")

source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # For TruncateLongEntries
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))

source(file.path(R_functions_directory, "24) Finding unintended targets of mutated gRNAs.R"))



# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(library_RData_directory, "1) General")
output_directory         <- file.path(s2r2_directory, "5) Output", "Tables", "Targets of mutated gRNAs")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2r2_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData"))

load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))






# Add gene information to the mutations_df data frame ---------------------

matches_vec <- match(mutations_df[, "Combined_ID"], library_df[, "Combined_ID"])

stopifnot(!(anyNA(matches_vec)))

mutations_df <- data.frame(
  mutations_df,
  library_df[matches_vec, c("Modality", "Target_ID" ,"TSS_number", "Entrez_ID", "Gene_symbol")],
  stringsAsFactors = FALSE,
  row.names = NULL
)




# Prepare the mutations data frames ---------------------------------------

TSS_mutations_df <- AnnotateMutations(mutations_df, c("CRISPRa", "CRISPRi"),
                                      FindNearbyTSSs, all_TSS_df
                                      )


CRISPRko_mutations_df <- AnnotateMutations(mutations_df, "CRISPRko",
                                           FindOverlappingGenes,
                                           CDS_or_exon_locations_df
                                           )




# Export data -------------------------------------------------------------

ExportMutatedDf(TSS_mutations_df, "Targets_of_mutated_gRNAs__CRISPRa_and_CRISPRi")
ExportMutatedDf(CRISPRko_mutations_df, "Targets_of_mutated_gRNAs__CRISPRko")





# Save data ---------------------------------------------------------------

save(list = c("TSS_mutations_df", "CRISPRko_mutations_df"),
     file = file.path(s2r2_R_objects_directory, "27) Annotate mutated sgRNAs with the genes they target.RData")
     )





