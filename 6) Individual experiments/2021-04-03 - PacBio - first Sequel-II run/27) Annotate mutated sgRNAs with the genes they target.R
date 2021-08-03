### 23rd July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts/1) R functions")
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(library_RData_directory, "1) General")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(sql2_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData"))


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

names(mutations_df)[names(mutations_df) == "Loci_0MM"] <- "Locations_0MM"

for (column_name in c("Chromosome", "Strand", "Start", "End")) {
  mutations_df[[column_name]] <- NA
}

are_TSS_based <- mutations_df[["Modality"]] %in% c("CRISPRa", "CRISPRi")
are_CRISPRko <- mutations_df[["Modality"]] == "CRISPRko"
have_0MM <- mutations_df[["Num_0MM"]] >= 1

TSS_mutations_df <- mutations_df[are_TSS_based & have_0MM, ]
CRISPRko_mutations_df <- mutations_df[are_CRISPRko & have_0MM, ]
row.names(TSS_mutations_df) <- NULL
row.names(CRISPRko_mutations_df) <- NULL




# Identify targeted genes -------------------------------------------------

nearby_list <- AlignSummaryDf(FindNearbyTSSs,
                              TSS_mutations_df,
                              all_TSS_df
                              )

CDS_or_exons_list <- AlignSummaryDf(FindOverlappingGenes,
                                    CRISPRko_mutations_df,
                                    CDS_or_exon_locations_df
                                    )

goo


# # Save data ---------------------------------------------------------------
#
# save(list = "mutations_df",
#      file = file.path(sql2_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData")
#      )
#

