### 11th January 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R")) # For Are4sg
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))




# Find overlapping genes for 4sg combinations (predicted deletions) -------

are_4sg <- Are4sg(merged_CRISPRko_df, sublibraries_all_entrezs_list)

deletions_CDS_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                CDS_or_exon_locations_df
                                                )
deletions_exon_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                 exon_locations_df
                                                 )

deletions_CDS_protein_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                        CDS_or_exon_locations_df,
                                                        only_protein_coding = TRUE
                                                        )
deletions_exon_protein_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                         exon_locations_df,
                                                         only_protein_coding = TRUE
                                                         )




# Make corrections to prevent downstream errors ---------------------------

are_to_replace <- (merged_CRISPRko_df[["Entrez_ID"]] %in% "7795") &
                  (merged_CRISPRko_df[["Gene_symbol"]] %in% "MEMO1")
merged_CRISPRko_df[["Entrez_ID"]][are_to_replace] <- "51072"





# Find overlapping genes for sgRNAs ---------------------------------------

guides_CDS_list <- AlignSummaryDf(FindOverlappingGenes,
                                  merged_CRISPRko_df,
                                  CDS_or_exon_locations_df
                                  )
guides_exon_list <- AlignSummaryDf(FindOverlappingGenes,
                                   merged_CRISPRko_df,
                                   exon_locations_df
                                   )

guides_CDS_protein_list <- AlignSummaryDf(FindOverlappingGenes,
                                          merged_CRISPRko_df,
                                          CDS_or_exon_locations_df,
                                          only_protein_coding = TRUE
                                          )
guides_exon_protein_list <- AlignSummaryDf(FindOverlappingGenes,
                                           merged_CRISPRko_df,
                                           exon_locations_df,
                                           only_protein_coding = TRUE
                                           )




# Prepare for saving data frames ------------------------------------------

deletion_names <- grep("^deletions_.+_list$", ls(), value = TRUE)
sg_names       <- grep("^guides_.+_list$", ls(), value = TRUE)

all_names        <- c(deletion_names, sg_names)
summary_df_names <- sub("_list$", "_df", all_names)
full_df_names    <- sub("_list$", "_full_df", all_names)

for (i in seq_along(all_names)) {
  assign(summary_df_names[[i]], get(all_names[[i]])[["summary_df"]])
  assign(full_df_names[[i]], get(all_names[[i]])[["full_df"]])
}




# Save data ---------------------------------------------------------------

save(list = summary_df_names,
     file = file.path(CRISPRko_RData_directory, "16) Find all genes targeted by each sgRNA - summary data.RData")
     )

save(list = full_df_names,
     file = file.path(CRISPRko_RData_directory, "16) Find all genes targeted by each sgRNA - full data.RData")
     )







