### 11th January 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R")) # For Are4sg
source(file.path(general_functions_directory, "15) Finding non-overlapping sgRNAs.R")) # For CheckThatFactorIsInOrder




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
                                                        CDS_or_exon_locations_df
                                                        )
deletions_exon_protein_list <- FindOverlapsWithDeletions(merged_CRISPRko_df[are_4sg, ],
                                                         exon_locations_df
                                                         )




# Determine which sgRNAs can be analyzed ----------------------------------

are_to_replace <- (merged_CRISPRko_df[["Entrez_ID"]] %in% "7795") &
                  (merged_CRISPRko_df[["Gene_symbol"]] %in% "MEMO1")
merged_CRISPRko_df[["Entrez_ID"]][are_to_replace] <- "51072"


have_single_entrez <- !(is.na(merged_CRISPRko_df[["Entrez_ID"]])) &
                      !(grepl(",", merged_CRISPRko_df[["Entrez_ID"]], fixed = TRUE))

indices_vec <- rep(NA, nrow(merged_CRISPRko_df))
indices_vec[have_single_entrez] <- seq_len(sum(have_single_entrez))




# Find overlapping genes for sgRNAs ---------------------------------------

guides_CDS_list <- AlignSummaryDf(FindOverlappingGenes,
                                  merged_CRISPRko_df[have_single_entrez, ],
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



# Prepare for saving summary data frames ----------------------------------

deletions_CDS_df          <- deletions_CDS_list[["summary_df"]]
deletions_exon_df         <- deletions_exon_list[["summary_df"]]
deletions_CDS_protein_df  <- deletions_CDS_list[["summary_df"]]
deletions_exon_protein_df <- deletions_exon_list[["summary_df"]]

guides_CDS_df             <- guides_CDS_list[["summary_df"]]
guides_exon_df            <- guides_CDS_list[["summary_df"]]
guides_CDS_protein_df     <- guides_CDS_list[["summary_df"]]
guides_exon_protein_df    <- guides_CDS_list[["summary_df"]]



# Save data ---------------------------------------------------------------

save(list = c("deletions_CDS_df", "deletions_exon_df",
              "deletions_CDS_protein_df", "deletions_exon_protein_df",
              "guides_CDS_df", "guides_exon_df",
              "guides_CDS_protein_df", "guides_exon_protein_df"
              ),
     file = file.path(CRISPRko_RData_directory, "16) Find all genes targeted by each sgRNA.RData")
     )







