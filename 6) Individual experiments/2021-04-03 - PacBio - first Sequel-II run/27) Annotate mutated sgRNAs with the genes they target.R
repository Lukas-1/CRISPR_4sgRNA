### 23rd July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts/1) R functions")
plate1_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory       <- file.path(plate1_directory, "1) R functions")

source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # For TruncateLongEntries
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))

source(file.path(R_functions_directory, "24) Finding unintended targets of mutated gRNAs.R"))



# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(library_RData_directory, "1) General")
output_directory         <- file.path(sql2_directory, "5) Output", "Tables", "Targets of mutated gRNAs")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(sql2_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))


load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))





# Identify contaminated ZMWs ----------------------------------------------

use_df <- ccs7_df_list[["individual_reads_df"]]
any_aligned_contam <- apply(as.matrix(use_df[, paste0("sg", 1:4, "_category")]) == "Contamination",
                            1,
                            any
                            )
any_search_contam <- use_df[, "Num_contaminating_guides"] > 0
are_contam_zmws <- use_df[["ZMW"]][any_search_contam | any_aligned_contam]






# Add gene information to the mutations_df data frame ---------------------

matches_vec <- match(mutations_df[, "Combined_ID"], library_df[, "Combined_ID"])

stopifnot(!(anyNA(matches_vec)))

mutations_df <- AddContamZMW(ccs7_df_list, mutations_df)

mutations_df <- data.frame(
  mutations_df,
  library_df[matches_vec, c("Modality", "Target_ID" ,"TSS_number", "Entrez_ID", "Gene_symbol")],
  stringsAsFactors = FALSE,
  row.names = NULL
)





# Prepare the mutations data frames ---------------------------------------

TSS_mut_list <- AnnotateMutations(mutations_df, c("CRISPRa", "CRISPRi"),
                                  FindNearbyTSSs, all_TSS_df
                                  )


CRISPRko_mut_list <- AnnotateMutations(mutations_df, "CRISPRko",
                                       FindOverlappingGenes,
                                       CDS_or_exon_locations_df
                                       )




# Export data on individual reads -----------------------------------------

ExportMutatedDf(TSS_mut_list[["annotated_df"]],
                "Targets_of_mutated_gRNAs__CRISPRa_and_CRISPRi"
                )
ExportMutatedDf(CRISPRko_mut_list[["annotated_df"]],
                "Targets_of_mutated_gRNAs__CRISPRko"
                )





# Draw doughnut/bar plots -------------------------------------------------

pdf(file = file.path(output_directory, "Donut charts - targets of mutated gRNAs.pdf"),
    width = 6, height = 4
    )

MutationsDonutBar(CRISPRko_mut_list[["gRNA_numbers"]],
                  "First run \u2013 CRISPRko (target: CDS or exon)"
                  )

MutationsDonutBar(TSS_mut_list[["gRNA_numbers"]],
                  "First run \u2013 CRISPRa/i (target: within 1000 bp of the TSS)"
                  )

dev.off()






# Save data ---------------------------------------------------------------

save(list = c("TSS_mut_list", "CRISPRko_mut_list"),
     file = file.path(sql2_R_objects_directory, "27) Annotate mutated sgRNAs with the genes they target.RData")
     )




