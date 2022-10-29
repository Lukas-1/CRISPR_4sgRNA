### 2022-09-14


# Load packages and source code -------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R"))



# Define paths ------------------------------------------------------------

first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
first_rdata_dir          <- file.path(first_illumina_trial_dir, "03_R_objects")
rdata_dir_2sg            <- file.path(experiments_directory, "2022-06-21 - Illumina paired-end 2sg - correct reference", "03_R_objects")
rdata_dir_4sg            <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing", "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir_2sg, "03_disambiguate_dJR072_CRISPRoff_library.RData"))
load(file.path(rdata_dir_4sg, "09_create_figures_from_count_data.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__2020Q2_gene_lists.RData"))
load(file.path(first_rdata_dir, "05_compile_data_on_essential_genes__essential_df.RData"))




# Identify genes common to both datasets ----------------------------------

essential_entrezs_2sg     <- GetAvailableGenes(essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)
non_essential_entrezs_2sg <- GetAvailableGenes(non_essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)

CRISPRoff_df <- sg_CRISPRoff_df
essential_entrezs_4sg     <- GetAvailableGenes(essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)
non_essential_entrezs_4sg <- GetAvailableGenes(non_essentials_2020Q2_df[, "Entrez_ID"], min_count = 0)

intersect_essential_entrezs <- intersect(essential_entrezs_2sg, essential_entrezs_4sg)
intersect_non_essential_entrezs <- intersect(non_essential_entrezs_2sg, non_essential_entrezs_4sg)

length(intersect_essential_entrezs)
length(intersect_non_essential_entrezs)




# Save data ---------------------------------------------------------------

save(list = c("intersect_essential_entrezs", "intersect_non_essential_entrezs"),
     file = file.path(rdata_dir_4sg, "10_identify_genes_present_in_both_libraries.RData")
     )





