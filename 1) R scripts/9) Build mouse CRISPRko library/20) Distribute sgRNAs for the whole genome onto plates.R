### 21st July 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "6) Mouse - General")
CRISPRko_RData_directory <- file.path(RData_directory, "8) Mouse - CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Mouse - CRISPRko")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "03) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "06) Read in gene lists.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))





# Create an empty SNP column ----------------------------------------------

merged_CRISPRko_df[[preferred_AF_max_column]] <- NA_real_




# Assign all guides to plates ---------------------------------------------

sg4_by_gene_df <- AllocateAllGuidesToPlates(merged_CRISPRko_df,
                                            sublibraries_entrezs_list = list("All" = collected_entrez_IDs),
                                            num_control_wells = 96,
                                            reorder_df = FALSE
                                            )
sg4_by_gene_df[["Old_order"]] <- seq_len(nrow(sg4_by_gene_df))

sg4_by_well_df <- ReorderPlates(sg4_by_gene_df)
sg4_by_well_df <- RenumberPlatesContinuously(sg4_by_well_df)
sg4_by_well_df <- AssignPlateStrings(sg4_by_well_df)

sg4_by_gene_df <- sg4_by_well_df[order(sg4_by_well_df[["Old_order"]]), ]

sg4_by_gene_df <- sg4_by_gene_df[, names(sg4_by_gene_df) != "Old_order"]
sg4_by_well_df <- sg4_by_well_df[, names(sg4_by_well_df) != "Old_order"]





# Export all shared / duplicated sgRNAs -----------------------------------

shared_sgRNAs_df <- SharedsgRNAsDf(sg4_by_gene_df)

use_sub_folder <- "Plate layout - all genes"

ExportSharedDf(shared_sgRNAs_df,
               file.path(use_sub_folder,
                         "Duplicated mouse CRISPRa sgRNAs (shared between genes)"
                         )
               )








# Export the whole-genome plate layouts -----------------------------------

ExportPlates(sg4_by_gene_df,
             "All_sublibraries_ordered_by_gene",
             use_sub_folder,
             add_primers = TRUE
             )
ExportPlates(sg4_by_well_df,
             "All_sublibraries_ordered_by_well",
             use_sub_folder,
             add_primers = TRUE
             )





# Remove the empty SNP column ---------------------------------------------

sg4_by_gene_df <- sg4_by_gene_df[, names(sg4_by_gene_df) != preferred_AF_max_column]
sg4_by_well_df <- sg4_by_well_df[, names(sg4_by_well_df) != preferred_AF_max_column]






# Save data ---------------------------------------------------------------

save(list = c("sg4_by_gene_df", "sg4_by_well_df"),
     file = file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData")
     )




