### 20 May 2020 ###



# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Analysis")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "17) Read in additional gene lists.RData"))

load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates - PD genes.RData"))
PD_4sg_CRISPRko_df <- PD_4sg_by_well_df
load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates - PD genes.RData"))
PD_4sg_CRISPRa_df <- PD_4sg_by_well_df

rm(PD_4sg_by_well_df)
rm(PD_4sg_by_gene_df)




# Check for the availability of PD genes ----------------------------------

are_present_CRISPRa  <- PD_mapped_df[["Entrez_ID"]] %in% PD_4sg_CRISPRa_df[["Entrez_ID"]]
are_present_CRISPRko <- PD_mapped_df[["Entrez_ID"]] %in% PD_4sg_CRISPRko_df[["Entrez_ID"]]

stopifnot(identical(are_present_CRISPRa, are_present_CRISPRko))




# Check for duplicated genes ----------------------------------------------

num_occurrences_vec <- table(PD_mapped_df[["Entrez_ID"]])[PD_mapped_df[["Entrez_ID"]]]
PD_mapped_df[(num_occurrences_vec > 1) %in% TRUE, ]




# Display unavailable genes -----------------------------------------------

PD_mapped_df[!(are_present_CRISPRa), ]




# Prepare the unavailable genes for export --------------------------------

export_df <- PD_mapped_df
export_df[["Gene_symbol"]] <- ifelse(is.na(export_df[["Gene_symbol"]]),
                                     export_df[["Original_symbol"]],
                                     export_df[["Gene_symbol"]]
                                     )
export_df <- export_df[, 1:2]




# Export unavailable genes ------------------------------------------------

write.table(export_df[!(are_present_CRISPRa), ],
            file = file.path(file_output_directory, "Unavailable_PD_genes.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )













