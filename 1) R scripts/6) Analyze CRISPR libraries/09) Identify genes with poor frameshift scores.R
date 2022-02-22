### 15th February 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))
source(file.path(general_functions_directory, "20) Randomly allocating sgRNAs to plate layouts.R")) # For CRISPRaAreTop4Mat
source(file.path(general_functions_directory, "32) Compiling data on essential genes.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
analysis_RData_directory <- file.path(RData_directory, "14) Analysis of CRISPR libraries")
annotation_intermediate_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Annotation")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Analysis", "Frameshift scores")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))

load(file.path(CRISPRko_RData_directory, "13) Summarize the human transcription factor sub-library - TF_overview_df.RData"))
load(file.path(CRISPRko_RData_directory, "14) Summarize the human secretome sub-library.RData"))
CRISPRko_TF_overview_df        <- TF_overview_df
CRISPRko_secretome_overview_df <- secretome_overview_df

load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))

load(file.path(analysis_RData_directory, "08) Annotate genes for which good antibodies exist.RData"))





# Define the set of included genes ----------------------------------------

TF_entrezs <- CRISPRko_TF_overview_df[["Entrez_ID"]][!(is.na(CRISPRko_TF_overview_df[["Num_total"]]))]
secretome_entrezs <- CRISPRko_secretome_overview_df[["Entrez_ID"]][!(is.na(CRISPRko_secretome_overview_df[["Num_total"]]))]
all_entrezs <- unique(c(collected_entrez_IDs, TF_entrezs, secretome_entrezs))





# Identify chosen sgRNAs --------------------------------------------------

all_CRISPRko_df <- merged_CRISPRko_df[merged_CRISPRko_df[["Entrez_ID"]] %in% all_entrezs, ]

CRISPRko_are_top4_mat <- CRISPRkoAreTop4Mat(all_CRISPRko_df)
are_core <- grepl("Brunello|TKOv3", all_CRISPRko_df[, "Source"])
are_4sg <- CRISPRko_are_top4_mat[, "Are_final_4sg"]
are_chosen <- (are_core | are_4sg) &
              (all_CRISPRko_df[, "Entrez_ID"] %in% all_CRISPRko_df[, "Entrez_ID"][are_4sg])

core_df <- all_CRISPRko_df[are_chosen, ]
core_df[, "Is_4sg"] <- CRISPRko_are_top4_mat[are_chosen, "Are_final_4sg"]

use_columns <- c("Entrez_ID", "Gene_symbol", "Source", "Rank", "Is_4sg",
                 "CRISPOR_out_of_frame", "CRISPOR_lindel_score",
                 "Rule_set_2_score", "GuideScan_efficiency", "CRISPOR_Doench_efficacy",
                 "CRISPOR_Moreno_Mateos",
                 "GuideScan_specificity", "CRISPOR_3MM_specificity"
                 )

core_df <- core_df[, use_columns]



# Rank sgRNA combinations by the likelihood of a frameshift ---------------

genes_fac <- factor(core_df[, "Entrez_ID"], levels = unique(core_df[, "Entrez_ID"]))
GetMax <- function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE)

max_all_oof <- tapply(core_df[, "CRISPOR_out_of_frame"], genes_fac, GetMax)
max_4sg_oof <- tapply(core_df[, "CRISPOR_out_of_frame"][core_df[, "Is_4sg"]], genes_fac[core_df[, "Is_4sg"]], GetMax)
max_all_lindel <- tapply(core_df[, "CRISPOR_lindel_score"], genes_fac, GetMax)
max_4sg_lindel <- tapply(core_df[, "CRISPOR_lindel_score"][core_df[, "Is_4sg"]], genes_fac[core_df[, "Is_4sg"]], GetMax)

core_df[, "Max_out_of_frame_4sg"] <- max_4sg_oof[core_df[, "Entrez_ID"]]
core_df[, "Max_out_of_frame_all"] <- max_all_oof[core_df[, "Entrez_ID"]]

core_df[, "Max_lindel_4sg"] <- max_4sg_lindel[core_df[, "Entrez_ID"]]
core_df[, "Max_lindel_all"] <- max_all_lindel[core_df[, "Entrez_ID"]]






# Integrate surfaceome data -----------------------------------------------

matches_vec <- match(core_df[, "Entrez_ID"], all_genes_df[, "Entrez_ID"])

exclude_columns <- c(names(all_genes_df)[1:4], "Antibody_RRIDs")

core_df <- data.frame(core_df,
                      all_genes_df[matches_vec, !(names(all_genes_df) %in% exclude_columns)],
                      stringsAsFactors = FALSE
                      )




# Re-order genes ----------------------------------------------------------

new_order <- order(!(core_df[, "CSC_HEK_2015"] %in% "Yes"),
                   core_df[, "Max_out_of_frame_all"]
                   )

core_df <- core_df[new_order, ]
row.names(core_df) <- NULL






# Create a data frame with one row per gene -------------------------------

are_gene_specific <- vapply(names(core_df), function(x) {
  all(tapply(core_df[, x], core_df[, "Entrez_ID"], function(y) length(unique(y)) == 1))
}, logical(1))

genes_df <- core_df[!(duplicated(core_df[, "Entrez_ID"])), names(core_df)[are_gene_specific]]




# Prepare for export ------------------------------------------------------

core_df[, "Color"] <- OnesAndZeros(core_df[["Entrez_ID"]]) + 1L




# Export the integrated data frame ----------------------------------------

write.table(ReplaceEssentialNAs(genes_df),
            file = file.path(file_output_directory, "frameshift_genes.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t"
            )
exclude_columns <- setdiff(names(all_genes_df), names(all_CRISPRko_df))
write.table(core_df[, !(names(core_df) %in% exclude_columns)],
            file = file.path(file_output_directory, "frameshift_sgRNAs.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t", na = ""
            )




