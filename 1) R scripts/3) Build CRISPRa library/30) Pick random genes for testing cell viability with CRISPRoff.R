### 8th March 2021 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRoff", "CRISPRa guides", "Viability test")
CRISPRoff_library_path   <- file.path(CRISPR_input_directory, "Gene lists", "CRISPRoff",
                                      "2021 - Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing - Table S3.xlsx"
                                      )



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "22) Compile data on essential genes - essential_df.RData"))
load(file.path(general_RData_directory, "23) Compile data on protein localization.RData"))
load(file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))



# Read in data ------------------------------------------------------------

Nunez_gRNAs_df <- data.frame(read_excel(CRISPRoff_library_path,
                                        sheet = "GW protospacers",
                                        skip = 3
                                        ),
                       stringsAsFactors = FALSE
                       )

Nunez_scores_df <- data.frame(read_excel(CRISPRoff_library_path,
                                         sheet = "Phenotype scores",
                                         skip = 3
                                         ),
                              stringsAsFactors = FALSE,
                              check.names = FALSE
                              )



# Explore the CRISPRoff libray of Nunez et al. ----------------------------

table(Nunez_gRNAs_df[, "source_A"])
table(Nunez_gRNAs_df[, "source_B"])
table(Nunez_gRNAs_df[, "representation"])

head(sort(table(Nunez_gRNAs_df[, "transcript_A"]), decreasing = TRUE), 10)
head(sort(table(Nunez_gRNAs_df[, "transcript_B"]), decreasing = TRUE), 10)

table(Nunez_gRNAs_df[, "source_A"], Nunez_gRNAs_df[, "representation"])

table(Nunez_scores_df[, "hit_and_run_hit"])
table(Nunez_scores_df[, "mutant_hit"])
table(Nunez_scores_df[, "hit_and_run_hit"], Nunez_scores_df[, "mutant_hit"])
table(Nunez_scores_df[, "CpG"])
table(Nunez_scores_df[, "source_A"])





# Categorize by essentiality ----------------------------------------------

essential_fac <- essential_df[["Four_categories"]]

are_NE <- essential_fac %in% "Non-essential"
are_E  <- essential_fac %in% "Essential"

essential_fac[are_NE] <- ifelse((essential_df[["CRISPR_mean_effect"]][are_NE] >= 0) &
                                (essential_df[["CRISPR_mean_effect"]][are_NE] < 0.05),
                                "Non-essential", "Other"
                                )

essential_fraction <- essential_df[["CRISPR_num_essential"]] /
                      essential_df[["CRISPR_num_cell_lines"]]

essential_fac[are_E] <- ifelse(essential_fraction[are_E] > 0.998,
                               "Essential",
                               "Other"
                               )




# Define eligible genes ---------------------------------------------------

essential_matches <- match(sgRNAs_overview_df[["Entrez_ID"]],
                           essential_df[["Entrez_ID"]]
                           )
are_eligible <- (sgRNAs_overview_df[["In_4sg_library"]] %in% "Yes") &
                !(is.na(essential_matches))




# Retrieve plate coordinates ----------------------------------------------

plate_matches <- match(sgRNAs_overview_df[["Entrez_ID"]][are_eligible],
                       full_4sg_by_well_df[["Entrez_ID"]]
                       )
plate_string <- sapply(strsplit(full_4sg_by_well_df[["Plate_string"]][plate_matches], "_sg", fixed = TRUE), "[", 1)




# Integrate the CRISPRa library with gene essentiality data ---------------

select_columns <- c("Entrez_ID", "Gene_symbol", "Num_transcripts")

include_essential_columns <- c(
  "CRISPR_num_essential", "CRISPR_num_cell_lines",
  "CRISPR_mean_effect",
  "DEMETER2_num_essential", "DEMETER2_num_cell_lines"
)

integrated_df <- data.frame(
  "Plate"                 = plate_string,
  "Well"                  = full_4sg_by_well_df[["Well_number"]][plate_matches],
  sgRNAs_overview_df[are_eligible, select_columns],
  "Essential_category"    = essential_fac[essential_matches][are_eligible],
  essential_df[essential_matches, include_essential_columns][are_eligible, ],
  stringsAsFactors        = FALSE
)

names(integrated_df)[names(integrated_df) == "Num_transcripts"] <- "Num_TSSs"




# Incorporate data from the Human Protein Atlas ---------------------------

stopifnot(!(anyNA(integrated_df[["Entrez_ID"]])))
matches_vec <- match(integrated_df[["Entrez_ID"]], HPA_df[["Entrez_ID"]])

select_columns <- "Subcellular_location" #"Reliability_IH", "Reliability_IF", "Antibody_RRIDs"
integrated_df <- data.frame(integrated_df,
                            HPA_df[matches_vec, select_columns, drop = FALSE],
                            stringsAsFactors = FALSE
                            )

new_order <- order(as.integer(integrated_df[["Entrez_ID"]]))

integrated_df <- integrated_df[new_order, ]
row.names(integrated_df) <- NULL





# Integrate the CRISPRoff data from Nunez et al. --------------------------

Nunez_scores_mapped_df <- MapToEntrezs(symbols_vec = Nunez_scores_df[, "gene"])

matches_vec <- match(integrated_df[["Entrez_ID"]], Nunez_scores_mapped_df[["Entrez_ID"]])

select_columns <- c("CRISPRoff_average", "hit_and_run_hit", "CpG")

integrated_df <- data.frame(
  integrated_df,
  Nunez_scores_mapped_df[matches_vec, "Original_symbol", drop = FALSE],
  Nunez_scores_df[matches_vec, select_columns],
  row.names = NULL
)

names(integrated_df)[names(integrated_df) == "hit_and_run_hit"] <- "CRISPRoff_hit"
names(integrated_df)[names(integrated_df) == "CpG"] <- "Has_CpG_island"




# Exclude ineligible genes (multiple TSSs / absent for CRISPRoff) ---------

are_eligible <- (integrated_df[["Num_TSSs"]] == 1) &
                (!(is.na(matches_vec)))

integrated_df <- integrated_df[are_eligible, ]
row.names(integrated_df) <- NULL

table(integrated_df[["Essential_category"]])




# Collect potential genes for a viability trial ---------------------------

viability_df <- integrated_df[integrated_df[["Essential_category"]] %in% c("Non-essential", "Essential"), ]
new_order <- order(viability_df[["Essential_category"]],
                   as.integer(viability_df[["Entrez_ID"]])
                   )
viability_df <- viability_df[new_order, ]
row.names(viability_df) <- NULL




# Select random genes (within viability categories) -----------------------

viability_df[["Randomized_rank"]] <- NA_integer_
set.seed(1)
for (use_level in levels(viability_df[["Essential_category"]])) {
  are_this_level <- viability_df[["Essential_category"]] == use_level
  random_ranks <- sample(seq_len(sum(are_this_level)))
  use_order <- order(random_ranks)
  randomized_rank <- order(use_order)
  viability_df[["Randomized_rank"]][are_this_level] <- randomized_rank
}




# Re-order randomized genes (within viability categories) -----------------

shuffled_order <- order(viability_df[["Essential_category"]],
                        viability_df[["Randomized_rank"]]
                        )
shuffled_viability_df <- viability_df[shuffled_order, ]
row.names(shuffled_viability_df) <- NULL

are_sel_75 <- shuffled_viability_df[["Randomized_rank"]] %in% 1:75

viability_150_genes_df <- shuffled_viability_df[are_sel_75, ]
row.names(viability_150_genes_df) <- NULL

table(viability_150_genes_df[["Essential_category"]],
      viability_150_genes_df[["CRISPRoff_hit"]]
      )




# Export selected genes (based on gene essentiality) ----------------------

ExportViabilityDf <- function(viab_df, file_name) {
  viab_df[["Color"]] <- as.integer(droplevels(viab_df[["Essential_category"]]))

  for (column_name in names(viab_df)) {
    viab_df[, column_name] <- ifelse(is.na(viab_df[, column_name]) |
                                     (viab_df[, column_name] %in% ""),
                                     " ",
                                     as.character(viab_df[, column_name])
                                     )
  }

  write.table(viab_df,
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )

}

ExportViabilityDf(shuffled_viability_df, "CRISPRoff viability test - all eligible genes")
ExportViabilityDf(viability_150_genes_df, "CRISPRoff viability test - 150 genes")





