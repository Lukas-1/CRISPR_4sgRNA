### 8th March 2021 ###




# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko", "Viability test")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "22) Compile data on essential genes.RData"))
load(file.path(general_RData_directory, "23) Compile data on protein localization.RData"))
load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))





# Categorize by deletion size ---------------------------------------------

deletion_categories <- c("Very small (< 1 kb)", "Small (1-10 kb)", "Large (10-100 kb)", "Very large (> 100 kb)")
deletion_vec <- ifelse(sgRNAs_overview_df[["Max_deletion_size"]] <= 999,
                       deletion_categories[[1]],
                       ifelse(sgRNAs_overview_df[["Max_deletion_size"]] <= 9999,
                              deletion_categories[[2]],
                              ifelse(sgRNAs_overview_df[["Max_deletion_size"]] <= 99999,
                                     deletion_categories[[3]],
                                     deletion_categories[[4]]
                                     )
                              )
                       )
deletion_fac <- factor(deletion_vec,
                       levels = deletion_categories,
                       ordered = TRUE
                       )




# Categorize by essentiality ----------------------------------------------

essential_fac <- essential_df[["Four_categories"]]

are_NE <- essential_fac %in% "Non-essential"

essential_fac[are_NE] <- ifelse((essential_df[["CRISPR_mean_effect"]][are_NE] >= 0) &
                                (essential_df[["CRISPR_mean_effect"]][are_NE] < 0.05),
                                "Non-essential", "Other"
                                )





# Define eligible genes ---------------------------------------------------

essential_matches <- match(sgRNAs_overview_df[["Entrez_ID"]],
                           essential_df[["Entrez_ID"]]
                           )
are_eligible <- (sgRNAs_overview_df[["In_4sg_library"]] %in% "Yes") &
                (essential_fac[essential_matches] %in% c("Non-essential", "Essential", "Intermediate")) &
                (!(is.na(sgRNAs_overview_df[["Max_deletion_size"]])))





# Retrieve plate coordinates ----------------------------------------------

plate_matches <- match(sgRNAs_overview_df[["Entrez_ID"]][are_eligible],
                       full_4sg_by_well_df[["Entrez_ID"]]
                       )
plate_string <- sapply(strsplit(full_4sg_by_well_df[["Plate_string"]][plate_matches], "_sg", fixed = TRUE), "[", 1)





# Collect potential genes for a viability trial ---------------------------

selected_columns <- c("Entrez_ID", "Gene_symbol", "Max_deletion_size",
                      "Min_deletion_size"
                      )

include_essential_columns <- c(
  "CRISPR_num_essential", "CRISPR_num_cell_lines",
  "CRISPR_mean_effect",
  "DEMETER2_num_essential", "DEMETER2_num_cell_lines"
)

viability_df <- data.frame(
  "Plate"                 = plate_string,
  "Well"                  = full_4sg_by_well_df[["Well_number"]][plate_matches],
  sgRNAs_overview_df[are_eligible, selected_columns],
  "Max_deletion_category" = droplevels(deletion_fac[are_eligible]),
  "Essential_category"    = droplevels(essential_fac[essential_matches][are_eligible]),
  essential_df[essential_matches, include_essential_columns][are_eligible, ],
  stringsAsFactors        = FALSE
)

comb_fac <- interaction(viability_df[["Max_deletion_category"]],
                        viability_df[["Essential_category"]],
                        sep = "__",
                        lex.order = TRUE
                        )
new_order <- order(comb_fac, as.integer(viability_df[["Entrez_ID"]]))

viability_df[["Combined_category"]] <- comb_fac

viability_df <- viability_df[new_order, ]
row.names(viability_df) <- NULL




# Incorporate data from the Human Protein Atlas ---------------------------

stopifnot(!(anyNA(viability_df[["Entrez_ID"]])))
matches_vec <- match(viability_df[["Entrez_ID"]], HPA_df[["Entrez_ID"]])

select_columns <- c("Reliability_IH", "Reliability_IF", "Subcellular_location", "Antibody_RRIDs")
viability_df <- data.frame(viability_df,
                           HPA_df[matches_vec, select_columns],
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )



# Select random genes -----------------------------------------------------

viability_df[["Randomized_rank"]] <- NA_integer_
set.seed(1)
for (use_level in levels(viability_df[["Combined_category"]])) {
  are_this_level <- viability_df[["Combined_category"]] == use_level
  random_ranks <- sample(seq_len(sum(are_this_level)))
  # use_order <- order(viability_df[["Reliability_IF"]][are_this_level],
  #                    viability_df[["Reliability_IH"]][are_this_level],
  #                    -(random_ranks),
  #                    decreasing = TRUE
  #                    )
  use_order <- order(random_ranks)
  randomized_rank <- order(use_order)
  viability_df[["Randomized_rank"]][are_this_level] <- randomized_rank
}




# Re-order randomized genes -----------------------------------------------

shuffled_order <- order(viability_df[["Combined_category"]],
                        viability_df[["Randomized_rank"]]
                        )
shuffled_viability_df <- viability_df[shuffled_order, ]
row.names(shuffled_viability_df) <- NULL


are_sel_28 <- ifelse(shuffled_viability_df[["Essential_category"]] %in% c("Essential", "Intermediate"),
                     shuffled_viability_df[["Randomized_rank"]] == 1,
                     ifelse(shuffled_viability_df[["Max_deletion_category"]] %in% deletion_categories[c(1, 4)],
                            shuffled_viability_df[["Randomized_rank"]] %in% 1:4,
                            shuffled_viability_df[["Randomized_rank"]] %in% 1:6
                            )
                     )
viability_28_genes_df <- shuffled_viability_df[are_sel_28, ]





# Export data -------------------------------------------------------------

ExportViabilityDf <- function(viab_df, file_name) {
  viab_df[["Color"]] <- as.integer(viab_df[["Combined_category"]])
  viab_df <- viab_df[, names(viab_df) != "Combined_category"]

  for (column_name in names(viab_df)) {
    viab_df[, column_name] <- ifelse(is.na(viab_df[, column_name]),
                                     " ",
                                     as.character(viab_df[, column_name])
                                     )
  }

  write.table(viab_df,
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )

}

ExportViabilityDf(shuffled_viability_df, "All eligible genes")
ExportViabilityDf(viability_28_genes_df, "28 random genes")






