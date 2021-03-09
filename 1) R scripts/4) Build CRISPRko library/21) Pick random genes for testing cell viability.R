### 8th March 2021 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRko", "Viability test")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "22) Compile data on essential genes.RData"))
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

full_4sg_by_well_df



# Categorize by essentiality ----------------------------------------------

essential_fraction <- essential_df[["Achilles_num_essential"]] /
                      essential_df[["Achilles_num_cell_lines"]]
essential_vec <- ifelse(essential_fraction < 0.1,
                        "Non-essential",
                        ifelse(essential_fraction > 0.9,
                               "Essential",
                               "Indeterminate"
                               )
                        )
essential_fac <- factor(essential_vec,
                        levels = c("Essential", "Indeterminate", "Non-essential"),
                        ordered = TRUE
                        )

table(essential_fac, essential_df[["Achilles_common"]])




# Define eligible genes ---------------------------------------------------

essential_matches <- match(sgRNAs_overview_df[["Entrez_ID"]],
                           essential_df[["Entrez_ID"]]
                           )
are_eligible <- (sgRNAs_overview_df[["In_4sg_library"]] %in% "Yes") &
                (essential_fac[essential_matches] %in% c("Non-essential", "Essential")) &
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

viability_df <- data.frame(
  "Plate"                 = plate_string,
  "Well"                  = full_4sg_by_well_df[["Well_number"]][plate_matches],
  sgRNAs_overview_df[are_eligible, selected_columns],
  "Max_deletion_category" = droplevels(deletion_fac[are_eligible]),
  "Essential_category"    = droplevels(essential_fac[essential_matches][are_eligible]),
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




# Select random genes -----------------------------------------------------

viability_df[["Randomized_rank"]] <- NA_integer_
set.seed(1)
for (use_level in levels(viability_df[["Combined_category"]])) {
  are_this_level <- viability_df[["Combined_category"]] == use_level
  viability_df[["Randomized_rank"]][are_this_level] <- sample(seq_len(sum(are_this_level)))
}




# Re-order randomized genes -----------------------------------------------

shuffled_order <- order(viability_df[["Combined_category"]],
                        viability_df[["Randomized_rank"]]
                        )
shuffled_viability_df <- viability_df[shuffled_order, ]
row.names(shuffled_viability_df) <- NULL

are_sel_20 <- ifelse(shuffled_viability_df[["Max_deletion_category"]] %in% deletion_categories[c(1, 4)],
                     shuffled_viability_df[["Randomized_rank"]] %in% 1:2,
                     shuffled_viability_df[["Randomized_rank"]] %in% 1:3
                     )
viability_20_genes_df <- shuffled_viability_df[are_sel_20, ]




# Export data -------------------------------------------------------------

ExportViabilityDf <- function(viab_df, file_name) {
  viab_df[["Color"]] <- as.integer(viab_df[["Combined_category"]])
  viab_df <- viab_df[, names(viab_df) != "Combined_category"]
  write.table(viab_df,
              file = file.path(file_output_directory, paste0(file_name, ".tsv")),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
              )

}

ExportViabilityDf(shuffled_viability_df, "All eligible genes")
ExportViabilityDf(viability_20_genes_df, "20 random genes")






