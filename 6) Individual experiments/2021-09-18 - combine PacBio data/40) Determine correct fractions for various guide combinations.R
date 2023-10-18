### 25th October 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"



# Define folder paths -----------------------------------------------------

experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(library_RData_directory, "2) CRISPRa")

file_output_directory    <- file.path(s2rC_directory, "5) Output", "Tables", "Unique combinations")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rC_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))

load(file.path(CRISPRa_RData_directory, "31) Look for non-unique genes and plasmids.RData"))




# Define functions --------------------------------------------------------

CountForCombo <- function(use_reads_df, use_wells_fac, combo_vec) {
  columns_vec <- paste0("sg", combo_vec, "_cr", combo_vec)
  numeric_mat <- as.matrix(use_reads_df[, columns_vec, drop = FALSE])
  stopifnot(!(anyNA(numeric_mat)))
  are_all_correct <- rowSums(numeric_mat == 1) == length(combo_vec)
  results_vec <- tapply(are_all_correct, use_wells_fac, sum)
  results_vec[is.na(results_vec)] <- 0L
  return(results_vec)
}




# Determine the fraction of correct reads for various sg combos -----------

reads_df <- ccs7_df_list[["individual_reads_df"]]

pass_bc   <- reads_df[["Passes_barcode_filters"]] == 1
pass_read <- reads_df[["Passes_read_quality"]] == 1

reads_df <- reads_df[pass_bc & pass_read, ]
row.names(reads_df) <- NULL

wells_fac <- factor(reads_df[["Combined_ID"]],
                    levels = unique(library_df[["Combined_ID"]])
                    )

all_4_vec <- CountForCombo(reads_df, wells_fac, 1:4)
stopifnot(identical(as.integer(all_4_vec),
                    ccs7_df_list[["filtered_summary_df"]][["Count_all_4"]]
                    )
          )

single_sg_mat <- do.call(cbind, lapply(1:4, function(x) CountForCombo(reads_df, wells_fac, x)))
colnames(single_sg_mat) <- paste0("sg", 1:4)

combn_2sg_mat <- combn(1:4, 2)
combos_2sg_mat <- apply(combn_2sg_mat, 2, function(x) CountForCombo(reads_df, wells_fac, x))
colnames(combos_2sg_mat) <- paste0("sg", combn_2sg_mat[1, ], "_sg", combn_2sg_mat[2, ])

wells_df <- data.frame(
  "Combined_ID" = levels(wells_fac),
  "Num_reads"   = tabulate(wells_fac),
  "All_4"       = all_4_vec,
  single_sg_mat,
  combos_2sg_mat,
  stringsAsFactors = FALSE
)

wells_fraction_mat <- apply(as.matrix(wells_df[, 3:ncol(wells_df)]), 2, function(x) x / wells_df[["Num_reads"]])




# Add summary statistics on % correct to unique_combos_df -----------------

unique_combos_df <- data.frame(
  unique_combos_df,
  "Fraction_correct" = colMeans(wells_fraction_mat, na.rm = TRUE)
)



# Export data -------------------------------------------------------------

write.table(unique_combos_df,
            file = file.path(file_output_directory, "Unique features.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
            )




