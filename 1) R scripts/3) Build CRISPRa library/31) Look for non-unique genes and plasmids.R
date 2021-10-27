### 11th October 2021 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "CRISPRa")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))




# Remove controls ---------------------------------------------------------

CRISPR_df <- sg4_by_gene_df[sg4_by_gene_df[["Is_control"]] %in% "No", ]
row.names(CRISPR_df) <- NULL




# Identify duplicated plasmids --------------------------------------------

IDs_fac <- factor(CRISPR_df[, "AltTSS_ID"],
                  levels = unique(CRISPR_df[, "AltTSS_ID"])
                  )
CheckThatFactorIsInOrder(IDs_fac)
stopifnot(as.integer(table(IDs_fac)) == 4)

stopifnot(length(unique(split(CRISPR_df[["Rank"]], IDs_fac))) == 1)

combos_list <- split(toupper(CRISPR_df[["sgRNA_sequence"]]), IDs_fac)

combos_vec <- vapply(combos_list, function(x) paste0(x, collapse = "_"), "")

combos_table <- table(combos_vec)
stopifnot(sum(table(combos_table)) == length(unique(combos_list)))

num_occurrences <- as.integer(combos_table[combos_vec])

plasmid_index <- match(combos_vec, combos_vec)

CRISPR_df[["Num_dupl_4sg"]] <- rep(num_occurrences, each = 4)
CRISPR_df[["Index_4sg"]] <- rep(plasmid_index, each = 4)




# Identify duplicated two-guide combinations ------------------------------

two_guides_combos <- unname(as.list(data.frame(combn(4, 2))))
two_guides_list_list <- lapply(two_guides_combos, function(x) {
  lapply(combos_list, function(y) y[x])
})

names(two_guides_list_list) <- vapply(two_guides_combos,
                                      function(x) paste0("sg", x[[1]], "_sg", x[[2]]),
                                      ""
                                      )

two_guides_mat <- do.call(cbind,
                          lapply(two_guides_list_list,
                                 function(x) vapply(x, function(y) paste0(y, collapse = "_"), "")
                                 )
                          )

two_guides_num_mat <- apply(two_guides_mat, 2, function(combos_vec) {
  combos_table <- table(combos_vec)
  num_occurrences <- as.integer(combos_table[combos_vec])
})

two_guides_summary_mat <- t(apply(two_guides_num_mat, 2, function(x) table(x == 1)))
colnames(two_guides_summary_mat) <- c("Unique", "Non_unique")

two_guides_summary_df <- data.frame(
  "sg_A" = sapply(two_guides_combos, "[[", 1),
  "sg_B" = sapply(two_guides_combos, "[[", 2),
  two_guides_summary_mat,
  stringsAsFactors = FALSE
)

two_guides_summary_df[["Description"]] <- paste0("sg",
                                                 two_guides_summary_df[["sg_A"]],
                                                 " and sg",
                                                 two_guides_summary_df[["sg_B"]]
                                                 )




# Examine duplicated plasmids ---------------------------------------------

new_order <- order(-(CRISPR_df[["Num_dupl_4sg"]]), CRISPR_df[["Index_4sg"]])

selected_columns <- c("Entrez_ID", "Other_target_Entrez_IDs",
                      "Gene_symbol", "Other_target_symbols",
                      "AltTSS_ID", "TSS_ID", "TSS_number", "Allocated_TSS",
                      "Num_TSSs", "Is_main_TSS",
                      "Rank", "Spacing", "Source",
                      "sgRNA_sequence", "PAM"
                      )

head(CRISPR_df[new_order, selected_columns])




# Identify duplicated sgRNAs ----------------------------------------------

sg_num_occurrences_list <- lapply(1:4, function(x) {
  are_this_sg <- CRISPR_df[["Rank"]] %in% x
  sg_vec <- toupper(CRISPR_df[["sgRNA_sequence"]][are_this_sg])
  sg_table <- table(sg_vec)
  message(paste0("Table for sg", x))
  table_table <- table(sg_table)
  print(table(sg_table, dnn = NULL))
  message("")
  num_occurrences <- sg_table[sg_vec]
  return(as.integer(num_occurrences))
})

sg_num_occurrences_mat <- do.call(cbind, sg_num_occurrences_list)
colnames(sg_num_occurrences_mat) <- paste0("sg", 1:4)

sg_summary_mat <- t(apply(sg_num_occurrences_mat, 2, function(x) table(x == 1)))




# Combine data ------------------------------------------------------------

tables_mat <- rbind(
  table(num_occurrences == 1),
  sg_summary_mat,
  two_guides_summary_mat
)[, c(2, 1)]

colnames(tables_mat) <- c("Unique", "Non_unique")

unique_combos_df <- data.frame(
  "Element" = c("4sg combination", paste0("sg", 1:4),
                two_guides_summary_df[["Description"]]
                ),
  tables_mat,
  stringsAsFactors = FALSE,
  row.names = NULL
)



# Export data -------------------------------------------------------------

write.table(unique_combos_df,
            file = file.path(file_output_directory, "Unique plasmids", "Unique features.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
            )



# Save data ---------------------------------------------------------------

save("unique_combos_df",
     file = file.path(CRISPRa_RData_directory, "31) Look for non-unique genes and plasmids.RData")
     )


